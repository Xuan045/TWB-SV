#!/usr/bin/env python3
"""
Script: 02_add_svlen_svtype.py
Description:
    Annotates a PLINK2 .pvar file with SVLEN and SVTYPE INFO fields by matching
    variants against a pre-built SV info TSV. Outputs an annotated VCF file.

    Matching strategy (in order):
      1. Exact match on (CHROM, POS, REF, ALT)
      2. Fuzzy match: position ±1 bp with compatible alleles
      3. Inference from REF/ALT lengths (for unmatched INS/DEL only)

    Variants with abs(SVLEN) < 50 are NOT tagged with SVLEN/SVTYPE,
    as they do not meet the conventional structural variant size threshold.

================================================================================
USAGE:
    python 02_add_svlen_svtype.py <sv_tsv> <pvar_in> <vcf_out> [log_file]

ARGUMENTS:
    sv_tsv    : Tab-separated SV info file (CHROM, POS, REF, ALT, SVLEN, SVTYPE)
    pvar_in   : Input PLINK2 .pvar file
    vcf_out   : Output annotated VCF file
    log_file  : (Optional) Path to log file; defaults to stderr

================================================================================
OUTPUT:
    - Annotated VCF with SVLEN and SVTYPE added to INFO field
    - Summary counts of exact matches, fuzzy matches, and unmatched variants
"""

import sys
import re

# ================================================================================
# Argument Parsing
# ================================================================================
if len(sys.argv) > 4:
    sv_tsv, pvar_in, vcf_out, log_file = sys.argv[1:]
    log_fh = open(log_file, "w")
else:
    sv_tsv, pvar_in, vcf_out = sys.argv[1:]
    log_fh = sys.stderr


def log_message(message):
    """Write a message to the log file (or stderr if no log file provided)."""
    print(message, file=log_fh)
    log_fh.flush()

# ================================================================================
# Step 1. Load SV info TSV into lookup structures
# ================================================================================
sv_dict = {}        # (chrom, pos, ref, alt) -> (svlen, svtype)
position_index = {} # chrom -> {pos_int -> key} for fuzzy matching

with open(sv_tsv) as f:
    for ln in f:
        chrom, pos, ref, alt, svlen, svtype = ln.rstrip("\n").split("\t")
        key = (chrom, pos, ref, alt)
        sv_dict[key] = (svlen, svtype)

        # Build position index for fuzzy matching (keyed by integer position)
        pos_int = int(pos)
        if chrom not in position_index:
            position_index[chrom] = {}
        position_index[chrom][pos_int] = key

# ================================================================================
# Step 2. Helper Functions
# ================================================================================
def find_fuzzy_match(chrom, pos_str, ref, alt):
    """
    Attempt to match a variant against sv_dict using:
      - Exact match first
      - Then position ±1 bp with allele compatibility checks

    Allele compatibility:
      - General: one REF is a substring of the other AND one ALT is a substring of the other
      - Insertions: compare inserted sequences (ALT[len(REF):]) for substring overlap

    Returns the matched key from sv_dict, or None if no match found.
    """
    if chrom not in position_index:
        return None

    pos_int = int(pos_str)

    # Exact match
    exact_key = (chrom, pos_str, ref, alt)
    if exact_key in sv_dict:
        return exact_key

    # Fuzzy match: try positions ±1 bp
    for offset in [-1, 0, 1]:
        test_pos = pos_int + offset
        if test_pos not in position_index[chrom]:
            continue

        stored_key = position_index[chrom][test_pos]
        _, _, stored_ref, stored_alt = stored_key

        # General allele compatibility: substring in either direction
        if (ref in stored_ref or stored_ref in ref) and \
           (alt in stored_alt or stored_alt in alt):
            return stored_key

        # Insertion-specific check: compare inserted sequences
        if len(alt) > len(ref) and len(stored_alt) > len(stored_ref):
            inserted_seq        = alt[len(ref):]
            stored_inserted_seq = stored_alt[len(stored_ref):]
            if inserted_seq == stored_inserted_seq or \
               inserted_seq in stored_inserted_seq or \
               stored_inserted_seq in inserted_seq:
                return stored_key

    return None


def infer_sv_info(ref, alt):
    """
    Infer SVLEN and SVTYPE from REF and ALT allele lengths.
    Only applies to insertions and deletions (not SNPs or MNPs).

    Returns:
        (svlen_str, svtype_str) if INS or DEL, else (None, None)

    SVLEN convention:
        - Positive for insertions
        - Negative for deletions
    """
    ref_len = len(ref)
    alt_len = len(alt)

    if alt_len > ref_len:
        return str(alt_len - ref_len), "INS"
    elif ref_len > alt_len:
        return str(-(ref_len - alt_len)), "DEL"
    else:
        # Same length: SNP or MNP — not a structural variant
        return None, None


def is_large_enough_sv(svlen_str, threshold=50):
    """
    Check whether a variant meets the minimum SV size threshold.
    Variants with abs(SVLEN) < threshold are excluded from SVLEN/SVTYPE tagging.

    Args:
        svlen_str : SVLEN as a string (may be negative for deletions)
        threshold : Minimum absolute SVLEN to be considered an SV (default: 50)

    Returns:
        True if abs(SVLEN) >= threshold, False otherwise.
    """
    try:
        return abs(int(svlen_str)) >= threshold
    except ValueError:
        return False

# ================================================================================
# Step 3. Read .pvar and write annotated VCF
# ================================================================================
# Regex to parse symbolic ALT alleles like <DEL:SVSIZE=500> -> <DEL>
alt_re = re.compile(r'^<([A-Z]+):SVSIZE=([0-9]+)>$')

# VCF INFO header lines to inject
info_svlen  = '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length of structural variant">\n'
info_svtype = '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">\n'

matched_count       = 0
fuzzy_matched_count = 0
unmatched_count     = 0
skipped_small_count = 0  # Variants with abs(SVLEN) < 50

with open(pvar_in, "r") as fin, open(vcf_out, "wt") as fout:
    for line in fin:
        # Pass through meta-information lines
        if line.startswith("##"):
            fout.write(line)
            continue

        # Inject SVLEN/SVTYPE INFO headers before the column header line
        if line.startswith("#CHROM"):
            fout.write(info_svlen)
            fout.write(info_svtype)
            fout.write(line)
            continue

        cols = line.rstrip("\n").split("\t")
        if len(cols) < 8:
            continue

        chrom, pos, _id, ref, alt = cols[:5]

        # Simplify symbolic ALT alleles: <DEL:SVSIZE=500> -> <DEL>
        m = alt_re.match(alt)
        if m:
            alt_clean = f"<{m.group(1)}>"
            cols[4]   = alt_clean
        else:
            alt_clean = alt

        # --- Attempt match and annotation ---
        matched_key = find_fuzzy_match(chrom, pos, ref, alt_clean)

        if matched_key:
            svlen, svtype = sv_dict[matched_key]

            if not is_large_enough_sv(svlen):
                # Skip tagging: variant is too small to be classified as SV
                skipped_small_count += 1
            else:
                cols[7] = f"SVLEN={svlen};SVTYPE={svtype}" if cols[7] == "." \
                          else cols[7] + f";SVLEN={svlen};SVTYPE={svtype}"

                if matched_key == (chrom, pos, ref, alt_clean):
                    matched_count += 1
                else:
                    fuzzy_matched_count += 1
                    log_message(f"[INFO] Fuzzy match: pvar({chrom}, {pos}, {ref}, {alt_clean}) -> sv_info{matched_key}")

        else:
            # No match: attempt inference from REF/ALT lengths
            svlen, svtype = infer_sv_info(ref, alt_clean)

            if svlen is not None:
                if not is_large_enough_sv(svlen):
                    skipped_small_count += 1
                else:
                    cols[7] = f"SVLEN={svlen};SVTYPE={svtype}" if cols[7] == "." \
                              else cols[7] + f";SVLEN={svlen};SVTYPE={svtype}"
                    log_message(f"[WARNING] Not in SV info, inferred: ({chrom}, {pos}, {ref}, {alt_clean}) -> SVLEN={svlen}, SVTYPE={svtype}")
            else:
                log_message(f"[WARNING] Not in SV info, not a structural variant: ({chrom}, {pos}, {ref}, {alt_clean})")

            unmatched_count += 1

        fout.write("\t".join(cols) + "\n")

# ================================================================================
# Step 4. Summary
# ================================================================================
log_message(f"\n[SUMMARY]")
log_message(f"  Exact matches  : {matched_count}")
log_message(f"  Fuzzy matches  : {fuzzy_matched_count}")
log_message(f"  Unmatched      : {unmatched_count}")
log_message(f"  Skipped (<50bp): {skipped_small_count}")
log_message(f"  Output VCF     : {vcf_out}")

if len(sys.argv) > 4:
    log_fh.close()
