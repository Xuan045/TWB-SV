#!/usr/bin/env python3
"""
Filter cross-cohort related samples from test set.

Usage:
    python 01_2_final_remove_list.py \
        --pairs twb2_rel_test_removed.2degree.txt \
        --test twb2_id_test_20.txt \
        --train twb2_id_train_80.txt \
        --output twb2_final

Logic:
    - Read related pairs (IID1, IID2)
    - For each pair, check if one ID is in test and the other is in train
    - If cross-cohort: remove the test-side ID from the test set
    - Same-cohort pairs (test-test or train-train) are ignored

Outputs:
    <o>.removed_test.txt   - IDs to remove from test set (1 field)
    <o>.removed_test.indv  - IDs to remove from test set (2 field, tab-delimited)
    <o>.removed_train.txt  - Related train-side IDs (1 field)
    <o>.removed_train.indv - Related train-side IDs (2 field, tab-delimited)
"""

import argparse
from pathlib import Path


def load_id_list(filepath: str) -> set:
    """Load a list of IDs from a file (one ID per line)."""
    ids = set()
    with open(filepath, "r") as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith("#"):
                ids.add(line)
    return ids


def load_pairs(filepath: str) -> list:
    """Load related pairs from a tab-separated file with header IID1 IID2."""
    pairs = []
    with open(filepath, "r") as f:
        for i, line in enumerate(f):
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                print(f"  [WARNING] Line {i+1} skipped (not enough columns): {line}")
                continue
            id1, id2 = parts[0], parts[1]
            # Skip header row
            if id1 == "IID1" and id2 == "IID2":
                continue
            pairs.append((id1, id2))
    return pairs


def classify_pair(id1, id2, test_set, train_set):
    """
    Classify a pair and return the test-side ID to remove (if cross-cohort).

    Returns:
        (label1, label2, remove_id or None)
    """
    in_test_1  = id1 in test_set
    in_train_1 = id1 in train_set
    in_test_2  = id2 in test_set
    in_train_2 = id2 in train_set

    label1 = "TEST"  if in_test_1  else ("TRAIN" if in_train_1 else "UNKNOWN")
    label2 = "TEST"  if in_test_2  else ("TRAIN" if in_train_2 else "UNKNOWN")

    # Cross-cohort: one in test, the other in train
    if in_test_1 and in_train_2:
        return label1, label2, id1
    if in_train_1 and in_test_2:
        return label1, label2, id2

    return label1, label2, None


def main():
    parser = argparse.ArgumentParser(
        description="Remove cross-cohort related samples from test set."
    )
    parser.add_argument("--pairs",  required=True, help="Related pairs file (TSV, header: IID1 IID2)")
    parser.add_argument("--test",   required=True, help="Test set ID list (one ID per line)")
    parser.add_argument("--train",  required=True, help="Train set ID list (one ID per line)")
    parser.add_argument("--output", required=True, help="Output prefix for removed ID lists")
    args = parser.parse_args()

    # ── Load inputs ──────────────────────────────────────────────────────────
    print("=" * 60)
    print("Loading files...")
    test_set  = load_id_list(args.test)
    train_set = load_id_list(args.train)
    pairs     = load_pairs(args.pairs)

    print(f"  Test  set : {len(test_set):>6} samples  ({args.test})")
    print(f"  Train set : {len(train_set):>6} samples  ({args.train})")
    print(f"  Pairs     : {len(pairs):>6} related pairs  ({args.pairs})")

    # ── Classify pairs ───────────────────────────────────────────────────────
    print("\n" + "=" * 60)
    print("Classifying pairs...")

    cross_cohort = []   # (id1, id2, label1, label2, remove_id)
    same_cohort  = []
    unknown      = []

    for id1, id2 in pairs:
        label1, label2, remove_id = classify_pair(id1, id2, test_set, train_set)

        if remove_id is not None:
            cross_cohort.append((id1, id2, label1, label2, remove_id))
        elif "UNKNOWN" in (label1, label2):
            unknown.append((id1, id2, label1, label2))
        else:
            same_cohort.append((id1, id2, label1, label2))

    # ── Report ───────────────────────────────────────────────────────────────
    print(f"\n  Cross-cohort pairs (test <-> train) : {len(cross_cohort)}")
    print(f"  Same-cohort pairs                   : {len(same_cohort)}")
    print(f"  Pairs with UNKNOWN IDs              : {len(unknown)}")

    if cross_cohort:
        print("\n" + "=" * 60)
        print("Cross-cohort pairs → removing from TEST:")
        for id1, id2, l1, l2, rm in cross_cohort:
            print(f"  {id1} [{l1}]  <->  {id2} [{l2}]  =>  REMOVE: {rm}")

    if same_cohort:
        print("\n" + "=" * 60)
        print("Same-cohort pairs (no action):")
        for id1, id2, l1, l2 in same_cohort:
            print(f"  {id1} [{l1}]  <->  {id2} [{l2}]")

    if unknown:
        print("\n" + "=" * 60)
        print("Pairs with IDs not found in test or train (no action):")
        for id1, id2, l1, l2 in unknown:
            print(f"  {id1} [{l1}]  <->  {id2} [{l2}]")

    # ── Collect removed IDs from each side ───────────────────────────────────
    test_to_remove  = set()
    train_to_remove = set()

    for id1, id2, l1, l2, rm in cross_cohort:
        test_to_remove.add(rm)
        train_side = id2 if rm == id1 else id1
        train_to_remove.add(train_side)

    # ── Summary ──────────────────────────────────────────────────────────────
    print("\n" + "=" * 60)
    print("Result summary:")
    print(f"  Original test  size : {len(test_set)}")
    print(f"  Test  removed       : {len(test_to_remove)}")
    print()
    print(f"  Original train size : {len(train_set)}")
    print(f"  Train removed       : {len(train_to_remove)}")

    if test_to_remove:
        print("\n  [TEST] Removed IDs:")
        for rid in sorted(test_to_remove):
            print(f"    - {rid}")

    if train_to_remove:
        print("\n  [TRAIN] Their related train-side IDs:")
        for rid in sorted(train_to_remove):
            print(f"    - {rid}")

    # ── Write output ─────────────────────────────────────────────────────────
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Removed list: test side
    with open(f"{output_path}.removed_test.txt", "w") as f:
        for sample_id in sorted(test_to_remove):
            f.write(sample_id + "\n")
    with open(f"{output_path}.removed_test.indv", "w") as f:
        for sample_id in sorted(test_to_remove):
            f.write(sample_id + "\t" + sample_id + "\n")

    # Removed list: train side
    with open(f"{output_path}.removed_train.txt", "w") as f:
        for sample_id in sorted(train_to_remove):
            f.write(sample_id + "\n")
    with open(f"{output_path}.removed_train.indv", "w") as f:
        for sample_id in sorted(train_to_remove):
            f.write(sample_id + "\t" + sample_id + "\n")

    print(f"\nOutputs:")
    print(f"  Removed from test  : {output_path}.removed_test.txt / .removed_test.indv")
    print(f"  Removed from train : {output_path}.removed_train.txt / .removed_train.indv")
    print("=" * 60)


if __name__ == "__main__":
    main()
