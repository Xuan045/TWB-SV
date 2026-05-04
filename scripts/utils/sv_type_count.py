#!/usr/bin/env python3
"""
Variant Analysis Script with Rare Variant Stratification
Analyzes genetic variants by chromosome and type, with granular MAF bins for rare variants.

MAF BINS:
- Ultra-rare: (0.0,0.001], (0.001,0.005]
- Rare: (0.005,0.01], (0.01,0.05]
- Low frequency: (0.05,0.1]
- Common: (0.1,0.2], (0.2,0.3], (0.3,0.4], (0.4,0.5]
"""

import pandas as pd
import numpy as np
import os
from glob import glob
from typing import List, Dict, Tuple


class VariantAnalyzer:
    """Class to handle variant analysis operations."""
    
    def __init__(self, data_path: str):
        self.data_path = data_path
        self.file_list = glob(os.path.join(data_path, "chr*.imp_results.txt"))
        self.quality_col = 'R2'
        self.maf_col = 'MAF'
        
        # Define ordering for categorical variables with granular rare variant bins
        self.maf_order = [
            '(0.0,0.001]',    # Ultra-rare
            '(0.001,0.005]',  # Ultra-rare
            '(0.005,0.01]',   # Rare
            '(0.01,0.05]',    # Rare
            '(0.05,0.1]',     # Low frequency
            '(0.1,0.2]',      # Common
            '(0.2,0.3]',      # Common
            '(0.3,0.4]',      # Common
            '(0.4,0.5]'       # Common
        ]
        self.var_type_order = ['DEL', 'INDEL', 'INS', 'SNV', 'Other']
        self.r2_order = ['(0.0,0.1]', '(0.1,0.2]', '(0.2,0.3]', '(0.3,0.4]', '(0.4,0.5]', 
                         '(0.5,0.6]', '(0.6,0.7]', '(0.7,0.8]', '(0.8,0.9]', '(0.9,1.0]']
        self.imp_status_order = ['Imputed', 'Genotyped', 'Both', 'Neither']
    
    @staticmethod
    def assign_r2_bin(r2: float) -> str:
        """
        Assign R2 to predefined quality bins (0.1 intervals).
        
        Args:
            r2: R2 quality score value
            
        Returns:
            str: R2 bin label
        """
        if pd.isna(r2):
            return None
        elif r2 <= 0.1:
            return '(0.0,0.1]'
        elif r2 <= 0.2:
            return '(0.1,0.2]'
        elif r2 <= 0.3:
            return '(0.2,0.3]'
        elif r2 <= 0.4:
            return '(0.3,0.4]'
        elif r2 <= 0.5:
            return '(0.4,0.5]'
        elif r2 <= 0.6:
            return '(0.5,0.6]'
        elif r2 <= 0.7:
            return '(0.6,0.7]'
        elif r2 <= 0.8:
            return '(0.7,0.8]'
        elif r2 <= 0.9:
            return '(0.8,0.9]'
        elif r2 <= 1.0:
            return '(0.9,1.0]'
        else:
            return 'Other'
    
    @staticmethod
    def classify_variant(row: pd.Series) -> str:
        """
        Classify variant type based on ID and REF/ALT lengths.
        
        Args:
            row: DataFrame row containing ID, REF, and ALT columns
            
        Returns:
            str: Variant type classification
        """
        variant_id = str(row['ID'])
        ref_len = len(row['REF'])
        alt_len = len(row['ALT'])
        
        if 'DEL' in variant_id or ref_len >= 50:
            return 'DEL'
        elif 'INS' in variant_id or alt_len >= 50:
            return 'INS'
        elif ref_len == 1 and alt_len == 1:
            return 'SNV'
        elif (1 < ref_len < 50) or (1 < alt_len < 50):
            return 'INDEL'
        else:
            return 'Other'
    
    @staticmethod
    def assign_maf_bin(maf: float) -> str:
        """
        Assign MAF to predefined bins with granular stratification for rare variants.
        
        Args:
            maf: Minor allele frequency value
            
        Returns:
            str: MAF bin label
        """
        if pd.isna(maf):
            return 'Other'
        elif maf <= 0.001:
            return '(0.0,0.001]'
        elif maf <= 0.005:
            return '(0.001,0.005]'
        elif maf <= 0.01:
            return '(0.005,0.01]'
        elif maf <= 0.05:
            return '(0.01,0.05]'
        elif maf <= 0.1:
            return '(0.05,0.1]'
        elif maf <= 0.2:
            return '(0.1,0.2]'
        elif maf <= 0.3:
            return '(0.2,0.3]'
        elif maf <= 0.4:
            return '(0.3,0.4]'
        elif maf <= 0.5:
            return '(0.4,0.5]'
        else:
            return 'Other'
    
    @staticmethod
    def classify_imputation_status(row: pd.Series) -> str:
        """
        Classify imputation status based on IMPUTED and TYPED columns.
        
        Args:
            row: DataFrame row containing IMPUTED and TYPED columns
            
        Returns:
            str: Imputation status classification
        """
        imputed = row['IMPUTED']
        typed = row['TYPED']
        
        imputed_flag = str(imputed) == '1' or imputed == 1
        typed_flag = str(typed) == '1' or typed == 1
        
        if imputed_flag and typed_flag:
            return 'Both'
        elif imputed_flag:
            return 'Imputed'
        elif typed_flag:
            return 'Genotyped'
        else:
            return 'Neither'
    
    def load_and_process_data(self) -> pd.DataFrame:
        """
        Load all files and process variant classifications.
        
        Returns:
            pd.DataFrame: Combined processed data
        """
        all_data = []
        
        for file_path in self.file_list:
            chr_id = os.path.basename(file_path).split('.')[0]
            
            try:
                df = pd.read_csv(file_path, sep='\t')
                
                required_cols = ['IMPUTED', 'TYPED', 'ID', 'REF', 'ALT', self.maf_col, self.quality_col]
                missing_cols = [col for col in required_cols if col not in df.columns]
                if missing_cols:
                    print(f"Warning: Missing columns in {file_path}: {missing_cols}")
                    continue
                
                df['var_type'] = df.apply(self.classify_variant, axis=1)
                df['MAF_bin'] = df[self.maf_col].apply(self.assign_maf_bin)
                df['R2_bin'] = df[self.quality_col].apply(self.assign_r2_bin)
                df['imp_status'] = df.apply(self.classify_imputation_status, axis=1)
                df['chr'] = chr_id
                
                processed_df = df[['var_type', 'MAF_bin', 'R2_bin', 'imp_status', 
                                 self.maf_col, 'chr', self.quality_col]].copy()
                all_data.append(processed_df)
                
                print(f"Processed {chr_id}: {len(df)} variants")
                
            except Exception as e:
                print(f"Error processing {file_path}: {e}")
                continue
        
        if not all_data:
            raise ValueError("No data files were successfully processed")
        
        return pd.concat(all_data, ignore_index=True)
    
    def create_chromosome_summary(self, combined_df: pd.DataFrame) -> pd.DataFrame:
        """Create summary table by chromosome and variant type."""
        chr_summary = (combined_df.groupby(['chr', 'var_type'])
                      .size()
                      .reset_index(name='count'))
        
        pivot_df = chr_summary.pivot_table(
            index='chr', 
            columns='var_type', 
            values='count', 
            fill_value=0
        ).reset_index()
        
        return pivot_df
    
    def create_maf_summary(self, combined_df: pd.DataFrame) -> pd.DataFrame:
        """Create summary statistics by MAF bin and variant type."""
        summary_stats = combined_df.groupby(['MAF_bin', 'var_type']).agg({
            self.maf_col: 'count',
            self.quality_col: ['mean', 'median']
        }).round(6)
        
        summary_stats.columns = ['count', 'mean_r2', 'median_r2']
        summary_stats = summary_stats.reset_index()
        
        summary_stats['MAF_bin'] = pd.Categorical(
            summary_stats['MAF_bin'], 
            categories=self.maf_order, 
            ordered=True
        )
        summary_stats['var_type'] = pd.Categorical(
            summary_stats['var_type'], 
            categories=self.var_type_order, 
            ordered=True
        )
        
        return summary_stats.sort_values(['MAF_bin', 'var_type']).reset_index(drop=True)
    
    def create_r2_summary(self, combined_df: pd.DataFrame) -> pd.DataFrame:
        """Create summary table by R2 bin and variant type."""
        filtered_df = combined_df.dropna(subset=['R2_bin'])
        
        if len(filtered_df) == 0:
            print("Warning: No valid R2 values found after filtering NA values")
            return pd.DataFrame()
        
        print(f"Filtered out {len(combined_df) - len(filtered_df):,} rows with NA R2 values")
        
        r2_summary = (filtered_df.groupby(['R2_bin', 'var_type'])
                     .size()
                     .reset_index(name='count'))
        
        r2_summary['R2_bin'] = pd.Categorical(
            r2_summary['R2_bin'], 
            categories=self.r2_order, 
            ordered=True
        )
        r2_summary['var_type'] = pd.Categorical(
            r2_summary['var_type'], 
            categories=self.var_type_order, 
            ordered=True
        )
        
        return r2_summary.sort_values(['R2_bin', 'var_type']).reset_index(drop=True)
    
    def create_imputation_summary(self, combined_df: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """Create summary tables for imputation status analysis."""
        imp_by_type = (combined_df.groupby(['var_type', 'imp_status'])
                      .size()
                      .reset_index(name='count'))
        
        imp_by_type['var_type'] = pd.Categorical(
            imp_by_type['var_type'], 
            categories=self.var_type_order, 
            ordered=True
        )
        imp_by_type['imp_status'] = pd.Categorical(
            imp_by_type['imp_status'], 
            categories=self.imp_status_order, 
            ordered=True
        )
        
        imp_by_type = imp_by_type.sort_values(['var_type', 'imp_status']).reset_index(drop=True)
        
        imp_by_maf = (combined_df.groupby(['MAF_bin', 'imp_status'])
                     .size()
                     .reset_index(name='count'))
        
        imp_by_maf['MAF_bin'] = pd.Categorical(
            imp_by_maf['MAF_bin'], 
            categories=self.maf_order, 
            ordered=True
        )
        imp_by_maf['imp_status'] = pd.Categorical(
            imp_by_maf['imp_status'], 
            categories=self.imp_status_order, 
            ordered=True
        )
        
        imp_by_maf = imp_by_maf.sort_values(['MAF_bin', 'imp_status']).reset_index(drop=True)
        
        return imp_by_type, imp_by_maf
    
    def create_comprehensive_summary(self, combined_df: pd.DataFrame) -> pd.DataFrame:
        """Create comprehensive summary table with MAF bins, variant types, and imputation status."""
        comprehensive = combined_df.groupby(['MAF_bin', 'var_type', 'imp_status']).agg({
            self.maf_col: 'count',
            self.quality_col: ['mean', 'median']
        }).round(6)
        
        comprehensive.columns = ['count', 'mean_r2', 'median_r2']
        comprehensive = comprehensive.reset_index()
        
        comprehensive['MAF_bin'] = pd.Categorical(
            comprehensive['MAF_bin'], 
            categories=self.maf_order, 
            ordered=True
        )
        comprehensive['var_type'] = pd.Categorical(
            comprehensive['var_type'], 
            categories=self.var_type_order, 
            ordered=True
        )
        comprehensive['imp_status'] = pd.Categorical(
            comprehensive['imp_status'], 
            categories=self.imp_status_order, 
            ordered=True
        )
        
        return comprehensive.sort_values(['MAF_bin', 'var_type', 'imp_status']).reset_index(drop=True)
    
    def run_analysis(self) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
        """Run complete analysis pipeline."""
        print("Loading and processing data...")
        combined_df = self.load_and_process_data()
        print(f"Total variants processed: {len(combined_df):,}")
        
        imp_overview = combined_df['imp_status'].value_counts()
        print(f"\nImputation Status Overview:")
        for status, count in imp_overview.items():
            print(f"  {status}: {count:,} ({count/len(combined_df)*100:.1f}%)")
        
        # Print MAF bin distribution
        maf_overview = combined_df['MAF_bin'].value_counts()
        print(f"\nMAF Bin Distribution:")
        for bin_label in self.maf_order:
            if bin_label in maf_overview.index:
                count = maf_overview[bin_label]
                print(f"  {bin_label}: {count:,} ({count/len(combined_df)*100:.1f}%)")
        
        print("\nCreating chromosome summary...")
        chr_summary = self.create_chromosome_summary(combined_df)
        
        print("Creating MAF summary...")
        maf_summary = self.create_maf_summary(combined_df)
        
        print("Creating R2 summary...")
        r2_summary = self.create_r2_summary(combined_df)
        
        print("Creating imputation summaries...")
        imp_by_type, imp_by_maf = self.create_imputation_summary(combined_df)
        
        print("Creating comprehensive summary...")
        comprehensive_summary = self.create_comprehensive_summary(combined_df)
        
        return chr_summary, maf_summary, r2_summary, imp_by_type, imp_by_maf, comprehensive_summary
    
    def save_results(self, chr_summary: pd.DataFrame, maf_summary: pd.DataFrame, 
                    r2_summary: pd.DataFrame, imp_by_type: pd.DataFrame, 
                    imp_by_maf: pd.DataFrame, comprehensive_summary: pd.DataFrame):
        """Save analysis results to files."""
        chr_file = "var_type_counts_by_chr.tsv"
        chr_summary.to_csv(chr_file, sep="\t", index=False)
        print(f"Chromosome summary saved to {chr_file}")
        
        maf_file = "var_maf_summary.tsv"
        maf_summary.to_csv(maf_file, sep="\t", index=False)
        print(f"MAF summary saved to {maf_file}")
        
        if not r2_summary.empty:
            r2_file = "var_r2_summary.tsv"
            r2_summary.to_csv(r2_file, sep="\t", index=False)
            print(f"R2 summary saved to {r2_file}")
        else:
            print("R2 summary not saved (no valid data)")
        
        imp_type_file = "var_imputation_by_type.tsv"
        imp_by_type.to_csv(imp_type_file, sep="\t", index=False)
        print(f"Imputation by type summary saved to {imp_type_file}")
        
        imp_maf_file = "var_imputation_by_maf.tsv"
        imp_by_maf.to_csv(imp_maf_file, sep="\t", index=False)
        print(f"Imputation by MAF summary saved to {imp_maf_file}")
        
        comp_file = "var_comprehensive_summary.tsv"
        comprehensive_summary.to_csv(comp_file, sep="\t", index=False)
        print(f"Comprehensive summary saved to {comp_file}")


def main():
    """Main execution function."""
    data_path = "/staging/biology/u4432941/sv/imputation/analysis/sv_snv_imp_info"
    
    try:
        analyzer = VariantAnalyzer(data_path)
        results = analyzer.run_analysis()
        chr_summary, maf_summary, r2_summary, imp_by_type, imp_by_maf, comprehensive_summary = results
        
        print("\n" + "="*50)
        print("CHROMOSOME SUMMARY")
        print("="*50)
        print(chr_summary)
        
        print("\n" + "="*50)
        print("MAF SUMMARY (First 30 rows)")
        print("="*50)
        print(maf_summary.head(30))
        if len(maf_summary) > 30:
            print(f"... and {len(maf_summary) - 30} more rows")
        
        print("\n" + "="*50)
        print("R2 SUMMARY")
        print("="*50)
        if not r2_summary.empty:
            print(r2_summary)
        else:
            print("No valid R2 data available")
        
        print("\n" + "="*50)
        print("IMPUTATION BY VARIANT TYPE")
        print("="*50)
        print(imp_by_type)
        
        print("\n" + "="*50)
        print("IMPUTATION BY MAF BIN")
        print("="*50)
        print(imp_by_maf)
        
        print("\n" + "="*50)
        print("COMPREHENSIVE SUMMARY (First 40 rows)")
        print("="*50)
        print(comprehensive_summary.head(40))
        if len(comprehensive_summary) > 40:
            print(f"... and {len(comprehensive_summary) - 40} more rows")
        
        print("\n" + "="*50)
        print("SAVING RESULTS")
        print("="*50)
        analyzer.save_results(chr_summary, maf_summary, r2_summary, 
                            imp_by_type, imp_by_maf, comprehensive_summary)
        
    except Exception as e:
        print(f"Analysis failed: {e}")
        raise


if __name__ == "__main__":
    main()
    