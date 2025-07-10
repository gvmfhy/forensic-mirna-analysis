#!/usr/bin/env python3
"""
GPR Parser for Forensic miRNA Analysis
Extracts and normalizes data from GenePix Result (GPR) files
"""

import pandas as pd
import numpy as np
import gzip
import os
from pathlib import Path
import logging

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class GPRParser:
    """Parse and process GenePix Result files for miRNA analysis"""
    
    def __init__(self, gpr_path):
        self.gpr_path = gpr_path
        self.metadata = {}
        self.data = None
        self.mirna_data = None
        
    def parse_file(self):
        """Parse GPR file and extract metadata and intensity data"""
        logger.info(f"Parsing {os.path.basename(self.gpr_path)}")
        
        with gzip.open(self.gpr_path, 'rt') as f:
            lines = f.readlines()
        
        # Parse metadata section
        data_start_idx = self._parse_metadata(lines)
        
        # Parse data section
        self._parse_data(lines, data_start_idx)
        
        # Filter for miRNA spots
        self._filter_mirna_spots()
        
        return self
    
    def _parse_metadata(self, lines):
        """Extract metadata from GPR header"""
        data_start_idx = 0
        
        for i, line in enumerate(lines):
            if line.startswith('"') and '=' in line:
                # Parse metadata lines
                key_value = line.strip().strip('"')
                if '=' in key_value:
                    key, value = key_value.split('=', 1)
                    self.metadata[key] = value.strip('"')
            elif line.startswith('"Block"'):
                # Found data section
                data_start_idx = i
                break
        
        logger.info(f"Platform: {self.metadata.get('Product', 'Unknown')}")
        logger.info(f"Scanner: {self.metadata.get('Scanner', 'Unknown')}")
        
        return data_start_idx
    
    def _parse_data(self, lines, start_idx):
        """Parse the data section of GPR file"""
        # Read column headers
        headers = lines[start_idx].strip().split('\t')
        headers = [h.strip('"') for h in headers]
        
        # Parse data rows
        data_rows = []
        for line in lines[start_idx + 1:]:
            if line.strip():
                fields = line.strip().split('\t')
                fields = [f.strip('"') for f in fields]
                if len(fields) == len(headers):
                    data_rows.append(fields)
        
        # Create DataFrame
        self.data = pd.DataFrame(data_rows, columns=headers)
        
        # Convert numeric columns
        numeric_cols = ['F635 Median', 'F635 Mean', 'B635', 'B635 Median',
                       'F532 Median', 'F532 Mean', 'B532', 'B532 Median',
                       'F635 SD', 'F532 SD', 'B635 SD', 'B532 SD']
        
        for col in numeric_cols:
            if col in self.data.columns:
                self.data[col] = pd.to_numeric(self.data[col], errors='coerce')
        
        logger.info(f"Loaded {len(self.data)} spots")
    
    def _filter_mirna_spots(self):
        """Filter data to keep only miRNA spots"""
        # Keep rows where Name contains 'miR'
        mirna_mask = self.data['Name'].str.contains('miR', na=False)
        self.mirna_data = self.data[mirna_mask].copy()
        
        # Remove empty miRNA names
        self.mirna_data = self.mirna_data[self.mirna_data['Name'] != '']
        
        logger.info(f"Found {len(self.mirna_data)} miRNA spots")
    
    def calculate_intensities(self, method='ratio'):
        """Calculate normalized intensities from raw data
        
        Args:
            method: 'ratio' for Cy5/Cy3, 'cy5' for Cy5 only, 'cy3' for Cy3 only
        """
        if self.mirna_data is None:
            raise ValueError("No miRNA data available. Run parse_file() first.")
        
        # Background correction
        cy5_corrected = self.mirna_data['F635 Median'] - self.mirna_data['B635 Median']
        cy3_corrected = self.mirna_data['F532 Median'] - self.mirna_data['B532 Median']
        
        # Handle negative values after background correction
        cy5_corrected = np.maximum(cy5_corrected, 1)
        cy3_corrected = np.maximum(cy3_corrected, 1)
        
        if method == 'ratio':
            # Calculate log2 ratio
            intensities = np.log2(cy5_corrected / cy3_corrected)
        elif method == 'cy5':
            # Use Cy5 channel only
            intensities = np.log2(cy5_corrected)
        elif method == 'cy3':
            # Use Cy3 channel only
            intensities = np.log2(cy3_corrected)
        else:
            raise ValueError(f"Unknown method: {method}")
        
        # Create result DataFrame
        result = pd.DataFrame({
            'miRNA': self.mirna_data['Name'].values,
            'intensity': intensities,
            'cy5_raw': self.mirna_data['F635 Median'].values,
            'cy3_raw': self.mirna_data['F532 Median'].values,
            'cy5_bg': self.mirna_data['B635 Median'].values,
            'cy3_bg': self.mirna_data['B532 Median'].values
        })
        
        # Remove infinite values
        result = result[~result['intensity'].isin([np.inf, -np.inf])]
        
        return result
    
    def quality_metrics(self):
        """Calculate quality metrics for the array"""
        if self.mirna_data is None:
            return None
        
        metrics = {
            'total_spots': len(self.data),
            'mirna_spots': len(self.mirna_data),
            'mirna_percentage': len(self.mirna_data) / len(self.data) * 100,
            'cy5_median': self.mirna_data['F635 Median'].median(),
            'cy3_median': self.mirna_data['F532 Median'].median(),
            'cy5_bg_median': self.mirna_data['B635 Median'].median(),
            'cy3_bg_median': self.mirna_data['B532 Median'].median()
        }
        
        return metrics


def process_all_gpr_files(input_dir, output_dir, method='ratio'):
    """Process all GPR files in a directory"""
    input_path = Path(input_dir)
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    gpr_files = list(input_path.glob("*.gpr.gz"))
    logger.info(f"Found {len(gpr_files)} GPR files")
    
    all_data = []
    quality_reports = []
    
    for gpr_file in gpr_files:
        try:
            # Parse file
            parser = GPRParser(gpr_file)
            parser.parse_file()
            
            # Calculate intensities
            intensities = parser.calculate_intensities(method=method)
            
            # Add sample name
            sample_name = os.path.basename(gpr_file).replace('.gpr.gz', '')
            intensities['sample'] = sample_name
            
            # Extract body fluid type from filename
            if 'saliva' in sample_name.lower():
                body_fluid = 'saliva'
            elif 'peripheral_blood' in sample_name.lower():
                body_fluid = 'peripheral_blood'
            elif 'menstrual_blood' in sample_name.lower():
                body_fluid = 'menstrual_blood'
            elif 'semen' in sample_name.lower():
                body_fluid = 'semen'
            elif 'vaginal' in sample_name.lower():
                body_fluid = 'vaginal_secretion'
            else:
                body_fluid = 'unknown'
            
            intensities['body_fluid'] = body_fluid
            
            all_data.append(intensities)
            
            # Quality metrics
            metrics = parser.quality_metrics()
            metrics['sample'] = sample_name
            metrics['body_fluid'] = body_fluid
            quality_reports.append(metrics)
            
        except Exception as e:
            logger.error(f"Error processing {gpr_file}: {e}")
    
    # Combine all data
    if all_data:
        combined_data = pd.concat(all_data, ignore_index=True)
        
        # Save raw combined data
        combined_data.to_csv(output_path / 'gpr_combined_raw.csv', index=False)
        logger.info(f"Saved combined data to {output_path / 'gpr_combined_raw.csv'}")
        
        # Create wide format (samples x miRNAs)
        # Handle duplicate miRNAs by taking the mean
        wide_data = combined_data.groupby(['sample', 'miRNA'])['intensity'].mean().unstack()
        wide_data.to_csv(output_path / 'gpr_expression_matrix.csv')
        logger.info(f"Saved expression matrix to {output_path / 'gpr_expression_matrix.csv'}")
        
        # Save quality report
        quality_df = pd.DataFrame(quality_reports)
        quality_df.to_csv(output_path / 'gpr_quality_report.csv', index=False)
        logger.info(f"Saved quality report to {output_path / 'gpr_quality_report.csv'}")
        
        return combined_data, quality_df
    
    return None, None


if __name__ == "__main__":
    # Process GPR files
    input_dir = "data/raw/GSE153135_GPR"
    output_dir = "data/processed/gpr"
    
    logger.info("Starting GPR processing pipeline")
    combined_data, quality_report = process_all_gpr_files(input_dir, output_dir, method='ratio')
    
    if combined_data is not None:
        # Summary statistics
        logger.info("\n=== Processing Summary ===")
        logger.info(f"Total samples: {combined_data['sample'].nunique()}")
        logger.info(f"Total unique miRNAs: {combined_data['miRNA'].nunique()}")
        logger.info(f"Body fluid distribution:")
        print(combined_data.groupby('body_fluid')['sample'].nunique())
        
        logger.info("\n=== Quality Summary ===")
        print(quality_report[['sample', 'mirna_spots', 'cy5_median', 'cy3_median']])