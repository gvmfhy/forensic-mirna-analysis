#!/usr/bin/env python3
"""
Check data formats and compatibility between CEL and GPR files
"""

import os
import gzip
import re
from pathlib import Path

def check_gpr_format(gpr_file):
    """Parse GPR file and extract key information"""
    print(f"\n=== Analyzing GPR file: {os.path.basename(gpr_file)} ===")
    
    with gzip.open(gpr_file, 'rt') as f:
        lines = f.readlines()
    
    # Extract metadata
    metadata = {}
    data_start_idx = 0
    
    for i, line in enumerate(lines):
        if line.startswith('"') and '=' in line:
            key_value = line.strip().strip('"')
            if '=' in key_value:
                key, value = key_value.split('=', 1)
                metadata[key] = value.strip('"')
        elif line.startswith('"Block"'):
            data_start_idx = i
            break
    
    # Parse column headers
    headers = lines[data_start_idx].strip().split('\t')
    headers = [h.strip('"') for h in headers]
    
    # Count miRNA entries
    mirna_count = 0
    mirna_examples = []
    
    for line in lines[data_start_idx + 1:]:
        fields = line.strip().split('\t')
        if len(fields) > 3:
            name = fields[3].strip('"')
            if 'miR' in name and name:
                mirna_count += 1
                if len(mirna_examples) < 5:
                    mirna_examples.append(name)
    
    print(f"Platform: {metadata.get('Product', 'Unknown')}")
    print(f"Scanner: {metadata.get('Scanner', 'Unknown')}")
    print(f"Array contents: {metadata.get('ArrayContents_DesignedFor', 'Unknown')}")
    print(f"Total spots: {len(lines) - data_start_idx - 1}")
    print(f"miRNA spots: {mirna_count}")
    print(f"Example miRNAs: {', '.join(mirna_examples)}")
    
    return {
        'platform': metadata.get('Product', 'Unknown'),
        'total_spots': len(lines) - data_start_idx - 1,
        'mirna_count': mirna_count,
        'channels': ['F635', 'F532'],  # Red and Green channels
        'file_type': 'GPR'
    }

def check_cel_annotation(annotation_file):
    """Parse CEL annotation file"""
    print("\n=== Analyzing CEL annotation file ===")
    
    mirna_probes = []
    total_probes = 0
    
    with gzip.open(annotation_file, 'rt', encoding='latin-1') as f:
        # Skip header lines
        for line in f:
            if not line.startswith('#'):
                break
        
        # Parse data
        for line in f:
            total_probes += 1
            fields = line.strip().split('","')
            if len(fields) > 1:
                probe_name = fields[1].strip('"')
                if 'miR' in probe_name:
                    mirna_probes.append(probe_name)
    
    print(f"Total probes: {total_probes}")
    print(f"miRNA probes: {len(mirna_probes)}")
    print(f"Example miRNAs: {', '.join(mirna_probes[:5])}")
    
    return {
        'platform': 'Affymetrix miRNA-3_0',
        'total_probes': total_probes,
        'mirna_count': len(mirna_probes),
        'channels': ['Single channel'],
        'file_type': 'CEL'
    }

def main():
    """Main analysis function"""
    print("=== Forensic miRNA Data Format Analysis ===")
    
    # Check GPR files
    gpr_dir = Path("data/raw/GSE153135_GPR")
    gpr_files = list(gpr_dir.glob("*.gpr.gz"))
    
    if gpr_files:
        gpr_info = check_gpr_format(gpr_files[0])
    
    # Check CEL annotation
    cel_annotation = Path("data/raw/GSE49630_CEL/GPL16384_miRNA-3_1-st-v1.annotations.20140513.csv.gz")
    if cel_annotation.exists():
        cel_info = check_cel_annotation(cel_annotation)
    
    # Summary
    print("\n=== Integration Challenges Summary ===")
    print("1. Channel differences: CEL (single) vs GPR (dual)")
    print("2. Different array designs with partial overlap")
    print("3. Need to map probes to common miRNA identifiers")
    print("4. Batch effects from different platforms")
    print("\nRecommended approach:")
    print("- Extract single intensity measure from each platform")
    print("- Map to miRBase IDs for common namespace")
    print("- Apply cross-platform normalization")

if __name__ == "__main__":
    main()