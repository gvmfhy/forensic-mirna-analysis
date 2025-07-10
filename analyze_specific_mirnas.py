import csv
import statistics

# Known body fluid specific miRNAs from literature
fluid_specific_mirnas = {
    'saliva': ['hsa-miR-205-5p', 'hsa-miR-203a-3p', 'hsa-miR-124-3p'],
    'blood': ['hsa-miR-451a', 'hsa-miR-16-5p', 'hsa-miR-486-5p'],
    'semen': ['hsa-miR-891a-5p', 'hsa-miR-888-5p', 'hsa-miR-10b-5p'],
    'vaginal': ['hsa-miR-124-3p', 'hsa-miR-372-3p', 'hsa-miR-617'],
    'menstrual': ['hsa-miR-412-3p', 'hsa-miR-451a', 'hsa-miR-205-5p']
}

# Read data
with open('data/processed/gpr/gpr_expression_matrix.csv', 'r') as f:
    reader = csv.reader(f)
    header = next(reader)
    
    # Find miRNA indices
    mirna_indices = {}
    for fluid, mirnas in fluid_specific_mirnas.items():
        mirna_indices[fluid] = {}
        for mirna in mirnas:
            if mirna in header:
                mirna_indices[fluid][mirna] = header.index(mirna)
    
    # Collect expression values
    expression_data = {}
    for row in reader:
        sample = row[0]
        expression_data[sample] = {}
        for fluid, mirna_dict in mirna_indices.items():
            for mirna, idx in mirna_dict.items():
                expression_data[sample][mirna] = float(row[idx])

# Analyze expression patterns
print("Expression of body fluid-specific miRNAs across samples:\n")

for sample, data in expression_data.items():
    # Determine body fluid type from sample name
    sample_type = None
    for fluid in ['saliva', 'blood', 'semen', 'vaginal', 'menstrual']:
        if fluid in sample.lower():
            sample_type = fluid
            break
    if 'peripheral' in sample.lower():
        sample_type = 'blood'
    
    print(f"\nSample: {sample} (Type: {sample_type})")
    
    # Check expression of each fluid's markers
    for fluid, mirnas in fluid_specific_mirnas.items():
        print(f"\n  {fluid.capitalize()} markers:")
        for mirna in mirnas:
            if mirna in data:
                value = data[mirna]
                status = "HIGH" if value > -2 else "MEDIUM" if value > -6 else "LOW"
                match = "✓" if fluid == sample_type and status != "LOW" else ""
                print(f"    {mirna}: {value:.3f} ({status}) {match}")

# Calculate detection thresholds
print("\n\nDETECTION THRESHOLD ANALYSIS:")
print("=" * 50)

# Collect values for each miRNA in its expected body fluid
expected_high = []
expected_low = []

for sample, data in expression_data.items():
    sample_type = None
    for fluid in ['saliva', 'blood', 'semen', 'vaginal', 'menstrual']:
        if fluid in sample.lower():
            sample_type = fluid
            break
    if 'peripheral' in sample.lower():
        sample_type = 'blood'
    
    if sample_type:
        # Expected high expression
        for mirna in fluid_specific_mirnas.get(sample_type, []):
            if mirna in data:
                expected_high.append(data[mirna])
        
        # Expected low expression (other fluids' markers)
        for other_fluid, mirnas in fluid_specific_mirnas.items():
            if other_fluid != sample_type:
                for mirna in mirnas:
                    if mirna in data:
                        expected_low.append(data[mirna])

if expected_high and expected_low:
    print(f"\nExpected HIGH expression (specific markers in correct fluid):")
    print(f"  Mean: {statistics.mean(expected_high):.3f}")
    print(f"  Median: {statistics.median(expected_high):.3f}")
    print(f"  Min: {min(expected_high):.3f}")
    print(f"  Max: {max(expected_high):.3f}")
    
    print(f"\nExpected LOW expression (specific markers in wrong fluid):")
    print(f"  Mean: {statistics.mean(expected_low):.3f}")
    print(f"  Median: {statistics.median(expected_low):.3f}")
    print(f"  Min: {min(expected_low):.3f}")
    print(f"  Max: {max(expected_low):.3f}")
    
    # Find optimal threshold
    high_sorted = sorted(expected_high)
    low_sorted = sorted(expected_low, reverse=True)
    
    print(f"\nLowest 10% of HIGH expression: {high_sorted[len(high_sorted)//10]:.3f}")
    print(f"Highest 10% of LOW expression: {low_sorted[len(low_sorted)//10]:.3f}")