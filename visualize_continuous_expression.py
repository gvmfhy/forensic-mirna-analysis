import csv
import statistics
from collections import defaultdict

def create_simple_heatmap():
    """Create a text-based heatmap showing continuous expression"""
    
    # Read specific miRNAs
    target_mirnas = ['hsa-miR-205-5p', 'hsa-miR-451a', 'hsa-miR-891a-5p', 
                     'hsa-miR-888-5p', 'hsa-miR-124-3p', 'hsa-miR-16-5p']
    
    with open('data/processed/gpr/gpr_expression_matrix.csv', 'r') as f:
        reader = csv.reader(f)
        header = next(reader)
        
        # Find indices
        mirna_indices = {}
        for mirna in target_mirnas:
            if mirna in header:
                mirna_indices[mirna] = header.index(mirna)
        
        # Collect data
        heatmap_data = []
        for row in reader:
            sample = row[0]
            # Simplify sample name
            if 'saliva' in sample:
                simple_name = 'Saliva_' + sample.split('-')[-1]
            elif 'blood' in sample and 'menstrual' not in sample:
                simple_name = 'Blood_' + sample.split('-')[-1]
            elif 'semen' in sample:
                simple_name = 'Semen_' + sample.split('-')[-1]
            elif 'vaginal' in sample:
                simple_name = 'Vaginal_' + sample.split('-')[-1]
            elif 'menstrual' in sample:
                simple_name = 'Menstrual_' + sample.split('-')[-1]
            else:
                simple_name = sample
                
            row_data = {'sample': simple_name}
            for mirna, idx in mirna_indices.items():
                row_data[mirna] = float(row[idx])
            heatmap_data.append(row_data)
    
    # Create visual heatmap
    print("\nCONTINUOUS EXPRESSION HEATMAP")
    print("=" * 80)
    print("\nExpression Scale:")
    print("█ = Very High (>-2)")
    print("▓ = High (-4 to -2)")  
    print("▒ = Medium (-6 to -4)")
    print("░ = Low (-8 to -6)")
    print("· = Very Low (-10 to -8)")
    print("  = Below Detection (<-10)")
    print()
    
    # Print header
    print(f"{'Sample':<20}", end='')
    for mirna in target_mirnas:
        short_name = mirna.replace('hsa-miR-', '')[:8]
        print(f"{short_name:<10}", end='')
    print("\n" + "-" * 80)
    
    # Print data with visual representation
    for data in heatmap_data:
        print(f"{data['sample']:<20}", end='')
        for mirna in target_mirnas:
            value = data.get(mirna, -99)
            # Convert to visual
            if value > -2:
                symbol = '█'
            elif value > -4:
                symbol = '▓'
            elif value > -6:
                symbol = '▒'
            elif value > -8:
                symbol = '░'
            elif value > -10:
                symbol = '·'
            else:
                symbol = ' '
            
            # Add exact value in parentheses
            print(f"{symbol} ({value:5.1f})  ", end='')
        print()
    
    print("\n" + "=" * 80)
    print("\nKEY OBSERVATIONS:")
    print("1. Expression exists on a continuous gradient, not binary presence/absence")
    print("2. Even 'specific' markers show variable expression across samples")
    print("3. No clear cutoff separates 'present' from 'absent'")
    print("4. Sample quality affects overall expression levels")

def create_expression_distribution():
    """Show the continuous distribution of a specific miRNA"""
    
    mirna = 'hsa-miR-205-5p'  # Known saliva marker
    
    with open('data/processed/gpr/gpr_expression_matrix.csv', 'r') as f:
        reader = csv.reader(f)
        header = next(reader)
        
        if mirna not in header:
            return
            
        idx = header.index(mirna)
        
        # Collect by body fluid
        fluid_values = defaultdict(list)
        for row in reader:
            sample = row[0]
            value = float(row[idx])
            
            if 'saliva' in sample:
                fluid_values['saliva'].append(value)
            elif 'blood' in sample and 'menstrual' not in sample:
                fluid_values['blood'].append(value)
            elif 'semen' in sample:
                fluid_values['semen'].append(value)
            elif 'vaginal' in sample:
                fluid_values['vaginal'].append(value)
            elif 'menstrual' in sample:
                fluid_values['menstrual'].append(value)
    
    print(f"\n\nDISTRIBUTION OF {mirna} ACROSS BODY FLUIDS")
    print("=" * 60)
    print("\nThis 'saliva-specific' marker shows continuous expression:")
    print()
    
    for fluid, values in fluid_values.items():
        if values:
            print(f"{fluid.upper()}:")
            print(f"  Range: {min(values):.2f} to {max(values):.2f}")
            print(f"  Mean: {statistics.mean(values):.2f}")
            print(f"  Visual: ", end='')
            
            # Create mini histogram
            for v in values:
                if v > -9:
                    print("█", end='')
                elif v > -10:
                    print("▓", end='')
                elif v > -11:
                    print("▒", end='')
                else:
                    print("░", end='')
            print(f" (n={len(values)})")
            print()
    
    print("\nIMPLICATIONS:")
    print("- Even in saliva, expression ranges from -10.7 to -12.6")
    print("- Overlaps with expression in other body fluids")
    print("- Continuous gradient, not binary presence/absence")

def suggest_quality_metrics():
    """Calculate and suggest quality control metrics"""
    
    print("\n\nQUALITY CONTROL METRICS FOR DEGRADED SAMPLES")
    print("=" * 60)
    
    with open('data/processed/gpr/gpr_expression_matrix.csv', 'r') as f:
        reader = csv.reader(f)
        header = next(reader)
        
        sample_metrics = []
        for row in reader:
            sample = row[0]
            values = [float(v) for v in row[1:]]
            
            # Calculate metrics
            mean_expr = statistics.mean(values)
            median_expr = statistics.median(values)
            stdev_expr = statistics.stdev(values)
            pct_detected = sum(1 for v in values if v > -10) / len(values) * 100
            
            sample_metrics.append({
                'sample': sample,
                'mean': mean_expr,
                'median': median_expr,
                'stdev': stdev_expr,
                'pct_detected': pct_detected
            })
    
    print("\nSample Quality Scores:")
    print(f"{'Sample':<30} {'Mean':<8} {'StDev':<8} {'%>-10':<8} {'Quality'}")
    print("-" * 70)
    
    for metrics in sample_metrics:
        # Calculate quality score
        quality = 'Good'
        if metrics['mean'] < -6.0 or metrics['stdev'] < 2.5:
            quality = 'Moderate'
        if metrics['mean'] < -7.0 or metrics['stdev'] < 2.0:
            quality = 'Poor'
            
        print(f"{metrics['sample']:<30} {metrics['mean']:>6.2f} {metrics['stdev']:>7.2f} "
              f"{metrics['pct_detected']:>6.1f}% {quality:>8}")
    
    print("\n\nRECOMMENDED QUALITY THRESHOLDS:")
    print("- Good quality: Mean > -6.0, StDev > 2.5")
    print("- Moderate quality: Mean > -7.0, StDev > 2.0")
    print("- Poor quality: Mean < -7.0 or StDev < 2.0")
    print("\nAdjust expression thresholds based on quality score")

if __name__ == "__main__":
    create_simple_heatmap()
    create_expression_distribution()
    suggest_quality_metrics()