import csv
import statistics
from collections import Counter

# Read expression values
values = []
sample_info = {}
with open('data/processed/gpr/gpr_expression_matrix.csv', 'r') as f:
    reader = csv.reader(f)
    header = next(reader)
    mirnas = header[1:]
    for row in reader:
        sample = row[0]
        sample_info[sample] = []
        for val in row[1:]:  # Skip sample name
            fval = float(val)
            values.append(fval)
            sample_info[sample].append(fval)

# Calculate statistics
print(f"Total values: {len(values)}")
print(f"Number of miRNAs: {len(mirnas)}")
print(f"Number of samples: {len(sample_info)}")
print(f"Min: {min(values):.3f}")
print(f"Max: {max(values):.3f}")
print(f"Mean: {statistics.mean(values):.3f}")
print(f"Median: {statistics.median(values):.3f}")
print(f"Stdev: {statistics.stdev(values):.3f}")

# Calculate percentiles
sorted_values = sorted(values)
n = len(sorted_values)
percentiles = [0, 1, 5, 10, 25, 50, 75, 90, 95, 99, 100]
print("\nPercentiles:")
for p in percentiles:
    idx = int(n * p / 100)
    if idx >= n:
        idx = n - 1
    print(f"{p}%: {sorted_values[idx]:.3f}")

# Count values in ranges
ranges = [
    (-20, -10),
    (-10, -8),
    (-8, -6),
    (-6, -4),
    (-4, -2),
    (-2, 0),
    (0, 1)
]
print("\nValue ranges:")
for low, high in ranges:
    count = sum(1 for v in values if low <= v < high)
    pct = 100 * count / len(values)
    print(f"[{low}, {high}): {count} ({pct:.1f}%)")

# Analyze by body fluid type
print("\nBody fluid types:")
fluid_types = {
    'saliva': [],
    'blood': [],
    'semen': [],
    'vaginal': [],
    'menstrual': []
}

for sample, vals in sample_info.items():
    for fluid in fluid_types:
        if fluid in sample.lower():
            fluid_types[fluid].extend(vals)
            break

print("\nStatistics by body fluid:")
for fluid, vals in fluid_types.items():
    if vals:
        print(f"\n{fluid.capitalize()} (n={len(vals)} values):")
        print(f"  Mean: {statistics.mean(vals):.3f}")
        print(f"  Median: {statistics.median(vals):.3f}")
        print(f"  Min: {min(vals):.3f}")
        print(f"  Max: {max(vals):.3f}")
        print(f"  Stdev: {statistics.stdev(vals):.3f}")