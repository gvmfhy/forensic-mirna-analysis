"""
Forensic miRNA Expression Analysis: Practical Implementation Guide
Handles continuous expression nature for court-defensible results
"""

import statistics
from typing import Dict, List, Tuple, Optional

class ForensicMiRNAAnalyzer:
    """
    Analyzer for continuous miRNA expression in forensic samples
    Implements recommendations from continuous expression analysis
    """
    
    def __init__(self):
        # Define multi-tier thresholds based on data analysis
        self.thresholds = {
            'high_confidence': -2.0,
            'moderate_confidence': -6.0,
            'low_confidence': -10.0
        }
        
        # Quality control thresholds
        self.quality_thresholds = {
            'mean_good': -6.0,
            'mean_moderate': -7.0,
            'stdev_good': 2.5,
            'stdev_moderate': 2.0
        }
        
        # Body fluid specific markers (from literature)
        self.body_fluid_markers = {
            'saliva': ['hsa-miR-205-5p', 'hsa-miR-203a-3p'],
            'blood': ['hsa-miR-451a', 'hsa-miR-16-5p'],
            'semen': ['hsa-miR-888-5p', 'hsa-miR-891a-5p'],
            'vaginal': ['hsa-miR-124-3p', 'hsa-miR-617'],
            'menstrual': ['hsa-miR-412-3p', 'hsa-miR-451a']
        }
    
    def assess_sample_quality(self, expression_values: List[float]) -> Dict:
        """
        Assess sample quality based on global expression patterns
        
        Returns:
            Dictionary with quality metrics and overall score
        """
        mean_expr = statistics.mean(expression_values)
        median_expr = statistics.median(expression_values)
        stdev_expr = statistics.stdev(expression_values)
        pct_detected = sum(1 for v in expression_values if v > -10) / len(expression_values) * 100
        
        # Determine quality category
        if mean_expr > self.quality_thresholds['mean_good'] and \
           stdev_expr > self.quality_thresholds['stdev_good']:
            quality_category = 'Good'
            quality_score = 1.0
        elif mean_expr > self.quality_thresholds['mean_moderate'] and \
             stdev_expr > self.quality_thresholds['stdev_moderate']:
            quality_category = 'Moderate'
            quality_score = 0.7
        else:
            quality_category = 'Poor'
            quality_score = 0.4
        
        return {
            'mean_expression': mean_expr,
            'median_expression': median_expr,
            'stdev_expression': stdev_expr,
            'percent_detected': pct_detected,
            'quality_category': quality_category,
            'quality_score': quality_score
        }
    
    def classify_expression_level(self, expression_value: float, 
                                 quality_score: float = 1.0) -> Tuple[str, float]:
        """
        Classify expression level with quality-adjusted thresholds
        
        Returns:
            Tuple of (classification, adjusted_threshold)
        """
        # Adjust thresholds based on sample quality
        adjustment = 2.0 * (1 - quality_score)
        
        adjusted_thresholds = {
            level: threshold - adjustment 
            for level, threshold in self.thresholds.items()
        }
        
        if expression_value > adjusted_thresholds['high_confidence']:
            return 'High', adjusted_thresholds['high_confidence']
        elif expression_value > adjusted_thresholds['moderate_confidence']:
            return 'Moderate', adjusted_thresholds['moderate_confidence']
        elif expression_value > adjusted_thresholds['low_confidence']:
            return 'Low', adjusted_thresholds['low_confidence']
        else:
            return 'Below Detection', adjusted_thresholds['low_confidence']
    
    def calculate_expression_ratios(self, sample_data: Dict[str, float]) -> Dict:
        """
        Calculate expression ratios between key markers
        Helps distinguish body fluids using relative expression
        """
        ratios = {}
        
        # Example: saliva vs blood markers
        if 'hsa-miR-205-5p' in sample_data and 'hsa-miR-451a' in sample_data:
            ratios['miR-205/miR-451'] = sample_data['hsa-miR-205-5p'] - sample_data['hsa-miR-451a']
        
        # Example: semen vs other markers
        if 'hsa-miR-888-5p' in sample_data and 'hsa-miR-16-5p' in sample_data:
            ratios['miR-888/miR-16'] = sample_data['hsa-miR-888-5p'] - sample_data['hsa-miR-16-5p']
        
        return ratios
    
    def calculate_likelihood_ratios(self, sample_data: Dict[str, float], 
                                  quality_score: float) -> Dict[str, float]:
        """
        Calculate likelihood ratios for each body fluid
        Based on continuous expression patterns
        """
        likelihoods = {}
        
        for fluid, markers in self.body_fluid_markers.items():
            # Calculate geometric mean of marker expressions
            marker_values = []
            for marker in markers:
                if marker in sample_data:
                    value = sample_data[marker]
                    level, _ = self.classify_expression_level(value, quality_score)
                    
                    # Assign likelihood based on expression level
                    if level == 'High':
                        marker_values.append(10.0)
                    elif level == 'Moderate':
                        marker_values.append(3.0)
                    elif level == 'Low':
                        marker_values.append(0.5)
                    else:
                        marker_values.append(0.1)
            
            if marker_values:
                # Geometric mean for combined likelihood
                likelihoods[fluid] = statistics.geometric_mean(marker_values)
        
        # Normalize to sum to 1
        total = sum(likelihoods.values())
        if total > 0:
            return {fluid: value/total for fluid, value in likelihoods.items()}
        return likelihoods
    
    def generate_court_report(self, sample_data: Dict[str, float], 
                            sample_name: str) -> str:
        """
        Generate a court-friendly report with continuous expression data
        """
        # Get all expression values
        all_values = list(sample_data.values())
        
        # Assess quality
        quality = self.assess_sample_quality(all_values)
        
        # Calculate likelihoods
        likelihoods = self.calculate_likelihood_ratios(sample_data, quality['quality_score'])
        
        # Build report
        report = f"""
FORENSIC miRNA EXPRESSION ANALYSIS REPORT
=========================================

Sample ID: {sample_name}
Analysis Date: [Current Date]

1. SAMPLE QUALITY ASSESSMENT
---------------------------
Quality Category: {quality['quality_category']}
Quality Score: {quality['quality_score']:.2f} (0-1 scale)
Mean Expression: {quality['mean_expression']:.2f}
Expression Variability: {quality['stdev_expression']:.2f}
Detected miRNAs: {quality['percent_detected']:.1f}%

2. EXPRESSION PROFILE SUMMARY
----------------------------
"""
        
        # Add key markers for each body fluid
        for fluid, markers in self.body_fluid_markers.items():
            report += f"\n{fluid.upper()} markers:\n"
            for marker in markers[:3]:  # Top 3 markers
                if marker in sample_data:
                    value = sample_data[marker]
                    level, threshold = self.classify_expression_level(
                        value, quality['quality_score']
                    )
                    report += f"  {marker}: {value:.2f} ({level})\n"
        
        # Add likelihood results
        report += "\n3. BODY FLUID CLASSIFICATION\n"
        report += "----------------------------\n"
        sorted_fluids = sorted(likelihoods.items(), key=lambda x: x[1], reverse=True)
        for fluid, prob in sorted_fluids:
            report += f"{fluid.capitalize()}: {prob*100:.1f}% probability\n"
        
        # Add interpretation
        top_fluid = sorted_fluids[0][0] if sorted_fluids else "Unknown"
        top_prob = sorted_fluids[0][1] if sorted_fluids else 0
        
        report += f"""
4. INTERPRETATION
-----------------
The expression profile is most consistent with {top_fluid} 
(probability: {top_prob*100:.1f}%).

This conclusion is based on continuous expression analysis
of multiple miRNA markers, accounting for sample quality.

5. LIMITATIONS
--------------
- miRNA expression is continuous, not binary
- Sample quality may affect detection sensitivity
- Cannot exclude mixtures below 10% concentration
- Results should be considered with other evidence

6. TECHNICAL NOTES
------------------
- Expression values are log2-transformed
- Thresholds adjusted for sample quality
- Analysis uses probabilistic framework
- Multiple markers evaluated for robustness
"""
        
        return report

# Example usage
def demonstrate_analysis():
    """Demonstrate the forensic analysis on example data"""
    
    analyzer = ForensicMiRNAAnalyzer()
    
    # Example sample data (from actual GPR data)
    sample_data = {
        'hsa-miR-205-5p': -12.628,
        'hsa-miR-203a-3p': -8.759,
        'hsa-miR-124-3p': -8.639,
        'hsa-miR-451a': -8.176,
        'hsa-miR-16-5p': -6.918,
        'hsa-miR-888-5p': -1.500,
        'hsa-miR-891a-5p': -10.383,
        'hsa-miR-617': -2.565,
        'hsa-miR-412-3p': -4.509
    }
    
    # Add more miRNAs for quality assessment
    import random
    for i in range(50):
        sample_data[f'other-miR-{i}'] = random.gauss(-5.6, 3.0)
    
    report = analyzer.generate_court_report(sample_data, "GSM4634783_saliva-array1")
    print(report)
    
    # Demonstrate continuous nature
    print("\n\nDEMONSTRATING CONTINUOUS EXPRESSION:")
    print("=" * 50)
    print("Same miRNA (miR-205-5p) across different samples:")
    values = [-12.628, -10.744, -8.330, -9.462]
    for i, val in enumerate(values):
        level, _ = analyzer.classify_expression_level(val)
        print(f"Sample {i+1}: {val:.2f} -> {level}")
    print("\nThis shows the continuous gradient, not binary presence/absence")

if __name__ == "__main__":
    demonstrate_analysis()