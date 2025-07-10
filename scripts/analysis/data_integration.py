#!/usr/bin/env python3
"""
Cross-Platform Data Integration for Forensic miRNA Analysis
Merges GPR and CEL data, performs normalization and batch correction
"""

import pandas as pd
import numpy as np
from pathlib import Path
import logging
from sklearn.preprocessing import StandardScaler, quantile_transform
import matplotlib.pyplot as plt
import seaborn as sns

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class ForensicMiRNAIntegrator:
    """Integrate miRNA data from multiple platforms for forensic analysis"""
    
    def __init__(self):
        self.gpr_data = None
        self.cel_data = None
        self.common_mirnas = None
        self.integrated_data = None
        
    def load_gpr_data(self, gpr_path):
        """Load processed GPR data"""
        logger.info("Loading GPR data...")
        
        # Load expression matrix
        expr_matrix = pd.read_csv(gpr_path / 'gpr_expression_matrix.csv', index_col=0)
        
        # Load metadata
        raw_data = pd.read_csv(gpr_path / 'gpr_combined_raw.csv')
        
        # Get sample metadata
        sample_metadata = raw_data[['sample', 'body_fluid']].drop_duplicates()
        
        self.gpr_data = {
            'expression': expr_matrix,
            'metadata': sample_metadata,
            'platform': 'GPR'
        }
        
        logger.info(f"Loaded GPR data: {expr_matrix.shape[0]} samples x {expr_matrix.shape[1]} miRNAs")
        
    def load_cel_data(self, cel_path):
        """Load processed CEL data from R output"""
        logger.info("Loading CEL data...")
        
        # Check if CEL data has been processed
        expr_file = cel_path / 'cel_expression_matrix.csv'
        if not expr_file.exists():
            logger.warning("CEL expression data not found. Please run R processing first.")
            logger.warning("Run: Rscript data/processed/cel/process_cel_files.R")
            return False
        
        # Load expression matrix
        expr_matrix = pd.read_csv(expr_file, index_col=0)
        
        # Load metadata
        metadata = pd.read_csv(cel_path / 'cel_metadata.csv')
        
        # Ensure sample order matches
        metadata = metadata.set_index('sample').loc[expr_matrix.index].reset_index()
        
        self.cel_data = {
            'expression': expr_matrix.T,  # Transpose to match GPR format (samples x miRNAs)
            'metadata': metadata,
            'platform': 'CEL'
        }
        
        logger.info(f"Loaded CEL data: {expr_matrix.shape[0]} samples x {expr_matrix.shape[1]} miRNAs")
        return True
    
    def find_common_mirnas(self):
        """Identify miRNAs present in both platforms"""
        gpr_mirnas = set(self.gpr_data['expression'].columns)
        cel_mirnas = set(self.cel_data['expression'].columns)
        
        self.common_mirnas = list(gpr_mirnas.intersection(cel_mirnas))
        
        logger.info(f"Found {len(self.common_mirnas)} common miRNAs")
        logger.info(f"GPR unique: {len(gpr_mirnas - cel_mirnas)}")
        logger.info(f"CEL unique: {len(cel_mirnas - gpr_mirnas)}")
        
        if len(self.common_mirnas) < 50:
            # If too few common miRNAs, use name matching
            logger.info("Few exact matches, trying fuzzy matching...")
            self.common_mirnas = self._fuzzy_match_mirnas(gpr_mirnas, cel_mirnas)
        
        return self.common_mirnas
    
    def _fuzzy_match_mirnas(self, gpr_mirnas, cel_mirnas):
        """Match miRNAs by core name (e.g., miR-21 matches hsa-miR-21-5p)"""
        matches = []
        
        for gpr_mir in gpr_mirnas:
            # Extract core miRNA number
            if 'miR-' in gpr_mir:
                core = gpr_mir.split('miR-')[1].split('-')[0]
                for cel_mir in cel_mirnas:
                    if f'miR-{core}' in cel_mir:
                        matches.append(gpr_mir)
                        break
        
        logger.info(f"Fuzzy matching found {len(matches)} common miRNAs")
        return matches[:200]  # Use top 200 for analysis
    
    def normalize_within_platform(self):
        """Apply quantile normalization within each platform"""
        logger.info("Applying within-platform normalization...")
        
        # GPR normalization
        gpr_expr = self.gpr_data['expression'][self.common_mirnas]
        gpr_norm = pd.DataFrame(
            quantile_transform(gpr_expr.T, n_quantiles=1000, random_state=42).T,
            index=gpr_expr.index,
            columns=gpr_expr.columns
        )
        
        # CEL normalization
        cel_expr = self.cel_data['expression'][self.common_mirnas]
        cel_norm = pd.DataFrame(
            quantile_transform(cel_expr.T, n_quantiles=1000, random_state=42).T,
            index=cel_expr.index,
            columns=cel_expr.columns
        )
        
        return gpr_norm, cel_norm
    
    def combat_batch_correction(self, gpr_norm, cel_norm):
        """Apply ComBat batch correction (simplified version)"""
        logger.info("Applying batch correction...")
        
        # Combine data
        all_samples = pd.concat([gpr_norm, cel_norm])
        
        # Create batch labels
        batch = ['GPR'] * len(gpr_norm) + ['CEL'] * len(cel_norm)
        
        # Create biological labels
        all_metadata = pd.concat([
            self.gpr_data['metadata'].set_index('sample'),
            self.cel_data['metadata'].set_index('sample')
        ])
        
        bio_labels = all_metadata.loc[all_samples.index, 'body_fluid']
        
        # Simplified batch correction (z-score normalization per batch)
        # In production, use pycombat or similar
        corrected_data = all_samples.copy()
        
        for platform in ['GPR', 'CEL']:
            mask = [b == platform for b in batch]
            platform_data = all_samples[mask]
            
            # Center and scale
            scaler = StandardScaler()
            corrected_data.loc[mask] = scaler.fit_transform(platform_data)
        
        # Preserve biological variation
        for fluid in bio_labels.unique():
            fluid_mask = bio_labels == fluid
            fluid_mean = corrected_data[fluid_mask].mean()
            
            # Add back biological signal
            corrected_data[fluid_mask] += fluid_mean * 0.3
        
        return corrected_data, batch, bio_labels
    
    def integrate_platforms(self):
        """Main integration pipeline"""
        logger.info("Starting platform integration...")
        
        # Find common miRNAs
        self.find_common_mirnas()
        
        if len(self.common_mirnas) == 0:
            logger.error("No common miRNAs found!")
            return None
        
        # Normalize within platforms
        gpr_norm, cel_norm = self.normalize_within_platform()
        
        # Batch correction
        integrated, batch, bio_labels = self.combat_batch_correction(gpr_norm, cel_norm)
        
        # Store results
        self.integrated_data = {
            'expression': integrated,
            'batch': batch,
            'body_fluid': bio_labels,
            'common_mirnas': self.common_mirnas
        }
        
        logger.info(f"Integration complete: {integrated.shape}")
        
        return integrated
    
    def quality_control_plots(self, output_dir):
        """Generate QC plots for integration"""
        output_path = Path(output_dir)
        output_path.mkdir(exist_ok=True)
        
        if self.integrated_data is None:
            logger.error("No integrated data available")
            return
        
        # PCA plot
        from sklearn.decomposition import PCA
        
        expr = self.integrated_data['expression']
        pca = PCA(n_components=2)
        pca_coords = pca.fit_transform(expr)
        
        # Create PCA plot colored by platform
        plt.figure(figsize=(10, 8))
        
        plt.subplot(2, 2, 1)
        for platform in ['GPR', 'CEL']:
            mask = [b == platform for b in self.integrated_data['batch']]
            plt.scatter(pca_coords[mask, 0], pca_coords[mask, 1], 
                       label=platform, alpha=0.7, s=100)
        plt.xlabel(f'PC1 ({pca.explained_variance_ratio_[0]:.1%})')
        plt.ylabel(f'PC2 ({pca.explained_variance_ratio_[1]:.1%})')
        plt.title('PCA by Platform')
        plt.legend()
        
        # PCA colored by body fluid
        plt.subplot(2, 2, 2)
        fluids = self.integrated_data['body_fluid'].unique()
        colors = plt.cm.Set3(np.linspace(0, 1, len(fluids)))
        
        for fluid, color in zip(fluids, colors):
            mask = self.integrated_data['body_fluid'] == fluid
            plt.scatter(pca_coords[mask, 0], pca_coords[mask, 1], 
                       label=fluid, alpha=0.7, s=100, color=color)
        plt.xlabel(f'PC1 ({pca.explained_variance_ratio_[0]:.1%})')
        plt.ylabel(f'PC2 ({pca.explained_variance_ratio_[1]:.1%})')
        plt.title('PCA by Body Fluid')
        plt.legend()
        
        # Expression distribution
        plt.subplot(2, 2, 3)
        for platform in ['GPR', 'CEL']:
            mask = [b == platform for b in self.integrated_data['batch']]
            platform_data = expr[mask].values.flatten()
            plt.hist(platform_data, bins=50, alpha=0.5, label=platform, density=True)
        plt.xlabel('Normalized Expression')
        plt.ylabel('Density')
        plt.title('Expression Distribution by Platform')
        plt.legend()
        
        # Heatmap of top variable miRNAs
        plt.subplot(2, 2, 4)
        var_mirnas = expr.var().nlargest(20).index
        subset_data = expr[var_mirnas].T
        
        # Create row colors for body fluids
        fluid_colors = {'blood': 'red', 'saliva': 'blue', 'semen': 'green',
                       'vaginal_secretion': 'purple', 'menstrual_blood': 'orange',
                       'peripheral_blood': 'darkred'}
        row_colors = [fluid_colors.get(f, 'gray') for f in self.integrated_data['body_fluid']]
        
        im = plt.imshow(subset_data, aspect='auto', cmap='RdBu_r')
        plt.colorbar(im, label='Expression')
        plt.title('Top 20 Variable miRNAs')
        plt.xlabel('Samples')
        plt.ylabel('miRNAs')
        
        plt.tight_layout()
        plt.savefig(output_path / 'integration_qc.png', dpi=300, bbox_inches='tight')
        plt.close()
        
        logger.info(f"QC plots saved to {output_path / 'integration_qc.png'}")
        
        # Save integrated data
        self.integrated_data['expression'].to_csv(output_path / 'integrated_expression.csv')
        
        # Save metadata
        metadata_df = pd.DataFrame({
            'sample': self.integrated_data['expression'].index,
            'platform': self.integrated_data['batch'],
            'body_fluid': self.integrated_data['body_fluid']
        })
        metadata_df.to_csv(output_path / 'integrated_metadata.csv', index=False)
        
        logger.info("Integration data saved")


def main():
    """Run the integration pipeline"""
    integrator = ForensicMiRNAIntegrator()
    
    # Load data
    integrator.load_gpr_data(Path('data/processed/gpr'))
    
    # Try to load CEL data
    cel_loaded = integrator.load_cel_data(Path('data/processed/cel'))
    
    if not cel_loaded:
        logger.error("CEL data not available. Cannot perform cross-platform integration.")
        logger.info("Proceeding with GPR-only analysis...")
        return
    
    # Integrate
    integrated_data = integrator.integrate_platforms()
    
    if integrated_data is not None:
        # Generate QC plots
        integrator.quality_control_plots('data/processed/integrated')
        
        # Summary statistics
        logger.info("\n=== Integration Summary ===")
        logger.info(f"Total samples: {integrated_data.shape[0]}")
        logger.info(f"Total miRNAs: {integrated_data.shape[1]}")
        logger.info(f"Platforms: GPR ({sum(b == 'GPR' for b in integrator.integrated_data['batch'])}), "
                   f"CEL ({sum(b == 'CEL' for b in integrator.integrated_data['batch'])})")
        
        logger.info("\nBody fluid distribution:")
        print(integrator.integrated_data['body_fluid'].value_counts())

if __name__ == "__main__":
    main()