#!/usr/bin/env python3
"""
GPR-Only Forensic miRNA Classifier
Analyzes only the real GPR data without CEL integration
"""

import pandas as pd
import numpy as np
from pathlib import Path
import logging
from sklearn.model_selection import StratifiedKFold, cross_val_score, train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import StandardScaler
from sklearn.feature_selection import SelectKBest, f_classif
import matplotlib.pyplot as plt
import seaborn as sns

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class GPRForensicClassifier:
    """Forensic classifier for GPR-only data"""
    
    def __init__(self):
        self.expression_data = None
        self.metadata = None
        self.selected_features = None
        self.model = None
        
    def load_gpr_data(self, gpr_path):
        """Load GPR expression data"""
        logger.info("Loading GPR data...")
        
        # Load expression matrix
        self.expression_data = pd.read_csv(gpr_path / 'gpr_expression_matrix.csv', index_col=0)
        
        # Load metadata
        raw_data = pd.read_csv(gpr_path / 'gpr_combined_raw.csv')
        self.metadata = raw_data[['sample', 'body_fluid']].drop_duplicates().set_index('sample')
        
        # Ensure sample order matches
        self.metadata = self.metadata.loc[self.expression_data.index]
        
        # Remove any NaN columns
        self.expression_data = self.expression_data.dropna(axis=1)
        
        logger.info(f"Loaded data: {self.expression_data.shape}")
        logger.info(f"Body fluid distribution:\n{self.metadata['body_fluid'].value_counts()}")
        
    def normalize_data(self):
        """Normalize expression data"""
        logger.info("Normalizing data...")
        
        # Z-score normalization
        scaler = StandardScaler()
        normalized_data = scaler.fit_transform(self.expression_data)
        
        self.expression_data = pd.DataFrame(
            normalized_data,
            index=self.expression_data.index,
            columns=self.expression_data.columns
        )
        
    def select_features(self, n_features=30):
        """Select top discriminative features"""
        logger.info(f"Selecting top {n_features} features...")
        
        X = self.expression_data.values
        y = self.metadata['body_fluid'].values
        
        # With only 2 samples per class, we'll use all data
        # This is a limitation of the small dataset
        valid_mask = [True] * len(y)
        
        X_valid = X[valid_mask]
        y_valid = np.array(y)[valid_mask]
        
        # ANOVA F-test
        selector = SelectKBest(f_classif, k=min(n_features, X_valid.shape[1]))
        selector.fit(X_valid, y_valid)
        
        # Get selected features
        feature_mask = selector.get_support()
        self.selected_features = self.expression_data.columns[feature_mask].tolist()
        
        logger.info(f"Selected {len(self.selected_features)} features")
        
        return self.selected_features
    
    def train_classifier(self):
        """Train Random Forest classifier"""
        logger.info("Training classifier...")
        
        # Use all data (small dataset limitation)
        y = self.metadata['body_fluid'].values
        valid_mask = [True] * len(y)
        
        X = self.expression_data[self.selected_features].values[valid_mask]
        y = y[valid_mask]
        
        # Train model
        self.model = RandomForestClassifier(
            n_estimators=100,
            max_depth=5,
            min_samples_split=2,
            random_state=42
        )
        
        # With only 2 samples per class, use leave-one-out cross-validation
        from sklearn.model_selection import LeaveOneOut
        cv = LeaveOneOut()
        scores = cross_val_score(self.model, X, y, cv=cv, scoring='accuracy')
        
        logger.info(f"Cross-validation accuracy: {scores.mean():.3f} (+/- {scores.std():.3f})")
        
        # Train on full data
        self.model.fit(X, y)
        
        # Feature importance
        importance_df = pd.DataFrame({
            'miRNA': self.selected_features,
            'importance': self.model.feature_importances_
        }).sort_values('importance', ascending=False)
        
        return importance_df
    
    def visualize_results(self, output_dir):
        """Create visualizations"""
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)
        
        # PCA visualization
        from sklearn.decomposition import PCA
        
        # Use all data
        y = self.metadata['body_fluid'].values
        valid_mask = [True] * len(y)
        
        X = self.expression_data[self.selected_features].values[valid_mask]
        y_valid = y[valid_mask]
        
        pca = PCA(n_components=2)
        pca_coords = pca.fit_transform(X)
        
        plt.figure(figsize=(10, 8))
        
        # PCA plot
        plt.subplot(2, 2, 1)
        for fluid in np.unique(y_valid):
            mask = y_valid == fluid
            plt.scatter(pca_coords[mask, 0], pca_coords[mask, 1], 
                       label=fluid, alpha=0.7, s=100)
        plt.xlabel(f'PC1 ({pca.explained_variance_ratio_[0]:.1%})')
        plt.ylabel(f'PC2 ({pca.explained_variance_ratio_[1]:.1%})')
        plt.title('PCA of GPR Data')
        plt.legend()
        
        # Feature importance
        plt.subplot(2, 2, 2)
        importance_df = pd.DataFrame({
            'miRNA': self.selected_features,
            'importance': self.model.feature_importances_
        }).sort_values('importance', ascending=False).head(10)
        
        plt.barh(range(len(importance_df)), importance_df['importance'])
        plt.yticks(range(len(importance_df)), importance_df['miRNA'])
        plt.xlabel('Importance')
        plt.title('Top 10 miRNA Markers')
        
        # Expression heatmap
        plt.subplot(2, 1, 2)
        
        # Get top features
        top_features = importance_df['miRNA'].head(15).tolist()
        
        # Create mean expression by body fluid
        mean_expr = []
        fluids = []
        for fluid in np.unique(y_valid):
            fluid_samples = self.expression_data.index[valid_mask][y_valid == fluid]
            mean_expr.append(self.expression_data.loc[fluid_samples, top_features].mean())
            fluids.append(fluid)
        
        mean_expr_df = pd.DataFrame(mean_expr, index=fluids)
        
        sns.heatmap(mean_expr_df.T, cmap='RdBu_r', center=0, 
                   xticklabels=True, yticklabels=True,
                   cbar_kws={'label': 'Normalized Expression'})
        plt.title('Mean Expression of Top 15 miRNAs by Body Fluid')
        plt.xlabel('Body Fluid')
        
        plt.tight_layout()
        plt.savefig(output_path / 'gpr_analysis_results.png', dpi=300, bbox_inches='tight')
        plt.close()
        
        # Save results
        importance_df.to_csv(output_path / 'gpr_feature_importance.csv', index=False)
        
        # Save predictions
        predictions = self.model.predict(X)
        pred_df = pd.DataFrame({
            'sample': self.expression_data.index[valid_mask],
            'true_label': y_valid,
            'predicted_label': predictions
        })
        pred_df.to_csv(output_path / 'gpr_predictions.csv', index=False)
        
        # Calculate accuracy by fluid
        from sklearn.metrics import classification_report
        report = classification_report(y_valid, predictions, output_dict=True)
        report_df = pd.DataFrame(report).transpose()
        report_df.to_csv(output_path / 'gpr_classification_report.csv')
        
        logger.info(f"Results saved to {output_path}")
        
        return importance_df

def main():
    """Run GPR-only analysis"""
    classifier = GPRForensicClassifier()
    
    # Load data
    classifier.load_gpr_data(Path('data/processed/gpr'))
    
    # Normalize
    classifier.normalize_data()
    
    # Feature selection
    classifier.select_features(n_features=30)
    
    # Train classifier
    importance_df = classifier.train_classifier()
    
    # Visualize results
    output_dir = Path('results/gpr_only_analysis')
    classifier.visualize_results(output_dir)
    
    logger.info("\n=== GPR-Only Analysis Complete ===")
    logger.info(f"Top 5 miRNA markers:")
    print(importance_df.head())
    
    logger.info("\nThis analysis uses only the REAL GPR data.")
    logger.info("No simulated data was used in this analysis.")

if __name__ == "__main__":
    main()