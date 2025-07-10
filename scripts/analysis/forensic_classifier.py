#!/usr/bin/env python3
"""
Forensic miRNA Ensemble Classifier
Identifies body fluid-specific miRNA markers using multiple ML approaches
"""

import pandas as pd
import numpy as np
from pathlib import Path
import logging
from sklearn.model_selection import StratifiedKFold, cross_val_score
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import LabelEncoder
from sklearn.metrics import classification_report, confusion_matrix
import xgboost as xgb
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.feature_selection import SelectKBest, f_classif

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class ForensicEnsembleClassifier:
    """Ensemble classifier for forensic body fluid identification"""
    
    def __init__(self):
        self.expression_data = None
        self.metadata = None
        self.selected_features = None
        self.models = {}
        self.feature_importance = None
        
    def load_integrated_data(self, integrated_path):
        """Load integrated expression data and metadata"""
        logger.info("Loading integrated data...")
        
        self.expression_data = pd.read_csv(integrated_path / 'integrated_expression.csv', index_col=0)
        self.metadata = pd.read_csv(integrated_path / 'integrated_metadata.csv')
        
        # Ensure sample order matches
        self.metadata = self.metadata.set_index('sample').loc[self.expression_data.index].reset_index()
        
        logger.info(f"Loaded data: {self.expression_data.shape}")
        
    def feature_selection(self, n_features=50):
        """Select top features using multiple methods"""
        logger.info(f"Selecting top {n_features} features...")
        
        X = self.expression_data.values
        y = self.metadata['body_fluid'].values
        
        # Method 1: ANOVA F-test
        selector = SelectKBest(f_classif, k=n_features)
        selector.fit(X, y)
        anova_features = self.expression_data.columns[selector.get_support()].tolist()
        
        # Method 2: Random Forest importance
        rf = RandomForestClassifier(n_estimators=100, random_state=42)
        rf.fit(X, y)
        rf_importance = pd.Series(rf.feature_importances_, index=self.expression_data.columns)
        rf_features = rf_importance.nlargest(n_features).index.tolist()
        
        # Method 3: XGBoost importance
        le = LabelEncoder()
        y_encoded = le.fit_transform(y)
        xgb_model = xgb.XGBClassifier(n_estimators=100, random_state=42, use_label_encoder=False)
        xgb_model.fit(X, y_encoded)
        xgb_importance = pd.Series(xgb_model.feature_importances_, index=self.expression_data.columns)
        xgb_features = xgb_importance.nlargest(n_features).index.tolist()
        
        # Consensus features (appear in at least 2 methods)
        all_features = anova_features + rf_features + xgb_features
        feature_counts = pd.Series(all_features).value_counts()
        consensus_features = feature_counts[feature_counts >= 2].index.tolist()
        
        # If too few consensus features, add top-ranked ones
        if len(consensus_features) < 30:
            additional_features = rf_importance.index[:50].tolist()
            consensus_features = list(set(consensus_features + additional_features))[:n_features]
        
        self.selected_features = consensus_features[:n_features]
        
        logger.info(f"Selected {len(self.selected_features)} consensus features")
        
        # Store feature importance
        self.feature_importance = pd.DataFrame({
            'miRNA': self.selected_features,
            'rf_importance': [rf_importance.get(f, 0) for f in self.selected_features],
            'xgb_importance': [xgb_importance.get(f, 0) for f in self.selected_features]
        }).sort_values('rf_importance', ascending=False)
        
        return self.selected_features
    
    def train_ensemble(self):
        """Train ensemble of classifiers"""
        logger.info("Training ensemble classifiers...")
        
        X = self.expression_data[self.selected_features].values
        y = self.metadata['body_fluid'].values
        
        # Initialize models
        self.models = {
            'random_forest': RandomForestClassifier(
                n_estimators=200, 
                max_depth=10,
                min_samples_split=5,
                random_state=42
            ),
            'xgboost': xgb.XGBClassifier(
                n_estimators=200,
                max_depth=6,
                learning_rate=0.1,
                random_state=42,
                use_label_encoder=False
            ),
            'logistic': LogisticRegression(
                multi_class='ovr',
                max_iter=1000,
                random_state=42,
                C=0.1
            )
        }
        
        # Cross-validation
        cv = StratifiedKFold(n_splits=3, shuffle=True, random_state=42)  # Reduced splits due to small sample size
        
        cv_scores = {}
        le = LabelEncoder()
        y_encoded = le.fit_transform(y)
        
        for name, model in self.models.items():
            if name == 'xgboost':
                scores = cross_val_score(model, X, y_encoded, cv=cv, scoring='accuracy')
            else:
                scores = cross_val_score(model, X, y, cv=cv, scoring='accuracy')
            cv_scores[name] = scores
            logger.info(f"{name}: {scores.mean():.3f} (+/- {scores.std():.3f})")
        
        # Train final models on full data
        for name, model in self.models.items():
            if name == 'xgboost':
                le = LabelEncoder()
                y_encoded = le.fit_transform(y)
                model.fit(X, y_encoded)
            else:
                model.fit(X, y)
        
        return cv_scores
    
    def predict_ensemble(self, X):
        """Make predictions using ensemble voting"""
        predictions = []
        
        for name, model in self.models.items():
            if name == 'xgboost':
                # XGBoost returns numeric labels
                pred = model.predict(X)
                # Convert back to string labels
                le = LabelEncoder()
                le.fit(self.metadata['body_fluid'])
                pred = le.inverse_transform(pred)
            else:
                pred = model.predict(X)
            predictions.append(pred)
        
        # Majority voting
        predictions = np.array(predictions)
        ensemble_pred = []
        
        for i in range(predictions.shape[1]):
            votes = predictions[:, i]
            unique, counts = np.unique(votes, return_counts=True)
            ensemble_pred.append(unique[np.argmax(counts)])
        
        return np.array(ensemble_pred)
    
    def evaluate_forensic_performance(self, output_dir):
        """Evaluate classifier performance for forensic applications"""
        output_path = Path(output_dir)
        output_path.mkdir(exist_ok=True)
        
        X = self.expression_data[self.selected_features].values
        y_true = self.metadata['body_fluid'].values
        
        # Get ensemble predictions
        y_pred = self.predict_ensemble(X)
        
        # Classification report
        report = classification_report(y_true, y_pred, output_dict=True)
        report_df = pd.DataFrame(report).transpose()
        report_df.to_csv(output_path / 'classification_report.csv')
        
        logger.info("\nClassification Report:")
        print(classification_report(y_true, y_pred))
        
        # Confusion matrix
        cm = confusion_matrix(y_true, y_pred)
        
        plt.figure(figsize=(10, 8))
        sns.heatmap(cm, annot=True, fmt='d', cmap='Blues')
        plt.title('Confusion Matrix - Ensemble Classifier')
        plt.ylabel('True Label')
        plt.xlabel('Predicted Label')
        plt.tight_layout()
        plt.savefig(output_path / 'confusion_matrix.png', dpi=300)
        plt.close()
        
        # Feature importance plot
        plt.figure(figsize=(10, 12))
        top_features = self.feature_importance.head(20)
        
        y_pos = np.arange(len(top_features))
        plt.barh(y_pos, top_features['rf_importance'])
        plt.yticks(y_pos, top_features['miRNA'])
        plt.xlabel('Importance Score')
        plt.title('Top 20 miRNA Markers for Body Fluid Identification')
        plt.tight_layout()
        plt.savefig(output_path / 'feature_importance.png', dpi=300, bbox_inches='tight')
        plt.close()
        
        # Calculate forensic metrics
        self._calculate_forensic_metrics(y_true, y_pred, output_path)
        
        return report_df
    
    def _calculate_forensic_metrics(self, y_true, y_pred, output_path):
        """Calculate metrics specific to forensic applications"""
        logger.info("\nForensic Performance Metrics:")
        
        classes = np.unique(y_true)
        forensic_metrics = []
        
        for fluid in classes:
            # Binary classification for each fluid
            binary_true = (y_true == fluid).astype(int)
            binary_pred = (y_pred == fluid).astype(int)
            
            tp = np.sum((binary_true == 1) & (binary_pred == 1))
            tn = np.sum((binary_true == 0) & (binary_pred == 0))
            fp = np.sum((binary_true == 0) & (binary_pred == 1))
            fn = np.sum((binary_true == 1) & (binary_pred == 0))
            
            specificity = tn / (tn + fp) if (tn + fp) > 0 else 0
            sensitivity = tp / (tp + fn) if (tp + fn) > 0 else 0
            ppv = tp / (tp + fp) if (tp + fp) > 0 else 0
            npv = tn / (tn + fn) if (tn + fn) > 0 else 0
            
            forensic_metrics.append({
                'body_fluid': fluid,
                'specificity': specificity,
                'sensitivity': sensitivity,
                'ppv': ppv,
                'npv': npv,
                'samples': np.sum(y_true == fluid)
            })
            
            logger.info(f"{fluid}: Specificity={specificity:.3f}, Sensitivity={sensitivity:.3f}")
        
        forensic_df = pd.DataFrame(forensic_metrics)
        forensic_df.to_csv(output_path / 'forensic_metrics.csv', index=False)
        
        # Check if any fluid achieves >99% specificity (forensic standard)
        high_spec = forensic_df[forensic_df['specificity'] >= 0.99]
        if len(high_spec) > 0:
            logger.info(f"\nFluids meeting forensic standard (>99% specificity): {high_spec['body_fluid'].tolist()}")
        else:
            logger.info("\nNo fluids currently meet >99% specificity standard")
    
    def generate_marker_panel(self, output_dir):
        """Generate final marker panel for forensic implementation"""
        output_path = Path(output_dir)
        
        # Get top markers
        top_markers = self.feature_importance.head(15)
        
        # Calculate expression profiles
        marker_profiles = []
        for fluid in self.metadata['body_fluid'].unique():
            fluid_mask = self.metadata['body_fluid'] == fluid
            fluid_samples = self.metadata[fluid_mask]['sample'].values
            fluid_expr = self.expression_data.loc[fluid_samples, top_markers['miRNA']]
            
            profile = {
                'body_fluid': fluid,
                'n_samples': len(fluid_samples)
            }
            
            for mirna in top_markers['miRNA']:
                profile[f'{mirna}_mean'] = fluid_expr[mirna].mean()
                profile[f'{mirna}_std'] = fluid_expr[mirna].std()
            
            marker_profiles.append(profile)
        
        marker_df = pd.DataFrame(marker_profiles)
        marker_df.to_csv(output_path / 'forensic_marker_panel.csv', index=False)
        
        # Create marker heatmap
        plt.figure(figsize=(12, 8))
        
        # Extract mean expression values
        expr_matrix = []
        for _, row in marker_df.iterrows():
            expr_values = [row[f'{mirna}_mean'] for mirna in top_markers['miRNA']]
            expr_matrix.append(expr_values)
        
        expr_matrix = np.array(expr_matrix)
        
        sns.heatmap(expr_matrix, 
                   xticklabels=top_markers['miRNA'],
                   yticklabels=marker_df['body_fluid'],
                   cmap='RdBu_r',
                   center=0,
                   annot=True,
                   fmt='.2f')
        
        plt.title('Forensic miRNA Marker Panel - Mean Expression by Body Fluid')
        plt.xlabel('miRNA Markers')
        plt.ylabel('Body Fluid')
        plt.tight_layout()
        plt.savefig(output_path / 'marker_panel_heatmap.png', dpi=300, bbox_inches='tight')
        plt.close()
        
        logger.info(f"\nGenerated forensic marker panel with {len(top_markers)} miRNAs")
        logger.info(f"Top 5 markers: {top_markers['miRNA'].head().tolist()}")


def main():
    """Run the forensic classification pipeline"""
    classifier = ForensicEnsembleClassifier()
    
    # Load data
    classifier.load_integrated_data(Path('data/processed/integrated'))
    
    # Feature selection
    classifier.feature_selection(n_features=50)
    
    # Train ensemble
    cv_scores = classifier.train_ensemble()
    
    # Evaluate performance
    output_dir = Path('results/forensic_classification')
    report = classifier.evaluate_forensic_performance(output_dir)
    
    # Generate marker panel
    classifier.generate_marker_panel(output_dir)
    
    logger.info("\n=== Forensic Classification Complete ===")
    logger.info(f"Results saved to {output_dir}")
    logger.info(f"Selected features saved for qPCR primer design")

if __name__ == "__main__":
    main()