#!/usr/bin/env python3
"""
GenePix Result (GPR) File Parser for Forensic miRNA Analysis

This module processes two-color microarray data from Agilent GenePix scanners,
specifically designed for GSE153135 dataset containing pooled body fluid samples.
The parser handles the unique challenges of forensic miRNA analysis including:

- Two-color array background correction (Cy5/Cy3 channels)
- Pooled sample processing (26 donors across 10 arrays)
- miRNA probe identification and filtering
- Intensity normalization for downstream analysis

Key Technical Details:
    - File Format: GenePix Result (.gpr) format from Agilent scanners
    - Array Type: Two-color miRNA arrays with Cy5 (sample) and Cy3 (reference)
    - Background Correction: Local background subtraction per spot
    - Normalization: Log2 ratio transformation for normal distribution
    - Quality Control: Spot flagging and intensity validation

Historical Context:
    This parser was developed after discovering that GSE153135 contains pooled 
    samples (not individual donors as initially assumed). The pooling strategy
    used 2-3 donors per array to enable robust statistical comparisons while
    managing sample availability constraints in forensic research.

WHY this approach:
    - Two-color arrays require specialized background correction
    - Log2 ratios normalize for probe-specific effects
    - Pooling increases biological validity but reduces statistical power
    - GPR format preserves all quality metrics needed for forensic validation

Author: Forensic miRNA Analysis Pipeline
Date: 2024
"""

import pandas as pd
import numpy as np
import gzip
import os
from pathlib import Path
import logging
from typing import Dict, List, Tuple, Optional, Union

# Configure logging with detailed format for debugging
logging.basicConfig(
    level=logging.INFO, 
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Analysis Constants
MIN_INTENSITY = 1.0  # Minimum intensity to prevent log(0) errors
DEFAULT_METHOD = 'ratio'  # Default calculation method for two-color arrays
MIRNA_PATTERN = 'miR'  # Pattern to identify miRNA probes
EXPECTED_CHANNELS = ['F635', 'F532', 'B635', 'B532']  # Required GPR columns

# Expected GPR file structure (for validation)
REQUIRED_COLUMNS = [
    'Block', 'Column', 'Row', 'Name', 'ID',  # Spot identification
    'F635 Median', 'F635 Mean', 'F635 SD',   # Cy5 foreground statistics
    'B635 Median', 'B635 Mean', 'B635 SD',   # Cy5 background statistics
    'F532 Median', 'F532 Mean', 'F532 SD',   # Cy3 foreground statistics
    'B532 Median', 'B532 Mean', 'B532 SD',   # Cy3 background statistics
]

# Body fluid identification patterns for GSE153135
FLUID_PATTERNS = {
    'saliva': ['saliva', 'oral'],
    'peripheral_blood': ['peripheral_blood', 'whole_blood'],
    'menstrual_blood': ['menstrual_blood', 'menstrual'],
    'semen': ['semen', 'seminal'],
    'vaginal': ['vaginal', 'vagina']
}


class GPRParser:
    """
    Parser for GenePix Result files from two-color miRNA microarray experiments.
    
    This class handles the complete processing pipeline for GPR files generated
    by Agilent GenePix scanners. It extracts both metadata and intensity data,
    performs background correction, and prepares normalized intensities for
    statistical analysis.
    
    The parser was specifically designed for GSE153135, which contains pooled
    body fluid samples analyzed on two-color miRNA arrays. Each array contains
    approximately 40,000 spots representing human miRNAs and controls.
    
    Key Features:
        - Handles gzip-compressed GPR files automatically
        - Extracts experimental metadata from file headers
        - Performs local background correction for each spot
        - Identifies and filters miRNA probes specifically
        - Calculates quality metrics for array assessment
        - Supports multiple intensity calculation methods
    
    Attributes:
        gpr_path (str): Path to the GPR file being processed
        metadata (dict): Experimental metadata extracted from file header
        data (pd.DataFrame): Complete spot data from the array
        mirna_data (pd.DataFrame): Filtered data containing only miRNA probes
        
    Example Usage:
        >>> parser = GPRParser('GSM123_sample.gpr.gz')
        >>> parser.parse_file()
        >>> intensities = parser.calculate_intensities(method='ratio')
        >>> quality = parser.quality_metrics()
        
    Note:
        The GSE153135 dataset uses a specific pooling strategy where 26 donors
        were distributed across 10 arrays (2-3 donors per array). This affects
        the interpretation of biological replicates and statistical power.
    """
    
    def __init__(self, gpr_path: Union[str, Path]):
        """
        Initialize the GPR parser with a file path.
        
        Args:
            gpr_path (Union[str, Path]): Path to the GPR file (.gpr or .gpr.gz)
                The file can be either plain text or gzip-compressed.
                
        Note:
            The parser automatically detects compressed files based on the
            .gz extension and handles them appropriately.
        """
        self.gpr_path = Path(gpr_path)
        self.metadata = {}  # Will store experimental metadata
        self.data = None    # Will store complete spot data
        self.mirna_data = None  # Will store filtered miRNA data
        
        # Validate file existence
        if not self.gpr_path.exists():
            raise FileNotFoundError(f"GPR file not found: {self.gpr_path}")
            
        logger.info(f"Initialized GPR parser for: {self.gpr_path.name}")
    
    def parse_file(self) -> 'GPRParser':
        """
        Parse the complete GPR file including metadata and spot data.
        
        This method orchestrates the complete parsing process:
        1. Reads the file (handling compression automatically)
        2. Extracts metadata from the header section
        3. Parses the tabular spot data
        4. Filters for miRNA probes specifically
        
        Returns:
            GPRParser: Self-reference for method chaining
            
        Raises:
            ValueError: If file format is not recognized as valid GPR
            UnicodeDecodeError: If file encoding is incompatible
            
        Note:
            GPR files have a specific structure:
            - Header lines with metadata (key=value pairs)
            - Column headers starting with "Block"
            - Tab-separated data rows
            
            This method is designed to handle variations in this format
            while maintaining compatibility with GenePix software output.
        """
        logger.info(f"Parsing GPR file: {self.gpr_path.name}")
        
        try:
            # Read file contents (auto-detect compression)
            if self.gpr_path.suffix == '.gz':
                with gzip.open(self.gpr_path, 'rt', encoding='utf-8') as f:
                    lines = f.readlines()
            else:
                with open(self.gpr_path, 'r', encoding='utf-8') as f:
                    lines = f.readlines()
                    
        except UnicodeDecodeError:
            # Try alternative encoding if UTF-8 fails
            logger.warning("UTF-8 decoding failed, trying latin-1 encoding")
            if self.gpr_path.suffix == '.gz':
                with gzip.open(self.gpr_path, 'rt', encoding='latin-1') as f:
                    lines = f.readlines()
            else:
                with open(self.gpr_path, 'r', encoding='latin-1') as f:
                    lines = f.readlines()
        
        # Validate minimum file structure
        if len(lines) < 10:
            raise ValueError(f"File appears too short to be a valid GPR file: {len(lines)} lines")
        
        # Parse metadata section and find data start
        data_start_idx = self._parse_metadata(lines)
        
        # Parse the tabular data section
        self._parse_data(lines, data_start_idx)
        
        # Filter to miRNA probes only
        self._filter_mirna_spots()
        
        logger.info(f"Successfully parsed {len(self.data)} total spots, "
                   f"{len(self.mirna_data)} miRNA spots")
        
        return self
    
    def _parse_metadata(self, lines: List[str]) -> int:
        """
        Extract experimental metadata from GPR file header.
        
        GPR files contain valuable metadata about the experimental conditions,
        scanner settings, and array information in the header section. This
        metadata is crucial for understanding data quality and experimental
        context.
        
        Args:
            lines (List[str]): All lines from the GPR file
            
        Returns:
            int: Index where tabular data begins (after metadata section)
            
        Extracted Metadata:
            - Product: Array type and version
            - Scanner: GenePix scanner model and version
            - Settings: Scanning parameters and PMT voltages
            - DateTime: When the scan was performed
            - Protocol: Experimental protocol information
            
        Note:
            Metadata format is typically key="value" pairs, but this can vary
            between GenePix software versions. The parser handles common
            variations while preserving all available information.
        """
        data_start_idx = 0
        metadata_lines_found = 0
        
        for i, line in enumerate(lines):
            line = line.strip()
            
            # Skip empty lines
            if not line:
                continue
                
            # Look for metadata lines (typically start with quotes)
            if line.startswith('"') and '=' in line:
                try:
                    # Parse key=value pairs, handling various quote formats
                    key_value = line.strip('"')
                    if '=' in key_value:
                        key, value = key_value.split('=', 1)
                        # Clean up key and value
                        key = key.strip().strip('"')
                        value = value.strip().strip('"')
                        self.metadata[key] = value
                        metadata_lines_found += 1
                except Exception as e:
                    logger.debug(f"Could not parse metadata line {i}: {line[:50]}... ({e})")
                    
            # Look for data section start (column headers)
            elif line.startswith('"Block"') or 'Block' in line:
                data_start_idx = i
                logger.info(f"Found data section at line {i}")
                break
        
        # Log important metadata for debugging
        logger.info(f"Extracted {metadata_lines_found} metadata fields")
        if 'Product' in self.metadata:
            logger.info(f"Array type: {self.metadata['Product']}")
        if 'Scanner' in self.metadata:
            logger.info(f"Scanner: {self.metadata['Scanner']}")
        if 'DateTime' in self.metadata:
            logger.info(f"Scan date: {self.metadata['DateTime']}")
            
        # Validate that we found the data section
        if data_start_idx == 0:
            logger.warning("Could not find data section start, assuming line 0")
            
        return data_start_idx
    
    def _parse_data(self, lines: List[str], start_idx: int) -> None:
        """
        Parse the tabular data section containing spot measurements.
        
        The data section contains the actual microarray measurements with one
        row per spot and columns for various intensity statistics. This includes
        foreground and background measurements for both Cy5 and Cy3 channels.
        
        Args:
            lines (List[str]): All lines from the GPR file
            start_idx (int): Line index where data section begins
            
        Raises:
            ValueError: If required columns are missing
            
        Data Columns (typical GPR format):
            - Spatial: Block, Column, Row (spot coordinates)
            - Identity: Name, ID (probe identifiers) 
            - Cy5 (635nm): F635 Median/Mean/SD, B635 Median/Mean/SD
            - Cy3 (532nm): F532 Median/Mean/SD, B532 Median/Mean/SD
            - Quality: Flags, Normalize, LogRatio, etc.
            
        Note:
            The parser focuses on median intensities as they are more robust
            to outliers than means, which is important for small spots that
            may have dust or other artifacts.
        """
        if start_idx >= len(lines):
            raise ValueError("Data section start index exceeds file length")
        
        # Parse column headers
        header_line = lines[start_idx].strip()
        headers = header_line.split('\t')
        # Clean headers (remove quotes and whitespace)
        headers = [h.strip().strip('"') for h in headers]
        
        logger.info(f"Found {len(headers)} data columns")
        logger.debug(f"Column headers: {headers[:10]}...")  # Log first 10 for debugging
        
        # Validate required columns are present
        missing_cols = [col for col in REQUIRED_COLUMNS if col not in headers]
        if missing_cols:
            logger.warning(f"Missing expected columns: {missing_cols}")
            # Don't fail completely, as some GPR variants may have different names
        
        # Parse data rows
        data_rows = []
        skipped_rows = 0
        
        for line_num, line in enumerate(lines[start_idx + 1:], start=start_idx + 1):
            line = line.strip()
            
            # Skip empty lines
            if not line:
                continue
                
            # Split and clean fields
            fields = line.split('\t')
            fields = [f.strip().strip('"') for f in fields]
            
            # Validate row length matches headers
            if len(fields) != len(headers):
                logger.debug(f"Row {line_num}: field count mismatch "
                           f"({len(fields)} vs {len(headers)})")
                skipped_rows += 1
                continue
                
            data_rows.append(fields)
        
        if skipped_rows > 0:
            logger.warning(f"Skipped {skipped_rows} malformed data rows")
        
        # Create DataFrame
        if not data_rows:
            raise ValueError("No valid data rows found in GPR file")
            
        self.data = pd.DataFrame(data_rows, columns=headers)
        logger.info(f"Loaded {len(self.data)} data rows")
        
        # Convert numeric columns to appropriate types
        self._convert_numeric_columns()
        
        # Basic data validation
        self._validate_data_quality()
    
    def _convert_numeric_columns(self) -> None:
        """
        Convert string columns to numeric types where appropriate.
        
        GPR files store all data as strings, but intensity measurements need
        to be converted to numeric types for mathematical operations. This
        method identifies and converts the relevant columns while handling
        missing or invalid values gracefully.
        """
        # Define columns that should be numeric
        numeric_patterns = [
            'F635', 'F532',  # Foreground intensities
            'B635', 'B532',  # Background intensities
            'Block', 'Column', 'Row',  # Spatial coordinates
            'X', 'Y',  # Physical coordinates
            'Dia',  # Spot diameter
            'LogRatio', 'LogPixSD'  # Calculated ratios
        ]
        
        # Find matching columns in our data
        numeric_cols = []
        for col in self.data.columns:
            if any(pattern in col for pattern in numeric_patterns):
                numeric_cols.append(col)
        
        logger.debug(f"Converting {len(numeric_cols)} columns to numeric")
        
        # Convert with error handling
        conversion_errors = 0
        for col in numeric_cols:
            try:
                # Convert to numeric, setting errors to NaN
                original_type = self.data[col].dtype
                self.data[col] = pd.to_numeric(self.data[col], errors='coerce')
                
                # Count conversion failures
                na_count = self.data[col].isna().sum()
                if na_count > 0:
                    logger.debug(f"Column {col}: {na_count} values could not be converted")
                    conversion_errors += na_count
                    
            except Exception as e:
                logger.warning(f"Failed to convert column {col}: {e}")
        
        if conversion_errors > 0:
            logger.info(f"Total numeric conversion errors: {conversion_errors}")
    
    def _validate_data_quality(self) -> None:
        """
        Perform basic validation of the loaded data quality.
        
        This method checks for common data quality issues that could affect
        downstream analysis, including missing values, extreme outliers,
        and suspicious patterns that might indicate scanning or processing errors.
        """
        # Check for completely missing intensity columns
        required_intensity_cols = ['F635 Median', 'F532 Median', 'B635 Median', 'B532 Median']
        missing_intensity_cols = [col for col in required_intensity_cols 
                                 if col not in self.data.columns]
        
        if missing_intensity_cols:
            raise ValueError(f"Missing required intensity columns: {missing_intensity_cols}")
        
        # Check for excessive missing values
        for col in required_intensity_cols:
            na_count = self.data[col].isna().sum()
            na_percent = (na_count / len(self.data)) * 100
            
            if na_percent > 50:
                logger.warning(f"Column {col} has {na_percent:.1f}% missing values")
            elif na_percent > 10:
                logger.info(f"Column {col} has {na_percent:.1f}% missing values")
        
        # Check for reasonable intensity ranges
        for col in ['F635 Median', 'F532 Median']:
            if col in self.data.columns:
                values = self.data[col].dropna()
                if len(values) > 0:
                    min_val, max_val = values.min(), values.max()
                    median_val = values.median()
                    
                    # Log intensity statistics
                    logger.debug(f"{col}: min={min_val}, median={median_val}, max={max_val}")
                    
                    # Check for suspicious patterns
                    if min_val == max_val:
                        logger.warning(f"{col}: All values are identical ({min_val})")
                    elif min_val < 0:
                        logger.warning(f"{col}: Contains negative values (min={min_val})")
    
    def _filter_mirna_spots(self) -> None:
        """
        Filter the data to retain only miRNA-specific probes.
        
        GPR files contain various types of probes including miRNAs, controls,
        blanks, and other non-coding RNAs. For forensic miRNA analysis, we
        need to focus specifically on miRNA probes while excluding controls
        and other probe types.
        
        Filtering Strategy:
            1. Identify probes with 'miR' in the name (primary filter)
            2. Exclude empty or missing probe names
            3. Exclude known control probe patterns
            4. Validate that we retain a reasonable number of miRNAs
            
        Note:
            The filtering is conservative to avoid excluding valid miRNAs
            while removing obvious non-miRNA probes. The exact number of
            miRNAs expected depends on the array version used.
        """
        if self.data is None:
            raise ValueError("No data loaded. Run parse_file() first.")
        
        original_count = len(self.data)
        
        # Primary filter: look for miRNA indicators in the Name column
        if 'Name' not in self.data.columns:
            raise ValueError("No 'Name' column found for miRNA identification")
        
        # Filter for miRNA probes using the pattern
        mirna_mask = self.data['Name'].str.contains(MIRNA_PATTERN, case=False, na=False)
        self.mirna_data = self.data[mirna_mask].copy()
        
        # Remove rows with empty names (additional quality filter)
        empty_name_mask = (self.mirna_data['Name'] == '') | self.mirna_data['Name'].isna()
        self.mirna_data = self.mirna_data[~empty_name_mask]
        
        # Optional: Remove known control patterns
        control_patterns = ['control', 'blank', 'empty', 'neg_']
        for pattern in control_patterns:
            control_mask = self.mirna_data['Name'].str.contains(pattern, case=False, na=False)
            removed_controls = control_mask.sum()
            if removed_controls > 0:
                logger.debug(f"Removing {removed_controls} {pattern} probes")
                self.mirna_data = self.mirna_data[~control_mask]
        
        final_count = len(self.mirna_data)
        filtered_count = original_count - final_count
        
        logger.info(f"miRNA filtering: {original_count} → {final_count} spots "
                   f"({filtered_count} filtered out)")
        
        # Validate results
        if final_count == 0:
            raise ValueError("No miRNA probes found. Check probe naming convention.")
        elif final_count < 100:
            logger.warning(f"Only {final_count} miRNA probes found. "
                          "This seems low for a miRNA array.")
        
        # Log some example miRNA names for verification
        example_names = self.mirna_data['Name'].head(5).tolist()
        logger.debug(f"Example miRNA names: {example_names}")
    
    def calculate_intensities(self, method: str = DEFAULT_METHOD) -> pd.DataFrame:
        """
        Calculate normalized intensities from raw fluorescence measurements.
        
        This is the core processing step that converts raw scanner measurements
        into normalized intensity values suitable for statistical analysis. The
        method includes background correction and appropriate transformations
        for two-color array data.
        
        Args:
            method (str): Calculation method for final intensities
                - 'ratio': Log2(Cy5/Cy3) after background correction (default)
                  This is standard for two-color arrays as it normalizes for
                  probe-specific binding effects and dye bias
                - 'cy5': Log2(Cy5) background-corrected, sample channel only
                - 'cy3': Log2(Cy3) background-corrected, reference channel only
                
        Returns:
            pd.DataFrame: Processed intensities with columns:
                - miRNA: miRNA probe identifier
                - intensity: Final calculated intensity value
                - cy5_raw: Raw Cy5 foreground intensity
                - cy3_raw: Raw Cy3 foreground intensity
                - cy5_bg: Cy5 background intensity
                - cy3_bg: Cy3 background intensity
                - cy5_corrected: Background-corrected Cy5
                - cy3_corrected: Background-corrected Cy3
                
        Raises:
            ValueError: If no miRNA data available or invalid method specified
            
        Background Correction:
            The local background correction subtracts the local background
            intensity (measured around each spot) from the foreground intensity.
            This removes non-specific signal and spatial artifacts on the array.
            
            Formula: corrected = foreground_median - background_median
            
        Log2 Transformation:
            Log2 transformation is applied because:
            1. It normalizes the data distribution
            2. Makes fold changes symmetric (2x up = -2x down)
            3. Is standard practice in microarray analysis
            4. Enables parametric statistical tests
            
        Note:
            The 'ratio' method is recommended for most analyses as it corrects
            for probe-specific effects and provides the most robust normalized
            measurements for comparative analysis.
        """
        if self.mirna_data is None:
            raise ValueError("No miRNA data available. Run parse_file() first.")
        
        # Validate method parameter
        valid_methods = ['ratio', 'cy5', 'cy3']
        if method not in valid_methods:
            raise ValueError(f"Method must be one of {valid_methods}, got '{method}'")
        
        logger.info(f"Calculating intensities using method: {method}")
        
        # Extract required columns
        required_cols = ['F635 Median', 'F532 Median', 'B635 Median', 'B532 Median']
        missing_cols = [col for col in required_cols if col not in self.mirna_data.columns]
        if missing_cols:
            raise ValueError(f"Missing required columns for intensity calculation: {missing_cols}")
        
        # Background correction
        # Subtract local background from foreground signal for each channel
        cy5_corrected = self.mirna_data['F635 Median'] - self.mirna_data['B635 Median']
        cy3_corrected = self.mirna_data['F532 Median'] - self.mirna_data['B532 Median']
        
        # Handle negative values after background correction
        # Negative values can occur when background is higher than foreground
        # This typically indicates very low or absent signal
        negative_cy5 = (cy5_corrected <= 0).sum()
        negative_cy3 = (cy3_corrected <= 0).sum()
        
        if negative_cy5 > 0:
            logger.debug(f"Found {negative_cy5} spots with negative/zero Cy5 after background correction")
        if negative_cy3 > 0:
            logger.debug(f"Found {negative_cy3} spots with negative/zero Cy3 after background correction")
        
        # Set minimum intensity to prevent log(0) errors
        # This is a standard approach in microarray analysis
        cy5_corrected = np.maximum(cy5_corrected, MIN_INTENSITY)
        cy3_corrected = np.maximum(cy3_corrected, MIN_INTENSITY)
        
        # Calculate final intensities based on selected method
        if method == 'ratio':
            # Standard two-color analysis: log2(Cy5/Cy3)
            # This normalizes for probe-specific binding differences
            intensities = np.log2(cy5_corrected / cy3_corrected)
            logger.debug(f"Ratio method: median log2(Cy5/Cy3) = {np.median(intensities):.2f}")
            
        elif method == 'cy5':
            # Sample channel only (useful for single-channel analysis)
            intensities = np.log2(cy5_corrected)
            logger.debug(f"Cy5 method: median log2(Cy5) = {np.median(intensities):.2f}")
            
        elif method == 'cy3':
            # Reference channel only (for reference analysis)
            intensities = np.log2(cy3_corrected)
            logger.debug(f"Cy3 method: median log2(Cy3) = {np.median(intensities):.2f}")
        
        # Create comprehensive result DataFrame
        result = pd.DataFrame({
            'miRNA': self.mirna_data['Name'].values,
            'intensity': intensities,
            'cy5_raw': self.mirna_data['F635 Median'].values,
            'cy3_raw': self.mirna_data['F532 Median'].values,
            'cy5_bg': self.mirna_data['B635 Median'].values,
            'cy3_bg': self.mirna_data['B532 Median'].values,
            'cy5_corrected': cy5_corrected.values,
            'cy3_corrected': cy3_corrected.values
        })
        
        # Remove infinite values (can occur if corrected intensity is 0)
        infinite_mask = np.isinf(result['intensity'])
        infinite_count = infinite_mask.sum()
        
        if infinite_count > 0:
            logger.warning(f"Removing {infinite_count} probes with infinite intensity values")
            result = result[~infinite_mask]
        
        # Basic quality checks on final intensities
        self._validate_intensity_results(result, method)
        
        logger.info(f"Successfully calculated {len(result)} miRNA intensities")
        
        return result
    
    def _validate_intensity_results(self, result: pd.DataFrame, method: str) -> None:
        """
        Validate the calculated intensity results for quality issues.
        
        Args:
            result (pd.DataFrame): Calculated intensity results
            method (str): Method used for calculation
        """
        intensities = result['intensity']
        
        # Basic statistics
        stats = {
            'count': len(intensities),
            'mean': intensities.mean(),
            'median': intensities.median(),
            'std': intensities.std(),
            'min': intensities.min(),
            'max': intensities.max()
        }
        
        logger.debug(f"Intensity statistics ({method}): {stats}")
        
        # Check for suspicious patterns
        if stats['std'] < 0.1:
            logger.warning("Very low intensity variance - check for processing issues")
        
        if method == 'ratio':
            # For ratios, expect values roughly centered around 0
            if abs(stats['median']) > 2:
                logger.warning(f"Ratio median far from 0: {stats['median']:.2f}")
        
        # Check dynamic range
        dynamic_range = stats['max'] - stats['min']
        if dynamic_range < 3:
            logger.warning(f"Low dynamic range: {dynamic_range:.2f}")
        elif dynamic_range > 20:
            logger.warning(f"Very high dynamic range: {dynamic_range:.2f}")
    
    def quality_metrics(self) -> Optional[Dict[str, Union[int, float]]]:
        """
        Calculate comprehensive quality metrics for the microarray data.
        
        These metrics help assess the technical quality of the array scan
        and can identify potential issues that might affect downstream analysis.
        The metrics are designed for forensic applications where data quality
        is crucial for reliable identification.
        
        Returns:
            Optional[Dict[str, Union[int, float]]]: Quality metrics including:
                - total_spots: Total number of spots on the array
                - mirna_spots: Number of miRNA-specific spots
                - mirna_percentage: Percentage of spots that are miRNAs
                - cy5_median: Median Cy5 foreground intensity
                - cy3_median: Median Cy3 foreground intensity
                - cy5_bg_median: Median Cy5 background intensity
                - cy3_bg_median: Median Cy3 background intensity
                - signal_to_noise_cy5: Ratio of foreground to background (Cy5)
                - signal_to_noise_cy3: Ratio of foreground to background (Cy3)
                - dynamic_range_cy5: Range of Cy5 intensities
                - dynamic_range_cy3: Range of Cy3 intensities
                
        Returns None if miRNA data is not available.
        
        Quality Thresholds (typical for miRNA arrays):
            - Signal-to-noise ratio: >3 is good, >10 is excellent
            - Background levels: Should be <1000 for most scanners
            - Dynamic range: Should span 3-4 orders of magnitude
            - miRNA percentage: Depends on array design (10-90% typical)
        """
        if self.mirna_data is None:
            logger.warning("No miRNA data available for quality metrics")
            return None
        
        logger.info("Calculating array quality metrics")
        
        # Basic count metrics
        total_spots = len(self.data) if self.data is not None else 0
        mirna_spots = len(self.mirna_data)
        mirna_percentage = (mirna_spots / total_spots * 100) if total_spots > 0 else 0
        
        # Intensity metrics (using miRNA spots only)
        intensity_cols = ['F635 Median', 'F532 Median', 'B635 Median', 'B532 Median']
        metrics = {
            'total_spots': total_spots,
            'mirna_spots': mirna_spots,
            'mirna_percentage': mirna_percentage
        }
        
        # Calculate intensity statistics
        for col in intensity_cols:
            if col in self.mirna_data.columns:
                values = self.mirna_data[col].dropna()
                if len(values) > 0:
                    col_name = col.lower().replace(' ', '_')
                    metrics[col_name] = values.median()
                    metrics[f'{col_name}_mean'] = values.mean()
                    metrics[f'{col_name}_std'] = values.std()
        
        # Calculate signal-to-noise ratios
        if all(col in metrics for col in ['f635_median', 'b635_median']):
            if metrics['b635_median'] > 0:
                metrics['signal_to_noise_cy5'] = metrics['f635_median'] / metrics['b635_median']
            else:
                metrics['signal_to_noise_cy5'] = float('inf')
        
        if all(col in metrics for col in ['f532_median', 'b532_median']):
            if metrics['b532_median'] > 0:
                metrics['signal_to_noise_cy3'] = metrics['f532_median'] / metrics['b532_median']
            else:
                metrics['signal_to_noise_cy3'] = float('inf')
        
        # Calculate dynamic ranges
        for channel in ['F635 Median', 'F532 Median']:
            if channel in self.mirna_data.columns:
                values = self.mirna_data[channel].dropna()
                if len(values) > 0:
                    channel_name = channel.lower().replace(' ', '_')
                    metrics[f'dynamic_range_{channel_name}'] = values.max() - values.min()
        
        # Log quality assessment
        self._assess_quality(metrics)
        
        return metrics
    
    def _assess_quality(self, metrics: Dict[str, Union[int, float]]) -> None:
        """
        Assess and log the quality of the array data based on calculated metrics.
        
        Args:
            metrics (Dict): Quality metrics dictionary
        """
        logger.info("=== Array Quality Assessment ===")
        
        # Assess signal-to-noise ratios
        for channel in ['cy5', 'cy3']:
            snr_key = f'signal_to_noise_{channel}'
            if snr_key in metrics:
                snr = metrics[snr_key]
                if snr > 10:
                    quality = "Excellent"
                elif snr > 3:
                    quality = "Good"
                elif snr > 1.5:
                    quality = "Acceptable"
                else:
                    quality = "Poor"
                
                logger.info(f"{channel.upper()} signal-to-noise: {snr:.2f} ({quality})")
        
        # Assess miRNA representation
        mirna_pct = metrics.get('mirna_percentage', 0)
        logger.info(f"miRNA probe percentage: {mirna_pct:.1f}%")
        
        # Assess background levels
        for channel in ['Cy5', 'Cy3']:
            bg_key = f'b{635 if channel == "Cy5" else 532}_median'
            if bg_key in metrics:
                bg_level = metrics[bg_key]
                if bg_level < 100:
                    bg_quality = "Low (good)"
                elif bg_level < 500:
                    bg_quality = "Moderate"
                else:
                    bg_quality = "High (concerning)"
                    
                logger.info(f"{channel} background level: {bg_level:.0f} ({bg_quality})")
        
        logger.info("=== End Quality Assessment ===")


def process_all_gpr_files(input_dir: Union[str, Path], 
                         output_dir: Union[str, Path], 
                         method: str = DEFAULT_METHOD) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Process all GPR files in a directory and create combined expression matrix.
    
    This function orchestrates the batch processing of multiple GPR files,
    combining them into a single expression matrix suitable for downstream
    statistical analysis. It handles the specific requirements of the GSE153135
    dataset including proper sample naming and body fluid identification.
    
    Args:
        input_dir (Union[str, Path]): Directory containing GPR files (.gpr.gz)
        output_dir (Union[str, Path]): Directory for output files
        method (str): Intensity calculation method (default: 'ratio')
            See calculate_intensities() for method descriptions
            
    Returns:
        Tuple[pd.DataFrame, pd.DataFrame]: 
            - Expression matrix (miRNAs x samples)
            - Sample metadata with body fluid annotations
            
    Creates:
        - gpr_expression_matrix.csv: Combined expression matrix
        - gpr_metadata.csv: Sample annotations and quality metrics
        - individual_samples/: Individual processed files (optional)
        
    Processing Steps:
        1. Discover all GPR files in input directory
        2. Parse each file and extract intensities
        3. Identify body fluid type from filename
        4. Calculate quality metrics for each array
        5. Combine into unified expression matrix
        6. Generate comprehensive metadata
        
    Sample Naming Convention (GSE153135):
        Files are named with patterns like:
        - "GSM*_saliva_*.gpr.gz"
        - "GSM*_peripheral_blood_*.gpr.gz"
        - "GSM*_menstrual_blood_*.gpr.gz"
        - "GSM*_semen_*.gpr.gz"
        - "GSM*_vaginal_*.gpr.gz"
        
    Note:
        This function assumes the GSE153135 naming convention. For other
        datasets, the body fluid identification logic may need modification.
        
    Quality Control:
        Each array is assessed for technical quality, and results are included
        in the metadata. Arrays with poor quality metrics are flagged but not
        automatically excluded to preserve statistical power.
    """
    input_path = Path(input_dir)
    output_path = Path(output_dir)
    
    # Validate inputs
    if not input_path.exists():
        raise FileNotFoundError(f"Input directory not found: {input_path}")
    
    # Create output directory
    output_path.mkdir(parents=True, exist_ok=True)
    
    # Find all GPR files
    gpr_files = list(input_path.glob("*.gpr.gz")) + list(input_path.glob("*.gpr"))
    if not gpr_files:
        raise ValueError(f"No GPR files found in {input_path}")
    
    logger.info(f"Found {len(gpr_files)} GPR files in {input_path}")
    
    # Process each file
    all_intensities = []
    sample_metadata = []
    processing_errors = []
    
    for i, gpr_file in enumerate(gpr_files, 1):
        logger.info(f"Processing file {i}/{len(gpr_files)}: {gpr_file.name}")
        
        try:
            # Parse the GPR file
            parser = GPRParser(gpr_file)
            parser.parse_file()
            
            # Calculate intensities
            intensities = parser.calculate_intensities(method=method)
            
            # Generate sample name from filename
            sample_name = gpr_file.stem
            if sample_name.endswith('.gpr'):
                sample_name = sample_name[:-4]  # Remove .gpr extension
                
            # Extract body fluid type from filename
            body_fluid = extract_body_fluid_from_filename(gpr_file.name)
            
            # Calculate quality metrics
            quality_metrics = parser.quality_metrics() or {}
            
            # Add sample information to intensities
            intensities['sample'] = sample_name
            all_intensities.append(intensities)
            
            # Collect sample metadata
            metadata_entry = {
                'sample': sample_name,
                'filename': gpr_file.name,
                'body_fluid': body_fluid,
                'method': method,
                'mirna_count': len(intensities),
                **quality_metrics  # Include all quality metrics
            }
            sample_metadata.append(metadata_entry)
            
            logger.info(f"  Processed: {len(intensities)} miRNAs, "
                       f"body_fluid={body_fluid}")
            
        except Exception as e:
            error_msg = f"Failed to process {gpr_file.name}: {str(e)}"
            logger.error(error_msg)
            processing_errors.append(error_msg)
            continue
    
    # Check if we have any successful processing
    if not all_intensities:
        error_summary = "\n".join(processing_errors)
        raise RuntimeError(f"No files could be processed successfully.\nErrors:\n{error_summary}")
    
    # Report processing summary
    logger.info(f"Successfully processed {len(all_intensities)}/{len(gpr_files)} files")
    if processing_errors:
        logger.warning(f"Failed to process {len(processing_errors)} files")
    
    # Combine intensities into expression matrix
    logger.info("Creating combined expression matrix...")
    expression_matrix = create_expression_matrix(all_intensities)
    
    # Create metadata DataFrame
    metadata_df = pd.DataFrame(sample_metadata)
    
    # Save results
    matrix_path = output_path / 'gpr_expression_matrix.csv'
    metadata_path = output_path / 'gpr_metadata.csv'
    
    expression_matrix.to_csv(matrix_path)
    metadata_df.to_csv(metadata_path, index=False)
    
    logger.info(f"Expression matrix saved: {matrix_path}")
    logger.info(f"Sample metadata saved: {metadata_path}")
    logger.info(f"Matrix dimensions: {expression_matrix.shape}")
    
    # Log summary statistics
    log_processing_summary(expression_matrix, metadata_df)
    
    return expression_matrix, metadata_df


def extract_body_fluid_from_filename(filename: str) -> str:
    """
    Extract body fluid type from GSE153135 filename patterns.
    
    Args:
        filename (str): GPR filename from GSE153135
        
    Returns:
        str: Body fluid type or 'unknown' if not identified
        
    GSE153135 Naming Patterns:
        The dataset uses descriptive filenames that include body fluid
        information. This function recognizes the standard patterns used
        in the dataset submission.
    """
    filename_lower = filename.lower()
    
    # Check each fluid pattern
    for fluid, patterns in FLUID_PATTERNS.items():
        if any(pattern in filename_lower for pattern in patterns):
            return fluid
    
    # Special handling for blood subtypes
    if 'blood' in filename_lower:
        if 'menstrual' in filename_lower:
            return 'menstrual_blood'
        else:
            return 'peripheral_blood'
    
    logger.warning(f"Could not identify body fluid from filename: {filename}")
    return 'unknown'


def create_expression_matrix(all_intensities: List[pd.DataFrame]) -> pd.DataFrame:
    """
    Combine individual sample intensities into a unified expression matrix.
    
    Args:
        all_intensities (List[pd.DataFrame]): List of intensity DataFrames
            from individual samples
            
    Returns:
        pd.DataFrame: Combined expression matrix with miRNAs as rows and
            samples as columns
            
    Note:
        This function handles missing miRNAs across samples by using an
        outer join strategy, filling missing values with NaN. This preserves
        the maximum amount of data while allowing downstream analyses to
        handle missing values appropriately.
    """
    if not all_intensities:
        raise ValueError("No intensity data provided")
    
    # Extract intensity data for matrix creation
    sample_matrices = []
    
    for intensity_df in all_intensities:
        # Create a sample-specific matrix
        sample_name = intensity_df['sample'].iloc[0]
        sample_matrix = intensity_df.set_index('miRNA')[['intensity']].copy()
        sample_matrix.columns = [sample_name]
        sample_matrices.append(sample_matrix)
    
    # Combine all samples using outer join to preserve all miRNAs
    logger.info("Combining intensity matrices...")
    expression_matrix = sample_matrices[0]
    
    for sample_matrix in sample_matrices[1:]:
        expression_matrix = expression_matrix.join(sample_matrix, how='outer')
    
    # Log matrix statistics
    logger.info(f"Combined matrix: {expression_matrix.shape[0]} miRNAs x "
               f"{expression_matrix.shape[1]} samples")
    
    missing_percent = (expression_matrix.isna().sum().sum() / 
                      (expression_matrix.shape[0] * expression_matrix.shape[1]) * 100)
    logger.info(f"Missing values: {missing_percent:.1f}%")
    
    return expression_matrix


def log_processing_summary(expression_matrix: pd.DataFrame, 
                          metadata_df: pd.DataFrame) -> None:
    """
    Log a comprehensive summary of the processing results.
    
    Args:
        expression_matrix (pd.DataFrame): Final expression matrix
        metadata_df (pd.DataFrame): Sample metadata
    """
    logger.info("=== GPR Processing Summary ===")
    
    # Sample distribution by body fluid
    if 'body_fluid' in metadata_df.columns:
        fluid_counts = metadata_df['body_fluid'].value_counts()
        logger.info("Sample distribution by body fluid:")
        for fluid, count in fluid_counts.items():
            logger.info(f"  {fluid}: {count} samples")
    
    # Quality summary
    quality_cols = [col for col in metadata_df.columns 
                   if any(term in col.lower() for term in ['signal', 'background', 'snr'])]
    
    if quality_cols:
        logger.info("Quality metrics summary:")
        for col in quality_cols[:3]:  # Show first 3 quality metrics
            if metadata_df[col].dtype in ['float64', 'int64']:
                mean_val = metadata_df[col].mean()
                logger.info(f"  {col}: mean = {mean_val:.2f}")
    
    # Expression matrix summary
    logger.info(f"Expression matrix: {expression_matrix.shape}")
    logger.info(f"Intensity range: {expression_matrix.min().min():.2f} to "
               f"{expression_matrix.max().max():.2f}")
    
    logger.info("=== End Processing Summary ===")


def main():
    """
    Main function for command-line usage of the GPR parser.
    
    Processes GPR files from the command line with configurable options.
    This provides a convenient interface for batch processing of GPR files
    in forensic miRNA analysis pipelines.
    """
    import argparse
    
    parser = argparse.ArgumentParser(
        description="Process GenePix Result (GPR) files for forensic miRNA analysis"
    )
    parser.add_argument('input_dir', help='Directory containing GPR files')
    parser.add_argument('output_dir', help='Directory for output files') 
    parser.add_argument('--method', choices=['ratio', 'cy5', 'cy3'], 
                       default='ratio', help='Intensity calculation method')
    parser.add_argument('--verbose', '-v', action='store_true',
                       help='Enable verbose logging')
    
    args = parser.parse_args()
    
    # Configure logging level
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    
    try:
        # Process all files
        expression_matrix, metadata_df = process_all_gpr_files(
            args.input_dir, args.output_dir, args.method
        )
        
        logger.info("GPR processing completed successfully")
        
    except Exception as e:
        logger.error(f"GPR processing failed: {e}")
        raise


if __name__ == "__main__":
    main()