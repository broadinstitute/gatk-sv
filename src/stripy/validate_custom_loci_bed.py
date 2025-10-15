#!/usr/bin/env python3
"""
Stripy Custom Loci BED Format Validator

Validates TSV files for use with Stripy's custom loci feature.
Provides detailed error messages for format violations.
"""

import sys
import csv
import argparse
from pathlib import Path
from typing import List, Dict, Tuple, Optional


class ValidationError:
    """Represents a validation error with context."""
    
    def __init__(self, line_num: int, severity: str, message: str, data: Optional[str] = None):
        self.line_num = line_num
        self.severity = severity  # 'ERROR', 'WARNING', 'INFO'
        self.message = message
        self.data = data
    
    def __str__(self):
        prefix = f"Line {self.line_num}: [{self.severity}]"
        if self.data:
            return f"{prefix} {self.message}\n  Data: {self.data}"
        return f"{prefix} {self.message}"


class StripyBEDValidator:
    """Validates Stripy custom loci BED files."""
    
    def __init__(self, strict_mode: bool = True):
        self.strict_mode = strict_mode
        self.errors: List[ValidationError] = []
        self.warnings: List[ValidationError] = []
        self.info: List[ValidationError] = []
        self.stats = {
            'total_lines': 0,
            'valid_lines': 0,
            'format_minimal': 0,
            'format_standard': 0,
            'format_extended': 0
        }
    
    def add_error(self, line_num: int, message: str, data: Optional[str] = None):
        """Add an error message."""
        self.errors.append(ValidationError(line_num, 'ERROR', message, data))
    
    def add_warning(self, line_num: int, message: str, data: Optional[str] = None):
        """Add a warning message."""
        self.warnings.append(ValidationError(line_num, 'WARNING', message, data))
    
    def add_info(self, line_num: int, message: str, data: Optional[str] = None):
        """Add an info message."""
        self.info.append(ValidationError(line_num, 'INFO', message, data))
    
    def validate_file(self, filepath: str) -> bool:
        """
        Validate a Stripy BED file.
        
        Returns:
            bool: True if file is valid, False otherwise
        """
        filepath = Path(filepath)
        
        # Check file exists
        if not filepath.exists():
            print(f"ERROR: File not found: {filepath}")
            return False
        
        # Check file is not empty
        if filepath.stat().st_size == 0:
            print(f"ERROR: File is empty: {filepath}")
            return False
        
        # Read and validate file
        try:
            with open(filepath, 'r', newline='', encoding='utf-8') as f:
                # Check for potential delimiter issues
                first_line = f.readline()
                f.seek(0)
                
                if ',' in first_line and '\t' not in first_line:
                    self.add_error(1, "File appears to be comma-delimited, but must be tab-delimited", 
                                 first_line.strip())
                
                # Detect potential header
                if self._looks_like_header(first_line):
                    self.add_warning(1, "First line looks like a header row. BED files should not have headers.",
                                   first_line.strip())
                
                # Parse as TSV
                reader = csv.reader(f, delimiter='\t')
                
                for line_num, row in enumerate(reader, start=1):
                    self.stats['total_lines'] += 1
                    
                    # Skip empty lines
                    if not row or all(cell.strip() == '' for cell in row):
                        self.add_info(line_num, "Empty line, skipping")
                        continue
                    
                    # Validate the row
                    if self._validate_row(row, line_num):
                        self.stats['valid_lines'] += 1
        
        except UnicodeDecodeError as e:
            print(f"ERROR: File encoding issue: {e}")
            return False
        except Exception as e:
            print(f"ERROR: Failed to read file: {e}")
            return False
        
        return len(self.errors) == 0
    
    def _looks_like_header(self, line: str) -> bool:
        """Check if line looks like a header."""
        line_lower = line.lower()
        header_terms = ['chromosome', 'chrom', 'start', 'end', 'motif', 
                       'locus', 'disease', 'inheritance']
        return any(term in line_lower for term in header_terms)
    
    def _validate_row(self, row: List[str], line_num: int) -> bool:
        """
        Validate a single row.
        
        Returns:
            bool: True if row is valid, False otherwise
        """
        row_valid = True
        raw_line = '\t'.join(row)
        
        # Check minimum columns
        if len(row) < 4:
            self.add_error(line_num, 
                         f"Insufficient columns: found {len(row)}, need at least 4 (chromosome, start, end, motif)",
                         raw_line)
            return False
        
        # Validate column count is sensible
        if len(row) > 9:
            self.add_warning(line_num,
                           f"Found {len(row)} columns, expected 4, 5, or 9. Extra columns will be ignored.",
                           raw_line)
        
        # Detect format type
        if len(row) == 4:
            self.stats['format_minimal'] += 1
        elif len(row) == 5:
            self.stats['format_standard'] += 1
        elif len(row) >= 9:
            self.stats['format_extended'] += 1
        elif 6 <= len(row) <= 8:
            self.add_error(line_num,
                         f"Partial extended format: found {len(row)} columns. "
                         "Extended format requires all 9 columns (chromosome, start, end, motif, locus_id, "
                         "disease, inheritance, normal_range, pathogenic_cutoff)",
                         raw_line)
            return False
        
        # Extract columns
        chromosome = row[0].strip()
        start_str = row[1].strip()
        end_str = row[2].strip()
        motif = row[3].strip()
        
        # Validate chromosome
        if not chromosome:
            self.add_error(line_num, "Chromosome (column 1) is empty", raw_line)
            row_valid = False
        
        # Validate start coordinate
        try:
            start = int(start_str)
            if start < 0:
                self.add_error(line_num, f"Start coordinate is negative: {start}", raw_line)
                row_valid = False
        except ValueError:
            self.add_error(line_num, f"Start coordinate (column 2) must be an integer, got: '{start_str}'", raw_line)
            return False
        
        # Validate end coordinate
        try:
            end = int(end_str)
            if end < 0:
                self.add_error(line_num, f"End coordinate is negative: {end}", raw_line)
                row_valid = False
        except ValueError:
            self.add_error(line_num, f"End coordinate (column 3) must be an integer, got: '{end_str}'", raw_line)
            return False
        
        # Validate start < end (Stripy requirement: end > start)
        if end <= start:
            self.add_error(line_num, 
                         f"End coordinate ({end}) must be greater than start coordinate ({start})",
                         raw_line)
            row_valid = False
        
        # Validate region length < 500 bp
        region_length = end - start
        if region_length >= 500:
            self.add_error(line_num,
                         f"Region length ({region_length} bp) must be less than 500 bp (start={start}, end={end})",
                         raw_line)
            row_valid = False
        
        # Validate motif
        if not motif:
            self.add_error(line_num, "Motif (column 4) is empty", raw_line)
            row_valid = False
        else:
            # Check if motif contains valid DNA characters
            valid_chars = set('ACGTNacgtn')
            invalid_chars = set(motif) - valid_chars
            if invalid_chars:
                self.add_warning(line_num,
                               f"Motif contains non-standard DNA characters: {invalid_chars}",
                               raw_line)
            
            # Check if lowercase
            if motif != motif.upper():
                self.add_info(line_num,
                            "Motif contains lowercase letters. Consider using uppercase for consistency.",
                            raw_line)
        
        # Validate optional columns
        if len(row) >= 5:
            locus_id = row[4].strip()
            if not locus_id:
                self.add_warning(line_num, "LocusID (column 5) is empty, will default to 'Custom'", raw_line)
        
        # Validate extended format columns
        if len(row) >= 9:
            disease = row[5].strip()
            inheritance = row[6].strip()
            normal_range = row[7].strip()
            pathogenic_cutoff = row[8].strip()
            
            if not disease:
                self.add_warning(line_num, "Disease (column 6) is empty", raw_line)
            
            if not inheritance:
                self.add_warning(line_num, "Inheritance (column 7) is empty", raw_line)
            else:
                # Check common inheritance patterns
                valid_patterns = ['AD', 'AR', 'XL', 'XLD', 'XLR', 'Digenic', 'Mitochondrial', 'NA']
                if inheritance not in valid_patterns:
                    self.add_info(line_num,
                                f"Inheritance pattern '{inheritance}' is non-standard. "
                                f"Common values: {', '.join(valid_patterns)}",
                                raw_line)
            
            # Validate normal range format
            if not normal_range:
                self.add_warning(line_num, "Normal range (column 8) is empty", raw_line)
            else:
                if '-' in normal_range:
                    parts = normal_range.split('-')
                    if len(parts) != 2:
                        self.add_error(line_num,
                                     f"Normal range format invalid: '{normal_range}'. "
                                     "Expected 'min-max' or single value",
                                     raw_line)
                        row_valid = False
                    else:
                        try:
                            min_val = int(parts[0])
                            max_val = int(parts[1])
                            if min_val > max_val:
                                self.add_error(line_num,
                                             f"Normal range minimum ({min_val}) is greater than maximum ({max_val})",
                                             raw_line)
                                row_valid = False
                        except ValueError:
                            self.add_error(line_num,
                                         f"Normal range values must be integers: '{normal_range}'",
                                         raw_line)
                            row_valid = False
                else:
                    # Single value
                    try:
                        int(normal_range)
                    except ValueError:
                        self.add_error(line_num,
                                     f"Normal range must be integer or 'min-max' format: '{normal_range}'",
                                     raw_line)
                        row_valid = False
            
            # Validate pathogenic cutoff
            if not pathogenic_cutoff:
                self.add_warning(line_num, "Pathogenic cutoff (column 9) is empty", raw_line)
            else:
                try:
                    cutoff = int(pathogenic_cutoff)
                    if cutoff < 0:
                        self.add_warning(line_num,
                                       f"Pathogenic cutoff is negative: {cutoff}",
                                       raw_line)
                except ValueError:
                    self.add_error(line_num,
                                 f"Pathogenic cutoff (column 9) must be an integer, got: '{pathogenic_cutoff}'",
                                 raw_line)
                    row_valid = False
        
        return row_valid
    
    def print_report(self):
        """Print validation report."""
        print("\n" + "="*80)
        print("STRIPY BED VALIDATION REPORT")
        print("="*80)
        
        # Statistics
        print(f"\nStatistics:")
        print(f"  Total lines processed: {self.stats['total_lines']}")
        print(f"  Valid lines: {self.stats['valid_lines']}")
        print(f"  Minimal format (4 cols): {self.stats['format_minimal']}")
        print(f"  Standard format (5 cols): {self.stats['format_standard']}")
        print(f"  Extended format (9 cols): {self.stats['format_extended']}")
        
        # Errors
        if self.errors:
            print(f"\n{'ERRORS':-^80}")
            print(f"Found {len(self.errors)} error(s):\n")
            for error in self.errors:
                print(error)
                print()
        
        # Warnings
        if self.warnings:
            print(f"\n{'WARNINGS':-^80}")
            print(f"Found {len(self.warnings)} warning(s):\n")
            for warning in self.warnings:
                print(warning)
                print()
        
        # Info messages
        if self.info and not self.strict_mode:
            print(f"\n{'INFORMATION':-^80}")
            print(f"Found {len(self.info)} info message(s):\n")
            for info in self.info:
                print(info)
                print()
        
        # Final result
        print("\n" + "="*80)
        if not self.errors:
            print("RESULT: VALID - File passes all validation checks")
            if self.warnings:
                print(f"        Note: {len(self.warnings)} warning(s) found but not critical")
        else:
            print("RESULT: INVALID - File contains errors that must be fixed")
        print("="*80 + "\n")


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description='Validate Stripy custom loci BED files',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s my_loci.bed
  %(prog)s --strict my_loci.bed
  %(prog)s --no-info my_loci.bed

Expected BED format (tab-delimited, no header):
  - Minimal (4 columns): chromosome, start, end, motif
  - Standard (5 columns): + locus_id
  - Extended (9 columns): + disease, inheritance, normal_range, pathogenic_cutoff
        """
    )
    
    parser.add_argument('bed_file', help='Path to BED file to validate')
    parser.add_argument('--strict', action='store_true',
                       help='Treat warnings as errors')
    parser.add_argument('--no-info', action='store_true',
                       help='Suppress informational messages')
    
    args = parser.parse_args()
    
    # Run validation
    validator = StripyBEDValidator(strict_mode=args.strict)
    
    print(f"Validating: {args.bed_file}")
    is_valid = validator.validate_file(args.bed_file)
    
    # Print report
    validator.print_report()
    
    # Exit with appropriate code
    if args.strict and validator.warnings:
        sys.exit(1)
    elif not is_valid:
        sys.exit(1)
    else:
        sys.exit(0)


if __name__ == '__main__':
    main()