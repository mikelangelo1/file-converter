#!/usr/bin/env python3
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import argparse
import os

class FastaToFastqConverter:
    def __init__(self, quality_score=40):
        """
        Initialize converter with default quality score
        
        Parameters:
        quality_score (int): Phred quality score to assign (default: 40)
        """
        self.quality_score = quality_score
        self.phred_quality = chr(33 + quality_score)  # Convert to ASCII
        
    def generate_quality_string(self, sequence_length):
        """Generate quality string of given length"""
        return self.phred_quality * sequence_length
    
    def convert_file(self, input_file, output_file):
        """
        Convert FASTA file to FASTQ format
        
        Parameters:
        input_file (str): Path to input FASTA file
        output_file (str): Path to output FASTQ file
        """
        try:
            # Read FASTA records
            with open(input_file, 'r') as fasta_handle:
                # Open output file
                with open(output_file, 'w') as fastq_handle:
                    # Process each sequence
                    for record in SeqIO.parse(fasta_handle, 'fasta'):
                        # Generate quality scores
                        quality = self.generate_quality_string(len(record.seq))
                        
                        # Write FASTQ format
                        fastq_handle.write(f"@{record.id}\n")  # Header
                        fastq_handle.write(f"{str(record.seq)}\n")  # Sequence
                        fastq_handle.write("+\n")  # Quality header
                        fastq_handle.write(f"{quality}\n")  # Quality scores
            
            return True
        except Exception as e:
            print(f"Error converting file: {str(e)}")
            return False
    
    def convert_string(self, fasta_string):
        """
        Convert FASTA string to FASTQ format
        
        Parameters:
        fasta_string (str): FASTA formatted string
        
        Returns:
        str: FASTQ formatted string
        """
        try:
            # Split input into lines
            lines = fasta_string.strip().split('\n')
            fastq_lines = []
            
            # Process each sequence
            current_header = ""
            current_sequence = []
            
            for line in lines:
                if line.startswith('>'):
                    # If we have a previous sequence, write it
                    if current_header and current_sequence:
                        sequence = ''.join(current_sequence)
                        quality = self.generate_quality_string(len(sequence))
                        
                        fastq_lines.append(f"@{current_header}")
                        fastq_lines.append(sequence)
                        fastq_lines.append("+")
                        fastq_lines.append(quality)
                    
                    # Start new sequence
                    current_header = line[1:]  # Remove '>'
                    current_sequence = []
                else:
                    current_sequence.append(line.strip())
            
            # Write last sequence
            if current_header and current_sequence:
                sequence = ''.join(current_sequence)
                quality = self.generate_quality_string(len(sequence))
                
                fastq_lines.append(f"@{current_header}")
                fastq_lines.append(sequence)
                fastq_lines.append("+")
                fastq_lines.append(quality)
            
            return '\n'.join(fastq_lines)
        except Exception as e:
            print(f"Error converting string: {str(e)}")
            return None

def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description='Convert FASTA to FASTQ format')
    parser.add_argument('input', help='Input FASTA file')
    parser.add_argument('output', help='Output FASTQ file')
    parser.add_argument('-q', '--quality', type=int, default=40,
                      help='Quality score to use (default: 40)')
    
    # Parse arguments
    args = parser.parse_args()
    
    # Create converter
    converter = FastaToFastqConverter(quality_score=args.quality)
    
    # Convert file
    if converter.convert_file(args.input, args.output):
        print(f"Successfully converted {args.input} to {args.output}")
        print(f"Used quality score: {args.quality}")
    else:
        print("Conversion failed")

if __name__ == "__main__":
    main()