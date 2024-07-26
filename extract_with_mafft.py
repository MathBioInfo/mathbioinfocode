import os
import subprocess
import sys
from Bio import AlignIO, SeqIO
from Bio.Align import MultipleSeqAlignment

def extract_region_with_reference(alignment_file, reference_file, output_file, temp_dir="temp", format="fasta"):
    # Ensure the temporary directory exists
    os.makedirs(temp_dir, exist_ok=True)
   
    # Read the reference sequence
    reference_seq = SeqIO.read(reference_file, format)
   
    # Read the alignment
    alignment = AlignIO.read(alignment_file, format)
   
    # Create a new alignment including the reference sequence
    combined_records = list(alignment) + [reference_seq]
    combined_alignment_file = os.path.join(temp_dir, "combined_alignment.fasta")
    SeqIO.write(combined_records, combined_alignment_file, format)
   
    # Run MAFFT to realign sequences including the reference
    aligned_output_file = os.path.join(temp_dir, "realigned_with_reference.fasta")
    mafft_cline = ["mafft", "--auto", combined_alignment_file]
    with open(aligned_output_file, "w") as outfile:
        subprocess.run(mafft_cline, stdout=outfile)
   
    # Read the realigned sequences
    realigned_alignment = AlignIO.read(aligned_output_file, format)
   
    # Find the index of the reference sequence in the alignment
    reference_index = None
    for i, record in enumerate(realigned_alignment):
        if record.id == reference_seq.id:
            reference_index = i
            break
   
    if reference_index is None:
        raise ValueError("Reference sequence not found in the alignment")
   
    # Extract the positions of the region from the reference sequence
    reference_record = realigned_alignment[reference_index]
    start, end = None, None
    for i, res in enumerate(reference_record.seq):
        if res != "-":
            if start is None:
                start = i
            end = i
   
    if start is None or end is None:
        raise ValueError("No non-gap positions found in the reference sequence")
   
    # Extract the regions from the aligned sequences
    spike_regions = []
    for record in realigned_alignment:
        spike_region = record[start:end + 1]
        spike_regions.append(spike_region)
   
    # Write the extracted regions to a new file
    AlignIO.write(MultipleSeqAlignment(spike_regions), output_file, format)

def print_help():
    help_text = """
Usage: python extract_with_mafft.py <alignment_file> <reference_file> <output_file>

Arguments:
    alignment_file    Path to the input alignment file (FASTA format).
    reference_file    Path to the reference spike sequence file (FASTA format).
    output_file       Path to the output file for the extracted spike regions (FASTA format).

Options:
    -h, --help        Show this help message and exit.

Example:
    python extract_spike_with_mafft.py aligned_sequences.fasta reference_spike.fasta spike_regions.fasta
    """
    print(help_text)

if __name__ == "__main__":
    if len(sys.argv) != 4:
        if len(sys.argv) == 2 and sys.argv[1] in ("-h", "--help"):
            print_help()
            sys.exit(0)
        else:
            print_help()
            sys.exit(1)
   
    alignment_file = sys.argv[1]
    reference_file = sys.argv[2]
    output_file = sys.argv[3]
   
    extract_region_with_reference(alignment_file, reference_file, output_file)
