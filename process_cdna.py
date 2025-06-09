import re
import random
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import argparse


class cdna_modification:
    """
        A class to handle cDNA modifications using guide RNAs.

        Attributes:
            cdna (str): The input cDNA sequence
            guide1 (str): First guide RNA sequence
            guide2 (str): Second guide RNA sequence
            modification_type (str): Type of modification to perform (delete/invert/none)
        """
    def __init__(self, cdna_fasta):
        """
        Initialize the cdna_modification object.
        Args:
            cdna: Input cDNA sequence
        """
        with open(cdna_fasta, 'r') as handle:
            record = next(SeqIO.parse(handle, "fasta"))
            self.cdna_id = str(record.id)
            self.cdna_seq = str(record.seq)


    def find_cutsites(self, guide):
        """
        Find the cut sites for the given guide RNA sequence in the cDNA.
        Args:
            guide(str): Guide RNA sequence

        Returns:
            int: Position of the cut site in the cDNA seq

        """
        # Search for forward sequence
        positions = re.search(guide, self.cdna_seq)

        if positions:
            guide_start, guide_end = positions.span()
            cut_site = guide_end - 3
            print(cut_site)
            return cut_site

        # Search for reverse complement
        guide_seq = Seq(guide)
        rc_guide = str(guide_seq.reverse_complement())
        positions = re.search(rc_guide, self.cdna_seq)

        if positions:
            guide_start, guide_end = positions.span()
            cut_site = guide_start + 3
            print(cut_site)
            return cut_site

        # If neither forward nor reverse complement is found
        raise ValueError("Guide sequence not found in cDNA")



    def delete_sequences(self, cut_site1, cut_site2):
        """
        Delete the sequences between the two cut sites in the cDNA.
        Args:
            cut_site1(int): Position of first cut site
            cut_site2(int): Position of second cut site

        Returns:
            str: Modified cDNA sequence with the region between cut sites deleted

        """
        # Ensure cut_site1 is smaller than cut_site2
        cut_site1, cut_site2 = sorted([cut_site1, cut_site2])
        cdna_del = self.cdna_seq[:cut_site1] + self.cdna_seq[cut_site2:]
        return cdna_del

    def invert_seq(self, cutsite1, cutsite2):
        """
        Invert the sequence between the two cut sites in the cDNA.
        Args:
            cutsite1 (int): Position of first cut site
            cutsite2 (int): Position of second cut site

        Returns:
            str: Modified cDNA sequence with the region between cut sites inverted

        """
        # Ensure cutsite1 is smaller than cutsite2
        cutsite1, cutsite2 = sorted([cutsite1, cutsite2])
        inv_seq = str(Seq(self.cdna_seq[cutsite1:cutsite2]).reverse_complement())
        cdna_inv = self.cdna_seq[:cutsite1] + inv_seq + self.cdna_seq[cutsite2:]
        return cdna_inv

    def modify_by_dual_guides(self, guide1, guide2):
        """
        Modify the cDNA sequence by deleting or inverting the region between the two guides.
        Returns:
            str: Modified cDNA sequence

        """
        self.guide1 = guide1
        self.guide2 = guide2

        cut_site1 = self.find_cutsites(self.guide1)
        cut_site2 = self.find_cutsites(self.guide2)

        # generate sequence with large deletion
        self.cdna_del = self.delete_sequences(cut_site1, cut_site2)
        # generate sequence with inversion
        self.cdna_inv = self.invert_seq(cut_site1, cut_site2)

    def generate_fastq(self, output1, output2, modification,num_reads,read_length=250):
        """
        Generate simulated FASTQ files for the modified cDNA sequence.
        Args:
            output1(str): Path for R1 fastq file
            output2(str): Path for R2 fastq file
            read_length(int): Length of each read
            num_reads(int): Number of reads to generate
            modification(str):modification type

        Returns:
            R1,R1 Fatsq files

        """

        self.modification = modification

        if self.modification == 'delete':
            cdna_seq = self.cdna_del
        elif self.modification== 'invert':
            cdna_seq = self.cdna_inv
        else:  # 'none' or default case
            cdna_seq = self.cdna_seq

        random.seed(10)
        insert_size_list = [random.randint(300, 400) for _ in range(num_reads)]

        output1, output2 = (f"target_{modification}_r1.fastq", f"target_{modification}_r2.fastq") if modification in [
            'delete', 'invert'] else (output1, output2)

        with open(output1, "w") as h1, \
                open(output2, "w") as h2:
            for i in range(num_reads):
                insert_size = insert_size_list[i]

                # R2 is sequenced from the polyA site (i.e., the tail of cDNA)
                r2_start = len(cdna_seq) - read_length
                r2 = cdna_seq[r2_start:len(cdna_seq)]
                r2_rc = Seq(r2).reverse_complement()

                # R1 is sequenced from a given distance from the R2
                r1_start = len(cdna_seq) - insert_size
                r1_end = r1_start + read_length
                r1 = cdna_seq[r1_start:r1_end]

                r1_scores = [random.randint(20, 40) for _ in range(read_length)]
                r1_record = SeqRecord(Seq(r1), id=f"read{i}_R1", description=f"modification={self.modification}",
                                      letter_annotations={"phred_quality": r1_scores})

                r2_scores = [random.randint(20, 40) for _ in range(read_length)]
                r2_record = SeqRecord(r2_rc, id=f"read{i}_R2",description=f"modification={self.modification}",
                                      letter_annotations={"phred_quality": r2_scores})

                SeqIO.write(r1_record, h1, "fastq")
                SeqIO.write(r2_record, h2, "fastq")


def parse_arguments():
    """
        Parse command line arguments for cDNA modification and FASTQ generation.

        Returns:
            argparse.Namespace: An object containing all the parsed arguments with the following attributes:
                - input_fasta: Path to input FASTA file
                - guide1: First guide RNA sequence
                - guide2: Second guide RNA sequence
                - modification: Type of modification (delete/invert/none)
                - output1: Path for first read FASTQ file
                - output2: Path for second read FASTQ file
        """

    parser = argparse.ArgumentParser(description='Process cDNA modifications with guide RNAs')

    parser.add_argument('-i', '--input_fasta',
                        help='Input FASTA file path',
                        required=True)

    parser.add_argument('-g1', '--guide1',
                        default="CGAGGTCTCAGGAAGGGTTC",
                        help='First guide RNA sequence')

    parser.add_argument('-g2', '--guide2',
                        default="AGGCAATAACCCCCTACACA",
                        help='Second guide RNA sequence')

    parser.add_argument('-n', '--num_reads',
                        type=int,
                        default=100,
                        help='Number of reads to generate (default: 100)')



    parser.add_argument('-m', '--modification',
                        choices=['delete', 'invert', 'none'],
                        default='none',
                        help='Type of modification to perform')

    parser.add_argument('-o1', '--output1',
                        default="target_wildtype_r1.fastq",
                        help='Output path for first read FASTQ file')

    parser.add_argument('-o2', '--output2',
                        default="target_wildtype_r2.fastq",
                        help='Output path for second read FASTQ file')

    return parser.parse_args()


def main():
    """
    Main function to process cDNA modifications and generate FASTQ files.
    """
    args = parse_arguments()


    # Create instance and process
    modified_cdna = cdna_modification(args.input_fasta)
    modified_cdna.modify_by_dual_guides(args.guide1, args.guide2)

    # write the three sequences into fasta file
    # Option 2: Using direct f-string with ternary operator
    with open(
            f"{'deletion' if args.modification == 'delete' else 'inversion' if args.modification == 'invert' else 'wildtype'}_cdna.fasta",
            "w") as handle:

        if args.modification == 'delete':
            SeqIO.write(SeqRecord(Seq(modified_cdna.cdna_del), id=f"{modified_cdna.cdna_id}_del",
                                   description=f"modification={args.modification}"), handle, "fasta")
        elif args.modification == 'invert':
            SeqIO.write(SeqRecord(Seq(modified_cdna.cdna_inv), id=f"{modified_cdna.cdna_id}_inv",
                                   description=f"modification={args.modification}"), handle, "fasta")
        else:
            SeqIO.write(SeqRecord(Seq(modified_cdna.cdna_seq), id=f"{modified_cdna.cdna_id}_ori",
                                   description=f"modification={args.modification}"), handle, "fasta")

       # Generate FASTQ files
    modified_cdna.generate_fastq(args.output1, args.output2,args.modification,args.num_reads)


if __name__ == "__main__":
    main()








