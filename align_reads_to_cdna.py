import os.path
import pysam
from Bio import SeqIO
from Bio.Seq import Seq
import re
import generate_fastq_from_bam
from Bio.SeqRecord import SeqRecord
import subprocess
import requests



class align_reads_to_cdna():

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

    def find_cutsites(self,guide):
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


    def delete_sequences(self,cut_site1, cut_site2):
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


    def invert_seq( self,cutsite1, cutsite2):
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

    def generate_and_write_variants(self, output_file, cut_site1, cut_site2, generate_deletion=True,
                                    generate_inversion=True):
        """
        Generate deletion and/or inversion variants and write sequences to a FASTA file.

        Args:
            output_file (str): Path to output FASTA file
            cut_site1 (int): First cut site position
            cut_site2 (int): Second cut site position
            generate_deletion (bool): Whether to generate deletion variant (default: True)
            generate_inversion (bool): Whether to generate inversion variant (default: True)

        Raises:
            RuntimeError: If there's an error generating variants or writing the FASTA file
        """




        try:
            # Initialize sequences list with wildtype
            sequences = [
                SeqRecord(
                    Seq(self.cdna_seq),
                    id=f"{self.cdna_id}_wildtype",
                    description="Original sequence"
                )
            ]

            # Generate deletion variant if requested
            if generate_deletion:
                cdna_del = self.delete_sequences(cut_site1, cut_site2)
                sequences.append(
                    SeqRecord(
                        Seq(cdna_del),
                        id=f"{self.cdna_id}_deletion",
                        description="Deletion variant"
                    )
                )

            # Generate inversion variant if requested
            if generate_inversion:
                cdna_inv = self.invert_seq(cut_site1, cut_site2)
                sequences.append(
                    SeqRecord(
                        Seq(cdna_inv),
                        id=f"{self.cdna_id}_inversion",
                        description="Inversion variant"
                    )
                )

            # Write sequences to FASTA file
            SeqIO.write(sequences, output_file, "fasta")
            print(f"Modified sequences written to {output_file}")

        except Exception as e:
            raise RuntimeError(f"Error generating variants or writing FASTA file: {str(e)}")

    def build_and_align_bowtie2(self, reference_fasta, r1_path, r2_path, output_sam, index_prefix):

        """
        Build bowtie2 index and perform paired-end alignment.

        Args:
            reference_fasta (str): Path to the reference FASTA file
            r1 (str): Path to first paired-end read file
            r2 (str): Path to second paired-end read file
            output_sam (str): Path for output SAM file
            index_prefix (str): Prefix for bowtie2 index files

        Raises:
            subprocess.CalledProcessError: If bowtie2 commands fail
        """
        # Conda environment name
        conda_env = "bioinfo"
        #created a directory for alignment to direct the index files and SAM files here
        alignment_dir = "alignment"
        os.makedirs(alignment_dir, exist_ok=True)

        try:
            # First run the build command separately
            build_cmd = f"conda run -n {conda_env} bowtie2-build {reference_fasta} ./{alignment_dir}/{index_prefix}"
            print(f"Running build command: {build_cmd}")
            subprocess.run(
                build_cmd,
                check=True,
                shell=True,
                executable='/bin/bash'
            )
            print("Index building completed successfully")


            align_cmd = f"conda run -n {conda_env} bowtie2 -x ./{alignment_dir}/{index_prefix} -1 {r1_path} -2 {r2_path} -S ./{alignment_dir}/{output_sam}"



            subprocess.run(
                align_cmd,
                check=True,
                shell=True,
                executable='/bin/bash',
                capture_output=True,
                text=True
            )

            print(f"Alignment completed. Results saved in {output_sam}")




        except Exception as e:
            print(f"Unexpected error: {str(e)}")
            raise



def main():

        #Using ENSEMBL API to fetch the cdna fasta sequences
        base_url = "https://rest.ensembl.org"
        transcript_id = "ENST00000558518"
        url = f"{base_url}/sequence/id/{transcript_id}?content-type=text/x-fasta;type=cdna"

        response = requests.get(url)

        with open("cdna_fasta.fasta", 'w') as handle:
            handle.write(response.text)

        #creating a instance of a class
        cdna_aligner = align_reads_to_cdna("cdna_fasta.fasta")

        #calling the main function from generate_fastq_from_bam script to generate r1 and r2 fastq from BAM
        generate_fastq_from_bam.main()

        #again processing the analysis sheet
        #doubtful in this step (check)
        analysis_sheet = generate_fastq_from_bam.download_analysis_sheet()
        samples = generate_fastq_from_bam.ProcessAnalysisSheet(analysis_sheet)


       #iterating over the samples dictionery to get the values
        for sample in samples:
            sample_name = sample["sample_name"]
            guide1 = sample["guide1"]
            guide2 = sample["guide2"]
            r1 = f"ldlr_{sample_name}_r1.fastq"
            r2 = f"ldlr_{sample_name}_r2.fastq"
            fastq_dir = "output_fastq"
            r1_path = os.path.join(fastq_dir, r1)
            r2_path = os.path.join(fastq_dir, r2)

            #wrote a if else block so that we can process the MOCK samples which dont have guides
            if guide1.lower() == "none" and guide2.lower() == "none":
                print(f"Sample {sample_name} has no guides - aligning directly to original cDNA")

                # aligning the reads to the orginal cdna fasta sequence
                # not sure if this is needed ( do check)

                # Align directly to the original cDNA sequence
                cdna_aligner.build_and_align_bowtie2(
                    "cdna_fasta.fasta",  # Use original cDNA file
                    r1_path,
                    r2_path,
                    f"alignment_{sample_name}.sam",
                    "cdna_index"
                )
            else:

               #goes into else block if the guide sequences are not none
                cut_site1 = cdna_aligner.find_cutsites(guide1)
                cut_site2 = cdna_aligner.find_cutsites(guide2)

                cdna_aligner.generate_and_write_variants(
                    output_file="modified_cdnas.fasta",
                    cut_site1=cut_site1,
                    cut_site2=cut_site2
                )





                # Build and align with bowtie2 using modified cDNAs
                cdna_aligner.build_and_align_bowtie2(
                    "modified_cdnas.fasta",
                    r1_path,
                    r2_path,
                    f"alignment_{sample_name}.sam",
                    "modified_cdna_index"
                )

                # Use the full path for the SAM file
                sam_path = os.path.join("alignment", f"alignment_{sample_name}.sam")
                samfile = pysam.AlignmentFile(sam_path)


                # Get the reference names from the BAM file header
                reference_names = samfile.references

                # Initialize a dictionary to keep track of reads aligned to each reference
                read_counts = {name: 0 for name in reference_names}
                total_reads = 0

                # Loop through each read in the BAM file
                for read in samfile.fetch():
                    total_reads += 1
                    if not read.is_unmapped:  # Make sure the read is mapped
                        reference_name = reference_names[read.reference_id]  # Get the reference name for the read
                        read_counts[reference_name] += 1

                # Calculate and print the percentage of reads aligned to each reference
                for reference_name, count in read_counts.items():
                    percentage = (count / total_reads) * 100
                    print(f"Percentage of reads aligned to {reference_name}: {percentage:.2f}%")

                # Close the BAM file
                samfile.close()


if __name__ == "__main__":
        main()

