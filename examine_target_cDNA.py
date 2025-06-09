import re
import pysam
import json
import subprocess
import requests
import pandas as pd

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from pathlib import Path
from pipe_utils import lims_utils


class CdnaEditGenerator():

    def __init__(self, cdna_id):
        """
        Initialize the cDNA object.
        Args:
            cdna: Input cDNA ENSEMBL ID
        """
        self.cdna_id = cdna_id.split('.')[0]
        base_url = "https://rest.ensembl.org"
        url = f"{base_url}/sequence/id/{cdna_id}?content-type=text/x-fasta;type=cdna"

        try:
            response = requests.get(url)
            raw_text = response.text.upper()
            self.cdna_seq = ''.join(raw_text.split('\n')[1:])
            response.raise_for_status()

        except requests.exceptions.JSONDecodeError:
            print("Failed to decode JSON response:")
            print(f"Response text: {response.text[:200]}...")  # Print first 200 chars

        except requests.exceptions.RequestException as e:
            print(f"Request failed: {e}")

    def get_guide_seq(self, guide_id):
        """
        Get the guide RNA sequence from Editas LIMS using lims_utils.
        Args:
            guide_id: Guide RNA id
        Returns:
            str: Guide RNA sequence
        """
        try:
            lims_obj = lims_utils.LimsProcess()
            my_entity = lims_obj.get_entity_by_registry_id(guide_id)
            guide_seq = my_entity.json()['bases']
            return guide_seq
        except:
            raise ValueError(f"Failed to retrieve guide sequence for {guide_id}")

    def find_cutsite(self, guide_seq):
        """
        Find the cut sites for the given guide RNA sequence in the cDNA.
        Args:
        guide_seq(str): Guide RNA sequence

        Returns:
            int: Position of the cut site in the cDNA seq

        """
        # Search for forward sequence
        positions = re.search(guide_seq, self.cdna_seq)

        if positions:
            guide_start, guide_end = positions.span()
            cut_site = guide_end - 3
            return cut_site

        # Search for reverse complement
        guide_seq = Seq(guide_seq)
        rc_guide = str(guide_seq.reverse_complement())
        positions = re.search(rc_guide, self.cdna_seq)

        if positions:
            guide_start, guide_end = positions.span()
            cut_site = guide_start + 3
            return cut_site

        # If neither forward nor reverse complement is found
        raise ValueError("Guide sequence not found in cDNA")

    def delete_sequence(self, cut_site1, cut_site2):
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

    def invert_sequence(self, cut_site1, cut_site2):
        """
        Invert the sequence between the two cut sites in the cDNA.
        Args:
        cut_site1 (int): Position of first cut site
        cut_site2 (int): Position of second cut site

        Returns:
            str: Modified cDNA sequence with the region between cut sites inverted

         """
        # Ensure cutsite1 is smaller than cutsite2
        cut_site1, cut_site2 = sorted([cut_site1, cut_site2])
        inv_seq = str(Seq(self.cdna_seq[cut_site1:cut_site2]).reverse_complement())
        cdna_inv = self.cdna_seq[:cut_site1] + inv_seq + self.cdna_seq[cut_site2:]
        return cdna_inv

    def generate_edited_seq_from_dual_guides(
            self, out_dir, guide1, guide2,
            generate_deletion=True, generate_inversion=True):
        """
        Generate deletion and/or inversion variants and write sequences to a FASTA file.

        Args:
            out_dir: the output directory
            cut_site1 (int): First cut site position
            cut_site2 (int): Second cut site position
            generate_deletion (bool): Whether to generate deletion variant (default: True)
            generate_inversion (bool): Whether to generate inversion variant (default: True)

        Raises:
            RuntimeError: If there's an error generating variants or writing the FASTA file
        """
        # specify output file
        outfile_name = f"{self.cdna_id}_{guide1}_{guide2}_editedSeq.fasta"
        output_file = out_dir / outfile_name

        if pd.isna(guide1) or pd.isna(guide2):
            generate_inversion = False
            generate_deletion = False
        else:
            # retrieve guide sequences
            guide1_seq = self.get_guide_seq(guide1)
            guide2_seq = self.get_guide_seq(guide2)
            # infer the cut sites
            cut_site1 = self.find_cutsite(guide1_seq)
            cut_site2 = self.find_cutsite(guide2_seq)

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
                cdna_del = self.delete_sequence(cut_site1, cut_site2)
                sequences.append(
                    SeqRecord(
                        Seq(cdna_del),
                        id=f"{self.cdna_id}_deletion",
                        description="Deletion variant"
                    )
                )

            # Generate inversion variant if requested
            if generate_inversion:
                cdna_inv = self.invert_sequence(cut_site1, cut_site2)
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
            return (output_file)

        except Exception as e:
            raise RuntimeError(f"Error generating variants or writing FASTA file: {str(e)}")


def ProcessBAM2FASTQ(s3_bam_path, output_dir):
    """

        Args:
            s3_bam_path:
            output_dir:

        Returns:
            reads_fq

        """
    filename = Path(s3_bam_path).name
    local_bam_path = str(Path(output_dir) / filename)
    out_prefix = '_'.join(filename.split('.')[:-1])

    try:
        # Add check=True to raise an exception if the command fails
        subprocess.run(['aws', 's3', 'cp', s3_bam_path, local_bam_path],
                       check=True,
                       stdout=subprocess.DEVNULL,
                       stderr=subprocess.DEVNULL
                       )

        # Create the full shell command
        # Fix the path formatting and add proper spacing
        reads_fq = str(Path(output_dir) / f'{out_prefix}_reads.fastq')
        command = ['samtools', 'fastq', '-0', '/dev/null', local_bam_path]

        # Execute samtools command
        with open(reads_fq, 'w') as fout:
            subprocess.run(command,
                           check=True,
                           stdout=fout,
                           stderr=subprocess.DEVNULL
                           )
        print(f"Converted {filename} to FASTQ")
        return (reads_fq)

    except Exception as e:
        print(f"Error processing {filename}: {e}")


def build_and_align_bowtie2(reference_fasta, reads_fq, out_dir):

    """
    Build bowtie2 index and perform paired-end alignment.

    Args:
        reference_fasta (str): Path to the reference FASTA file
        reads_fq (str): Path to single-end read file
        out_dir (str): Path to the output directory
    Raises:
        subprocess.CalledProcessError: If bowtie2 commands fail
    """
    try:
        ref_dir = str(Path(reference_fasta).parent)
        ref_prefix = (Path(reference_fasta).name).split('.fasta')[0]
        build_cmd = ['bowtie2-build', reference_fasta, f'{ref_dir}/{ref_prefix}']
        print(f"Running build command: {build_cmd}")
        subprocess.run(
            build_cmd,
            check=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL
        )
        print("Index building completed successfully")

    except Exception as e:
        raise RuntimeError(f"Unexpected error during bowtie2-build: {str(e)}")

    try:
        sample_prefix = (Path(reads_fq).name).split('.fastq')[0]
        out_sam = f"{out_dir}/{sample_prefix}.sam"
        align_cmd = ['bowtie2', '-x', f'{ref_dir}/{ref_prefix}',
                     '--local', '-U', reads_fq, '-S', out_sam]
        print(align_cmd)
        subprocess.run(
            align_cmd,
            check=True,
            stdout = subprocess.DEVNULL,
            stderr = subprocess.DEVNULL
        )
        print(f"Alignment completed. Results saved in {out_sam}")
        return (out_sam)

    except Exception as e:
        raise RuntimeError(f"Unexpected error during bowtie2 alignment: {str(e)}")


def process_cdna_alignment(samfile_path, min_mapq=20):
    """
    Process the SAM file to extract read counts. Keeps reads if:
    - Both reads map to same reference and at least one read passes MAPQ
    - Only one read is mapped and passes MAPQ
    Filters out chimeric pairs (reads mapping to different references)

    Args:
        samfile_path (str): Path to the SAM file
        min_mapq (int): Minimum mapping quality score (default: 20)

    Returns:
        dict: A dictionary with reference names as keys and read counts as values
    """
    try:
        with pysam.AlignmentFile(samfile_path, "r") as samfile:
            reference_names = samfile.references
            read_counts = {name: 0 for name in reference_names}
            aligned_reads = {}

            # load in each read mapping info
            for read in samfile.fetch():
                query_name = read.query_name

                # Skip unmapped reads
                if read.is_unmapped:
                    continue

                if query_name not in aligned_reads:
                    aligned_reads[query_name] = {
                        'aligned_ref': read.reference_name,
                        'read_mapq': read.mapq
                    }

        # count reads passing the mapping quality and filter chimeric reads
        for read in aligned_reads:
            aligned_ref = aligned_reads[read]['aligned_ref']
            read_mapq = aligned_reads[read]['read_mapq']
            # Check if read is aligned to one reference and passes MAPQ
            if aligned_ref is not None and read_mapq >= min_mapq:
                read_counts[aligned_ref] += 1

        return read_counts

    except Exception as e:
        raise RuntimeError(f"Error processing SAM file: {str(e)}")



def main():
    input_json = "input.json"
    with open(input_json, 'r') as file:
        config = json.load(file)

    req_id = config["req_id"]
    gene_name = config["gene_name"]
    transcript_id = config['transcript_id']
    s3_bam_key = config["s3_bam_key"]
    sample_df = pd.read_excel(config["sample_info_file"], engine='openpyxl')

    data_dir = Path(__file__).parent / f"data/{req_id}"
    data_dir.mkdir(parents=True, exist_ok=True)
    result_dir = Path(__file__).parent / f"result/{req_id}"
    result_dir.mkdir(parents=True, exist_ok=True)

    results_list = []
    for idx, row in sample_df.iterrows():
        sample_id = row["sample_id"]
        guide1 = row["guide_1"]
        guide2 = row["guide_2"]

        # generate cDNA reference sequence
        target_cdna = CdnaEditGenerator(transcript_id)
        reference_fasta = target_cdna.generate_edited_seq_from_dual_guides(data_dir, guide1, guide2)

        #  Download BAM file from S3
        bam_file = f"{sample_id}.{gene_name}.bam"
        s3_bam_path = f"{s3_bam_key}/{req_id}/{sample_id}/alignment/{bam_file}"
        reads_fq = ProcessBAM2FASTQ(s3_bam_path, data_dir)

        cdna_aligned_sam = build_and_align_bowtie2(
            reference_fasta,
            reads_fq,
            result_dir
        )

        # Process the alignment file to get count dictionary
        count_dict = process_cdna_alignment(cdna_aligned_sam, min_mapq=10)

        # from count dictionary to fraction dictionary
        total_count = sum(count_dict.values())
        if total_count > 0:
            fraction_dict = {key: round(value/total_count, 2) for key, value in count_dict.items()}
        else:
            fraction_dict = {key: 0 for key, value in count_dict.items()}

        # Create a result dictionary
        result = {'sample_id': sample_id, 'total_count': total_count}
        result.update({f'count_{ref.split("_")[1]}': count for ref, count in count_dict.items()})
        result.update({f'fraction_{ref.split("_")[1]}': frac for ref, frac in fraction_dict.items()})
        results_list.append(result)

    # Generate the result DataFrame
    columns = ['sample_id', 'total_count',
               'count_wildtype', 'count_deletion', 'count_inversion',
               'fraction_wildtype', 'fraction_deletion', 'fraction_inversion']
    result_df = pd.DataFrame(results_list, columns=columns)

    merged_df = pd.merge(
        sample_df,
        result_df,
        on='sample_id',
        how='left'  # keeps all samples from sample_df
    )

    # Save the output to file and upload to S3
    result_file = result_dir / f"{req_id}_{gene_name}_cDNA_edit_result.xlsx"
    merged_df.to_excel(result_file, index=False)

    # Upload the result file to S3
    s3_remote_path = f"{s3_bam_key}/{req_id}/custom_analysis/cDNA_edit_analysis"
    upload_cmd = ['aws', 's3', 'cp', result_dir, s3_remote_path, "--recursive"]
    subprocess.run(upload_cmd,
                   check=True,
                   stdout=subprocess.DEVNULL,
                   stderr=subprocess.DEVNULL
                   )


if __name__ == "__main__":
        main()

