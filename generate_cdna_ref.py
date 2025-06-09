from Bio import SeqIO
from Bio.Seq import Seq
import re


def find_cutsites(self, guide):
    """
    Find the cut sites for the given guide RNA sequence in the cDNA.
    Args:
        guide(str): Guide RNA sequence

    Returns:
        int: Position of the cut site in the cDNA seq

    """
    # Search for forward sequence
    positions = re.search(guide, cdna_seq)

    if positions:
        guide_start, guide_end = positions.span()
        cut_site = guide_end - 3
        print(cut_site)
        return cut_site

    # Search for reverse complement
    guide_seq = Seq(guide)
    rc_guide = str(guide_seq.reverse_complement())
    positions = re.search(rc_guide, cdna_seq)

    if positions:
        guide_start, guide_end = positions.span()
        cut_site = guide_start + 3
        print(cut_site)
        return cut_site

    # If neither forward nor reverse complement is found
    raise ValueError("Guide sequence not found in cDNA")


def delete_sequences(cut_site1, cut_site2):
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
    cdna_del = cdna_seq[:cut_site1] + cdna_seq[cut_site2:]
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
    inv_seq = str(Seq(cdna_seq[cutsite1:cutsite2]).reverse_complement())
    cdna_inv = cdna_seq[:cutsite1] + inv_seq + cdna_seq[cutsite2:]
    return cdna_inv


def write_to_fasta(cdna_id, cdna_seq, cdna_del, cdna_inv, output_file):
    """
    Write wildtype, deletion, and inversion sequences to a single FASTA file.

    Args:
        cdna_id (str): ID of the cDNA sequence
        cdna_seq (str): Original/wildtype cDNA sequence
        cdna_del (str): Deletion modified sequence
        cdna_inv (str): Inversion modified sequence
        output_file (str): Path to output FASTA file
    """
    sequences = [
        SeqRecord(Seq(cdna_seq), id=f"{cdna_id}_wildtype", description="Original sequence"),
        SeqRecord(Seq(cdna_del), id=f"{cdna_id}_deletion", description="Deletion variant"),
        SeqRecord(Seq(cdna_inv), id=f"{cdna_id}_inversion", description="Inversion variant")
    ]

    SeqIO.write(sequences, output_file, "fasta")
    print(f"Modified sequences written to {output_file}")


if __name__ == "__main__":
    try:
        # Read input FASTA
        with open(cdna_fasta, 'r') as handle:
            record = next(SeqIO.parse(handle, "fasta"))
            cdna_id = str(record.id)
            cdna_seq = str(record.seq)

        # Find cut sites
        guide1 = "GUIDE1SEQUENCE"
        guide2 = "GUIDE2SEQUENCE"
        cut_site1 = find_cutsites(guide1, cdna_seq)
        cut_site2 = find_cutsites(guide2, cdna_seq)

        # Generate modified sequences
        cdna_del = delete_sequences(cdna_seq, cut_site1, cut_site2)
        cdna_inv = invert_seq(cdna_seq, cut_site1, cut_site2)

        # Write to FASTA file
        write_to_fasta(cdna_id, cdna_seq, cdna_del, cdna_inv, "modified_cdnas.fasta")

    except Exception as e:
        print(f"Error: {e}")




mv pipe_utils-latest.prod-py3-none-any.whl pipe_utils-1.0.0-py3-none-any.whl
pip install pipe_utils-1.0.0-py3-none-any.whl


from pipe_utils import lims_utils
my_lims = lims_utils.LimsProcess()
my_entity = my_lims.get_entity_by_registry_id('RSQ73391')
my_entity.json()['bases']



import pysam

# Open the BAM file
samfile = pysam.AlignmentFile("/Users/harshini.muthukumar/PycharmProjects/harshini-coop/scripts/alignment_SNA219633.sam", "r")

# Get the reference names from the BAM file header
reference_names = samfile.references
name = reference_names[0]

# Print the reference names
print(name)

