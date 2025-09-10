# 🧬 Pipeline Scripts  

This folder contains the core scripts for building and running the CRISPR editing analysis workflow.  

- **`examine_target_cDNA.py`** → The **main pipeline script**, which integrates all components to quantify transcript-level CRISPR editing outcomes.  
- **`align_reads_to_cdna.py`** → Performs read alignment to the cDNA reference.  
- **`generate_cdna_ref.py`** → Generates reference files for cDNA-based alignment.  
- **`process_cdna.py`** → Generates cDNA reference sequences used in the alignment process.  
- **`generate_fastq_from_bam.py`** → Converts BAM files into FASTQ format for downstream processing.  

### 📌 Summary  
Together, these scripts form a **modular and reusable pipeline** for RNA-seq data processing and CRISPR editing quantification. The design allows each step to be run independently or as part of the full pipeline, making it scalable, flexible, and well-suited for real-world bioinformatics workflows.  
