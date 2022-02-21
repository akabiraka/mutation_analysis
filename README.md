# Mutation Analysis

* **Download and clean**: `python generators/download_pdb_and_gen_fasta.py`
* **Evolutionary features** (dimension=16):

  * Based on wild and mutant sequences.
  * PSI-Blast setup:

    * Download the latest version ncbi-blast executables from https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
    * File name: ncbi-blast-2.12.0+-x64-linux.tar.gz
    * ncbi-blast executables contain blastp, blastn, psiblast and so on.
    * Quick start on BLAST: https://www.ncbi.nlm.nih.gov/books/NBK569856/
    * After installation the path should be: `"3rd_party_items/ncbi-blast-2.12.0+/bin/psiblast"`
  * Database setup:

    * Download rp-seq-75 (or rp-seq-55) from this link: https://proteininformationresource.org/rps/
    * Make BlastDB from fasta following this link: https://www.ncbi.nlm.nih.gov/books/NBK569856/
    * The particular command should be looks like this: `3rd_party_items/ncbi-blast-2.12.0+/bin/makeblastdb -dbtype prot -in path/to/rp-seq-75.fasta -input_type fasta -out 3rd_party_items/rp_seq_75/rp_seq_75`
  * Since PSSM feature generation takes longer, therefore generated distributedly in GMU argo cluster:

    * `sbatch jobs/distributed_pssm_generator.sh`
  * To generate the **evolutionary features**:

    * `python data_generators/feature_pssm.py`
  * The following items are optional (practice purposes only):

    * Download the database from https://ftp.ncbi.nlm.nih.gov/blast/db/*

      * File name: swissprot.tar.gz (because it is small to setup and test)
    * Biopython provides `NcbipsiblastCommandline` to use psi-blast software.

      * `from Bio.Blast.Applications import NcbipsiblastCommandline`
    * Sample example of BLAST usage from command line:

      * To download swissprot: `3rd_party_items/ncbi-blast-2.12.0+/bin/update_blastdb.pl --decompress swissprot`
      * To download all "nr" databases: "nr" is the non-redundent version of the sequence database.
        * `3rd_party_items/ncbi-blast-2.12.0+/bin/update_blastdb.pl --decompress nr [*]`
        * `ncbi-blast-2.12.0+/bin/update_blastdb.pl --decompress nr.02`
      * To create blast-db from fasta sequences:
        * `3rd_party_items/ncbi-blast-2.12.0+/bin/makeblastdb -dbtype prot -in path/4eiuA.fasta -input_type fasta -out path_to_save/db_name`
        * `3rd_party_items/ncbi-blast-2.12.0+/bin/makeblastdb -dbtype prot -in 3rd_party_items/rp-seqs-15.fasta -input_type fasta -out 3rd_party_items/rp_req_15/rp_req_15`

#### Analyzers
* To see the distribution of any mutation data file:
    * `python analyzers/data.py`