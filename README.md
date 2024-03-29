# Mutation Analysis

#### Feature generation

* **Download and clean**: `python generators/download_pdb_and_gen_fasta.py`
* **Static amino acid descriptor** (finish later):
  | Descriptor                   | Source |
  | ---------------------------- | ------ |
  | AA_FORMAL_CHARGE             | URL    |
  | NORMALIZED_VAN_DER_WAALS_VOL | URL    |
* **Distance**:
  * Sample usage:
    * Remove the commented line at the end of the file.
    * Run: `python features/Distance.py`
* **Orientation**:
  * The angles are in [-pi, pi].
  * Dihedral angles (phi, psi, omega):
    * Use the two links ([url1](https://biopython.org/docs/dev/api/Bio.PDB.internal_coords.html#Bio.PDB.internal_coords.IC_Residue.pick_angle), [url2](https://biopython.org/docs/latest/api/Bio.PDB.vectors.html?highlight=calc_dihedral#Bio.PDB.vectors.calc_dihedral)) to compute phi, psi and omega angels.
      * (sC, nN, nCA, nC)   # phi i+1
      * (sN, sCA, sC, nN)   # psi
      * (sCA, sC, nN, nCA)  # omega i+1
  * Planer angle:
    * Use SCONES 2021 Figure 4: (sCA, sCB, nCB, nCA)
  * Sample usage:
    * Remove the commented line at the end of the file.
    * Run: `python features/Orientation.py`
* **MSA using hhblits**:
  * hh-suite3.0 setup:
    * Homepage: [URL](https://github.com/soedinglab/hh-suite).
    * User guide: [URL](https://github.com/soedinglab/hh-suite/wiki#generating-a-multiple-sequence-alignment-using-hhblits).
    * Download precompiled: [URL](https://mmseqs.com/hhsuite/) (hhsuite-linux-sse2.tar.gz).
    * An example:
      * Download a small database (i.e. scop40_01Mar17.tgz (381M)) from [URL](https://wwwuser.gwdg.de/~compbiol/data/hhsuite/databases/hhsuite_dbs/).
      * To generate MSA given a fasta file: `./3rd_party_items/hhsuite3.0/bin/hhblits -i path/to/filename.fasta -o path/to/filename.msa -d ./3rd_party_items/scop40_01Mar17/scop40`.
    * To search against a compiled large database:
      * Download Uniprot (18G): [URL](http://wwwuser.gwdg.de/~compbiol/data/hhsuite/benchmark/uniprot20_2015_06_bench.tgz).

    
* **PSSM generation**:
  * Based on wild and mutant sequences.
  * PSI-Blast setup:
    * Download the latest version ncbi-blast executables from [URL](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)
    * File name: ncbi-blast-2.12.0+-x64-linux.tar.gz
    * ncbi-blast executables contain blastp, blastn, psiblast and so on.
    * Quick start on BLAST: [URL](https://www.ncbi.nlm.nih.gov/books/NBK569856/)
    * After installation the path should be: `"3rd_party_items/ncbi-blast-2.12.0+/bin/psiblast"`
  * Database setup:
    * Download rp-seq-75 (or rp-seq-55) from this link: [URL](https://proteininformationresource.org/rps/)
    * Make BlastDB from fasta following this link: [URL](https://www.ncbi.nlm.nih.gov/books/NBK569856/)
    * The particular command should be looks like this: `3rd_party_items/ncbi-blast-2.12.0+/bin/makeblastdb -dbtype prot -in path/to/rp-seq-75.fasta -input_type fasta -out 3rd_party_items/rp_seq_75/rp_seq_75`
  * Since PSSM feature generation takes longer, therefore generated distributedly in GMU argo cluster:
    * `sbatch jobs/distributed_pssm_generator.sh`
  * The following items are optional (practice purposes only):
    * Download the database from [URL](https://ftp.ncbi.nlm.nih.gov/blast/db/*)
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
* **Secondary structure (SS), solvent accessibility (SA), relative accessible surface area (RASA) using DSSP module**
  * Bases on wild-type protein structure.
  * Install DSSP: `sudo apt install dssp`
    * from: [URL](https://ssbio.readthedocs.io/en/latest/instructions/dssp.html)
    * This installs 3.0.0 but current version is 4.0
    * 4.0 is not compatible with Biopython DSSP module. It can be installed by following:
      * Home: [URL](https://swift.cmbi.umcn.nl/gv/dssp/)
      * Github: [URL](https://github.com/PDB-REDO/dssp)
  * Implementation:
    * `features/SS_SA_RASA.py`

#### Analyzers

* To see the distribution of any mutation data file:
  * `python analyzers/data.py`
* To see the distribution of the features:
  * `python analyzers/features_statistics.py`
* To analyze the WaveEncoding which is same as positional encoding:
  * `python analyzers/encoding.py`
