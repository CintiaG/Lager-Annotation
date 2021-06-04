# MAKER Crash Course

Author: Cintia Gómez-Muñoz

Created: November 7, 2018

Updated: June 3, 2021

---

## Introduction
An annotation is a collection of [gene models](http://www.informatics.jax.org/glossary/gene_model). **MAKER** is a pipeline for *de novo* genome annotation of newly sequenced genomes (also to update annotations). It combines evidence alignment with Open Reading Frames (ORFs) prediction. The annotation process is divided in two parts:

* **Structural:** where the features are (start and end of exons, introns, UTRs, etc).
* **Functional:** what the features do (biological processes assignation).

**MAKER** automatically can perform the structural annotation by a main command, and the functional annotation through a series of other scripts.

### How does **MAKER** works?

| Stage | Description | Program |
| --- | --- | --- |
| 1. Identify and mask repeats | Repeats provide little information and give many useless hits| RepeatMasker and RepeatRunner|
| 2. *Ab initio* prediction| Gene prediction with mathematical models | SNAP, Augustus, GeneMark, and FGENESH |
| 3. Align Expressed Sequence Tags (ESTs) | Actively transcribed genes evidence | BLASTN |
| 4. Align proteins | Homology to known proteins | BLASTX |
| 5. Alignments polishing | Improves BLAST alignments | Exonerate |
| 6. Sythesis of information | Collection of gene models | **MAKER** |
| 7. Evidence quality values | For curation after the process | **MAKER** |

MAKER is an annotation tool, not a predictor.

> Prediction != Annotation

## Installation

Unfortunately, [MAKER](https://bioconda.github.io/recipes/maker/README.html) has a lot of prerequisites; one option would be to work in a server that already has the proper environment. Another option would be to use **Conda** to install **MAKER**.

In case you don't have **Conda**, go to the troubleshooting section at the end of this document.

```bash
conda install maker
#solving environment
```

It is likely that you will encounter many problems during the set up and running of **MAKER**, but for the brave, I appended a '**Troubleshooting installing MAKER**' at the end of this document.

## Dummy example

### Structural annotation

Now that you have successfully installed **MAKER**, feel free to read the manual.

```bash
maker -h
```

As stated above, **MAKER** uses a combination of tools to generate an annotation. Therefore, running **MAKER** in a command line with all options is not feasible. Instead, options are read from three control or '**.ctl**' files:

1. Executables and algorithms location (**maker_exe.ctl**).
2. Genome and evidence files, and other options (**maker_opts.ctl**).
3. BLAST and Exonerate Statistics Thresholds (**maker_bopts.ctl**).

To automatically create control files with default parameters, go to your desired working directory and type:

```bash
maker -CTL
```

You can open and explore each file to see which options you can change.

#### Executables and algorithms location

The **maker_exe.ctl** is the file that **MAKER** uses to know where the programs are. Make sure to have the routes of at least the following programs:
* Executables
  - makeblastdb
  - blastn
  - blastx
  - tblastx
  - RepeatMasker
  - exonerate


* Algorithms
  - snap
  - augustus

<!--Although I am not sure of this-->
Having extra programs such as tRNAscan-SE or snoscan cannot hurt, as long as they are turned off.

#### General options file

The **maker_opts.ctl** file includes most of the indications of how **MAKER** should behave and it has three main aspects: **genome** and **evidence** files, **RepeatMasker** and **gene prediction** options.

##### Genome and evidence files settings

To run **MAKER** you will need the direction of the three following files:
1. Genome file we want to annotate.
2. ESTs or RNA/Transcripts evidence file.
3. Proteins evidence file.

All the files should be in [FASTA format](http://www.bioinformatics.nl/tools/crab_fasta.html). Files from different organisms can be concatenated in a single file. Be aware that the longer your files are, the longer it would take alignments to be completed; however the more accurate your annotation can be. Nevertheless, the choose of these files is the tricky part.

How to choose the right organisms and evidence files for genome annotation is an open question. In brief, you should retrieve files from highly curated databases and from the same or a closely related organism. The first is to avoid propagating annotation errors. The second could be somewhat complicated if there is no close model organism; but this is not the case for yeasts of the *Saccharomyces* genus. Alternatively, you could use experimental evidence obtained from your study organism (e.g. RNA-Seq).

For the sake of a quick example, we are going to run **MAKER** without giving it too much thought to the evidence files (we would have to check the origin of each file). Our goal is to annotate the genome of a *Saccharomyces pastorianus* strain C-1082. *S. pastorianus* is a natural hybrid between *Saccharomyces cerevisiae* and *Saccharomyces eubayanus*, therefore we will use the concatenated files of transcripts and proteins (separately) provided in the parents genome's assemblies.

```bash
cat path_to_file/S288C_genome/orf_coding_all_R64-2-1_20150113.fasta path_to_file/S_eubayanus_genome/GCF_001298625.1_SEUB3.0_cds_from_genomic.fna > both_genomes_rna.fasta

cat path_to_file/S288C_genome/orf_trans_all_R64-2-1_20150113.fasta path_to_file/S_eubayanus_genome/GCF_001298625.1_SEUB3.0_protein.fasta > both_genomes_proteins.fasta
```

Be careful that each **FASTA** identifier is unique. Notice that **MAKER** does not allow the following characters in the EST files: 'R, Y, K, M, S, W, B, D, H and V', so edit your files with this command if necessary.

```bash
sed -i '/^>/!s/[R,Y,K,M,S,W,B,D,H,V]/N/g' both_genomes_rna.fasta
#use -i ony if you want to edit your files permanently
```

Edit the **maker_opts.ctl** file with the directions of the genome and evidence files. You can use a text editor such as `nano`. For example:

```bash
nano maker_opts.ctl
#do it manually
```

Enter manually (or copy and paste) the direction to the genome, ESTs and protein files after `genome=`, `est=`, and `protein=`, respectively. Be careful of not putting any space between the '=' sign and the directory file route. Then `Ctrl+O` and `Enter` to overwrite changes.

Or use `sed` to avoid modifying any other parts of the file accidentally.

```bash
sed -i 's/^genome=/genome=dir_to_genome_file/' maker_opts.ctl
sed -i 's/^est=/est=dir_to_est_file/' maker_opts.ctl
sed -i 's/^protein=/protein=dir_to_protein_file/' maker_opts.ctl
#can be piped to automatically generate multiple maker_opts.ctl files
#do not use -i option if you do not want to permanently modify your files
```

At the end, your **maker_opts.ctl** file should look something like this:

```bash
#-----Genome (these are always required)
genome=/data/scaffolds/820_v2_scaffolds.fasta #genome sequence (fasta file or fasta embeded in GFF3 file)
organism_type=eukaryotic #eukaryotic or prokaryotic. Default is eukaryotic


#-----EST Evidence (for best results provide a file for at least one)
est=both_genomes_rna.fasta #set of ESTs or assembled mRNA-seq in fasta format
altest= #EST/cDNA sequence file in fasta format from an alternate organism
est_gff= #aligned ESTs or mRNA-seq from an external GFF3 file
altest_gff= #aligned ESTs from a closly relate species in GFF3 format


#-----Protein Homology Evidence (for best results provide a file for at least one)
protein=both_genomes_proteins.fasta #protein sequence file in fasta format (i.e. from mutiple oransisms)
protein_gff=  #aligned protein homology evidence from an external GFF3 file
```

I omitted the 'Re-annotation' part, but just leave it unchanged in your original file.

Take not that predictions without protein evidence file does not work at all, so make sure to always provide one.

##### RepeatMasker settings

[RepeatMasker](http://www.repeatmasker.org/ "RepeatMasker web page") is a program that searches for repeats and replaces them with N's (hard masking) or lower case characters (soft masking). For this purpose, **RepeatMasker** needs a database in **FASTA** format. **RepeatMasker** already has three [databases](http://www.repeatmasker.org/libraries/ "RepeatMasker libraries") called Dfam, Dfam_consensus and RepeatPeps. Another option would be to use a [RepBase](https://www.girinst.org/ "RepBase web page") database. RepBase is a service of the Genetic Information Research Institute (GIRI) of representative repetitive elements in eukaryotic species, but soon it will be available only under a subscription. The **RepeatMasker** databases are stored in the Libraries's directory in '.embl' (detailed description) or '.lib' (actually, a **FASTA** file). The current Release of RepBase (dc20170127) is also included in the same directory.

```bash
cd /home/cintia/programs/anaconda3/envs/maker_py/share/RepeatMasker/Libraries
#directory direction in my case
```

You can specify the direction of one the **RepeatMasker** databases in the `rmlib=` option or you can choose a RepBase model organism in `model_org=`. We chose to use the RepBase specific *Saccharomyces* library only as follows:

```bash
#-----Repeat Masking (leave values blank to skip repeat masking)
model_org=saccharomyces
repeat_protein=/data/opt/bin/maker/data/te_proteins.fasta
```
<!---I do not know from where Marcel downloaded the te proteins file--->

You can also choose between soft or hard masking `softmask=1` (1 on or 0 off). For the rest, leave the default parameters.

##### Gene Prediction settings

This step basically consist of detecting Open Reading Frames (ORFs) using mainly a probabilistic model. For the former you can:

* Use an already existing probabilistic model from a well studied organism (*e.g. S. cerevisiae*).
* Create an organism specific probabilistic model from the experimental evidence alignments or from a protein database of closely related organisms, by training the gene predictor in a iterative process (the more information, the better).

**MAKER** supports several gene predictors but we will only discuss **SNAP** and **Augustus**. Both programs use [Hidden Markov Models](https://www.nature.com/articles/nbt1004-1315) (HMMs). HMMs are probabilistic models representing the *emission and transition probabilities* of a gene structure, in this case (e.g. start codon, exon, intron, etc.).

* [SNAP](https://github.com/KorfLab/SNAP/blob/master/README.md 'SNAP GitHub') (Semi-HMM-based Nucleic Acid Parser) uses a six state HMM ([see original article](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-5-59)) to generate a gene model in a '.hmm' file that can be input to MAKER.
* [Augustus](https://github.com/Gaius-Augustus/Augustus/blob/master/README.md 'Augustus GitHub') ([see original article](https://academic.oup.com/bioinformatics/article/19/suppl_2/ii215/180603?searchresult=1)) uses a HMM with a more complete structure, and it is pre-loaded with trained annotation files for several species including *S. cerevisiae*.

<!---como encuentran el archivo hmm de Sc?--->

```bash
#-----Gene prediction
snaphmm=saccharomyces_cerevisiae_S288C.length.hmm #SNAP HMM file
augustus_species=saccharomyces_cerevisiae_S288C
```

<!---For this test, I did not use the SNAP file, but I included it in the manual for reference puroposes--->

In case you need to generate an organism specific gene model from the ESTs or protein files, you can turn on these options:

```bash
est2genome=1
protein2genome=1
#default is 0 (off)
```

And then you train each gene predictor accordingly.

#### BLAST alignment options

The **maker_bopts.ctl** file is where **MAKER** reads the BLAST options. You can choose the type of BLAST among **ncbi+**, **ncbi**, or **wublast**. The three following sections are the percent coverage, percent identity, e-value, bit score and depth thresholds and cutoffs for the blastn, blastx and tblastx programs. Last, there are options for Transposable Element Masking and Exonerate. You can leave all of this file unchanged, unless you have special needs.

#### Running MAKER

After we have set up our control files, we can run **MAKER** with this simple command.
```bash
maker
```

This might take a while depending on your computer capacities. **MAKER** is very verbose, so if you wish to only see handy messages run **MAKER** with `maker -q` (quiet) or `maker -qq` (super quiet).

MAKER creates an output directory using the base name of your genome file and the '**.maker.output**' suffix. There it stores the resulting:

* Raw files (**theVoid**)
* Databases (**mpi_blastdb**)
* Log files (**used control files**)
* Results (**_datastore**)
  - Actual results are in a complex data structure storage. The '**_master_datastore_index.log**' is a guide to see if each contig was processed correctly and to know the directory of the following contig results:
    - The **proteins and transcripts** in **FASTA** files
    - A '**GFF3**' file that contains all annotations and alignments

If no probabilistic models were specified (empty gene prediction variables in maker_opts.ctl file), then no **FASTA** files are going to be created.

At this point, all the information is scattered if different directories. To obtain a single file with all the contigs information you can do the following:

```bash
fasta_merge -d genome_base_name_master_datastore_index.log
gff3_merge -d genome_base_name_master_datastore_index.log
#genome base name should be the name of your original file
```

To get only the gene models and no the *ab initio* predictions and evidence alignments, you use the `-g` flag.

```bash
gff3_merge -d "$OUTPUT_DIR"/"$OUTPUT_PREFIX"_master_datastore_index.log -g -o "$OUTPUT_FILE"
```

You can inspect the 'GFF3' file but it contains too much information to be understood by the human eye; therefore, it would be better to use a visualization tool. I used the Integrative Genomics Viewer (IGV), but you can also use JBrowser or Artemis.

<!---check JBrowser and Artemis--->

### Functional annotation

This step can be performed using inferred homology to a protein database using blastp. For that matter, I used the UniProt [*S. cerevisiae* database](https://www.uniprot.org/proteomes/UP000002311), and I created a new directory for the alignments results.

```bash
mkdir blast_homology
```

First, you need to construct a **BLAST** database if there is not one available with the following command:

```bash
makeblastdb -dbtype prot -in UP000002311_559292.fasta
```

Second, you can align the protein database with your **MAKER** proteins file:

```bash
blastp -db UP000002311_559292.fasta -query 820_v2_scaffolds.maker.output/820_v2_scaffolds.all.maker.proteins.fasta -out blast_homology/protein.blastp -evalue 1e-6 -outfmt 6 -max_target_seqs 1 -num_threads 2
```

Third, you can parse the information of the blastp output file to the GFF3 and the **FASTA** files with the following command.

```bash
#For GFF3 file
maker_functional_gff UP000002311_559292.fasta blast_homology/protein.blastp original.gff > new_file.gff
#For proteins fasta file
maker_functional_fasta UP000002311_559292.fasta blast_homology/protein.blastp original_protein.fasta > new_protein_file.fasta
#For transcripts fasta file
maker_functional_fasta UP000002311_559292.fasta blast_homology/protein.blastp original_transcripts.gff > new_transcripts_file.fasta
```

<!---Marcel made a script to perform the above automatically, but I don't understand why he uses cd .. and wait--->

You can check your results using the `less` command.

Notice that this only works with [UniProt fasta headers](https://www.uniprot.org/help/fasta-headers), so make sure that your database follows this structure:

```none
>db|UniqueIdentifier|EntryName ProteinName OS=OrganismName OX=OrganismIdentifier [GN=GeneName ]PE=ProteinExistence SV=SequenceVersion
```

Now, you have finished your first attempt of annotating a genome!

### Additional useful scripts

You may have noticed that the resulting ORFs IDs are long and difficult to understand, so perhaps you would like to have a shorter and systematic name. For that, **MAKER** has additional built in scripts that allow you to perform the former.

Firts, you need to create a map of the current and new IDs as follows:

```bash
maker_map_ids --prefix Sp820_ --justify 5 annotated_820.all.gff > annotated_820.map
```

The following scripts overwrite the original files, so if you want to keep them, make a copy first with `cp`. Afterwards, you can apply the changes to all files:

```bash
map_gff_ids annotated_820.map annotated_820.renamed.gff
map_fasta_ids annotated_820.map annotated_820.maker.proteins.renamed.fasta
map_fasta_ids annotated_820.map annotated_820.maker.transcripts.renamed.fasta
```

Take into consideration that this function does not work well if the reference genome sequences ID's contains the pipe symbol `|`, so it would be helpful to remove all of them prior this step.

---

## Troubleshooting installing Maker

This is a summary of the problems I found with my Linux machine with Ubuntu 18.04.1 LTS during **MAKER** installation.

### Conda

[Conda](https://conda.io/docs/) is a package, dependency and environment manager for several languages. In other words, it makes a little bit easier to check if you have all the required tools or to manage programs that use different versions of Python. To use it, you should already have **Anaconda**.

There are two versions:
* Anaconda (all packages).
* Miniconda (no packages).

You can download the appropriate version for your OS [here](https://www.anaconda.com/download/ "Anaconda download"). For Linux users you can follow the instructions of this [video](https://www.youtube.com/watch?v=3i0YNcgmGZ8). For the impatient, here is the terminal code.

```bash
#execute installation
bash Anaconda3-5.3.0-Linux-x86_64.sh
```
Review the license agreement and type 'yes' if you agree. Accept default directory to store or specify your own. I chose a special directory where I store my programs. Open a new terminal and check if it is already installed.

```bash
conda --version
#outputs conda 4.5.11
```

To install **MAKER**, additionally you need to set up this [channels](https://bioconda.github.io/index.html#set-up-channels) in that order.

```bash
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```

In case any requirements are not met, a warning will appear. In my case, I had Python 3 version instead of the Python 2.7 that is required by **MAKER**.

To view all available **Python** versions, type:

```bash
conda search python
#outputs all available versions
```

To create a **MAKER** suitable environment you can type:

```bash
conda create -n maker_py python=2.7
#maker_py in the name I assigned to the new environment
```

To see all available **Conda** environments:

```bash
conda info --envs

# conda environments:
#
#base                  *  /home/cintia/programs/anaconda3
#maker_py                 /home/cintia/programs/anaconda3/envs/maker_py
```

All available environments outputs; '**base**' is your default environment. To set the new environment:

```bash
conda activate maker_py
```

You can type `conda deactivate` to go back to **base** environment.

Re-run the `conda install maker` command and accept the suggested required packages. When it is done, type `which maker` to see if it is already installed. The direction where it is stored should appear.

### RepeatMasker troubleshooting

Another problem I encountered is that **MAKER** could not find **RepeatMasker** because it was not included in my PATH. This can be easily solved by modifying this variable.

```bash
cd /home/cintia/
nano .bashrc
```

Append to the end of the file the following:

```bash
export PATH=$PATH:/home/cintia/programs/anaconda3/envs/maker_py/bin
#or whatever direction the program is located
```

And then type the following to update the changes:

```bash
source .bashrc
```

Further problems with **RepeatMasker**, like 'Could not determine if RepBase is installed', were due to misconfiguration. To fix this, you can go to the **RepeatMasker** directory, run the configuration file and follow the instructions.

```bash
cd /home/cintia/programs/anaconda3/envs/maker_py/share/RepeatMasker
perl ./configure
```

Check and provide the right information. For the search engine, I chose 'RMBlast - NCBI Blast with RepeatMasker extensions'.

Finally, I had to install one additional algorithm called 'Text::Soundex' with the following command:

```bash
sudo cpan Text::Soundex
```

And with all this, I was able to run **MAKER** in my Linux machine.

## Reading material

* [MAKER Tutorial for WGS Assembly and Annotation Winter School 2018](http://weatherby.genetics.utah.edu/MAKER/wiki/index.php/MAKER_Tutorial_for_WGS_Assembly_and_Annotation_Winter_School_2018 "MAKER Tutotial").
