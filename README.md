# MetaCompare

MetaCompare is a computational pipeline for prioritizing resistome risk by estimating the potential for ARGs to be disseminated into human pathogens from a given environmental sample based on metagenomic sequencing data.

## Getting Started

These instructions will get you a copy of the project up and running on your local machine or server.

### System Requirements (*tested on linux Ubuntu 14.04*)

* git installed
* Python3 with pandas package installed
* Blast 2.2.8 or higher version installed

### Installing

**Step 1:** Change the current working directory to the location where you want the cloned `MetaCompare` directory to be made.
**Step 2:** Clone the repository using git command
```
~$ git clone https://github.com/minoh0201/MetaCompare
```

**Step 3:** Make directory `BlastDB` and change woring directory to it.

```
~$ mkdir BlastDB
~$ cd BlastDB
```

**Step 4:** Download the compressed Blast Database file from the web server (25 GB) and uncompress it.

```
~/BlastDB$ wget http://bench.cs.vt.edu/ftp/data/metacomp/BlastDB.tar.gz
~/BlastDB$ tar -zxvf BlastDB.tar.gz
```

**Step 5:** Get back to working directory `MetaCompare` and run `metacmp.py`

```
~/BlastDB$ cd ..
~$ ./metacmp.py
```

## Running MetaCompare

### Preparing input files

MetaCompare requires two input FASTA files, one for the assembled contigs and the other for predicted gene list derived from the assembled contigs using prodigal. 

If you have raw reads you can submit them to MetaStorm (http://bench.cs.vt.edu/MetaStorm/) and run assembly pipeline to get assembled contigs and predicted gene list. 

These files can be downloaded from MetaStorm by clicking `Scaffolds` button and `Genes` button in the assembled sample page.

### Running

Suppose you have the assembled contigs file, `S1.fa`, and predicted gene list `S1_genes.fa`. 

The following command runs MetaCompare with 128 threads.

```
~$ ./metacmp.py -c S1.fa -g S1_genes.fa -t 128
```
The output should be look like as follows:
```
Running blastn on ACLAME
Running blastx on CARD
Running blastn on PATRIC
Reading files...
Computing resistome risk score..
Resistome risk score: 38.64014990951873
```

You can see detailed description for command line options by using `-h` option.
```
~$ ./metacmp.py -h
```
## Citation

TBA
