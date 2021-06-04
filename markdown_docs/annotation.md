# Annotating Bacterial Genomes, with a focus on PGAP

Annotating a bacterial genome involves two major steps: structural annotation and functional annotation. Structural annotation is to identify all relevant genomic sequences for protein coding genes, structural RNA genes, as well as other types of genomic features. The predicted protein coding gene sequences will be used for functional annotation to predict the functions of the proteins, and/or relationship to known pathways.

![Pipeline](./annotation_figures/annopipeline.png)


There are a few annotation pipelines designed for annotating bacterial genomes.

1. NCBI's Prokaryotic Genome Annotation Pipeline (PGAP)[](https://github.com/ncbi/pgap)
2. Prokka [prokarytoic annotation](https://github.com/tseemann/prokka)
3. RAST [Rapid Annotations using Subsystem Technology](https://rast.nmpdr.org/)
4. DRAM [Distilled and Refined Annotation of Metabolism](https://github.com/shafferm/DRAM)

## Introduction to PGAP
Running PGAP annotation for newly assembled bacterial genomes is easy to setup and the requirement in computing resources is reasonable (8 CPUs with 16GB memory or higher). The advantages of using PGAP is not only that it produces NCBI/GenBank recognized file formats, but also the extremely well organized and curated databases used in PGAP. 

![PGAP](./annotation_figures/PGAP_flowchart.png)


## Running PGAP
The computing resources required by PGAP is a machine has 8 CPUs and minimum 2GB memory per CPU. Though I have successfully run PGAP on a laptop that has 4 Cores and 32GB memory. So, the number of CPUs could be lower than 8, but the memory per CPU should be sufficient. Otherwise, one important step in the pipeline (blastp) will fail. PGAP can be run quite easily using docker. The first step in setting up the environment for running PGAP is to install docker on your laptop, or on an AWS instance. Simply download the proper version of [docker](https://docs.docker.com/get-docker/). The second package required for PGAP is python3.6 or higher. This could be easily installed by installing [Anaconda](https://www.anaconda.com/products/individual#Downloads) on your machine.

Once the above two components are properly installed, one could follow the instructions below to obtain PGAP, as well as the most up-to-date databases.

    curl -OL https://github.com/ncbi/pgap/raw/prod/scripts/pgap.py
    ./pgap.py --update

The input to PGAP is the final assembled genome fasta file. One could download the file using the following command ** if the previous assembly step has not finished **.

    wget https://bioshare.bioinformatics.ucdavis.edu/bioshare/download/ua4gjm2jyj12oxl/baccanu.contigs.fasta

For those who have Windowns machines, the file could be found [here](https://bioshare.bioinformatics.ucdavis.edu/bioshare/view/2021-ASM-genome-assembly/). Please right click on the file *baccanu.contigs.fasta* and save it to your working directory.


The _Bacillus thuringiensis_ RZ2MS9 is a plant growth-promoting bacteria, so it harbors several genes related with plant growth-promoting traits, such as the production of indole acetic acid, solubilization of phosphate, and more.

![Article](../fig_bact_tax/Batista_2018.png)

Read the paper: [Batista et al. 2018](../https://www.sciencedirect.com/science/article/pii/S0944501317309229).

## Contacts

## Computing needs

The quantity of data and the amount of processing needed should be sufficient on most laptop systems; however, many of the application used in the workshop requires the ability to compile code on a command line. As such we expect you to have these tools available on your system in order to fully participate.

### The applications that need to be first installed are:

    1. [samtools](http://www.htslib.org/)

     
## Workshop Materials

* This Bacterial Genome Assembly Workshop

   [https://ucdavis-bioinformatics-training.github.io/2021-ASM-genome-assembly/](https://ucdavis-bioinformatics-training.github.io/2021-ASM-genome-assembly/)
