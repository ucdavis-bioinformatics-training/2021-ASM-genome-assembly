# Annotating Bacterial Genomes, with a focus on PGAP

Annotating a bacterial genome involves two major steps: structural annotation and functional annotation. Structural annotation is to identify all relevant genomic sequences for protein coding genes, structural RNA genes, as well as other types of genomic features. The predicted protein coding gene sequences will be used for functional annotation to predict the functions of the proteins, and/or relationship to known pathways.

![Pipeline](./annotation_figures/annopipeline.png)


There are a few annotation pipelines designed for annotating bacterial genomes.

1. NCBI's Prokaryotic Genome Annotation Pipeline (PGAP)[](https://github.com/ncbi/pgap)
2. Prokka [prokarytoic annotation](https://github.com/tseemann/prokka)
3. RAST [Rapid Annotations using Subsystem Technology](https://rast.nmpdr.org/)
4. DRAM [Distilled and Refined Annotation of Metabolism](https://github.com/shafferm/DRAM)

## PGAP
Running PGAP annotation for newly assembled bacterial genomes is easy to setup and the requirement in computing resources is reasonable (8 CPUs with 16GB memory or higher). The advantages of using PGAP is not only that it produces NCBI/GenBank recognized file formats, but also the extremely well organized and curated databases used in PGAP. 

![PGAP](./annotation_figures/PGAP_flowchart.png)

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