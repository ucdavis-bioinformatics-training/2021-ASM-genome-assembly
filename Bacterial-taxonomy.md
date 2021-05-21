# Bacterial Taxonomy with Whole-genome Data
Dr. Maria Bonatelli; Univ. of São Paulo, Piracicaba, Brazil.

## Introduction to Bacterial Taxonomy

### What is bacterial taxonomy?

Bacterial taxonomy consists of three things: classification, nomenclature, and identification of bacteria.

Definition of [Cowan et al. 1964](https://www.microbiologyresearch.org/content/journal/micro/10.1099/00221287-39-1-143?crawler=true):

>(1) Classification, the orderly arrangement of units into groups of larger units. A simple analogy can be found in a pack of cards; the individual cards can first be sorted by color, then into suits. Within each suit the cards can be arranged in a numerical sequence, and the face cards placed in some order of seniority.
>
>(2) Nomenclature, the naming of the units defined and delineated by the classification. In the example of cards, the face cards are given names and more than one name, for example, jack or knave, may be given to the same card.
>
>(3) Identification of unknown units with known units of the classification developed in (1) and bearing names given in (2).

In the beginning, phenotypic characteristics were used to classify and identify bacteria (e.g., morphological and physiological). Now, we use a polyphasic taxonomy that integrates phenotypic, genotypic, and phylogenetic data to that end.

## 16S rRNA gene

### What is 16S rRNA gene?
The 16S rRNA gene encodes for the structural element of the small subunit of the ribosome. This gene is ubiquitous in prokaryotes and due to the highly conserved nature of the 16S rRNA, it became an important phylogenetic marker for this group.

Interestingly, the 16S rRNA contains several conserved and variable regions – that are connect with its structure.

![16S RNAr gene](fig_bact_tax/16S_domain_structure.png) 

Stacy Yudina, CC BY-SA 4.0 <https://creativecommons.org/licenses/by-sa/4.0>, via Wikimedia Commons.

Such regions are target in phylogenetic studies to access different regions of the 16S rRNA, depending on the objective of the study.

For our purpose here today, we want to recover from the genome the full length rRNA gene, that has ~ 1,500 bp. Typically, we will consider for species circumscription of two bacterial isolates if they have 98.7% of rRNA 16S gene identity.

>Tip: If you have a new bacterial isolate that you want to identify it using 16S rRNA gene, you should use primers that target different regions of the gene to get the full-length 16S rRNA gene. In this publication of (Thompson et al. 2001)[ https://www.sciencedirect.com/science/article/pii/S0723202004700673] you can find a set of primers to that end.

### Using Barrnap to recover 16S rRNA gene from your genome sequence
You have a high-quality complete genome, and you want to find all ribosomal RNA gene. To that, you can use [Barrnap](https://github.com/tseemann/barrnap).

>Barrnap predicts the location of ribosomal RNA genes in genomes. It supports bacteria (5S,23S,16S), archaea (5S,5.8S,23S,16S), metazoan mitochondria (12S,16S) and eukaryotes (5S,5.8S,28S,18S).
>
>It takes FASTA DNA sequence as input, and write GFF3 as output. It uses the new ‘nhmmer’ tool that comes with HMMER 3.1 for HMM searching in RNA:DNA style.

To use barrnap, you can install it via Conda or Homebrew/Linuxbrew, depending on your operational system. If you are working in a cluster, you can look for ‘prokka’ software that also runs barrnap. For more information about prokka, a cool software that performs whole genome annotation, [see here]( https://github.com/tseemann/prokka).

Once you have barrnap ready to run, you will only need your genome in .fasta file. In this example, I am using barrnap 0.9:

`barrnap --help`

![barrnap_Help](fig_bact_tax/Figures_barrnap1.png)

`barrnap -o 16S_bacillus.fasta < bacillus.fasta > 16S_bacillus.gff3`

`head 16S_bacillus.gff3`

![barrnap_gff3](fig_bact_tax/Figures_barrnap2.png)

More information about .gff3 file [here](https://m.ensembl.org/info/website/upload/gff3.html).

Bacteria can have more than one copy of the 16S gene. Our isolate here have 14 copies of this gene.

`head -n2 16S_bacillus.fasta`

![barrnap_seq](fig_bact_tax/Figures_barrnap3.png)

Once you have one or more 16S rRNA genes from your genome, you can use different databases to identify your bacteria. I recommend three: 
1.	Ribosomal Database Project [Classifier Tool](https://rdp.cme.msu.edu/help/CL.jsp) allows classification of both bacterial and archaeal 16S rRNA sequences to the taxonomical hierarchy.

2.	Silva rRNA database project [ACT: Alignment, Classification and Tree Service](https://www.arb-silva.de/aligner/).

>**SINA Aligner**
>SINA (v1.2.11) will align your rRNA gene sequences accoding to the global SILVA alignment for rRNA genes. The results can be combined with any other sequences aligned by SINA or taken from the SILVA databases by concatenation of FASTA files or using the ARB MERGE tool.
>
>**SINA Search and Classify**
>Enabling "Search and classify" will force SINA to additionally classify your sequences with the least common ancestor (LCA) method based on the taxonomies hosted by SILVA.

3.	Blastn of NCBI

## Average Nucleotide Identity


