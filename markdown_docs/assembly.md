# Assembly, Assembly Metrics, and Quality Control

*Hannah Lyman, PhD, University of California, Davis*

## Assembling a bacterial genome with PacBio sequence data

The data and pipeline in this tutorial are adapted from [**Complete Genome Sequence of *Bacillus* sp. Strain Rz2MS9, a Multitrait Plant Growth Promoter**](https://journals.asm.org/doi/full/10.1128/MRA.00623-20 When referring to this pipeline, please cite the publication.

### Explore the sequence data

We will be working with two PacBio sequencing runs from a *Bacillus* species. Before we begin the assembly, we can first take a look at the reads, beginning with the demultiplexed output from Lima.

#### How many reads are there?

The bam stats tool from [bamutil](https://github.com/statgen/bamUtil/) allows us to quickly count the reads in each run.

    mkdir 01-bamStats
    bam stats --in 00-RawData/isi_run_01/lima_output.bc1002_BAK8A_OA--bc1002_BAK8A_OA.bam \
              --basic \
              2> 01-bamStats/isi_run_01.txt # redirect stderr to output file
    bam stats --in 00-RawData/isi_run_02/lima_output.bc1002_BAK8A_OA--bc1002_BAK8A_OA.bam \
              --basic \
              2> 01-bamStats/isi_run_02.txt

#### What does the ZMW count distribution look like?

First, we will need to count the occurrence of each ZMW. Using [samtools](https://github.com/samtools/samtools/releases/tag/1.9) and a few basic text manipulation tools on the command line, we can derive a comma separated counts table of ZMWs from the PacBio read data. The multi-line command below will:

1. Create a CSV file containing a header
2. Convert from SAM to BAM
3. Extract the read name field from the SAM formatted read file
4. Extract the ZMW from the read name
5. Sort ZMWs as a prerequisite to counting occurences
6. Display counts for each ZMW
7. Remove excess whitespace
8. Reorder columns, separating them with a comma
9. Append output to the CSV

While steps 2 through 9 could be run on a single line, we have split the command over multiple lines to make it easier to see what is going on.

    mkdir 02-zmwCounts
    echo "zmw,counts" > 02-zmwCounts/isi_run_01.csv
    samtools view 00-RawData/isi_run_01/lima_output.bc1002_BAK8A_OA--bc1002_BAK8A_OA.bam | \
    cut -f1 | \
    cut -d '/' -f2 | \
    sort | \
    uniq -c | \
    tr -s ' ' | \
    awk -F ' ' '{print $2 "," $1}' \
    >> 02-zmwCounts/isi_run_01.csv
    echo "zmw,counts" > 02-zmwCounts/isi_run_02.csv
    samtools view 00-RawData/isi_run_02/lima_output.bc1002_BAK8A_OA--bc1002_BAK8A_OA.bam | \
    cut -f1 | \
    cut -d '/' -f2 | \
    sort | \
    uniq -c | \
    tr -s ' ' | \
    awk -F ' ' '{print $2 "," $1}' \
    >> 02-zmwCounts/isi_run_02.csv

The ZMW counts CSV files can be imported into R for summary statistics and visualizations.

    run_01 <- read.csv("02-zmwCounts/isi_run_01.csv")
    dim(run_01)
    [1] 8274    2
    sum(run_01$counts)
    [1] 66292
    summary(run_01$counts)
    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
    1.000   2.000   7.000   8.012  12.000  97.000
    run_02 <- read.csv("02-zmwCounts/isi_run_02.csv")
    dim(run_02)
    [1] 15686     2
    sum(run_02$counts)
    [1] 122450
    summary(run_02$counts)
    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
    1.000   2.000   6.000   7.806  11.000 183.000

![Histogram of ZMW counts generated in R](assembly_figures/zmw_hist.pdf)

While the second run has approximately twice as many ZMW as the first, the distribution of ZMW counts is similar between the two runs.

### Convert CLR to CCS

The reads generated from the two PacBio runs are consenusus long reads (CLR). However, circular consensus sequencing (CCS) produces a better a more accurate, higher quality sequence by creating consenus sequence from sub-reads. It's possible to convert CLR to CCS, and if enough read coverage remains after conversion, these reads will generate a higher quality assembly. In order to convert our CLR to CCS, we will run PacBio's [ccs/4.2.0](https://github.com/PacificBiosciences/ccs/).

    ccs --reportFile=03-CCS/isi_run_01.report \
         00-RawData/isi_run_01/lima_output.bc1002_BAK8A_OA--bc1002_BAK8A_OA.bam \
         03-CCS/isi_run_01.ccs.bam

    ccs --reportFile=03-CCS/isi_run_02.report \
        00-RawData/isi_run_02/lima_output.bc1002_BAK8A_OA--bc1002_BAK8A_OA.bam \
        03-CCS/isi_run_02.ccs.bam

### Assemble the genome

The CCS reads are ready to assemble. [Canu](https://canu.readthedocs.io/en/latest/index.html) takes FASTA formatted reads as input, so we will need to convert our BAM formatted CCS reads to FASTA using samtools. The Canu pipeline consists of three tasks: read correction, read trimming, and unitig construction. Because CCS reads are already corrected, we will only run the final assembly step of the pipeline, using the "hifi" flag to let Canu know that we are providing corrected reads. This is referred to as "HiCanu."

    samtools fastq 03-CCS/isi_run_01.ccs.bam > 03-CCS/isi_run_01.ccs.fq
    samtools fastq 03-CCS/isi_run_02.ccs.bam > 03-CCS/isi_run_02.ccs.fq
    mkdir 04-HiCanu
    canu -assemble \
         useGrid=False \
         -p bacillus \
         -d 04-HiCanu \
         genomeSize=5.5m \
         -pacbio-hifi 03-CCS/isi_run_0?.ccs.fq

The HiCanu assembly produces a single contig of 5364948 bp in length. We will take a closer look at it in later steps.

If there are not sufficient reads to assemble after converting to CCS, we can assemble the CLR reads instead. In this case, we run the entire Canu pipeline, rather than the assembly step by itself.

    samtools fastq 00-RawData/isi_run_01/lima_output.bc1002_BAK8A_OA--bc1002_BAK8A_OA.bam > clr_assembly/isi_run_01.clr.fq
    samtools fastq 00-RawData/isi_run_02/lima_output.bc1002_BAK8A_OA--bc1002_BAK8A_OA.bam > clr_assembly/isi_run_02.clr.fq
    canu useGrid=False -p bacillus.clr -d clr_assembly genomeSize=5.5m -pacbio clr_assembly/isi_run_0?.clr.fq

The CLR assembly is composed of 4 contigs, ranging in length from 16328 bp to 5367862 bp.

### Polish the assembly using PacBio reads

PacBio's [GenomicConsensus](https://github.com/PacificBiosciences/GenomicConsensus) package includes tools for variant calling and polishing. Here we use the arrow algorithm, version 2.3.3. GenomicConsensus's variantCaller tool uses a file of file names ([FOFN](https://pb-falcon.readthedocs.io/en/latest/tutorial.html#create-fofn)), which is simply a text file containing the path to each read file. Before running variantCaller, we need to create the FOFN, align the PacBio reads to the assembly, and index the assembly FASTA file.


    mkdir 05-GenomicConsensus
    ls 00-RawData/isi_run_0*/*.bam > raw.fofn
    pbalign raw.fofn 04-HiCanu/bacillus.ccs.contigs.fasta 05-GenomicConsensus/raw_canu_align.bam
    samtools flagstat 05-GenomicConsensus/raw_canu_align.bam -O tsv > 05-GenomicConsensus/samFlags.tsv
    samtools faidx 04-HiCanu/bacillus.ccs.contigs.fasta
    variantCaller 05-GenomicConsensus/raw_canu_align.bam \
                  --algorithm=arrow \
                  -r 04-HiCanu/bacillus.ccs.contigs.fasta \
                  -o 05-GenomicConsensus/variants-raw.gff \
                  -o 05-GenomicConsensus/bacillus.arrow.fasta

The read-correction process is extremely effective. GenomicConsensus identifies 10 variants: 5 single basepair insertions and 5 small deletions in the raw reads relative to the reference.

    tail 05-GenomicConsensus/variants-raw.gff
    tig00000001	.	insertion	1598	1598	.	.	.	reference=.;variantSeq=A;coverage=22;confidence=41
    tig00000001	.	insertion	1770	1770	.	.	.	reference=.;variantSeq=T;coverage=22;confidence=68
    tig00000001	.	deletion	2382	2382	.	.	.	reference=G;variantSeq=.;coverage=23;confidence=93
    tig00000001	.	deletion	2438	2450	.	.	.	reference=CCCCTCCTCCCCC;variantSeq=.;coverage=23;confidence=93
    tig00000001	.	deletion	2452	2453	.	.	.	reference=CC;variantSeq=.;coverage=23;confidence=93
    tig00000001	.	insertion	297977	297977	.	.	.	reference=.;variantSeq=A;coverage=60;confidence=93
    tig00000001	.	insertion	301329	301329	.	.	.	reference=.;variantSeq=T;coverage=53;confidence=93
    tig00000001	.	insertion	302211	302211	.	.	.	reference=.;variantSeq=A;coverage=44;confidence=93
    tig00000001	.	deletion	3061047	3061048	.	.	.	reference=AT;variantSeq=.;coverage=100;confidence=93
    tig00000001	.	deletion	5363706	5363706	.	.	.	reference=A;variantSeq=.;coverage=30;confidence=71


### Incorporate short read data

In addition to the two PacBio runs, we have a paired-end Illumina library. Before the development of CCS, it was common practice to perform error correction on long read assemblies using more accurate short read technology. CCS reads are generally so accurate that this is not necessary. However, since we have the data, we will use it.

#### Preprocess reads

Before aligning the Illumina reads to the PacBio contigs, we will preprocess the reads to remove contaminants, duplicates, and low quality sequence using [HTStream](https://s4hts.github.io/HTStream/).

    mkdir 06-HTStream
    hts_Stats          --stats-file 06-HTStream/genomeStats.log \
                       --notes "initial statistics" \
                       --read1-input 00-RawData/RZ2M29_S1_L001_R1_001.fastq.gz \
                       --read2-input 00-RawData/RZ2M29_S1_L001_R2_001.fastq.gz | \
    hts_SeqScreener    --append-stats-file 06-HTStream/genomeStats.log \
                       --notes "remove PhiX" | \
    hts_SuperDeduper   --append-stats-file 06-HTStream/genomeStats.log \
                       --notes "deduplicate using 10bp key at 10th position" | \
    hts_Overlapper     --append-stats-file 06-HTStream/genomeStats.log \
                       --notes "merge pairs with >= 8bp overlap, with < 100 mismatches and >= 75% identity in overlapping region" | \
    hts_QWindowTrim    --append-stats-file 06-HTStream/genomeStats.log \
                       --notes "remove low quality bases from ends of reads. minimum window quality 20" | \
    hts_NTrimmer       --append-stats-file 06-HTStream/genomeStats.log \
                       --notes "return longest subsequence containing no Ns" | \
    hts_LengthFilter   --append-stats-file 06-HTStream/genomeStats.log \
                       --notes "remove all reads < 200bp, and all orphaned reads" \
                       --min-length 200 \
                       --no-orphans | \
    hts_Stats          --append-stats-file 06-HTStream/genomeStats.log \
                       --notes "final statistics" \
                       --fastq-output 06-HTStream/illumina.htstream

#### Align reads to assembly

The preprocessed reads can now be aligned to the assembly. We will use [bwa-mem2](https://github.com/bwa-mem2/bwa-mem2), a very fast aligner, followed by samtools to reformat, sort, and index the reads, and collect alignment statistics.

    mkdir 07-Pilon
    bwa-mem2 index 05-GenomicConsensus/bacillus.arrow.fasta
    bwa-mem2 mem 05-GenomicConsensus/bacillus.arrow.fasta 06-HTStream/illumina.htstream_SE.fastq.gz | \
    samtools view -bS | samtools sort > 07-Pilon/bacillus.arrow.SE.bam
    bwa mem 05-GenomicConsensus/bacillus.arrow.fasta 06-HTStream/illumina.htstream_R1.fastq.gz 06-HTStream/illumina.htstream_R2.fastq.gz | \
    samtools view -bS | samtools sort > 07-Pilon/bacillus.arrow.PE.bam
    samtools flagstat 07-Pilon/bacillus.arrow.SE.bam -O tsv > 07-Pilon/samFlags_SE.tsv
    samtools flagstat 07-Pilon/bacillus.arrow.PE.bam -O tsv > 07-Pilon/samFlags_PE.tsv
    samtools index 07-Pilon/bacillus.arrow.SE.bam
    samtools index 07-Pilon/bacillus.arrow.PE.bam

#### Polish assembly using short read data

[Pilon](https://github.com/broadinstitute/pilon/wiki) is a tool from the Broad Institute designed to identify small and large indels, block substitutions, local misassemblies, and single-base differences between the input genome and the provided reads.

    pilon --genome 05-GenomicConsensus/bacillus.arrow.fasta \
          --frags 07-Pilon/bacillus.arrow.PE.bam \
          --unpaired 07-Pilon/bacillus.arrow.SE.bam \
          --output 07-Pilon/bacillus.pilon \
          --changes \
          --vcf \
          --tracks

Running Pilon generates a series of genome browser tracks, a .changes file, and a corrected assembly FASTA. In this case, Pilon has made only one change to the genome sequence: correcting a 2 bp indel.

Further examination of the output files reveals that Pilon has identified a number of other regions where the short read alignments appeared "suspicious". This may mean that the coverage in the region dips or peaks abruptly, or that there are an unusually high number of improperly mapped read pairs in the region. Pilon will break the assembly at these locations, and attempt to re-assemble. In this case, Pilon was unable to find any solution for these regions, and no changes were made to the genome sequence as a result.

### Make final adjustments to assembly

#### Look for overlap between beginning and end of genome

Because the species we're assembling has a single circular chromosome, we may see overlap between the beginning and end of our genome sequence. Using samtools and BLAST on the command line, we can check the ends of our sequence to see if they overlap one another.

    mkdir 08-Final
    sed 's/|arrow|pilon//' 07-Pilon/bacillus.pilon.fasta > 08-Final/pilon.complete.fasta
    samtools faidx 08-Final/pilon.complete.fasta tig00000001:1-50000 --output 08-Final/pilon.beg.fasta
    samtools faidx 08-Final/pilon.complete.fasta tig00000001:5314937-5364936 --output 08-Final/pilon.end.fasta
    blastn -query 08-Final/pilon.end.fasta -subject 08-Final/pilon.beg.fasta -outfmt 6

The last 7744 bp of our contig are 100% identical to the first 7744 bp.

    tig00000001:5314937-5364936	tig00000001:1-50000	100.000	7744	0	0	42257	50000	1	7744	0.0	14301

We can resolve this by removing the redundant bases.

    samtools faidx 08-Final/pilon.complete.fasta tig00000001:1-5357192 --output 08-Final/pilon.trim.fasta


#### Verify strand orientation of assembly

The orientation of the assembly should be consistent with those of other genome sequences in closely related taxa. Because this assembly is of the genome of an unidentified *Bacillus* species, the conventional representation of the genome places the first base of the sequence at the beginning of the gene encoding DnaA, the chromosomal replication initiator protein. To verify the directionality of the strand, we can do a BLAST search against the protein sequence of dnaA in [another species](https://www.ncbi.nlm.nih.gov/protein/16077069).

    tblastn -query 08-Final/B.subtilis.subtilis.168.DnaA.fasta -subject 08-Final/pilon.trim.fasta -outfmt 6

This search returns a handful of results, only one of which is a full-length alignment.

    NP_387882.1	tig00000001:1-5357192	83.857	446	72	0	1	446	3526746	3528083	0.0	776
    NP_387882.1	tig00000001:1-5357192	34.737	95	54	3	122	216	2655266	2655006	3.44e-07	50.4
    NP_387882.1	tig00000001:1-5357192	33.333	48	31	1	70	116	937226	937369	0.45	30.4
    NP_387882.1	tig00000001:1-5357192	31.915	47	32	0	330	376	2507658	2507798	3.6	27.7
    NP_387882.1	tig00000001:1-5357192	35.897	39	25	0	220	258	3356476	3356360	3.9	27.3
    NP_387882.1	tig00000001:1-5357192	28.571	91	56	4	108	194	4979649	4979906	4.4	27.3
    NP_387882.1	tig00000001:1-5357192	28.333	60	38	2	222	279	2601964	2601794	5.8	26.9
    NP_387882.1	tig00000001:1-5357192	47.368	19	10	0	142	160	2009262	2009206	7.0	26.6

The start coordinate of the alignment is smaller than the end coordinate for both the query and the subject, which indicates that our assembly has the same orientation as that of the *Bacillus subtilis* strain we used as a reference. If, on the other hand, the query start and end coordinates had been reversed, we would need to reverse complement our genome sequence.

#### Set origin of replication

Finally, we will set the beginning of the chromosome to the first base of dnaA, as identified in our BLAST search. First we need to remove the coordinate information from the FASTA header  in order to allow samtools to correctly process the chromosme name.

    cut -d ":" -f1 08-Final/pilon.trim.fasta > 08-Final/pilon.trim.renamed.fasta

Then, we divide the contig in two at the first base of dnaA. Omitting one coordinate tells samtools to use the start or end coordinate by default. The length option ensures that our FASTA entries are not split over multiple lines. This will make the file easier to manually join in the next step.

    echo "tig00000001:3526746-" > 08-Final/regions.txt
    echo "tig00000001:-3526745" >> 08-Final/regions.txt
    faidx 08-Final/pilon.trim.renamed.fasta --region-file 08-Final/regions.txt --length 5357192 --output 08-Final/pilon.reorder.fasta

Finally, we can manually merge our split contig to form the correctly oriented chromosome, and rename the entry. To check that the sequence has the expected length, we can count the characters in the file.

    grep -v ">" 08-Final/genome.fasta | wc
          1       1 5357193

Excluding the header, there are 5357193 characters in the file, one of which is a newline character, which leaves us with a total of 5357192, the expected length of our genome sequence after Pilon removed 2 bp.

## Assembly Metrics and Quality Control

How well assembled is the genome? What are its characteristics? We can look at a few different metrics to get an idea of whether or not the assembly meets our expectations for the genome in terms of contiguity, completeness, and correctness.

### Contiguity

The size of an assembled genome tells us only how many base pairs of assembled sequence were generated. Looking only at assembly size, we have no way of knowing how those base pairs are distributed. Is this 5.5 MB assembly in 1 KB pieces? Ideally, an assembly is highly contiguous, with a single assembled sequence representing each chromosome or plasmid. We will look at two measures of assembly contiguity: N50 and L50.

#### N50

The N50 of an assembly is the length of the shortest contig in 50% of the total genome length. This calculated by ordering contigs from longest to shortest, summing lengths until reaching 50% of the total genome size, and reporting the length of the shortest contig used.

In a perfectly assembled single chromosome genome, the N50 and the size of the assembled genome will be identical. In a highly fragmented assembly, the N50 may be very small, as many small contigs will need to be included to reach 50% of the genome size. Here is a graphical representation of a hypothetical assembly:

![Graphical example of N50](assembly_figures/N50_example.pdf)

Our assembly is a single contig, so the N%0 is the same as the genome sequence length. Another measure that is used to assess contiguity is **NG50**, which uses the genome size where the N50 uses the assembly size, making the metric comparable between assemblies of different sizes. Since we only have one contig, the N50 and the NG50 are identical.

#### L50

The L50 is the smallest number of contigs whose length total at least 50% of the genome size. This is the number of contigs used in the N50. A lower number represents higher contiguity.

![Graphical example of L50](assembly_figures/L50_example.pdf)

Because we were able to assemble the sequence into a single contig, our genome has an L50 of 1.

### Completeness and correctness

The assembly size is very similar to the expected genome size of 5.5 MB. The variantCaller and Pilon tools identified and corrected a small number of assembly errors. We can be reasonably confident that the genome is complete and does not contain any major errors. One way to assess the completeness and correctness of the assembled sequence is with Benchmarking Universal Single-Copy Orthologs (BUSCO), which will be discussed in the next section.
