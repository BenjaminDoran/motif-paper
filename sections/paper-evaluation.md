
# Evaluation

In the previous section, we showed UniDip can, in generated examples, isolate a conserved motif from symbolic data. And, in more complex simulations still isolate regions of conservation from aligned samples. In this section, we will evaluate the UniDip algorithm by presenting a case study showcasing UniDip's abilities when working with real genomic sequences, and describe how UniDip will help in future motif discovery projects.

## Case Study with FOXK1 Transcription Factor Binding Site

We will evaluate UniDip on real genomic sequences by performing a sequence analysis case study on MA0852.2, the binding site to the FOXK1 transcription factor. As categorized in the JASPER database [@mathelier_a._fornes_o._arenillas_d.j._chen_c._denay_g._lee_j._shi_w._shyr_c._tan_g._worsley-hunt_r._et_al._jaspar_2015], FOXK1 is a transcription factor in the class of Fork head / winged helix factors, and is known to regulate antiviral gene expression [@panda_transcription_2015].

The JASPER database contains 1056 instances of this TFBS, which compile into the below logo. 

![Fig. 17: JASPER logo for MA0852.2](./imgs/MA0852.2.svg)

The instances in the JASPER database are not useful for evaluating UniDip on real sequences. While the JASPER instances do have some of the surrounding sequence, the extensions are short and centered around the motif instance. JASPER does provide a `.bed` file which contains the genomic coordinates. For the latest version of MA0852.2 all sequences are located in the GRCh38.p12 human genome assembly [@ncbiresourcecoordinators_database_2017]. Our preferred method for extracting these motif instances would be to track the gene associated with each instance and pull the 1000bp upstream sequence from the transcription start site. Unfortunately, JASPER does not track the associated genes. So instead, we emulate how the motif instances would likely occur throughout the 1000bp sequences, by randomly selecting 1000 bases up and down stream of the motif instance such that the overall sequence remains 1000bp long. 

We wrote a Python script that downloads the `.bed` file and the human genome's individual `.fasta` files for each chromosome. The same script then maps the motif locations to the corresponding chromosomal positions before extracting and saving the sequence instances as a single `.fasta` file. JASPER's record of MA0852.2 contains instances across all chromosomes except the Y chromosome.

To test our sequence extraction process and have a control without any preprocessing, we ran MEME on the raw MA0852.2 sequences [@bailey_meme_2009]. The process took 16 hours of real time on a single CPU, and resulted in finding the MA0852.2 motif in the top three results nearly exactly to how the logo is represented in the JASPER database.

![Fig. 18: MEME logo of MA0852.2 discovered from raw sequences](./imgs/MA0852RAW.png)

Having confirmed our extraction method, we test the preprocessing steps by aligning the sequences with MUSCLE and trimming the gaps after calculating SNE by position. *Fig 19* shows the sequence conservation from the original sequences, after alignment, and after gap trimming. 

![](./imgs/RealOriConservation.png)
![](./imgs/RealAlignedConservation.png)
![](./imgs/RealTrimmedConservation.png)

*Fig. 19: Raw Sequences (top), Aligned Sequences (middle), Aligned and Trimmed Sequences (bottom)*

Among the first things to notice as different between the raw sequences and the prior simulated examples is the lack strong conservation in specific regions. The highest level of conservation is 0.012 from a maximum of 2 bits of information. Running UniDip on across SNE calculated from the raw sequences would be meaningless as there are not any significant changes in conservation level drawn out from the metric.

However, after alignment and trimming, the maximum information content as measured by SNE has increased from 0.012 to approximately 0.68. From alignment and trimming, the sequence lengths have shortened to less than 900bp long, and there are clear peaks in conservation along the sequences. Running the UniDip algorithm does not isolate the smaller peaks, but it does isolate the regions where these peaks reside. With the default significance threshold of $\alpha = 0.05$, UniDip detects two regions of higher conservation. Raising this threshold would allow us to isolate more and smaller clusters, but at a higher risk of the clusters merely being a phenomenon of the sampling and processing steps. MEME and its fellow algorithms are designed to function on sequences with some background noise, which gives another reason to not make UniDip too specific.

![Fig. 20: UniDip isolates regions of higher conservation](./imgs/RealTrimmedCluster.png)

There may not be a clear indication of a single motif like in the simulated examples, but the alignment tool and UniDip working in conjunction are able to locate enhanced concentrations of motif instances in the aligned sequences. The left cluster covers 33% of the sequence but overlaps 40% of the motif instances. Likewise, the right cluster covers 20% of the sequence length but overlaps 25% of the motif instances. This concentration of motif instances enhances MEME's performance. We noted above that MEME took over 16 hours to find the motif from the raw sequences. This is mostly due to how many possible combinations there are in 1056 sequences 1000bp long. Limiting the search space to the isolated clusters shrinks that search space considerably. Indeed, even with the larger cluster being 330bp long, the search space is reduced by an order of magnitude. 

It is hard to separate how much of a role MUSCLE plays compared to UniDip. MUSCLE is aligning the sequences, allowing effective use of SNE. But, after initial trimming, sequences are still 90% of their initial length. It is only after running UniDip that specific regions are isolated to 33% and 20% of the raw lengths. Further examinations should study the difference in role, and compare UniDip's selections to random selections. Such a study is dependent on the sequences given, but we are confident that in datasets where alignment levels are variable across the sequence's length UniDip is able to isolate those regions of highest conservation and concentration.

We can rank the relevance of the two isolated regions by their average conservation. This ranking allows prioritization of the regions, such that regions with more prominent motifs are searched first. Ranking the regions, we see that the left region has an average SNE of 0.31 and the right region has 0.28. They are very close, so we should search both, but the left region seems to have slightly more conserved peaks. We use an average both to control for the length of the isolated region and because the average is sensitive to outliers. We want to search those regions with particularly high peaks and exclude those regions with particularly low gaps.

Running MEME on the isolated regions, we find the MA0852.2 motif in the top three matches of both regions 70% faster than discovery took on the raw sequences. Including preprocessing, MEME took a total time of 280 minutes or 4 hours and 40 minutes compared to the original 16 hours. Breaking down the total time, MUSCLE alignment took 30 minutes; UniDip clustering took 10 minutes; MEME on the left isolated cluster took 150 minutes; and MEME on the right cluster took 90 minutes. 

![Fig. 21: MEME logo of MA0852.2 discovered from left-most isolated cluster](./imgs/MA0852TRIMMED.png) 

From the total time of 280 minutes, UniDip took just 10 of those minutes. UniDip, in this implementation, is still largely unoptimized and has ample room for improvements in speed. This speed is a result of aggregated metrics' power in condensing information. Transforming the genomic data into a univariate measure of conservation allows for a rapid filtering based on that variable. Due to the detail loss from aggregation, UniDip may not be able to always function as a standalone motif finding algorithm. But, it can play a strong role in preprocessing genomic sequences, greatly reducing the time it takes to run modern MFAs like MEME. This case study shows that UniDip is beneficial to the search for biological motifs in real genomic data. 
