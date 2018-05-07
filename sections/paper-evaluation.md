
# Evaluation

In the previous section we showed how we are able to, in simple examples, isolate a conserved motif from symbolic data. And in more complex examples still isolate regions of conservation from aligned samples. In this section we will evaluate the UniDip algorithm by presenting a case study showcasing UniDip's abilities when working with real genomic sequences, and describe how UniDip will help in future motif discovery projects.

## Case Study with FOXK1 Transcription Factor Binding Site

We will evaluate UniDip on real genomic sequences by performing a sequence analysis case study on MA0852.2 the binding site to the FOXK1 transcription factor. As categorized in the JASPER database [@mathelier_a._fornes_o._arenillas_d.j._chen_c._denay_g._lee_j._shi_w._shyr_c._tan_g._worsley-hunt_r._et_al._jaspar_2015], FOXK1 is a transcription factor in the class of Fork head / winged helix factors, and is known to regulate antiviral gene expression [@panda_transcription_2015].

The JASPER database contains 1056 instances of this TFBS, which compile into the below logo. 

![](./imgs/MA0852.2.svg)

The instances themselves are not terribly useful for testing. While they do have some of the surrounding sequence, the extensions are short and centered around the motif instance. JASPER does provide a `.bed` file which contains the genomic coordinates. For the latest version of MA0852.2 all sequences are located in the GRCh38.p12 human genome assembly [@ncbiresourcecoordinators_database_2017]. Our preferred method for extracting these motif instances would be to track the gene associated with each instance and pull the 1000bp upstream sequence from the transcription start site. Unfortunately, JASPER does not track the associated genes. So instead, we emulate how the motif instances would likely occur throughout the 1000bp sequences, by randomly selecting 1000 bases up and down stream of the motif instance such that the overall sequence remains 1000bp long. 

We wrote a python script that downloads the `.bed` file and the human genome's individual `.fasta` files for each chromosome. The same script then maps the motif locations to the corresponding chromosomal positions before extracting and saving the sequence instances as a single `.fasta` file. JASPER's record of MA0852.2 contains instances across all chromosomes except the Y chromosome.

To test our sequence extraction process and have a control without any preprocessing, we ran MEME on the raw MA0852.2 sequences [@bailey_meme_2009]. The process took 16 hours of real time on a single CPU, and resulted in finding the MA0852.2 motif in the top three results nearly exactly to how the logo is represented in the JASPER database.

![](./imgs/MA0852RAW.png)

Having confirmed our extraction method, we test the preprocessing steps by aligning the sequences with MUSCLE and trimming the gaps after calculating the SNE by position. *Fig X-X* show the sequence conservation from the original sequences, after alignment, and after gap trimming. 

![](./imgs/RealOriConservation.png)

![](./imgs/RealAlignedConservation.png)

![](./imgs/RealTrimmedConservation.png)

Among the first things to notice as different between the real sequence and the prior simulated examples is the length of the sequences. An 8bp section of higher conservation, like the MA0852.2 motif, might be easy to dismiss as an outlier even if it were aligned. And indeed, the sequences are not aligned; the highest level of conservation is 0.012 from a maximum of 2 bits of information. Running UniDip on the raw sequences would be largely meaningless as there does not seem to be any significant changes in conservation level along the sequence.

After alignment and trimming, the max information content as measured by SNE has increased from 0.012 to approximately 0.68. From alignment and trimming the sequence lengths have also shortened to less than 900bp long. Running the UniDip algorithm isolates more specific regions. Running the algorithm with the default significance threshold of $\alpha = 0.05$ we find 2 regions of higher conservation. Raising this threshold would allow us to isolate more and smaller clusters, but at a higher risk of the clusters merely being a phenomenon of the sampling and processing steps. MEME and its fellow algorithms are designed to function on sequences with some background noise, which gives another reason to not make UniDip too specific.

![](./imgs/RealTrimmedCluster.png)

UniDip may not be able to directly isolate the motif, but the alignment tool and UniDip working in conjunction are able to locate enhanced concentrations of motif instances in the aligned sequences. UniDip's left-most cluster covers 33% of the sequence but overlaps 40% of the motif instances, the second cluster likewise covers 20% of the sequence length but overlaps 25% of the motif instances. This concentration of motif instances enhances MEME's performance. We noted above that MEME took over 16 hours to find the motif from the raw sequences. This is mostly due to how many possible combinations there are in 1056 sequences 1000bp long. Limiting the search space to the isolated clusters shrinks that search space considerably. Indeed, even with the large cluster 330bp long the search space has been reduced by an order of magnitude. 

It is hard to separate how much of a role MUSCLE plays compared to UniDip. MUSCLE is aligning the sequences, allowing effective use of SNE. But, after initial trimming, sequences are still 90% of their initial length. It is only after running UniDip that specific regions are isolated to 40% and 20% of the raw lengths. Further examinations should study the difference in role, and compare UniDip's selections to random selections. Such a study is dependent on the sequences given, but we are confident that in datasets where alignment levels are variable across the sequence's length UniDip is able to isolate those regions of highest conservation and concentration.

Having isolated 2 regions of higher conservation, we are posed the question of which region is the better to search first. We can rank these regions by there average conservation in the cluster. Ranking the regions we see that the left region has an average SNE of 0.31 and the right region has 0.28. They are very close, so we should search both, but the left regions seems to have slightly more conserved peaks. We use an average both to control for the length of the isolated region and because the average is sensitive to outliers. We want to search those regions with particularly high peaks and exclude those regions with particularly low gaps.

Running MEME on the original sequences as selected by UniDip we again find the MA0852.2 motif in the top three matches of both regions, albeit in 30% of the time. MUSCLE alignment taking 30 minutes and UniDip clustering taking 10 minutes. MEME on the left isolated cluster taking 150 minutes and the right cluster 90 minutes. For a total time of 280 minutes or 4 hours and 40 minutes.

![](./imgs/MA0852TRIMMED.png) 

From the total time of 280 minutes, UniDip took just 10 of those minutes. UniDip, in this implementation, is still largely unoptimized and has ample room for improvements in speed. This speed is the power of aggregated metrics. Transforming the genomic data into a measure of conservations allows for a fast understanding of that variable. Due to the detail loss from aggregation, UniDip may not be able to always function as a standalone motif finding algorithm, but it can play a strong role in preprocessing genomic sequences. Drastically reducing the time it takes to run current MFAs like MEME. What this case study shows is that UniDip is beneficial to the search for biological motifs in real genomic data. 
