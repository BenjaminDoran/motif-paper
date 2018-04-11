
# Evaluation

In the previous section we showed how we are able to, in simple examples, isolate a conserved motif from symbolic data. In this section we will evaluate these methods by determining UniDip's ability to detect motifs and conserved regions of variable length.

After testing with generated data, we will present a case study showcasing UniDips abilities when working with real genomic sequences, and how UniDip will help in future motif discovery projects.

## Motif Length Evaluation
<!-- 4, 8, 16, 24, 32 -->

Motifs come in many sizes anywhere from 8-30bp, to evaluate UniDip's ability to find these motifs we will test against motifs of lengths: 4, 8, 16, 24, and 32. For the purposes of this test we will hold our other variables constant such that we are testing on generated datasets of 20 sample sequences 500bp long. We will not include degenerate mutations in this test, but we will be allowing the motif instance to fall in any position along the sequence instance. We will be measuring performance by calculating the ratio between how many motif instances fall within a cluster and the percent coverage of the cluster across the aligned sequence. A cluster that covers the entire sequence and captures every motif instance would have a capture ratio of 1, cluster that covers half the sequence and captures every motif instance will have a capture ratio of 2. We choose this metric because it shows us how concentrated the motifs are in a set region. Allowing us to measure how much the global alignment is boosting the motif instance conservation. As well as UniDip's specificity in finding those regions.

Visualizing this capture ratio we find unsurprisingly that larger motifs are easier to detect. Comparing an 8bp and a 30bp motif in a sequence of 1000bp, the chance of perfectly aligned pair of motif instances changes from `1/992` to `1/970`(assuming the motif does not fall of the end of the sequence). However the chance that the pair will overlap changes from `16/992` to `60/970`, neither are hugely likely but one becomes inevitable much sooner. There are only 33 unique non-overlapping mappings onto the 1000bp sequence for a 30bp motif, compared to the 125 unique mappings for the 8bp motif. What this means in practice is that it is much easier for MUSCLE to get larger motifs into alignment then smaller motifs.

![](./imgs/EvalMotifLength.png)

Due to the prevalence of positive ratios, the above plot shows that MUSCLE and UniDir are boosting and locating the regions of higher conservation, even in the trials with short motifs. Once the motifs become long enough between the 8-16bp range we see a slight increase in the ability of MUSCLE to align the motif instances and UniDip to locate these aligned regions. 

UniDip does seem to work best on longer motifs and conserved regions. And is reliant on well aligned sequences to directly isolate motifs. We can enhance UniDip's abilities by preparing our data, such as with alignment and gap trimming. And, even when not able to directly locate a motif, we have shown that MUSCLE and UniDip can work together to locate higher conserved subregions in a dataset, which could assist in narrowing the focus of other motif discovery algorithms.

## Case Study with FOXK1 Transcription Factor Binding Site

We will evaluate UniDip on real genomic sequences by performing a sequence analysis case study on MA0852.2 the binding site to the FOXK1 transcription factor, as categorized in the JASPER database [@mathelier_a._fornes_o._arenillas_d.j._chen_c._denay_g._lee_j._shi_w._shyr_c._tan_g._worsley-hunt_r._et_al._jaspar_2015]. FOXK1 is a transcription factor in the class of Fork head / winged helix factors, and is known to regulate antiviral gene expression [@noauthor_transcription_nodate] among many other predicted functions [@noauthor_foxk1_nodate].

The JASPER database contains 1056 instances of this TFBS, which it compiles into the below logo. 

![](./imgs/MA0852.2.svg)

The instances themselves are not terribly useful for testing. While they do have some of the surrounding sequence, the extensions are short and centered around the motif instance. JASPER does provide a `.bed` file which contains the genomic coordinates. For the latest version of MA0852 all sequences are located in the GRCh38.p12 human genome assembly [@ncbiresourcecoordinators_database_2017]. Our preferred method for extracting these motif instances would be to track the gene associated with each instance and pull the 1000bp upstream sequence from the transcription start site. Unfortunately, JASPER does note track the associated genes. So, instead we emulate how the motif instances would likely occur throughout the 1000bp sequences, by randomly selecting 1000 bases up and down stream of the motif instance such that the overall sequence remains 1000bp long. 

We wrote a python script that downloads the `.bed` file and the human genome individual `.fasta` files for each chromosome. The same script then maps the motif locations to the corresponding chromosomal positions before extracting as saving the sequence instances as a single `.fasta` file. JASPER's record of MA0852 contains instances across all chromosomes except the Y chromosome.

To test our sequence extraction process and have a comparison to a production read motif finding algorithm, we ran MEME on the MA0852 sequences.[@bailey_meme_2009] The process took 6 hours of real time, and resulted in finding the MA0852 motif although reversed due to a slight issue with handling the two strands of DNA [Will be fixed in next draft].

![](./imgs/MA0852RAW.png)

Having confirmed our process, we continue to align the sequences with MUSCLE and trim the gaps after calculating the SNE by position. For reference we see the sequence conservation from the original sequences, after alignment, and after gap trimming. 

![](./imgs/RealOriConservation.png)
![](./imgs/RealAlignedConservation.png)
![](./imgs/RealTrimmedConservation.png)

Running UniDip on this trimmed data we can see that we locate 2 regions of higher conservation, where MUSCLE has made better alignments. 

![](./imgs/RealTrimmedCluster.png)

Due to the motif length, it is unsurprising that we are not able to directly isolate the motif. As we saw in our simulated evaluations, a motif of 8bp had 992 unique non-aligned mappings onto a sequence of 1000bp. In our unaligned example where we have 1056 sample sequences, we can expect overall 65 motif instance pairs to be perfectly aligned. These will still be drowned out by the other 1054 background nucleotides at that position. The global alignment will help bring nearby matches together, but the alignment tool does have its limits. 

We may not be able to directly isolate the motif, but the alignment tool and UniDip are able to locate enhanced the concentrations of motif instances in the aligned sequences. UniDip left-most cluster covers 33% of the sequence but overlaps 40% of the motif instance, the second cluster likewise covers 20% of the sequence length but overlaps 25% of the motif instances. We can use this feature to enhance MEME's performance. We noted above that MEME took over 6 hours to find our motif from the raw sequences. This is mostly due to how many possible combinations there are in 1056 sequences 1000bp long. Limiting our search space to our isolated clusters shrinks that search space considerably. Indeed, even with our large cluster of (330bp) long we are shrinking the search space by an order of magnitude. 

Running MEME on the original sequences as selected by UniDip we find find the MA0852 motif in the top three matches, albeit in half the time. MUSCLE alignment taking half an hour, UniDip cluster taking 10 minutes,and MEME on the isolated cluster sequences taking a little over two and a half hours. 

![](./imgs/MA0852TRIMMED.png) 

This real biological is not the best candidate for UniDips search capabilities, it is smaller than it needs to be to give a strong conservation signal in alignment. UniDip may not be able to always function as a standalone motif finding algorithm, but few motif finding algorithms do. What this case study shows is that even in this challenging circumstance UniDip is beneficial to the search for biological motifs. 

