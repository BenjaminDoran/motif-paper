
# Introduction

## Biologic Motifs

### Definition
Biological motifs are short conserved degenerate sequences of similar or identical sequence instances as found in DNA, RNA, or Protein sequences. To clarify, genomic and proteomic sequences are determined in a lab from various sequencing technologies [@heather_sequence_2016], and are as accurate as possible to what one would find in a specific individual.  Sub-sections from similar locations of these individual sequences are then compared to determine similarities between these sub-sections. This comparison takes many forms, but it is this process that is responsible for how we discover biological motifs. Because individual sequences may vary slightly due to non-harmful mutations, when looking for motifs we are not looking for exact matches or sequences. Instead, we model an idealized sequence that takes into account probable mutations. This model is the biologic motif.

### Role
Biologic Motifs have a rather general definition and much like the grouping of genes may serve many functions. As a general rule, motifs are short often 8-30bp (base pairs) long [@das_survey_2007]. In the genome, among other roles, they may serve as transcription factor binding sites (TFBS), such as in the lac operon or the TATA box.[@hannenhalli_eukaryotic_2008]. Finding these biologic motifs is useful for categorizing genes and proteins, and assists in furthering our understanding of how to regulate gene expression and protein formation. Research into motif discovery aids in the development of new treatments and cures for a host of human diseases. 

### Representation
Motif models are most commonly represented as a "point frequency matrix" (PFM) [cite pfm]. In PFMs each row is labeled by a base nucleotide, in the case of DNA from the set (A, C, G, T), and the columns indicate the aligned sequences positions. Each cell in the matrix is then calculated as the number of bases of its row type at that position. For instance, given these generated sequences:

```
1 CCCCT TAACC 10
1 ATCCT TAGGT 10
1 GACCT AACTA 10
1 CGCCT ACAGA 10
1 AGACT ACAAG 10
1 GAACC ACAGA 10
```
We can generate the PFM:

```
   1 2 3 4 5 6 7 8 9 10
--------------------------
A: 2 2 2 0 0 4 3 4 1 3
C: 2 1 4 6 1 0 3 1 1 1
G: 2 2 0 0 0 0 0 1 3 1
T: 0 1 0 0 5 2 0 0 1 1
``` 
This representation has multiple benefits. By shrinking the number of data points we need to look at, we can, even by eye, start to see areas of interest. (this example only has 8 sequences, but it is common to work with hundreds of sequences thousands of bases long in the best cases when searching for motifs.) Note the prevalence of 0s and counts greater than 4 in the index range 3-8 of the above PFM. This lessening of variability in which bases occur at a position in the sequences is called "conservation" and is what biologists watch for when searching for motifs. Another important advantage of this representation is it exchanges the previous symbolic representation of the sequences as letters, for a numeric representation. This expands the types of algorithms we can apply to this data in our search for biologic motifs. 

By normalizing the PFM we can get the probabilities for each base at each position. And with those probabilities we can also generate a score for how well a new sequence matches to the motif. We score by the simple equation 

$$
S = \prod_{i=1}^{n}P(x_i)
$$

where $x$ is the new sequence and $P(x_i)$ is the probability of that base at that position. This results in new sequences that match well having a higher score and dissimilar sequences having very low scores. This is a powerful technique that is useful in many motif finding algorithms as a metric for improving their motif model.


### Visualization
To aid in visualizing this PFM, we can plot this PFM as sequence logo:

![Counts logo](./imgs/CountsLogo.png)

Just like the PFM in this image the height of each letter is determined by the number of times it show up in that position. But, notice that the graphic is visually busy, and it is difficult to differentiate between levels of importance for each position. However, earlier biologists have improved on this visualization by instead using Shannon's Information Theory and equations for Information Content and Expected Information [@schneider_sequence_1990]

**Information Content:** 
$$ I(x_i) = \log_2 \Bigg( \frac{1}{P(x_i)} \Bigg) $$  
**Expected Information:** 
$$ H = -\sum_{i = 1}^{n} P(x_i) \log_2(P(x_i)) $$  

These equations measure the level of conservation at each position individually by comparing the probabilities we would expect to see if the nucleotide bases were evenly distributed (1/4 per base) to what is actually seen in the data. There is far more we could discuss about these equations [@nalbantoglu_data_2009], but all that is important for our purposes is to know that information content is measured in "bits of information" and that a higher measure of information content indicates more conservation of nucleotides at that position. These equation allows us to better weigh each position by how well conserved it is, as seen in the below image.

![Bits Logo](./imgs/BitsLogo.png)

From this image we can see that positions with low variability, such as position 4 where of all the possibilities only the nucleotide Cytosine is represented, correspond to high visibility in the sequence logo. Whereas positions of high variability, such as position 1, disappear.

Sequence logos are the primary method of visualizing biologic motifs [@hung_motif_2017], because they effectively communicate that a motif is an idealized sequence, a compilation of instances. Logos show how each position is likely to behave, position 3 in the above example is likely to be a Cytosine, but about a third of the time it might be an Arginine. 

## Current Motif Discovery Methods

### Profile Analysis

So far we have explored what motifs are, but we have yet to explain how they are discovered. In the simplest cases this can be done mostly manually, in a process called profile analysis [@hannenhalli_eukaryotic_2008]. If searching for TFBSs, The biologist gathers samples of the upstream region of genes known to be co-regulated, and performs a local multi-alignment using BLAST [@ncbiresourcecoordinators_database_2017] or similar tool. The produced alignments can be filtered based on their information content scores, and PFMs are generated from the highest scoring results.

Profile analysis requires large amounts of prior knowledge and human judgment, which is inconsistent at the best of times. This analysis can also struggle with instances that are less than reasonably similar. This analysis is also slow, running at human speed, clearly more robust and batch process-like tools are needed. 

### Combinatoric algorithms

The combinatoric family of motif discovery algorithms uses enumeration to solve for the problem of how likely it is to find a specific pattern `p` of length `l` length with `m` mutations among sequences of length `n`. These methods are related to the brute force approach of checking every possible sub-sting in the data, but use various tricks to reduce runtime, for example equating this brute force example to the median string problem [@jones_introduction_2004]. These algorithms, such as WEEDER and PSMILE [@pavesi_algorithm_2001; @carvalho_efficient_2004], are exhaustive, meaning that they will find a globally optimal solution given the data. This does not mean that the globally optimal solution is an actual motif, that can only be determined by _in-vito_ testing. It is common practice, for all motif finding algorithms, to proceed with the top few results rather than any single finding, greatly increasing the chances of discovering biologically functional motifs. Combinatoric algorithms in particular are known to struggle with longer motifs with high variability, and must deal with an exponentially growing search-space, which is to say that they are extremely effective at finding shorter motifs, so long as the underlying sequence does not become too long. 

### Probabilistic Algorithms

The other is common grouping of motif finding algorithms, are the probabilistic algorithms, of which the expectation maximization algorithms take up a large portion [@das_survey_2007]. Expectation maximization follows a repetitive pattern. We start setting a length of motif `m` we would like to find. We then choose randomly a starting point in each sequence we have. Once we generate a PFM from the sub-sequences starting at these random points, we can score against every sub-string in the larger sequences. We take the best scoring sub-strings and update our starting points to these new positions, and update our PFM. We do this over and over until the scores stop improving. 

In pseudo-code:

```
decide length of motif we would like to find
choose random starting points in sequences
do:
	generate PFM from starting points + decided length
	score every sub-string in sequence against PFM
	from top scoring sub-strings:
		update starting points
until scores stop increasing
```

This process is remarkably effective. And with improvements to deciding where to initially set the starting points such as Gibb's sampling [@das_survey_2007] among others. It should be noted that these methods implement a large amount of randomness, meaning that these probabilistic algorithms are apt to find a local optimum in their search. Biologists often run these algorithms multiple times, to increase the chances of finding the global optimum, but this takes increasing time with the length of the sequences. Despite these drawbacks, probabilistic algorithms are especially good at finding longer motifs with high variability.


## The Limits of Current Algorithms

### Benchmark Evaluation

For all algorithmic motif discovery methods, there is no known way to determine the accuracy of their search results completely computationally. The only sure way to determine if a discovered motif is truly biologically functional is to perform the correct _in-vito_ experiments, which in the case of TFBS means CHIP-seq or CHIP-chip testing [@hannenhalli_eukaryotic_2008]. Benchmark datasets have been compiled [@tompa_assessing_2005] using the previously experimentally determined motifs in the Transfac and Jasper databases [@matys_transfac_2006;  @mathelier_a._fornes_o._arenillas_d.j._chen_c._denay_g._lee_j._shi_w._shyr_c._tan_g._worsley-hunt_r._et_al._jaspar_2015], but these benchmarks do have their own limitations [@sandve_improved_2007]. For one, it is challenging to use real sequences as a benchmark as there may be unknown true motifs that an algorithm might be penalized for finding. However, generated sequences are not any better, because it is still unknown how to modal an accurate background sequence. Attempts with base distributions or Markov models are still inadequate [@tompa_assessing_2005; @sandve_improved_2007]. As a general rule, it seems sensible to test against both real and generated sequences with the understanding that any metric obtained will have a large degree of uncertainty. 

### Current Limits

From various surveys and comparisons, it appears that current motif finding algorithms have an accuracy of 15-20% at the nucleotide level (making an accurate motif model), and a slightly higher accuracy of 30-35% at the binding site level (finding the correct locations of motif instances) [@das_survey_2007; @tompa_assessing_2005; @sandve_improved_2007; @simcha_limits_2012; @hu_limitations_2005]. These finding may be slightly lower than would occur in practice, as many of these measurements were looking only at the top result whereas in practice biologists would usually be using the top few results. Several of these examinations also specify the largest challenges facing motif discovery today [@simcha_limits_2012; @hu_limitations_2005], 

In no particular order these challenges include:

- **Robustness to noise:** As mentioned in regards to our lack of understanding of the background sequences surrounding these motifs, it is often difficult to separate the true motifs from spurious motifs that arise in the data from random chance. This also relates to the extreme difference in signal to noise ratio, trying to rind a 8-30bp motif in many sequences over 1000bp long.
- **Ability to handle to different sized motifs:** As seen in both the combinatoric and probabilistic algorithms, we often need to decide what length of motif we would like to search for prior to starting our search. This is a significant issue since we usually don't know how long the motif should be. This requires that we run the algorithm multiple times, increasing both the runtime and the number of results we must sift through at the end. 
- **Efficiency with increasingly large sequences:** A continual issue, most motif finding algorithms grow exponentially with an increase in the length of sequences, and with increasing affordability to sequence large sections of the genome [@heather_sequence_2016] algorithmic efficiency will be a pressing challenge.
- **Simplicity of use:** Many motif finding algorithms possess numerous tuning parameters. These parameters can greatly increase effectiveness, but require an experts knowledge to actually use. This actually becomes a real problem in testing and comparing these algorithms [@tompa_assessing_2005; @sandve_improved_2007; @simcha_limits_2012; @hu_limitations_2005] as no individual or small team can be an expert in all of the available algorithms. Simplifying existing algorithms or developing a new simplified algorithm could reduce some of the uncertainty in algorithm comparisons as well as assist with motif finding in practice.

Solutions or improvements to any of these issues would count as a success to motif discovery research. 

