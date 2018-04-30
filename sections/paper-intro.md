# Introduction

This section describes biologic motif's common representations, including point frequency matrices and motif logos. It details the current state of research into motif finding algorithms (MFAs), touching on profile analysis, combinatoric algorithms, and probabilistic algorithms. It will reference the current challenges faced in the field of motif discovery [@das_survey_2007; @tompa_assessing_2005; @sandve_improved_2007; @simcha_limits_2012; @hu_limitations_2005]. And finishes by introducing the recently developed UniDip algorithm that has the potential to address many of these challenges.

## Biologic Motifs

### Representation

<!-- TODO Break up sentences -->
Biological motifs are short conserved patterns of similar or identical sequence instances as found in DNA. An example motif is the ALX3 transcription factor binding site (TFBS), MA0634.1. By controlling the expression of key genes, MA0634.1 serves a role regulating skeletal and mesoderm morphogenesis [@mathelier_a._fornes_o._arenillas_d.j._chen_c._denay_g._lee_j._shi_w._shyr_c._tan_g._worsley-hunt_r._et_al._jaspar_2015]. The ALX3 TFBS is a good example of a TFBS motif. Falling in the range of 6-30bp long, the TFBS shows has nearly 8000 reviewed samples, providing extensive supporting evidence of its role. 

![](./imgs/MA0634.1.svg)

The best visualization of MA0634.1, is with a motif logo, as seen above. At the core of this image is a numeric representation known as a "point frequency matrix" (PFM) [@hung_motif_2017]. In PFMs each row is labeled by a base nucleotide, and the columns indicate the tally of bases at the column's position in the sequences. 

The PFM for MA0634.1:

```PFM
A [ 1251   987   794  7877  7877    76   697  7877  2597  2759 ]
C [ 2174  3402  3367   506   317    72   729   449  1671  2446 ]
G [ 1062  1150   345   885   141     8   239    87  2205  1503 ]
T [ 3390  2337  4510   629   231  7877  7877   977  1403  1168 ]
```

A basic count representation shows all positions as equally important. However, in DNA, some positions do hold more information, i.e. containing one nucleotide far more often than any other. This departure from a uniform distribution of bases at any position is called conservation and motif logos are able to visualize conservation levels with the use of entropy and Shannon's information content [@schneider_sequence_1990]. 

Using Shannon's Information Theory and equations for Information Content and Expected Information, 

**Information Content:**
$$I(x_i) = -\log_2 \Big(P(x_i) \Big)$$

**Expected Information:**
$$H = \sum_{i = 1}^{n} P(x_i) I(x_i)$$

we can measure the level of conservation at each position individually by comparing the probabilities we would expect to an evenly distributed baseline ($1/4$ per base). These equations allow us to better weigh each position by how well conserved it is. We see in the logo above that positions with higher conservation can be visualized larger, making it much easier to locate the start and end sites of the motif at 3 and 8 respectively.

It should be noted that by default entropy and information content measure disorder, with the highest scores corresponding to when each possible options is equally represented. Logos show the increase of conservation by subtracting the entropy measurement from the max possible level of uncertainty. In genomic data with 4 nucleotides the maximum amount of entropy is 2 bits. Thus, for each base $b$, the level of conservation is calculated by the expression $2 - (H(b) + e(n))$ where $e(n)$ is a correction factor required when one only has a few sample sequences.

Sequence logos are the primary method of visualizing biologic motifs [@hung_motif_2017], because they effectively communicate the idealized sequence. Motif instances may differ, but they will follow the probabilities from the motif model.

## Current Motif Discovery Methods

Generating motif logos and PFMs make a large assumption that the sequences will be aligned. MFAs cannot make this assumption because in genomic sequences the instances of a motif will rarely be perfectly aligned [@hannenhalli_eukaryotic_2008]. Biologists must then use a range of sampling and alignment of the raw sequences to find motifs. Current motif finding approaches can generally be split into three branches as listed below:

**Profile Analysis**  
TFBS's are generally searched for by selecting a 1000bp upstream sequence from co-regulated genes. These may be different genes all within the same individual, or the same gene across species. For individual cases, analysis can be performed mostly manually in a process called profile analysis [@hannenhalli_eukaryotic_2008]. The biologist performs a multi-sequence local alignment using BLAST [@ncbiresourcecoordinators_database_2017] or similar tools. The produced alignments can be filtered based on their information content scores, and PFMs are generated from the highest scoring results. Profile analysis requires human observation at nearly every stage of the process, thus for large batches of cases, more self-contained MFAs are used. 

**Combinatoric algorithms**  
The combinatoric family of motif discovery algorithms uses enumeration to solve for the problem of how likely it is to find a specific pattern $p$ of length $l$ with $m$ mutations among sequences $n$ long. These methods are related to the brute force approach of checking every possible sub-string, but use various tricks to reduce runtime. For example, many equate this brute force example to the median string problem [@jones_introduction_2004]. These algorithms, such as WEEDER and PSMILE [@pavesi_algorithm_2001; @carvalho_efficient_2004], will find a globally optimal solution given the data. But, Combinatoric algorithms in particular are also known to struggle with longer motifs with high variability, and must deal with an exponentially growing search-space.

**Probabilistic Algorithms**  
The other is common grouping of MFAs, are the probabilistic algorithms, of which the expectation maximization algorithms are the most common [@das_survey_2007]. In pseudo-code we can see the general steps of the expectation maximization algorithm

```pseudocode
decide length of motif we would like to find
choose random starting points in sequences
do:
    generate PFM from starting points + decided length
    score every sub-string in sequence against PFM
    from top scoring sub-strings:
        update starting points
until scores stop increasing
```

This process is remarkably effective, particularly with improvements to deciding where to initially set the starting points such as with Gibb's sampling [@das_survey_2007]. It should be noted that these methods implement a large amount of randomness, meaning that these probabilistic algorithms are apt to find a local optimum in their search. Biologists often run these algorithms multiple times, to increase the chances of finding the global optimum, but this takes increasing time with the length of the sequences. Despite these drawbacks, probabilistic algorithms are especially good at finding longer motifs with high variability.

By working on the raw sequences these three branches encounter certain common limitations in terms of memory use and speed. Various studies have benchmarked and tested the current MFAs and specified several common challenges [@das_survey_2007; @tompa_assessing_2005; @sandve_improved_2007; @simcha_limits_2012; @hu_limitations_2005].

These challenges include:

**Robustness to noise:**  
As mentioned in regards to our lack of understanding of the background sequences surrounding these motifs, it is often difficult to separate the true motifs from spurious motifs that arise in the data from random chance. This also relates to the extreme difference in signal to noise ratio, trying to rind a 6-30bp motif in sequences over 1000bp long.

**Ability to handle to different sized motifs:**  
As seen in both the combinatoric and probabilistic algorithms, we often need to decide what length of motif we would like to search for prior to starting our search. This is a significant issue since we usually don't know how long the motif should be. This requires that we run the algorithm multiple times, increasing both the runtime and the number of results we must sift through at the end.

**Efficiency with increasingly large sequences:**  
A continual issue, most MFAs grow exponentially with an increase in the length of sequences, and with the lowering cost to sequence large sections of the genome [@heather_sequence_2016] algorithmic efficiency is a pressing challenge.

## Potential Application of UniDip

Some of these challenges can be mitigated by combining algorithms in ensemble methods. Different types of MFAs, are able to find different types of motifs [@das_survey_2007; @simcha_limits_2012]. Combinatoric algorithms can be better at finding shorter motifs, probabilistic algorithms can be better at longer motifs. This specialization of MFAs, means that ensemble methods can work very well [@simcha_limits_2012; @hu_limitations_2005]. But some other issues can be exacerbated -- in particular, the time needed to fit each individual model. What is needed is a way to reduce the time taken by existing algorithms. Time reduction can be achieved most flexibly by introducing data preprocessing steps. Improvements to hardware are expensive, and improvements to a core algorithm only affect that algorithm. Preprocessing steps can be performed on any machine, and the output is useful to many downstream algorithms.

One of the core difficulties of motif discovery is the search space involved. We are attempting to find a 6-30bp degenerate (non-exact) pattern in DNA samples 1000bp long or longer. Most MFAs asymptotically grow at an exponential rate [@das_survey_2007, @simcha_limits_2012; @hu_limitations_2005], meaning that just doubling the sample lengths to 2000bp presents a real challenge to motif discovery. Yet, the challenges of signal extraction from data with high noise to signal ratios are not unique to motif discovery. We should take advantage of this prior work and apply it to the problem of motif finding. The purposefully noise-robust UniDip algorithm could prove well suited to the motif discoveries particular set of challenges.

The UniDip algorithm is a core component of the SkinnyDip algorithm, developed by Maurus and Plant in 2016, and as Maurus and Plant describe "Practically, [SkinnyDip's] run-time grows linearly with the data" [@maurus_skinny-dip:_2016]. The big difference between UniDip and other clustering algorithms is that UniDip starts from an assumption that most of the data is noise. In centroid based methods all data is considered important, and thus all data points will be assigned to a cluster. Even density and expectation maximization based methods struggle with excluding noise in the data. UniDip and SkinnyDip are able to effectively isolate clusters, ignoring noise. In the case of motif discovery we will be able to isolate areas of high conservation in genomic sequences, trimming away large sections of the sample sequences that only clog the current MFAs.

UniDip is potentially powerful, but with these features also come some limitations. UniDip was developed to cluster continuous numeric distributions, data significantly removed from the symbolic letter-based format common to genomic datasets. We can transform the symbolic data to numeric data via aggregated representations taking inspiration from motif logos and PFMs. These representations make UniDip heavily reliant on sequence alignment when isolating clusters. But, the aggregations have the benefit of reducing the dimensionality of our data, greatly increasing speed. We also have an easier task, as we are not required to have exact motif instance extraction, we are only attempting to find the concentrations of conserved bases, such that we can run current MFAs, like MEME, on reduced data sets and still locate the embedded motifs. UniDip has the potential to make current MFAs scalable to the much larger sample sequences that are becoming more common with the reduced price of genomic sequencing.
