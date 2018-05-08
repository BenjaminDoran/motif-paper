# Methodology

If we would like to apply UniDip to genomic data we will need to make modifications, both by transforming the data and the algorithm itself. In this section we will detail what steps we took to apply the SkinnyDip family of algorithms to symbolic data.

## Applying to Histogram Data

To start, we will make some changes to the UniDip algorithm so that it can work with more discrete forms of data. The core of UniDip, Hartigan's dip-test, only requires access to the ECDF of the data. Thus, any data must be reducible to an ECDF. An ECDF can be calculated either directly from the sample of a continuous random variable, or an approximation can be made from the Y coordinates of the bins in a histogram.

*Eq. 3: Calculate ECDF from continuous random variable as X*
$$\text{ECDF } = X_1 ... X_n / \max(X)\text{; Where }X_i \leq X_{i+1}$$

*Eq. 4: Approximate ECDF From histogram bin heights as H:*
$$\text{ECDF } \approx \frac{H_k}{n} \in H \text{;  Where:} H_k = \sum_{i=1}^{k}H_i$$

In case of a random continuous variable, grouping points inevitably loses detail, but increasing the number of bins reduces this loss. 

![Fig. 4: Histogram and ECDF calculated from random continuous variable and histogram bin heights](./imgs/histogramConversion.png)

This image plots a random sample of 500 points from a normal distribution, and shows that the cumulative sum of the histogram bin heights is a reasonable estimate of the ECDF.

The reason it is useful to generate an ECDF from histogram heights, is because thinking back to the sequence logos, there are a lot of similarities to a histogram. Each position can be compared to a bin, with the height determined by the frequency or level of conservation at that position.

To apply UniDip to histogram data; however, we must change a core assumption of the algorithm. Looking back to *Eq. 3*, we see that the ECDF is calculated on a vector sorted lowest to highest by value. This means that, by default, UniDip sorts its inputted data. This works with a random sample where order is only determined by value and not position. But, when working with histogram data, sorting by value would destroy the ordering that indicates the density at each position.

This problem with sorting springs up in even more obscure parts of the algorithm. When presented with unimodal data Hartigan's dip-test will tend to return an extremely narrow model interval. Often this occurs when recursing into modal intervals, but in that case UniDip can simply return the end points of the parent modal interval as it is significantly certain that all of the data points are a part of the same peak.

*Code 4: UniDip returning modal intervals*
```pseudocode
if unimodal
    return (if isMod) ? (X[0], X[-1]) : (Mod[0], Mod[-1])
```

However, an issue becomes apparent when recursing to the left or right of the modal interval. If there is only one modal peak, following the pseudo-code from above, UniDip will return an interval cutting off a large portion of the peak's tails. The solution that Maurus and Plant developed is to "mirror" the data, such that if the data is `[1, 2, 3]` the mirrored dataset is `[-2, -1, 0, 1, 2]`. What this means is that when UniDip is recursing to the left or right and it has found a single peak, UniDip will perform the dip test on a bimodal mirror-set of the data to extract the whole peak. 

Data mirroring becomes a problem when applied to histogram data because the specific mirroring algorithm Maurus and Plant use sorts by value. To allow histogram data, we replaced this with a mirror function that flips the data by index rather than value. 

*Code 5: Mirroring data*
```pseudocode
if flip_left:
    mirror_data = concatenate((flip(X[1:]), X))
else:
    mirror_data = concatenate((X[:-1], flip(X)))
```

Making these modifications we are able to apply the UniDip algorithm to histogram data and isolate the same clusters as we would using the raw random samples.

![Fig. 5: UniDip clustering on continuous and histogram data](./imgs/histogramClustered.png)

## Applying to Symbolic Data

Taking the next step, we can apply UniDip to symbolic data. We will start on the simplest data we can generate. Our motif will be a 15bp long poly-A sequence with no mutations. We will surround this sequence on each side with a 100bp background sequence of 1/4 uniform sampling of all four nucleotides. We will generate 20 sequences of this type and visualize them with point frequency and individual information content.

![Fig. 6: Aggregated numeric representations of genomic data](./imgs/polyAMetricsSingle.png)

The first item to notice about both these metrics is that they are calculated individually by base. There are separate measurements for Adenine, Cytosine, Guanine, and Thymine. This can be important later on in representing motifs as it is important to know if position `i` is always "A" or perhaps can be either "A" or "C". At this stage though, it is more important to measure overall conservation so that UniDip can isolate the motif regardless of its specific pattern. Later analysis can uncover the individual nucleotide differences.

An overall metric is possible using expected information content, also called entropy. Entropy can be calculated at each position as the sum of each base's information weighted by its percentage. 

*Eq 5: Entropy by position*
$$H= -\sum_{i=a}^{t} P(x_i) \log_2(P(x_i))$$ 

Applying an entropy calculation column-wise along the sequence will allows measurement of general conservation non-specific to any particular nucleotide. But, entropy by itself is not well suited to UniDip. Looking at the individual information content in *Fig. 6*, we see that higher scores actually correspond to less frequency. This makes sense because a sequence that is always a single symbol cannot convey any information. 

However, to measure conservation, we would like to invert the Y axis and scale to above 0. This transformation is important for two reasons. First, inverting and scaling makes thinking about conservation easier, making higher scores indicate higher conservation. Second, inverting and scaling makes the metric conform to the format of histogram data that we have already shown interfaces well with UniDip. Following the example set by motif logos, we invert the Y axis and scale in one step by subtracting the positional entropy from the maximum score possible with four possible nucleotides $2-H(X)$. After this transformation, we can visualize the simple sequence as scaled negative entropy (SNE).

![Fig 7: Scaled Negative Entropy (SNE)](./imgs/polyANegEnt.png)

We can see the  metric shows a high peak in the region of the poly-A motif, indicating its conservation. The surrounding regions vary with each position, but are nowhere near the height of the generated motif.

By default, with unimodal data, UniDip would return all the data as a single cluster. We can make another change to the algorithm, so that, in the event that we return a single cluster, we subsequently perform data mirroring to isolate the single region of highest conservation.

![Fig. 8: UniDip clustering on symbolic data](./imgs/polyACluster.png)

In this example, UniDip returns the interval spanning positions 100-116. The actual motif instances span positions 100-115, meaning UniDip returned an interval one position too wide. This is still a good match, and for such simple examples with large difference between conserved and non-conserved regions, UniDip could be used on its own for motif discovery. Real genomic sequences are much more complex, so this result is not indicative of real world performance. However, this example does show that we are able to isolate a region of high conservation in symbolic data.

### Applying to Multi-Symbol Motifs

Because SNE is an overall metric, a function of position not nucleotide, UniDip works with motifs that are not all the same base. Our motif will be the 15bp long sequence `ACTGTGCACGTGACG` with no mutations. Other parameters to the background sequence remain the same. Visualizing with SNE, it is hard to see much of a difference.

![Fig. 9: SNE with multi-symbol motif](./imgs/MultiLet02.png)

It is only when looking at the individual metrics that the difference is clear.

![Fig. 10: Individual metrics with multi-symbol motif](./imgs/MultiLet01.png)

This consistency of SNE shows why it is important to have this grouped univariate metric. While the conservation of each nucleotide varies, the univariate metric, SNE, shows clearly where the motif is positioned. Utilizing univariate measures of conservation, UniDip can detect multi-symbol motifs in genomic data.

![Fig. 11: Univariate metric (SNE) allows UniDip clustering](./imgs/MultiLet03.png)

### Applying to Degenerate Data

Introducing mutations does increase the difficulty in finding the regions of conservation because conservation is lowered closer to the background sequence. We introduced 6 mutations into each instance of the motif across the 20 samples, with few effects.

![Fig. 12: SNE with 6 mutations per motif instance](./imgs/degenerateSeq03.png)

But higher levels of degeneracy do have an effect making UniDip return a wider interval. In this case we increase degeneracy to 10 potential mutations per instance.  

![Fig. 13: SNE with 10 mutations per motif instance](./imgs/degenerateSeq04.png)

This level of degeneracy is high for a motif, but a similar effect will be seen where the background sequence becomes more conserved. With less of a threshold between the signal and background, UniDip is apt to lose specificity. UniDip is still useful even when it is unable to exactly match motif sites. Note that the cluster above still contains the motif instances and has by exclusion marked out a large portion of the sequence to rule out. Despite the lessened difference in signal to noise, UniDip has still narrowed the search space considerably.

## Applying to Differently Aligned Motifs

Due to how SNE is calculated, UniDip is heavily reliant on the alignment of sequences to be able to measure nucleotide conservation. Even a misalignment of few nucleotides can obfuscate the entropy calculations. Thus, misaligning motifs is the largest challenge in applying UniDip to motif discovery. For reference, compare the below SNE plots that show data sets of 20 sequences. In one we have added a random +/-5bp misalignment and the other has perfect alignment.

![](./imgs/maskedMotifMisaligned.png)
![](./imgs/maskedMotifAligned.png)

*Fig. 14: Misaligned (top) and aligned (bottom) motif instances*

If the amount of misalignment is minor, this problem with misalignment could possibly be alleviated by increasing the amount of data. Once we have enough sequences, even placed randomly, motif instances will overlap such that we can detect the conservation. For 15bp motifs with mis-alignments of +/-5 adjacent positions we could be assured of a perfect alignment with 10 or more samples. But, this also means that the conservation would grow not just at a single motif site but at all overlapped motif sites. For the 15bp motif, conservation would increase in a 25bp region. 

However, in real sequences there are very few assumptions on where the motif instances will fall [@hannenhalli_eukaryotic_2008]. The general assumption for TFBS's is that, both with single-species co-regulated genes and across species orthologs, the TFBS may fall anywhere within 1000bp upstream of the transcription start site. With misalignments of 1000bp, we would only expect there to be a perfect alignment by chance after collecting 500 sequences. This requirement for data samples is multiple orders of magnitude larger than other motif finding algorithms.

Global alignment tools can reduce this number of required samples. The MUSCLE alignment tool [@edgar_muscle:_2004], performs a multi-sequence alignment minimizing the number of mutations, insertions, and deletions as much as possible. Other alignment tools exist and an in-depth comparison of their merits is warranted for further research, but we will only be using MUSCLE for this project. Global alignment tools move motif instances to overlap more and boost conservation in select regions. Of course, performing the alignment does introduce its own issues. First, is by introducing gaps where previously sequences were contiguous. Second, because MUSCLE is actively forcing the sequences to align better, SNE is no longer directly measuring sequence conservation. However, a few adjustments will handle these issues. 

Regarding the introduction of gaps, complications arise from needing to handle a new symbol, "-", indicating an insertion. It is not immediately clear how to count this insertion. If "-" is counted like any letter there will be large sections of perfect conservation, which will mask the signals of true conservation. However, there are equally large problems if "-" is automatically set to $0$. Remembering that entropy calculations divide by the counts of each base, having sections of $0$ could lead to division errors. Instead, for any insertion we add a tally to all other bases. Consequently, positions with many insertions have their conservation level lowered as the counts of bases at that position become more uniform. 

So, is this tallying method enough to be able to effectively handle misalignments? Generating 20 sample sequences with a motif at a random misalignment of +/-10bp, we see a high peak where the motifs have been aligned. However, after running UniDip on this data, UniDip is now clustering not just the motif but all alignments as well. This is because there are now three levels of conservation: the gaps, the background sequence, and the motif. UniDip is unable to find nested clusters, where there are multiple steps of density. 

![Fig. 15: UniDip clustering on aligned sequences](./imgs/GappedAlignedClustered.png)

Trimming the gaps from the data overcomes this limitation. The gaps are at a relatively even level less than 0.1 SNE. By filtering and concatenating only regions that are greater than that threshold, UniDip is able to correctly isolate just the conserved motif once again.

![Fig. 16: UniDip clustering on aligned and trimmed sequences](./imgs/UnGappedAlignedClustered.png)

Removing these gaps might present a problem for being able to map back to the original sequence instances, but by keeping track of the aligned indices while performing the filtering we are able to maintain the ability to map back to the original instances.

## Methodology Summary

In this section we have shown the methods that make it possible to isolate regions of conservation from symbolic genomic data using the UniDip algorithm. We have found that scaled negative entropy is easily utilized by UniDip and well represents positional conservation. SNE mitigates the challenges of multi-symbol degenerate motifs very well, though it does falter when facing misaligned motif instances. To mitigate this factor, we used the global alignment tool MUSCLE, which allows us to concentrate motif instances in select regions. And, by trimming gaps, we are able to isolate regions of conservation from the background sequence.