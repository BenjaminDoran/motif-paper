# The UniDip Algorithm

The UniDip algorithm is able to find clusters in univariate continuous numeric distributions. It can perform this feat by layering itself over the Hartigan Test of Unimodality [@hartigan_dip_1985]. This formal statistical test is not terribly well known outside of its statistical circles, but does provide some useful features. These features include linear time complexity and non-parametric compatibility to any unimodal distribution e.g. Gaussian, Uniform, or Beta distributions. SkinnyDip is a layer on top of UniDip that recursively runs UniDip across each dimension.

## Hartigan's Dip Test of Unimodality

Hartigan's dip-test makes use of a distribution's empirical cumulative distribution function (ECDF). As can be seen from the plot below, this function's gradient increases approaching a peak in the histogram, and decreases after. In unimodal data, this creates a stretched S shape.

![](./imgs/1PeakWithECDF.png)

![](./imgs/3PeakWithECDF.png)

The dip-test performs a best fit to this shape, finding an minimal width path such that, from left to right on the `x` axis, the gradient only increases until it reaches a point where-after the gradient only decreases. The dip statistic is defined as the width of this path divided by two, and does not vary with shifting or scaling.

Upon return of this dip statistic we may compare against a suitable unimodal null distribution, Hartigan suggests that, to obtain a p-value, the Uniform distribution is preferred [@hartigan_dip_1985]. If the p-value is greater than `alpha` we accept the null hypothesis that the distribution is unimodal, otherwise we accept the alternative that the distribution is _at least_ bi-modal.

The dip statistic and p-value are not enough in themselves to help us locate peaks. However, the dip-test also provides a modal interval, which specifies a lower and upper index to the isolated cluster. From this four pieces of information we can begin to search a univariate dataset.

## UniDip, Recursive Application of the Dip Test

UniDip takes us the next step by allowing us to recursively search a data sample to find peaks of density. We start by dipping along the entire single dimensional set of data. If the data is unimodal than we should return the set of data as the modal interval. otherwise we should perform UniDip within our located modal interval. We then try to take the same steps to the left and right of the modal interval. Since the dip-test can only differentiate between unimodal and multi-modal distributions, to determine if there are indeed peaks of interest to the left or right we need to include the left-most or right-most detected peak. By including left-most or right-most peak, the dip-test will only show evidence for multiple peaks if there are indeed additional peaks to the left or right.

UniDip in condensed pseudo-code:

```pseudocode
INPUT: 
    X, (1d sorted vector)
    alpha, (significance level)
    isMod, (is modal interval always true at start)

OUTPUT: set of modal intervals
    where a modal interval := (lower index, upper index)

UniDip(X, alpha=0.05, isMod=True)
    dip, pval, li, ui = DipTest(X)

    if pval > alpha
        return (if isMod) ? (X[0], X[-1]) : (li, ui)

    // recurse into interval
    Mm = UniDip(X[li, ui], alpha, True)

    // find left and right most intervals
    U = min(Mm, key=>(t) t[-1]); L = max(li, key=>(t) t[0])

    // check if left side is at least bi-modal when including
    // the left most mode then do same check but to the right
    pL = DipTest(X[X <= U]); pU = DipTest(X[X >= L])

    // recurse left if at least bi-modal
    Ml = (if pL <= alpha) ? UniDip(X[X < li], alpha, False) : ()

    // recurse right if at least bi-modal
    Mu = (if pU <= alpha) ? UniDip(X[X > ui], alpha, False) : ()

    return Mm & Ml & Mu
```

At the end of running this algorithm we collect the union of the all our recursive steps. We may need to merge any touching intervals, but none will overlap.

## SkinnyDip Recursive Application of UniDip

The multidimensional variant of UniDip, SkinnyDip, is beyond the scope of this project. But, is briefly described as another layer on top of UniDip.

```pseudocode
for each dimension
     get UniDip intervals
     within each interval
         get SkinnyDip hyperinterval from subsequent dimension
return hyperintervals
```
We recursively run UniDip on each dimension, and, in the event that we find clusters, we check the subsequent dimensions to determine intersections. In 2 dimensions our intersections or "hyper-intervals" will be squares, in 3 dimensions they will be cubes. 

While motif discovery starts with data in the form of many 1000bp long genomic sequences, our aim is to measure the level of nucleotide conservation, an inherently univariate metric. In the motif logos and PFMs, We have seen representations such as frequency and information content measure conservation as a function of each nucleotide. Distilling those representations down to a single function of position will allow us to isolate regions of conservation with the univariate UniDip algorithm.

## Implementing the Basic UniDip Algorithm

Before attempting to introduce any symbolic letter-based data we implemented our version of the UniDip algorithm in Python to the same specifications as the original UniDip source code [@samhelmholtz_skinny-dip:_2017]. 
The only real challenge we faced in implementing this algorithm is that there are very few libraries for calculating the Hartigan dip statistic. We settled on an implementation by Johannes Bauer [@bauer_dip_test:_2018] because its test suite proves it gives the same results, to within 8 decimal places, as those of the R package used by the original UniDip algorithm.

Bauer's dip implementation is slightly limited in that while it does calculate the dip statistic, it does not calculate a p-value for that statistic. To run a dip test, we calculate the dip in our data and then compare against a null distribution of dip statistics generated from a random samples of the standard uniform distribution as recommended by Hartigan [@hartigan_dip_1985]. Since, for the most part we would like to find mode clusters with a significance "alpha" of 0.05, we can run a dip test to this significance with `N = 1000` samples. While our approach is effective for this prototype stage, it is not optimal for two reasons. One, we are introducing randomness into an otherwise deterministic algorithm. And two, generating 1000 samples every time we need a p-value is a drain of computational resources. Martin Mächler has already studied these issues and found that p-value estimates can be interpolated from generated tables for much faster and deterministic performance [@machler_dip_nodate]. If we move forward with this algorithm, one necessary step will be implementing this type of p-value look up table. But, currently basic testing against the uniform distribution is sufficient.

After implementing the basic UniDip algorithm and Hartigan dip test we are able to isolate peaks in univariate numeric samples. We successfully tested our implementation with concatenated random samples from normal distributions, with 1, 2, 5, and 10 peaks among 80% noise. 

![](./imgs/plots-from-random-normal.png)

As we can see from this image, UniDip is able to isolate the regions of higher density. Returning the lowest and highest index of the clusters from the sorted dataset.
