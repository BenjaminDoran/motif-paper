# Methodology

We begin applying unimodality clustering to motif discovery with a discussion on implementing the UniDip algorithm, because algorithms in practice are never as simple as they are in theory.

The original SkinnyDip algorithm was implemented in R, We decided to port this example to Python for because other primary bioinformatics tools [cite BioPython] are already implemented in Python, We would need potentially be needing to make changes to the algorithm, and the original source code included extra algorithms we would not need for our tests. Overall porting this software went 

P-vals

mirroring

> ...Need to Write...

## Applying to Histograms

> ...Need to Write...

Sorted

modified mirroring

losing detail

## Applying to Symbolic Data

> ...Need to Write...

simplest case insert mono-symbolic exact motif into uniformly distributed background noise

counts to scaled negative entropy.

## Applying to Multi-Letter Motifs

> ...Need to Write...

scaled negative entropy works in this case as well

## Applying to Degenerate Data

> TODO

adding in randomness to the motif (degeneracy)

## Single Motifs

use mirroring

## Applying to Differently Aligned Motifs

> TODO

This is the biggest limitation

It can be moderately alleviated by using a global alignment. But more extensive workarounds are beyond the scale of this project.

## Filtering Out Spurious Motifs

> TODO
	
Ranking

## Summary