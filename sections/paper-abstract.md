# Abstract

## Motivation

Common motif finding algorithms work directly on raw sequences. This focus has advantages and disadvantages. Working with raw sequences does give more detail and is more grounded to the core of genomic data, but it limits the number of available algorithms. Current motif finding algorithms struggle with speed, underlying randomness, and high noise environments. The UniDip algorithm, developed outside the field of biology, is fast, deterministic, and noise robust. Levering these strengths, UniDip will be a powerful addition to sequence analysis of symbolic genomic data.

## Results

Inspired by the representation of biologic motifs in motif logos, we present a method to measure the conservation level of aligned sequences providing a numerical representation accessible to the UniDip algorithm. This metric is based on Shannon's information content and entropy formulas. We show that UniDip is able to take this numeric representation to isolate the regions of high conservation in simulated sequences, working on degenerate motifs with up to 55% mutation. We also show a case study isolating the transcription factor binding site of FOXK1. UniDip serves as a powerful processing tool that is able to trim out low conservation regions, shrinking the search space for conventional motif finding algorithms. With MEME, we are able to find the FOXK1 transcription factor binding site 70% faster preprocessing with UniDip versus running MEME directly on raw sequences.

