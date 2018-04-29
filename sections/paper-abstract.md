# Abstract

## Motivation

Common motif finding algorithms work directly on raw sequences. This focus has advantages and disadvantages. Working with raw sequences does give more detail and is more grounded to the core of genomic data. But it does limit the number of available algorithms. Current motif finding algorithms struggle with speed, underlying randomness, and high noise environments. The UniDip algorithm, developed outside the field of biology, is fast, deterministic, and noise robust. If applied well, the UniDip algorithm will be a powerful addition to sequence analysis.

## Results

We present a method to measure the conservation level of aligned sequences providing a numerical representation accessible to the UniDip algorithm. This metric is based on Shannon's information content and entropy formulas as popularized though motif logos. We show that UniDip is able to take this numeric representation to isolate the regions of high conservation in simulated sequences, working on degenerate motifs with up to $50%$ mutation. We also show a case study isolating finding the transcription factor binding site to FOXK1. UniDip serves as a powerful processing tool that is able to trim out low conservation regions, shrinking the search space for conventional motif finding algorithms. With MEME, We are able to find the FOXK1 transcription factor binding site in $1/5$ the time using our preprocessing steps versus running MEME directly on the raw sequences.

