
# Abstract

## Introduction and Background

In the _Introduction_ we will define what we mean by biologic motifs, and discuss their role in biology, including as transcription factor binding sites, as well as their common representations, such as point frequency matrices. We will continue on to describe the current state of research into motif finding algorithms, touching on profile analysis, combinatoric algorithms, and probabilistic algorithms. And finish by referencing current challenges faced in the field of motif discovery as reported by numerous comparison projects and survey papers [@das_survey_2007; @tompa_assessing_2005; @sandve_improved_2007; @simcha_limits_2012; @hu_limitations_2005]. These challenges include the complexity of finding the right tuning parameters, required repetition to find motifs of different sizes, lack of efficiency in the face of expanding sequence lengths, and sensitivity to noise.

Starting in the _Background_ section we will attempt to address some of these issues by proposing to apply the SkinnyDip algorithm and its derivatives as developed by [@maurus_skinny-dip:_2016]. We will explain SkinnyDip's function as well as its underlying components, and detail what steps we will need to take to apply this family of algorithms to the problem of motif discovery.

## Methodology

In this section we will describe our process of applying the SkinnyDip family of algorithms to the problem of motif discovery. Starting with the ability to find peaks from the bin heights of a histogram and how that will allow us to cluster symbolic data, we will incrementally apply a modified SkinnyDip to increasingly complex data. 

> TODO

## Evaluation

> TODO
