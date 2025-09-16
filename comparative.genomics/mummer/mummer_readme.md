# Description

MUMmer is a system for rapidly aligning large DNA sequences to one another. It can align:



whole genomes to other genomes

large genome assemblies to one another

partial (draft) genomes sequences to one another

or (with release 4) a set of reads to a genome.



MUMmer is very fast and easy to run. The current version, release 4.x, can find all 20-bp maximal exact matches between two bacterial genomes in just a few seconds on a typical desktop or laptop computer. MUMmer handles the 100s or 1000s of contigs from a draft genome with ease, and will align them to another set of contigs using the nucmer utility included with the system. The promer utility takes this a step further by generating alignments based upon the six-frame translations of both input sequences. Most of the improvements in release 4.x are to nucmer and MUMmer. Promer is unchanged from release 3.0. See the nucmer and promer readme files in the "docs/" subdirectory for more details.

MUMmer is free and open source, so all we ask is that you cite our most recent paper in any publications that use it. Here are the relevant publications:



##### Citation

Version 4.0+

MUMmer4 and nucmer4 are described in "MUMmer4: A fast and versatile genome alignment system" G. Marçais , A.L. Delcher, A.M. Phillippy, R. Coston, S.L. Salzberg, A. Zimin, PLoS computational biology (2018), 14(1): e1005944.

