**hc-OTU: A Fast and Accurate Method for Clustering Operational Taxonomic Units based on Homopolymer Compaction**

Seunghyun Park, Hyun-Soo Choi, Byunghan Lee, Jongsik Chun, Joong-Ho Won, and Sungroh Yoon, IEEE/ACM Transactions on Computational Biology and Bioinformatics, vol. 15, no. 2, pp. 441-451, March/April 2018.

[paper link](https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=7420667&tag=1)

Abtract - To assess the genetic diversity of an environmental sample in metagenomics studies, the amplicon sequences of 16s rRNA genes need to be clustered into operational taxonomic units (OTUs). Many existing tools for OTU clustering trade off between accuracy and computational efficiency. We propose a novel OTU clustering algorithm, hc-OTU, which achieves high accuracy and fast runtime by exploiting homopolymer compaction and k-mer profiling to significantly reduce the computing time for pairwise distances of amplicon sequences. We compare the proposed method with other widely used methods, including UCLUST, CD-HIT, MOTHUR, ESPRIT, ESPRIT-TREE, and CLUSTOM, comprehensively, using nine different experimental datasets and many evaluation metrics, such as normalized mutual information, adjusted Rand index, measure of concordance, and F-score. Our evaluation reveals that the proposed method achieves a level of accuracy comparable to the respective accuracy levels of MOTHUR and ESPRIT-TREE, two widely used OTU clustering methods, while delivering orders-of-magnitude speedups.


- C++
- OpenMP

```
hc-OTU assigns similar sequence reads to operational taxonomic units (OTUs) by clustering 
sequences based on a user-defined similarity threshold. Sequences which are similar at or 
above the threshold level are taken to represent the presence of a taxonomic unit (e.g., 
a species, when the similarity threshold is set at 0.97) in the sequence collection.
ver.1.0 Jan. 21, 2018
Usage : hc-OTU [file name] [similarity cut off ex. 0.97] [number of threads]
DS&AIL and ChunLab, Inc. in Seoul National University (c) 2009-2018
```
