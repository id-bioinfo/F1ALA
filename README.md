# An ultrafast and memory-efficient lineage annotator for huge SARS-CoV-2 phylogeny using F1-score

The unprecedented scale of the global SARS-CoV-2 phylogeny overwhelmed most common tree annotation methods (such as matUtils and PastML) when adding the PANGO lineage information to the tree. Furthermore, the accuracy of these lineage annotation methods has not been clearly elucidated yet. To resolve these challenges, we develop an efficient annotation method using F1-score by transversing all potential annotated nodes, namely TIPars2. Our approach achieved a high classification accuracy with ultrafast speed (less than 13 minutes) with significant less memory usage (3.6 GB) for the annotation of 2277 Pango lineages in a phylogeny with 5.26 million taxa, allowing to complete the real-time lineage tracking in a laptop computer. We also demonstrate this efficient annotation could help to optimize the tree topology by pruning all inconsistent labeled taxa against their closest annotated ancestors and re-inserting them back. Overall, TIPars2 can provide real-time notable information of the transmission lineages along the global spread of SARS-CoV-2.

If there is any question about using TIPars2, please send email to tipars@d24h.hk.

# How It Works 

<img src="/img/illustration.png" width="600">
Given a tree with 5 taxa (Node5-9) and 4 internal nodes (Node1-4) where Node5-7 are labelled lineage A and Node8-9 are labelled lineage B, F1-score is computed in three steps. Step1: Extract ancestor nodes for any taxon of lineage A (Node1-3 and 5-7) and B (Node1, 3-4 and 8-9). Step2: Determine the order of lineages to assign the annotation according to their largest F1-scores, i.e., B = 1 and A = 4/5 (marked by underlines in the top two tables). Step3: Assign the annotation for B at Node4 first (middle table), then for A at Node1 because the F1-score recalculation excludes the counting of Node8 and Node9 (bottom table).

# Installation

A precompiled executable program is available as TIPars2.jar (required Java 11 or above). 

For users to compile TIPars2 source code from GitHub, 
```bash
git clone https://github.com/id-bioinfo/TIPars2.git
cd TIPars2
make
```

# Quick Usage

## Pango lineage annotation

```bash
cd /home/ytye/tipars2_github/Benchmark_datasets/100k
/home/ytye/tipars2_github/tipars2 --annotation -t 100k_tree_InnodeNameAdded.nwk --label 100k_pangolin.tsv  --output 1248_in_100k_annotation.tsv -T 8 
```

## Annotation statistics and visualization 

```bash
cd /home/ytye/tipars2_github/Benchmark_datasets/100k
/home/ytye/tipars2_github/tipars2 --annotation_details -t 100k_tree_InnodeNameAdded.nwk --label 100k_pangolin.tsv --assignment 1248_in_100k_annotation.tsv --output 1248_in_100k_annotation_details.tsv -T 8
```

+ Write the annotation details to the output file, including annotated_node, annotated_node_precedor, distance_to_root, pangolineage, F1score and samples.
+ Write the annotation visualization to file graph-data-generated.js that should be moved to the provided 'visual' folder and open the 'graph.html' in a browser.
+ Collapse the tree by lineages and write this collapsed tree to file [#lineages]_collapsedTree.nwk.
+ Remove inconsistent taxa and write this pruned tree to file [#consistent_taxa]_removedTree.nwk (used for tree refinement using other phylogenetic insertion methods, e.g., UShER).
+ Write the inconsistent taxa names and their lineages to file [#inconsistent_taxa]_unKeepSamples_[#consistent_taxa]_tree.tsv (used for tree refinement using other phylogenetic insertion methods, e.g., UShER).




# How to Cite



# Acknowledgements

This project is supported by the Hong Kong Research Grants Council General Research Fund (17150816), the NSFC Excellent Young Scientists Fund (Hong Kong and Macau) (31922087),
the Health and Medical Research Fund (COVID1903011-549 WP1) and the Innovation and Technology Commission’s InnoHK funding (D<sup>2</sup>4H).


