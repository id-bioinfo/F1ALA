# F1ALA: ultrafast and memory-efficient Ancestral Lineage Annotation for huge SARS-CoV-2 phylogeny using F1-score

The unprecedented scale of the global SARS-CoV-2 phylogeny overwhelmed most common ancestral lineage annotation methods (such as matUtils and PastML) when annotating the PANGO lineage information to unlabeled nodes in a rooted tree. Furthermore, the accuracy of these annotation methods has not been clearly elucidated. To resolve these challenges, we developed an efficient and accurate ancestral lineage annotation method (F1ALA). It utilizes F1-score to evaluate the confidence for assigning a lineage annotation at a specific ancestor node given lineage labels of taxa in the tree. F1ALA achieved ultrafast speed (less than 13 minutes) with significant less memory usage (3.6 GB) for annotating 2,277 PANGO lineages in a phylogeny with 5.26 million taxa, allowing real-time lineage tracking to be performed on a laptop computer. Benchmarking on three phylogenies with 100K, 660K and 5.26M taxa, F1ALA significantly outperformed matUtils and was comparable to PastML on the annotation accuracy in both empirical and simulated tests. The high efficiency of F1ALA enables the refinement of a huge SARS-CoV-2 phylogeny by pruning all taxa with inconsistent label compared to their closest annotation nodes and re-inserting them back. We demonstrated that this refinement was able to optimize the SARS-CoV-2 phylogenetic tree topology achieving a larger tree log-likelihood and a smaller parsimony score.

If there is any question about using F1ALA, please send email to tipars@d24h.hk.

# How It Works 

<img src="/img/illustration.png" width="600">
Illustration of the algorithm for ancestral lineage annotation. Given a tree with 5 taxa (Node5-9) and 4 internal nodes (Node1-4) where Node5-7 are labeled lineage A and Node8-9 are labeled lineage B, ancestral lineage annotation is computed in three steps. Step1: Extract potential annotated ancestor nodes for lineage A (Node1-3 and 5-7) and B (Node1, 3-4 and 8-9) (shown in the headers (black background) of top two tables). Step2: Determine the order of lineages to assign the annotation based on the annotation confidence score (the largest F1-score in each lineage, i.e., A = 4/5 and B = 1, marked by underlines in top two tables). It is first to assign lineage B and then A. The order is shown as ① and ② in bottom two tables. Step3: Assign the annotation for B at Node4 firstly (middle table), then for A at Node1. When recalculating F1-scores for potential annotated ancestor nodes of lineage A, taxa Node8-9 are excluded due to taxa Node8-9 have been assigned to the confirmed annotation of Node4 in the previous processing (bottom table). Note that F1-score tables for lineage B are the same in Step1 and Step3 shown in the middle.

# Installation

**A precompiled executable program is available as F1ALA.jar (required Java 11 or above).**
```bash
git clone https://github.com/id-bioinfo/F1ALA.git
cd F1ALA
chmod a+x f1ala
# If users want to compile F1ALA from source code, 
make
```

For conda installation, 
```bash
conda create -n f1ala
conda activate f1ala
conda install f1ala::f1ala
```

# Quick Usage

## Ancestral lineage annotation
Infer the lineage information at the ancestor nodes in a given rooted tree with labeled taxa.

```bash
cd /home/ytye/f1ala_github/Benchmark_datasets/100k
/home/ytye/f1ala_github/f1ala --annotation -t 100k_tree_InnodeNameAdded.nwk --label 100k_pangolin.tsv  --output 1248_in_100k_annotation.tsv -T 8 
```

## Annotation statistics and visualization 

+ Write the annotation details to the output file, including annotation_node, annotation_node_precedor, distance_to_root, pangolineage, F1score and samples.
+ Write the annotation visualization to file graph-data-generated.js that should be moved to the provided 'visual' folder and open the 'graph.html' in a browser.
+ Collapse the tree by lineages and write this collapsed tree to file [#lineages]_collapsedTree.nwk.
+ Remove inconsistent taxa and write this pruned tree to file [#consistent_taxa]_removedTree.nwk (used for tree refinement using other phylogenetic insertion methods, e.g., UShER).
+ Write the inconsistent taxa names and their lineages to file [#inconsistent_taxa]_unKeepSamples_[#consistent_taxa]_tree.tsv.

```bash
cd /home/ytye/f1ala_github/Benchmark_datasets/100k
/home/ytye/f1ala_github/f1ala --annotation_details -t 100k_tree_InnodeNameAdded.nwk --label 100k_pangolin.tsv --assignment 1248_in_100k_annotation.tsv --output 1248_in_100k_annotation_details.tsv -T 8
```

## Tree refinement

Refine of a phylogeny by pruning all taxa with inconsistent label compared to their closest annotated ancestors and re-inserting them back using phylogenetic insertion methods such as TIPars and UShER.

### Phylogenetic insertion using TIPars
+ include processing ancestral lineage annotation
```bash
cd /home/ytye/f1ala_github/Benchmark_datasets/100k
/home/ytye/f1ala_github/f1ala --refinement -t 100k_tree_InnodeNameAdded.nwk -s 100k_taxa.fas -a 100k_anc.fas --label 100k_pangolin.tsv --output refined_tree.nwk -T 8 -x 8G
```

+  exclude processing ancestral lineage annotation that could be computed by other methods, e.g. PastML and matUtils
```bash
cd /home/ytye/f1ala_github/Benchmark_datasets/100k
/home/ytye/f1ala_github/f1ala --refinement_from_annotation -t 100k_tree_InnodeNameAdded.nwk -s 100k_taxa.fas -a 100k_anc.fas --label 100k_pangolin.tsv --assignment 1248_in_100k_annotation.tsv --output refined_tree.nwk -T 8 -x 8G
```

### Phylogenetic insertion using UShER 
+  use the [#consistent_taxa]_removedTree.nwk after ancestral lineage annotation by TIPars2
```bash
cd /home/ytye/f1ala_github/Benchmark_datasets/100k
/home/ytye/f1ala_github/f1ala --annotation -t 100k_tree_InnodeNameAdded.nwk --label 100k_pangolin.tsv  --output 1248_in_100k_annotation.tsv -T 8
/home/ytye/f1ala_github/f1ala --annotation_details -t 100k_tree_InnodeNameAdded.nwk --label 100k_pangolin.tsv --assignment 1248_in_100k_annotation.tsv --output 1248_in_100k_annotation_details.tsv -T 8
usher -v taxa.vcf -t 81784_removedTree.tree -d ./usher -o ./usher/81784_AddTo_100k.pb
```
+ The refined tree by UShER is ./usher/final-tree.nh.

## Tree bubbling (to be done)
Collapse the tree into multiple clusters based on the ancestal lineage annotation.
Large clusters (>exploreTreeNodeLimit) will further to stratify into multple bubbles by BFS search.
Small clusters (<smallClusterLimit) and bubbles (<smallBubbleLimit) will be merged.
Clusters will link to bubbles.

```bash
cd /home/ytye/f1ala_github/Benchmark_datasets/100k
/home/ytye/f1ala_github/f1ala --tree_BFS -t 100k_tree_InnodeNameAdded.nwk --label 100k_pangolin.tsv  --output 100k_tree_bfs.tsv --exploreTreeNodeLimit 2000 --smallBubbleLimit 5 --smallClusterLimit 5 -T 8 
```
+ Output tsv file includes 8 items.
1. bubble_type : 1 is cluster and 2 is bubble
2. annotation_node : root of the subtree for cluster or bubble
3. annotation_node_precedor : precedor of this annotation_node where precedor is also a cluster or bubble
4. dist_to_precedor : total branch length from annotation_node_precedor to annotated_node
5. parent_node : parent node of this annotation_node in the input tree
6. pangolineage : lineage label annotated by 'Ancestral lineage annotation'
7. num_nodes : number of nodes in the cluster or bubble
8. nodes : a list of nodes in the cluster or bubble (separated by comma)

# How to Cite



# Acknowledgements

This project is supported by the Theme Based Research Scheme (T11-705/21-N), 
the Health and Medical Research Fund (COVID1903011-549 WP1) and the Innovation and Technology Commission’s InnoHK funding (D<sup>2</sup>4H).


