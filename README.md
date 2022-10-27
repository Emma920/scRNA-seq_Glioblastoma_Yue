# scRNA-seq_Glioblastoma_Yue
1. under the /bash folder, there are pre-processing script, using STAR, trimmomatic, fastqc, etc.
2. Under the /R folder, go directly to main_prom1_md_Rmd.Rmd and main_plaur_rank_md.Rmd. The former project was meant to design a quantitative marker selection system for glioblastoma stem-like cells, using multiple reality-relevant parameters, for their potential use as a biomedical marker in radionucleotide therapy, phage display, in-vitro diagnosis, etc. The latter project investigated a specific protein as a biomedical marker regarding its subtype specificity and protein networks, this study is to support a phase-II clinical trial on this marker. 
3. Some other folders contain different methods used for different stages of pre-processing, Like different ways of clusting (leiden, louvain, with defalut or manually-set resolution and k values), different ways to find highly-viriable-genes(HVG, I ended up not doing any dimension reduction with HVG to maximize the findings), different normalization methods (pre-clusting, library_size, SC_norm, etc... I found normalization to be the most important part of any pre-processing, Normalization is easily over-done with a "fancy" method, that could also delete the real biological signal, I ended up using the most simple normalization methods), Different combinations of pre-processing, normalization, and clusting methods were tested.
Under the folder /correlation_between_features, contains a project to test correlation and distance correlation between two biological features of glioblastoma, I used bootstrapping, distance correlation between the nearest neighbors. This project did not find positive results. But it was fun:)
I also put my on-going projects in a folder  
I arranged my code to different folders before I upload them, it might not match the directories in the code  

# Tools
scanpy
SingleCellExperiment
igraph
edgeR
monocle
ggplot
GSEA 
Cytoscape
g:Profiler
scater
Combat
...
