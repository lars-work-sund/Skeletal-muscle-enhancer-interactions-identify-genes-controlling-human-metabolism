# Skeletal-muscle-enhancer-interactions-identify-novel-genes-controlling-whole-body-metabolism-in-huma
Code for the paper: Skeletal muscle enhancer interactions identify novel genes controlling whole body metabolism in humans

This repository contains the needed scripts to regenerate the results used in the paper: Skeletal muscle enhancer interactions identify novel genes controlling whole body metabolism in humans

It includes the raw counts for both the RNA and ChIP experiments. In order to save space the raw counts for the cHiC experment are not included. Instead the processed counts linking enhancers to promoters are included.

# Installation instructions
All scripts are R scripts. A working installation of R is necessary.
The scripts have been tested with R version 3.6.1.

Clone the repository, and install the necessary packages, listed below.
Note that many packages are distributed through Bioconductor (https://bioconductor.org/).

Once Bioconductor has been installed, the following piece of code will install the necessary packages.
```
BiocManager::install(c("data.table", "magrittr", "edgeR", "stringr", "clusterProfiler", "openxlsx", "TxDb.Hsapiens.UCSC.hg38.knownGene", "ggplot2", "pheatmap", "RColorBrewer", "PoiClaClu", "org.Hs.eg.db", "SummarizedExperiment", "reactome.db", "GO.db", "biomaRt", "ChIPseeker", "rtracklayer", "liftOver", "cbmR", "scales", "factoextra", "reshape2", "gridExtra", "vegan", "dplyr", "zoo", "gtable", "grid", "diffHic", "BSgenome.Hsapiens.UCSC.hg38", "csaw", "GenomicRanges", "clipr", "GenomicInteractions", "xlsx", "InteractionSet"))
```

The script runAnalysis.R runs the individual scripts in the correct order. Install time depends on internet speed, but should be less than 30 minutes for most users. Runtime is approx. 15 minutes on the computer used for the analysis.

After running, excel sheets with summary data for the analysis as well as figures included in the publication will be generated.

# Packages used in the analysis
The following packages are loaded during the analysis:
	data.table
	magrittr
	edgeR
	stringr
	clusterProfiler
	openxlsx
	TxDb.Hsapiens.UCSC.hg38.knownGene
	ggplot2
	pheatmap
	RColorBrewer
	PoiClaClu
	org.Hs.eg.db
	SummarizedExperiment
	reactome.db
	GO.db
	biomaRt
	ChIPseeker
	rtracklayer
	liftOver
	scales
	factoextra
	reshape2
	gridExtra
	vegan
	dplyr
	zoo
	gtable
	grid
	diffHic
	BSgenome.Hsapiens.UCSC.hg38
	csaw
	GenomicRanges
	clipr
	GenomicInteractions
	xlsx
	InteractionSet
	AnnotationDbi

# Complete list of packages used
The scripts have been tested with the following sessionInfo():
```
R version 3.6.1 (2019-07-05)
Platform: x86_64-apple-darwin15.6.0 (64-bit)
Running under: macOS Mojave 10.14.6

Matrix products: default
BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
LAPACK: /Library/Frameworks/R.framework/Versions/3.6/Resources/lib/libRlapack.dylib

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:
 [1] grid      parallel  stats4    stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] GenomicInteractions_1.18.1              clipr_0.7.0                            
 [3] csaw_1.18.0                             BSgenome.Hsapiens.UCSC.hg38_1.4.1      
 [5] BSgenome_1.52.0                         Biostrings_2.52.0                      
 [7] XVector_0.24.0                          diffHic_1.16.0                         
 [9] InteractionSet_1.12.0                   gtable_0.3.0                           
[11] zoo_1.8-6                               dplyr_0.8.3                            
[13] vegan_2.5-6                             lattice_0.20-38                        
[15] permute_0.9-5                           gridExtra_2.3                          
[17] reshape2_1.4.3                          factoextra_1.0.5                       
[19] scales_1.0.0                            biomaRt_2.40.4                         
[21] reactome.db_1.68.0                      clusterProfiler_3.12.0                 
[23] ChIPseeker_1.20.0                       PoiClaClu_1.0.2.1                      
[25] RColorBrewer_1.1-2                      pheatmap_1.0.12                        
[27] TxDb.Hsapiens.UCSC.hg38.knownGene_3.4.6 edgeR_3.26.8                           
[29] limma_3.40.6                            SummarizedExperiment_1.14.1            
[31] DelayedArray_0.10.0                     BiocParallel_1.18.1                    
[33] matrixStats_0.55.0                      cbmR_0.1.0                             
[35] stringr_1.4.0                           magrittr_1.5                           
[37] ggplot2_3.2.1                           data.table_1.12.3                      
[39] openxlsx_4.1.0.1                        liftOver_1.8.0                         
[41] rtracklayer_1.44.4                      gwascat_2.16.0                         
[43] Homo.sapiens_1.3.1                      TxDb.Hsapiens.UCSC.hg19.knownGene_3.2.2
[45] org.Hs.eg.db_3.8.2                      GO.db_3.8.2                            
[47] OrganismDbi_1.26.0                      GenomicFeatures_1.36.4                 
[49] GenomicRanges_1.36.1                    GenomeInfoDb_1.20.0                    
[51] AnnotationDbi_1.46.1                    IRanges_2.18.3                         
[53] S4Vectors_0.22.1                        Biobase_2.44.0                         
[55] BiocGenerics_0.30.0                    

loaded via a namespace (and not attached):
  [1] R.utils_2.9.0            tidyselect_0.2.5         htmlwidgets_1.3          RSQLite_2.1.2           
  [5] munsell_0.5.0            statmod_1.4.32           withr_2.1.2              colorspace_1.4-1        
  [9] GOSemSim_2.10.0          knitr_1.25               rstudioapi_0.10          DOSE_3.10.2             
 [13] labeling_0.3             urltools_1.7.3           GenomeInfoDbData_1.2.1   polyclip_1.10-0         
 [17] bit64_0.9-7              farver_1.1.0             rhdf5_2.28.0             vctrs_0.2.0             
 [21] xfun_0.9                 biovizBase_1.32.0        R6_2.4.0                 graphlayouts_0.5.0      
 [25] locfit_1.5-9.1           AnnotationFilter_1.8.0   bitops_1.0-6             fgsea_1.10.1            
 [29] gridGraphics_0.4-1       assertthat_0.2.1         nnet_7.3-12              ggraph_2.0.0            
 [33] enrichplot_1.4.0         ensembldb_2.8.0          tidygraph_1.1.2          rlang_0.4.0             
 [37] zeallot_0.1.0            splines_3.6.1            lazyeval_0.2.2           acepack_1.4.1           
 [41] dichromat_2.0-0          checkmate_1.9.4          europepmc_0.3            BiocManager_1.30.4      
 [45] backports_1.1.4          qvalue_2.16.0            Hmisc_4.2-0              RBGL_1.60.0             
 [49] tools_3.6.1              RMariaDB_1.0.6           ggplotify_0.0.4          gridBase_0.4-7          
 [53] gplots_3.0.1.1           ggridges_0.5.1           Rcpp_1.0.2               plyr_1.8.4              
 [57] base64enc_0.1-3          progress_1.2.2           zlibbioc_1.30.0          purrr_0.3.2             
 [61] RCurl_1.95-4.12          prettyunits_1.0.2        rpart_4.1-15             viridis_0.5.1           
 [65] cowplot_1.0.0            ggrepel_0.8.1            cluster_2.1.0            DO.db_2.9               
 [69] triebeard_0.3.0          ProtGenerics_1.16.0      hms_0.5.1                XML_3.98-1.20           
 [73] compiler_3.6.1           tibble_2.1.3             KernSmooth_2.23-15       crayon_1.3.4            
 [77] htmltools_0.3.6          R.oo_1.22.0              mgcv_1.8-29              Formula_1.2-3           
 [81] tidyr_1.0.0              DBI_1.0.0                tweenr_1.0.1             MASS_7.3-51.4           
 [85] boot_1.3-23              Matrix_1.2-17            R.methodsS3_1.7.1        gdata_2.18.0            
 [89] Gviz_1.28.3              igraph_1.2.4.1           pkgconfig_2.0.3          rvcheck_0.1.3           
 [93] GenomicAlignments_1.20.1 foreign_0.8-72           xml2_1.2.2               VariantAnnotation_1.30.1
 [97] digest_0.6.21            graph_1.62.0             fastmatch_1.1-0          htmlTable_1.13.2        
[101] curl_4.2                 Rsamtools_2.0.1          gtools_3.8.1             lifecycle_0.1.0         
[105] nlme_3.1-141             jsonlite_1.6             Rhdf5lib_1.6.1           viridisLite_0.3.0       
[109] pillar_1.4.2             Rhtslib_1.16.2           httr_1.4.1               plotrix_3.7-6           
[113] survival_2.44-1.1        glue_1.3.1               zip_2.0.4                UpSetR_1.4.0            
[117] bit_1.1-14               ggforce_0.3.1            stringi_1.4.3            blob_1.2.0              
[121] latticeExtra_0.6-28      caTools_1.17.1.2         memoise_1.1.0       
```