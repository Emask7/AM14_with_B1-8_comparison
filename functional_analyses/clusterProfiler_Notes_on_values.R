### clusterProfiler
## https://guangchuangyu.github.io/software/clusterProfiler/
## https://yulab-smu.top/biomedical-knowledge-mining-book/index.html

## Notes on ORA (over-representation analysis) results 
# geneRatio = input DEGs in a term / total # DEGs
# BgRatio = background genes in that term / total # background genes
# RichFactor = # DEGs annotated in a term / # background genes annotated to that term
# FoldEnrichment = geneRatio / BgRatio
  #              = (input DEGs in a term/total # DEGs) / (background genes in that term/total # background genes)
  #              = (input DEGs in a term/background genes in a term) / (total # DEGS/total # background genes)

## Notes on enrichplot
# enrichplot is built on ggplot2, so you can use ggplot2 syntax for enrichplot





### Resources for functional analysis

# g:Profiler - http://biit.cs.ut.ee/gprofiler/index.cgi 
# DAVID - http://david.abcc.ncifcrf.gov/tools.jsp 
# clusterProfiler - http://bioconductor.org/packages/release/bioc/html/clusterProfiler.html
# GeneMANIA - http://www.genemania.org/
# GenePattern -  http://www.broadinstitute.org/cancer/software/genepattern/ (need to register)
# WebGestalt - http://bioinfo.vanderbilt.edu/webgestalt/ (need to register)
# AmiGO - http://amigo.geneontology.org/amigo
# ReviGO (visualizing GO analysis, input is GO terms) - http://revigo.irb.hr/ 
# WGCNA - http://www.genetics.ucla.edu/labs/horvath/CoexpressionNetwork
# GSEA - http://software.broadinstitute.org/gsea/index.jsp
# SPIA - https://www.bioconductor.org/packages/release/bioc/html/SPIA.html
# GAGE/Pathview - http://www.bioconductor.org/packages/release/bioc/html/gage.html
# Reactome Pathway Analysis - https://reactome.org/
# Metascape - https://metascape.org/gp/index.html#/main/step1
# Cytoscape - http://cytoscape.org
# Chea3 - https://maayanlab.cloud/chea3/ (transcription factor enrichment analysis)