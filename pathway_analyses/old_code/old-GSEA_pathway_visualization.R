# pathview(gene.data = PL23_ORA_entrez_list, gene.idtype = "entrez",
#          pathway.id = 'mmu05322',
#          species = "mmu",
#          limit = list(gene=max(abs(PL23_ORA_entrez_list)), cpd=1))



res_PL23_full
res_R848_full
res_PL23vR848_full

pathview(gene.data = PL23_ORA_entrez_list, gene.idtype = "entrez",
                  pathway.id = 'mmu05322',
                  species = "mmu",
                  limit = list(gene=max(abs(PL23_ORA_entrez_list)), cpd=1))

PL23_pathways <- res_PL23_full[, c(6)]
names(PL23_pathways) <- res_PL23_full$entrezgene
# colnames(PL23_pathways) <- c("ENTREZID", "LFC", "padj")
length(PL23_pathways)
PL23_pathways <- subset(PL23_pathways, !duplicated(PL23_pathways))
head(PL23_pathways)

# REACTOME	SLC_MEDIATED_TRANSMEMBRANE_TRANSPORT
viewPathway("SLC-mediated transmembrane transport", 
            organism = "mouse", foldChange = PL23_pathways)
