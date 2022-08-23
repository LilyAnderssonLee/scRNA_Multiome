convertRatGeneList<-function(geneList){
    require("biomaRt")
    rat = useEnsembl("ensembl", dataset = "rnorvegicus_gene_ensembl",host = "https://dec2021.archive.ensembl.org/")
    mouse = useEnsembl(biomart = "ensembl", dataset = "mmusculus_gene_ensembl",host = "https://dec2021.archive.ensembl.org/")
    
    genesV2 = getLDS(attributes = c("rgd_symbol"), 
                 filters = "rgd_symbol", 
                 values = geneList , 
                 mart = rat, 
                 attributesL = c("mgi_symbol"), 
                 martL = mouse, 
                 uniqueRows=T)
    return (genesV2)

}

df<-read.table('../../../../../bulkRNA/DE/DE_markers/P10vsP60.txt',header = T,sep = '\t')
ConvertedGeneList<-convertRatGeneList(rownames(df))
write.table(ConvertedGeneList,'RatGeneConvertedMouseGenes_list.txt',row.names = T,col.names = T,sep = '\t')

convertMouseGeneList<-function(geneList){
    require("biomaRt")
    mouse = useEnsembl(biomart = "ensembl", dataset = "mmusculus_gene_ensembl",host = "https://dec2021.archive.ensembl.org/")
    human = useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl",host = "https://dec2021.archive.ensembl.org/")
    
    genesV2 = getLDS(attributes = c("mgi_symbol"), 
                 filters = "mgi_symbol", 
                 values = geneList , 
                 mart = mouse, 
                 attributesL = c("hgnc_symbol"), 
                 martL = mouse, 
                 uniqueRows=T)
    return (genesV2)

}

convertHumanGeneList<-function(geneList){
    require("biomaRt")
    human = useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl",host = "https://dec2021.archive.ensembl.org/")
    mouse = useEnsembl(biomart = "ensembl", dataset = "mmusculus_gene_ensembl",host = "https://dec2021.archive.ensembl.org/")
    
    genesV2 = getLDS(attributes = c("hgnc_symbol"), 
                 filters = "hgnc_symbol", 
                 values = geneList , 
                 mart = human, 
                 attributesL = c("mgi_symbol"), 
                 martL = mouse, 
                 uniqueRows=T)
    return (genesV2)

}


