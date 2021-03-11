# JSON has funky structure
library(jsonlite)
library(data.table)
file="~/Downloads/hsa00001.json"
document <- fromJSON(txt=file,flatten = T)

# from hgnc... to define ens to ncbi mapping. otherwise it is a mess!
# https://biomart.genenames.org/martform/#!/default/HGNC?datasets=hgnc_gene_mart&attributes=hgnc_gene__approved_symbol_1010%2Chgnc_gene__approved_name_1010%2Chgnc_gene__ensembl_gene__ensembl_gene_id_104%2Chgnc_gene__ncbi_gene__gene_id_1026&hgnc_gene__locus_group_1010=protein-coding+gene&hgnc_gene__status_1010=Approved

hgnc <- fread("data-raw/hgnc.biomart.txt",select = c("Approved symbol","NCBI gene ID","Ensembl gene ID"),col.names = c("geneidA","ncbiid","ens"),header=T)
hgnc[duplicated(ens)]
# ENSG00000255374
hgnc[ens==hgnc[duplicated(ens)]$ens,geneidA]
# [1] "TAS2R43" "TAS2R45"
# already defind go genes
genes <- fread("data-raw/hg19_nonredundant_noORS.gene.ens.txt",col.names = c("geneidA","ens"),header=F)

genes[ens==hgnc[duplicated(ens)]$ens]$geneidA
# [1] "TAS2R43"
hgnc <- hgnc[!ens=="TAS2R45"]
genes <- genes[hgnc,on=.(ens)][,.(geneidA,ens,ncbiid)]
# so now have a consistent 1-1 mapping
# write out ens ids from HGNC. bo to ensembl biomart hg10=9 and get gene coordinates for these!
write.table(hgnc[,ens],file="data-raw/hgnc.ens.txt",sep="\n",col.names = F,row.names = F,quote = F)

# got ens gene defs....
hgnc.ens.defs <- fread("data-raw/hgnc.biomart.ens.chr.s.e.txt",col.names = c("ens","chr","start","end"))
# get them all?
hgnc[!ens %in% hgnc.ens.defs$ens]
# 287 missing...
# this time use HGNC symbols and try again....
write.table(hgnc[!ens %in% hgnc.ens.defs$ens][,geneidA],file="data-raw/hgnc.ens.intitally.missingtxt",sep="\n",col.names = F,row.names = F,quote = F)

hgnc.ens.defs2 <- fread("data-raw/hgnc.biomart.ens.chr.s.e.initially.missing.txt",col.names = c("ens","chr","start","end","geneidA"))
hgnc.ens.defs2 <- hgnc.ens.defs2[geneidA %in% hgnc$geneidA]

hgnc.ens.defs <- rbindlist(list(hgnc.ens.defs,hgnc.ens.defs2[,.(ens,chr,start,end)]))
# filter for autosomes
hgnc.ens.defs <- hgnc.ens.defs[chr %in% 1:22]
hgnc.ens.defs <- hgnc[hgnc.ens.defs,on=.(ens)][!geneidA %in% NA][!ncbiid %in% NA][order(geneidA)]





all.kegg.genes.map <- data.table()
for (i in seq_along(document$children$children)){
  tmp <- document$children$children[[i]]
  for(j in seq_along(tmp$children)){
    for(k in seq_along(tmp$children[[j]]$name)){
      description <- strsplit(tmp$children[[j]]$name[k]," ")[[1]]
      keggid <- strsplit(description[length(description)],split = ":|]")[[1]][2]
      name <- paste(description[-length(description)][-1],collapse = "_")
      for(g in seq_along(tmp$children[[j]]$children[[k]]$name)){
        geneinfo <- strsplit(tmp$children[[j]]$children[[k]]$name[g], " |;")[[1]]
        dt <- data.table(keggid, name,ncbiid=as.integer(geneinfo[1]),geneid=geneinfo[2])
        all.kegg.genes.map <- rbindlist(list(all.kegg.genes.map,dt))
      }
    }
  }
}
# genes wiht more than one ncbi id?
unique(all.kegg.genes.map[,.(ncbiid,geneid)])[,.N,by=geneid][N > 1]
# geneid N
# 1:    CYP2D6 2
# 2:  putative 4
# 3:      ICOS 2
# 4: olfactory 4

# olfactory and putative are nonsense
all.kegg.genes.map <- all.kegg.genes.map[!geneid %in% c("olfactory", "putative") ]

# for CYP2D6 and ICOS, correct ids are 1565 and 29851
unique(all.kegg.genes.map[,.(ncbiid,geneid)])[grep("CYP2D6",geneid)]
# ncbiid geneid
# 1:      1565 CYP2D6
# 2: 107987479 CYP2D6
unique(all.kegg.genes.map[,.(ncbiid,geneid)])[grep("^ICOS$",geneid)]
# ncbiid geneid
# 1: 102723996   ICOS
# 2:     29851   ICOS

all.kegg.genes.map <- all.kegg.genes.map[!ncbiid %in% c("102723996", "107987479") ]



unique(all.kegg.genes.map[,.(ncbiid,geneid)])[,.N,by=geneid][N > 1]
# empty
# so a unique gene name - ncbi mapping
# how many genes in
# total rows
# 46705
all.kegg.genes.map <- all.kegg.genes.map[genes,on=.(ncbiid)]
all.kegg.genes.map <- all.kegg.genes.map[!keggid %in% NA]

all.kegg.genes.map <- unique(all.kegg.genes.map[!gene %in% NA,.(id=keggid,name,gene=geneidA)])
setkey(all.kegg.genes.map,gene)

