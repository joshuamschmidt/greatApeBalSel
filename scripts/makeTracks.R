## generate tracks of ncd scores for hg19 upload to ucsc
library(data.table)
datasuffix="_all_great_apes_ncd_scan_results_data.table.rds"
outgroup="human"
basedir="/Users/joshua_schmidt/Downloads/"
speciesNCDTablesDT <- readRDS(paste(basedir,outgroup,datasuffix,sep=""))
## need to add hg19 coordiantes....
setwd("/Users/joshua_schmidt/Documents/EVA_work/joao_NCD_greatApes/repeat_scan/enrichments/")
id.window.coords <- readRDS("human.window.coords.hg19.rds")
setkey(id.window.coords,id)
setwd("/Users/joshua_schmidt/Documents/EVA_work/joao_NCD_greatApes/repeat_scan/ucsc_tracks")
colours <- c("30,144,255","70,130,180","123,104,238","72,61,139","0,0,128", "255,69,0","255,127,80","139,0,0","220,20,60")
names(colours) <- c("troglodytes","schweinfurthii","ellioti","verus","bonobo","gorilla","graueri","abelii","pygmaeus")
# for (sp in speciesNCDTablesDT[,unique(species)]) {
#     col <- colours[sp]
#     for (c in 1:22) {
#         outfile <- paste(sp,"_chr.",c,"_ncd.log10p_ucsc.track.bedGraph",sep="")
#         headerString <- paste('track type=bedGraph name=',sp,
#             ' description="',sp,
#             ' log10 ncd pvalues" visibility=full color=',col,
#             ' viewLimits=0:10 maxHeightPixels=30:30:30 graphType=bar yLineMark=2 yLineOnOff=on',sep="")
#         writeLines(headerString, outfile)
#         tmp <- speciesNCDTablesDT[species==sp]
#         tmp <- tmp[!ncd2Zp %in% NA]
#         setkey(tmp,id)
#         tmp <- id.window.coords[tmp,on=.(id)][,.(chr,as.integer(hg19start+1449),as.integer(hg19end-1450), -log(ncd2Zp,base=10))]
#         tmp <- tmp[chr == c]
#         tmp[,chr:=paste("chr",chr,sep="")]
#         setkey(tmp,chr,V2,V3)
#         fwrite(tmp,outfile,sep="\t",quote=FALSE,col.names=FALSE,row.names=FALSE,append=TRUE)
#     }
# }

for (sp in speciesNCDTablesDT[,unique(species)]) {
    col <- colours[sp]
        outfile <- paste(sp,"_ncd.log10p_ucsc.track.bedGraph",sep="")
        headerString <- paste('track type=bedGraph name=',sp,
            ' description="',sp,
            ' log10 ncd pvalues" visibility=full color=',col,
            ' viewLimits=0:10 maxHeightPixels=30:30:30 graphType=bar yLineMark=2 yLineOnOff=on',sep="")
        writeLines(headerString, outfile)
        tmp <- speciesNCDTablesDT[species==sp]
        tmp <- tmp[!ncd2Zp %in% NA]
        setkey(tmp,id)
        tmp <- id.window.coords[tmp,on=.(id)][,.(chr,as.integer(hg19start+1449),as.integer(hg19end-1450), -log(ncd2Zp,base=10))]
        tmp[,chr:=paste("chr",chr,sep="")]
        tmp <- tmp[V3-V2==101]
        setkey(tmp,chr,V2,V3)
        fwrite(tmp,outfile,sep="\t",quote=FALSE,col.names=FALSE,row.names=FALSE,append=TRUE)
}
for f in $(ls *.bedGraph | sed 's/\.bedGraph//g'); do sort -k1,1 -k2,2n ${f}.bedGraph |awk '$3-$2==101' > ${f}.bedGraph.sorted;./bedGraphToBigWig ${f}.bedGraph.sorted hg19.chrom.sizes ${f}.bw; rm ${f}.bedGraph.sorted; rm ${f}.bedGraph; done
## OAS1 region.... chr12:113,275,278-113,455,017


https://github.com/joshuamschmidt/ncd_scores/raw/master/bonobo_ncd.log10p_ucsc.track.bw
https://github.com/joshuamschmidt/ncd_scores/raw/master/ellioti_ncd.log10p_ucsc.track.bw
https://github.com/joshuamschmidt/ncd_scores/raw/master/verus_ncd.log10p_ucsc.track.bw
https://github.com/joshuamschmidt/ncd_scores/raw/master/schweinfurthii_ncd.log10p_ucsc.track.bw
https://github.com/joshuamschmidt/ncd_scores/raw/master/troglodytes_ncd.log10p_ucsc.track.bw
https://github.com/joshuamschmidt/ncd_scores/raw/master/gorilla_ncd.log10p_ucsc.track.bw
https://github.com/joshuamschmidt/ncd_scores/raw/master/graueri_ncd.log10p_ucsc.track.bw
https://github.com/joshuamschmidt/ncd_scores/raw/master/abelii_ncd.log10p_ucsc.track.bw
https://github.com/joshuamschmidt/ncd_scores/raw/master/pygmaeus_ncd.log10p_ucsc.track.bw