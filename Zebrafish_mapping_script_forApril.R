data.dir<-"/gpfs/ysm/project/hoffman/yl883/RNAseq/HumanBrainRef/sc_dev_cortex_geschwind/"
setwd(data.dir)
old_map<-read.csv("/gpfs/ysm/project/hoffman/yl883/old_humanGeneMapping.csv")
tpm<-as.data.frame(read_tsv('/gpfs/ysm/project/hoffman/yl883/RNAseq/FromDejian/snc1_TPM.tsv'))

old_map[,1]<-as.character(old_map[,1])
old_map[,2]<-as.character(old_map[,2])
old_map[,3]<-as.character(old_map[,3])
tpm_geneid<-gsub(":.*","",tpm[,1])




old_map[grep("^[0-9]",old_map$Human.associated.gene.name),"Human.associated.gene.name"]
library(ensembldb)
library(EnsDb.Hsapiens.v86)
edb <- EnsDb.Hsapiens.v86
matchid_df<-select(edb, keys = old_map[grep("^[0-9]",old_map$Human.associated.gene.name),2], keytype = "GENEID",columns = c("GENENAME"))
rownames(matchid_df)<-matchid_df[,1]
old_map[grep("^[0-9]",old_map$Human.associated.gene.name),2]<-matchid_df[old_map[grep("^[0-9]",old_map$Human.associated.gene.name),2],2]
db_mapping<-select(edb, keys = old_map[,2], keytype = "GENEID",columns = c("GENENAME")) #this is to map the gene ID you have (old_map[,2]) to the gene ID in ensembldb to make sure gene names are correct. If gene names in old_map don't have false names (i.e. gene names turned into dates) then this is not necessary 
colnames(db_mapping)[1]<-"Human.Ensembl.Gene.ID"
old_mapa<-merge(old_map,db_mapping,by="Human.Ensembl.Gene.ID")
old_mapb<-merge(old_mapa,tpm,by="Ensembl.Gene.ID")
sort(table(old_mapb[,1]),decreasing = T)[1:10]
old_mapb[which(old_mapb[,1]=="ENSDARG00000058212"),c(4,5)]
sort(table(old_mapb[,3]),decreasing = T)[1:10]
old_mapb[which(old_mapb[,3]=="TRIM47"),c(4,5)]
sum(old_mapb[,4]%in%rownames(raw_counts_mat))
length(intersect(unique(old_mapb[,4]),rownames(raw_counts_mat)))
old_mapc<-old_mapb[which(old_mapb[,4]%in%rownames(raw_counts_mat)),]
sort(table(old_mapc[,1]),decreasing = T)[1:10]
old_mapc[which(old_mapc[,1]=="ENSDARG00000086494"),c(4,5)]
sum(table(old_mapc[,1])>1)
sort(table(old_mapc[,4]),decreasing = T)[1:10]
old_mapc[which(old_mapc[,4]=="TRIM47"),c(4,5)]
sum(table(old_mapc[,4])>1)
sum(table(old_mapc[,4])[table(old_mapc[,4])>1])

# for(i in names(table(old_mapc[,4]))){
#   print(i)
#   tmp.df<-old_mapc[which(old_mapc[,4]==i),]
#   if(nrow(tmp.df)>1){
#     tmp.df<-tmp.df[which.max(rowMeans(tmp.df[,-c(1:5)])),]
#   }
#
# }

old_mapd<-Reduce(rbind,lapply(names(table(old_mapc[,4])),function(i){
  tmp.df<-old_mapc[which(old_mapc[,4]==i),]
  if(nrow(tmp.df)>1){
    tmp.df<-tmp.df[which.max(rowMeans(tmp.df[,-c(1:5)])),]
  }
  return(tmp.df)
}))
old_mapd[is.na(old_mapd)]<-0
saveRDS(old_mapd,"./MappedBulkMutants_6mutants.rds")