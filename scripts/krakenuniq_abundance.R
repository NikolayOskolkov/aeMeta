library("data.table")

Path<-"./results/"; setwd(Path)
my_fastq<-list.dirs(path = ".", recursive = FALSE, full.names = FALSE)
annot<-fread(paste0(Path,"../helping_files/taxDB.orig"),header=FALSE,sep="\t")

my_abund<-list()
for(i in 1:length(my_fastq))
{
print(paste0("Working with sample ",my_fastq[i]," ######################################################################################################"))

df_mammals<-fread(paste0(Path,my_fastq[i],"/krakenuniq/","sequences_",my_fastq[i],".krakenuniq_MAMMALIAN"),header=FALSE,sep="\t")
df_plants<-fread(paste0(Path,my_fastq[i],"/krakenuniq/","sequences_",my_fastq[i],".krakenuniq_PLANT"),header=FALSE,sep="\t")
df_vertebrates<-fread(paste0(Path,my_fastq[i],"/krakenuniq/","sequences_",my_fastq[i],".krakenuniq_OTHER"),header=FALSE,sep="\t")
df_invertebrates<-fread(paste0(Path,my_fastq[i],"/krakenuniq/","sequences_",my_fastq[i],".krakenuniq_INVERTEBRATE"),header=FALSE,sep="\t")
df_microbes_ncbi<-fread(paste0(Path,my_fastq[i],"/krakenuniq/","sequences_",my_fastq[i],".krakenuniq_MICROBE"),header=FALSE,sep="\t")
#df_microbes_gtdb<-fread(paste0(Path,my_fastq[i],"/krakenuniq/","sequences_",my_fastq[i],".krakenuniq_MICROBE_GTDB"),header=FALSE,sep="\t")
df_microbes_gtdb <- data.table()

c_mammals<-as.character(df_mammals[as.character(df_mammals$V1)=="C",]$V2); c_plants<-as.character(df_plants[as.character(df_plants$V1)=="C",]$V2)
c_vertebrates<-as.character(df_vertebrates[as.character(df_vertebrates$V1)=="C",]$V2); c_invertebrates<-as.character(df_invertebrates[as.character(df_invertebrates$V1)=="C",]$V2)
c_microbes_ncbi<-as.character(df_microbes_ncbi[as.character(df_microbes_ncbi$V1)=="C",]$V2); c_microbes_gtdb<-as.character(df_microbes_gtdb[as.character(df_microbes_gtdb$V1)=="C",]$V2)

# mammalian reads
unique_mammalian_reads<-setdiff(c_mammals, c(c_vertebrates, c_plants, c_invertebrates, c_microbes_ncbi, c_microbes_gtdb))
abund_mammals<-as.data.frame(sort(table(df_mammals[as.character(df_mammals$V2)%in%unique_mammalian_reads,]$V3),TRUE))
colnames(abund_mammals)<-c("TAXID","READS")
abund_mammals$NAME<-annot$V3[match(as.character(abund_mammals$TAXID),as.character(annot$V1))]; abund_mammals$RANK<-annot$V4[match(as.character(abund_mammals$TAXID),as.character(annot$V1))]
abund_mammals<-na.omit(abund_mammals)
abund_mammals<-abund_mammals[as.character(abund_mammals$RANK)%in%c("species","subspecies","genus"),]; abund_mammals<-abund_mammals[abund_mammals$READS>5,]

# plant reads
unique_plant_reads<-setdiff(c_plants, c(c_vertebrates, c_mammals, c_invertebrates, c_microbes_ncbi, c_microbes_gtdb))
abund_plants<-as.data.frame(sort(table(df_plants[as.character(df_plants$V2)%in%unique_plant_reads,]$V3),TRUE))
colnames(abund_plants)<-c("TAXID","READS")
abund_plants$NAME<-annot$V3[match(as.character(abund_plants$TAXID),as.character(annot$V1))]; abund_plants$RANK<-annot$V4[match(as.character(abund_plants$TAXID),as.character(annot$V1))]
abund_plants<-na.omit(abund_plants)
abund_plants<-abund_plants[as.character(abund_plants$RANK)%in%c("species","subspecies","genus"),]; abund_plants<-abund_plants[abund_plants$READS>5,]

# vertebrate reads
unique_vertebrate_reads<-setdiff(c_vertebrates, c(c_plants, c_mammals, c_invertebrates, c_microbes_ncbi, c_microbes_gtdb))
abund_vertebrates<-as.data.frame(sort(table(df_vertebrates[as.character(df_vertebrates$V2)%in%unique_vertebrate_reads,]$V3),TRUE))
colnames(abund_vertebrates)<-c("TAXID","READS")
abund_vertebrates$NAME<-annot$V3[match(as.character(abund_vertebrates$TAXID),as.character(annot$V1))]; abund_vertebrates$RANK<-annot$V4[match(as.character(abund_vertebrates$TAXID),as.character(annot$V1))]
abund_vertebrates<-na.omit(abund_vertebrates)
abund_vertebrates<-abund_vertebrates[as.character(abund_vertebrates$RANK)%in%c("species","subspecies","genus"),]; abund_vertebrates<-abund_vertebrates[abund_vertebrates$READS>5,]

# invertebrate reads
unique_invertebrate_reads<-setdiff(c_invertebrates, c(c_plants, c_mammals, c_vertebrates, c_microbes_ncbi, c_microbes_gtdb))
abund_invertebrates<-as.data.frame(sort(table(df_invertebrates[as.character(df_invertebrates$V2)%in%unique_invertebrate_reads,]$V3),TRUE))
colnames(abund_invertebrates)<-c("TAXID","READS")
abund_invertebrates$NAME<-annot$V3[match(as.character(abund_invertebrates$TAXID),as.character(annot$V1))]; abund_invertebrates$RANK<-annot$V4[match(as.character(abund_invertebrates$TAXID),as.character(annot$V1))]
abund_invertebrates<-na.omit(abund_invertebrates)
abund_invertebrates<-abund_invertebrates[as.character(abund_invertebrates$RANK)%in%c("species","subspecies","genus"),]; abund_invertebrates<-abund_invertebrates[abund_invertebrates$READS>5,]

abund_total<-rbind(abund_invertebrates,rbind(abund_vertebrates,rbind(abund_mammals,abund_plants)))
abund_total<-abund_total[order(-abund_total$READS),]
abund_total$SAMPLE<-my_fastq[i]

my_abund[[i]]<-abund_total
}
merged<-Reduce(rbind,my_abund)
merged$NAME<-paste0(merged$NAME," ",merged$RANK)
merged<-merged[merged$READS>=10,]
print(head(merged,20))

unique_samples<-unique(merged$SAMPLE)
unique_species<-unique(merged$NAME)
unique_taxids<-unique(merged$TAXID)

#Computing abundance matrix
abundance_matrix<-matrix(NA,ncol=length(unique_samples),nrow=length(unique_species))
for(i in 1:length(unique_species))
{
 for(j in 1:length(unique_samples))
 {
  if(length(merged[merged$NAME==unique_species[i] & merged$SAMPLE==unique_samples[j],]$READS)==0)
  {
   abundance_matrix[i,j]<-0
  }
  else
  {
   abundance_matrix[i,j]<-merged[merged$NAME==unique_species[i] & merged$SAMPLE==unique_samples[j],]$READS
  }
 }
}
rownames(abundance_matrix)<-unique_species
colnames(abundance_matrix)<-unique_samples
abundance_matrix<-abundance_matrix[order(rownames(abundance_matrix)),]
print(head(abundance_matrix))
write.table(abundance_matrix,file="abundance_matrix_eukaryotes.txt",col.names=TRUE,row.names=TRUE,quote=FALSE,sep="\t")

