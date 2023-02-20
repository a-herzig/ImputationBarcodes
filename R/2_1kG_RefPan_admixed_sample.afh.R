R
library(gaston)
library(gaston.utils)
library(Mozza)
options(scipen=999)



for(chr in 1:22){

directory_sample <- paste("/WORKING_DIRECTORY/aherzig/tdekNOV22/CHR_",chr,sep="")
directory_refpan <- paste("/WORKING_DIRECTORY/aherzig/tdekNOV22/CHR_",chr,sep="")

KG_haplo_tagged <- read.bed.matrix(paste("/WORKING_DIRECTORY/aherzig/tdekNOV22/haplo_tagged_chr",chr,sep=""))


## Study sample: population from "ASW","ACB" and "MXL" of 1kG.

system(paste("ls ",directory_sample," 2>/dev/null||mkdir ",directory_sample,sep=""))

# creating virtual position to format file with shapeit.
position_correspondence <- cbind(1:length(KG_haplo_tagged@snps$pos),KG_haplo_tagged@snps$pos)
position_for_vcf <- cbind(rep(chr,length(KG_haplo_tagged@snps$pos)), KG_haplo_tagged@snps$pos)
write.table(position_for_vcf,paste(directory_sample,"/real_pos.vcf",sep=""),col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t")


studyS <- KG_haplo_tagged[which(KG_haplo_tagged@ped$population%in%c("ACB","ASW","MXL")),]
studyS@snps$pos <- position_correspondence[,1]

write.hap.file(studyS,paste(directory_sample,"/studyS.hap",sep=""))

leg <- cbind(chr,studyS@snps$id,studyS@snps$pos,studyS@snps$A1,studyS@snps$A2)
write.table(leg,paste(directory_sample,"/studyS.legend",sep=""),col.names=F,row.names=FALSE,quote=FALSE,sep=" ")

samp<-cbind(substr(studyS@ped$id,1,7)[c(TRUE,FALSE)],substr(studyS@ped$id,1,7)[c(TRUE,FALSE)],0,0,0,sample(c(1,2),length(studyS@ped$id)/2,replace=TRUE),-9)
samp<-rbind(c("0","0","0","D","D","D","B"),samp)
write.table(samp,paste(directory_sample,"/studyS.sample",sep=""),col.names=c("ID_1","ID_2","missing","father","mother","sex","plink_pheno"),row.names=FALSE,quote=FALSE,sep=" ")

pasteL<-paste(directory_sample,"/studyS.legend ",directory_sample,"/studyS.hap",sep="")
system(paste("paste ",pasteL," > ",directory_sample,"/studyS.haps",sep=""))

system(paste("rm ",directory_sample,"/studyS.hap",sep=""))
system(paste("rm ",directory_sample,"/studyS.legend",sep=""))

hap1<-paste(directory_sample,"/studyS",sep="")
out1<-paste(directory_sample,"/studyS.phased.vcf",sep="")

system(paste("shapeit -convert --thread 10 --input-haps ",hap1," --output-vcf ",out1,sep=""))

system(paste("rm ",directory_sample,"/studyS.haps",sep=""))
system(paste("rm ",directory_sample,"/studyS.sample",sep=""))

# Genotype array positions:
ukbb_array<-read.table("/WORKING_DIRECTORY/share/Imputation_UKBB_array.txt",header=F,as.is=T)

ukbb_array_chr <- ukbb_array[ukbb_array[,1] == chr,]
ukbb_array_chr <- ukbb_array_chr[order(ukbb_array_chr[,2]),]
ukbb_array_chr <- ukbb_array_chr[!(duplicated(ukbb_array_chr[,2])),]

m <- match(ukbb_array_chr[,2],position_correspondence[,2])
ukbb_array_chr[,3] <- position_correspondence[m,1]

leg2 <- leg[substr(leg[,2],1,2)=="rs",]
w<-which(paste(leg2[,1],leg2[,3],sep=":")%in%paste(ukbb_array_chr[,1],ukbb_array_chr[,3],sep=":"))

write.table(cbind(leg2[w,1],leg2[w,3]),paste(directory_sample,"/ukbbStudyS.txt",sep=""),row.names=F,col.names=F,quote=F,sep="\t")

system(paste("vcftools --vcf ",directory_sample,"/studyS.phased.vcf --out ",directory_sample,"/studyS.phased.ukbb --recode --positions ",directory_sample,"/ukbbStudyS.txt",sep=""))
system(paste("mv ",directory_sample,"/studyS.phased.ukbb.recode.vcf ",directory_sample,"/studyS.phased.ukbb.vcf",sep=""))
system(paste("rm ",directory_sample,"/studyS.phased.ukbb.log",sep=""))

system(paste("head -6 ",directory_sample,"/studyS.phased.vcf > ",directory_sample,"/studyS_temporary.vcf",sep=""))
system(paste("sed 1,6d ",directory_sample,"/studyS.phased.vcf | paste ",directory_sample,"/real_pos.vcf - | cut -f -2,5- >> ",directory_sample,"/studyS_temporary.vcf",sep=""))
system(paste("mv ",directory_sample,"/studyS_temporary.vcf ",directory_sample,"/studyS.phased.vcf",sep=""))

ukbbStudyS <- read.table(paste(directory_sample,"/ukbbStudyS.txt",sep=""))
ukbbStudyS[,2] <- position_correspondence[ukbbStudyS[,2],2]
write.table(ukbbStudyS,paste(directory_sample,"/ukbbStudyS.txt",sep=""),col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t")

system(paste("head -6 ",directory_sample,"/studyS.phased.ukbb.vcf > ",directory_sample,"/studyS_temporary.phased.ukbb.vcf",sep=""))
system(paste("sed 1,6d ",directory_sample,"/studyS.phased.ukbb.vcf | paste ",directory_sample,"/ukbbStudyS.txt - | cut -f -2,5- >> ",directory_sample,"/studyS_temporary.phased.ukbb.vcf",sep=""))
system(paste("mv ",directory_sample,"/studyS_temporary.phased.ukbb.vcf ",directory_sample,"/studyS.phased.ukbb.vcf",sep=""))

# separate Haplotype and create 2n pseudo individuals to ensure haplotye-by-haplotype imputation

ukbb_samp_simple <- read.table(paste(directory_sample,"/studyS.phased.ukbb.vcf",sep=""),header=F)

lineChrom<-strsplit(system(paste("grep -n -m 1 \"CHROM\" ",directory_sample,"/studyS.phased.ukbb.vcf",sep=""),intern=T),":")[[1]][1]
lineChrom<-as.numeric(lineChrom)

namesChrom<-strsplit(system(paste("grep -n -m 1 \"CHROM\" ",directory_sample,"/studyS.phased.ukbb.vcf",sep=""),intern=T),"::")[[1]]
namesChrom<-unlist(strsplit(namesChrom,"\t"));namesChrom[1]<-"#CHROM"

namesChromN<-namesChrom[1:9]

ukbb_samp_double<-ukbb_samp_simple[,1:9]
for (k in 10:ncol(ukbb_samp_simple)){
  h1<-substr(ukbb_samp_simple[,k],1,1)
  h2<-substr(ukbb_samp_simple[,k],3,3)
  ukbb_samp_double<-cbind(ukbb_samp_double,paste(h1,h1,sep="|"),paste(h2,h2,sep="|"))
  namesChromN<-c(namesChromN,paste(namesChrom[k],"h1",sep="_"),paste(namesChrom[k],"h2",sep="_"))
}

system(paste("head -",lineChrom-1," ",directory_sample,"/studyS.phased.ukbb.vcf > ",directory_sample,"/studyS.phased.ukbb.h1h2.vcf",sep=""))
write.table(ukbb_samp_double,paste(directory_sample,"/studyS.phased.ukbb.h1h2.vcf",sep=""),sep="\t",quote=F,col.names=namesChromN,row.names=F,append=T)

system(paste("bgzip -f ",directory_sample,"/studyS.phased.ukbb.h1h2.vcf",sep=""))
system(paste("tabix -f ",directory_sample,"/studyS.phased.ukbb.h1h2.vcf.gz",sep=""))

system(paste("bgzip -f ",directory_sample,"/studyS.phased.ukbb.vcf",sep=""))
system(paste("tabix -f ",directory_sample,"/studyS.phased.ukbb.vcf.gz",sep=""))

system(paste("bgzip -f ",directory_sample,"/studyS.phased.vcf",sep=""))
system(paste("tabix -f ",directory_sample,"/studyS.phased.vcf.gz",sep=""))
system(paste("rm ",directory_sample,"/real_pos.vcf",sep=""))


position_correspondence <- cbind(1:length(KG_haplo_tagged@snps$pos),KG_haplo_tagged@snps$pos)
position_for_vcf <- cbind(rep(chr,length(KG_haplo_tagged@snps$pos)), KG_haplo_tagged@snps$pos)
write.table(position_for_vcf,paste(directory_refpan,"/real_pos.vcf",sep=""),col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t")

### version of Ref Panel:
vrp_1 <- c("GBR","FIN","CHS","PUR","CDX","CLM","IBS","PEL","PJL","KHV","GWD","ESN","BEB","MSL","STU","ITU","CEU","YRI","CHB","JPT","LWK","TSI","GIH") # == All !%in% c(ACB,ASW,MXL)


RefPanel <- KG_haplo_tagged[which(KG_haplo_tagged@ped$population %in% vrp_1),]
RefPanel@snps$pos <- position_correspondence[,1]

write.hap.file(RefPanel,paste(directory_refpan,"/RefPanel.hap",sep=""))

leg <- cbind(chr,RefPanel@snps$id,RefPanel@snps$pos,RefPanel@snps$A1,RefPanel@snps$A2)
write.table(leg,paste(directory_refpan,"/RefPanel.legend",sep=""),col.names=F,row.names=FALSE,quote=FALSE,sep=" ")

samp <- cbind(substr(RefPanel@ped$id,1,7)[c(TRUE,FALSE)],substr(RefPanel@ped$id,1,7)[c(TRUE,FALSE)],0,0,0,1,-9)
samp <- rbind(c("0","0","0","D","D","D","B"),samp)
write.table(samp,paste(directory_refpan,"/RefPanel.sample",sep=""),col.names=c("ID_1","ID_2","missing","father","mother","sex","plink_pheno"),row.names=FALSE,quote=FALSE,sep=" ")

pasteL <- paste(directory_refpan,"/RefPanel.legend ",directory_refpan,"/RefPanel.hap",sep="")
system(paste("paste ",pasteL," > ",directory_refpan,"/RefPanel.haps",sep=""))

system(paste("rm ",directory_refpan,"/RefPanel.hap",sep=""))
system(paste("rm ",directory_refpan,"/RefPanel.legend",sep=""))

hap1 <- paste(directory_refpan,"/RefPanel",sep="")
out1 <- paste(directory_refpan,"/RefPanel.phased.vcf",sep="")

system(paste("shapeit -convert --input-haps ",hap1," --output-vcf ",out1,sep=""))

system(paste("rm ",directory_refpan,"/RefPanel.haps",sep=""))
system(paste("rm ",directory_refpan,"/RefPanel.sample",sep=""))

system(paste("head -6 ",directory_refpan,"/RefPanel.phased.vcf > ",directory_refpan,"/RefPanel_temporary.vcf",sep=""))
system(paste("sed 1,6d ",directory_refpan,"/RefPanel.phased.vcf | paste ",directory_refpan,"/real_pos.vcf - | cut -f -2,5- >> ",directory_refpan,"/RefPanel_temporary.vcf",sep=""))
system(paste("mv ",directory_refpan,"/RefPanel_temporary.vcf ",directory_refpan,"/RefPanel.phased.vcf",sep=""))

system(paste("bgzip -f ",directory_refpan,"/RefPanel.phased.vcf",sep=""))
system(paste("tabix -f ",directory_refpan,"/RefPanel.phased.vcf.gz",sep=""))

system(paste("rm ",directory_refpan,"/real_pos.vcf",sep=""))






input_file_refpan <- paste(directory_refpan,"/RefPanel.phased.vcf.gz",sep="")
input_file_sample <- paste(directory_sample,"/studyS.phased.ukbb.h1h2.vcf.gz",sep="")

###Coordinates
system(paste("/PROGS/EXTERN/IMPUTE5/imp5Chunker_1.1.5_static --h ",input_file_refpan," --g ",input_file_sample," --r ",chr," --o ",directory_refpan,"/coordinates.txt",sep=""))
system(paste("/PROGS/EXTERN/IMPUTE5/imp5Converter_1.1.5_static --h ",input_file_refpan," --r ",chr," --o ",directory_refpan,"/RefPanel.phased.imp5",sep=""))
coordinates <- read.table(paste(directory_refpan,"/coordinates.txt",sep=""),header=F,as.is=T)

for (j in 1:nrow(coordinates)){

system(paste("/PROGS/EXTERN/IMPUTE5/impute5_1.1.5_static --m /PUBLIC_DATA/ReferencePanels/1kG/beagle_vcf_ref/plink.chr",chr,".GRCh37.map --h ",directory_refpan,"/RefPanel.phased.imp5 --g ",input_file_sample," --r ",coordinates[j,4]," --buffer-region ",coordinates[j,3]," --o ",directory_refpan,"/studyS.phased.ukbb.h1h2_chunk",coordinates[j,1],".vcf.gz --l ",directory_refpan,"/studyS.phased.ukbb.h1h2_chunk",coordinates[j,1],".log",sep=""))

if(j==1){
write(paste(directory_refpan,"/studyS.phased.ukbb.h1h2_chunk",coordinates[j,1],".vcf.gz",sep=""),paste(directory_refpan,"/ligate.txt",sep=""))
} else {
write(paste(directory_refpan,"/studyS.phased.ukbb.h1h2_chunk",coordinates[j,1],".vcf.gz",sep=""),paste(directory_refpan,"/ligate.txt",sep=""),append=T)
}
}

system(paste("bcftools concat -n -f ",directory_refpan,"/ligate.txt -Oz -o ",directory_refpan,"/studyS.phased.ukbb.h1h2_chunkAll.vcf.gz",sep=""))


for (j in 1:nrow(coordinates)){
system(paste("rm ",directory_refpan,"/studyS.phased.ukbb.h1h2_chunk",coordinates[j,1],".vcf.*",sep=""))
}

system(paste("rm ",directory_refpan,"/coordinates.txt",sep=""))
system(paste("rm ",directory_refpan,"/RefPanel.phased.imp5*",sep=""))
system(paste("rm ",directory_refpan,"/ligate.txt",sep=""))
system(paste("rm ",directory_refpan,"/studyS.phased.ukbb.h1h2_chunk*.log",sep=""))

}



gc(reset = TRUE)
rm(list=ls())
q()
n

