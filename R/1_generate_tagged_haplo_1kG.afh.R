options(scipen=999)

# chargement des fonctions set.stats_haplotype() et select.snps_haplotype()
select.snps_haplotype <-function (x, condition)
{
    if (!is(x, "bed.matrix"))
        stop("x is not a bed.matrix")
    w <- eval(substitute(condition), x@snps, parent.frame())
    miss <- is.na(w)
    if (sum(miss) > 0) {
        warning(paste(sum(miss), "SNP(s) with undefined condition are removed from bed.matrix"))
        w <- w & !miss
    }
    x <- x[, w]
	drops <- c("N2","N2.x","N2.y","N2.mt","hz","hz.x","hz.y","hz.mt","N2.f")
	x <-list(x@ped[,!(names(x@ped) %in% drops)],x@snps,x@bed,x@p,x@mu,x@sigma,x@standardize_p,x@standardize_mu_sigma)
	x <- new("bed.matrix", bed = x[[3]], snps = x[[2]], ped = x[[1]], p = x[[4]], mu = x[[5]], sigma = x[[6]], standardize_p = x[[7]], standardize_mu_sigma = x[[8]])
	x
}

set.stats_haplotype <- function (x, set.p = TRUE, set.mu_sigma = TRUE, verbose = getOption("gaston.verbose",
    TRUE)){
{
    if (is(x) != "bed.matrix")
        stop("x must be an object of class bed.matrix")
    if (!is.logical(set.p) | !is.logical(set.mu_sigma))
        stop("set.* arguments must be logical")
    w.a <- is.autosome(x@snps$chr)
    w.x <- is.chr.x(x@snps$chr)
    w.y <- is.chr.y(x@snps$chr)
    w.mt <- is.chr.mt(x@snps$chr)
    w.f <- x@ped$sex == 2
    w.f[is.na(w.f)] <- FALSE
    st <- .Call("gg_geno_stats", PACKAGE = "gaston", x@bed, w.x,
        w.y, w.mt, w.f)
	drops <- c("N2","N2.x","N2.y","N2.mt","hz","hz.x","hz.y","hz.mt","N2.f")
	st <- list(st[["snps"]][,!(names(st[[1]]) %in% drops)],st[["inds"]][,!(names(st[[2]]) %in% drops)])
	names(st) <- c("snps","inds")
	nb.f <- sum(x@ped$sex == 2)
    nb.h <- nrow(x) - nb.f
    st$snps$callrate <- 1 - st$snps$NAs/nrow(x)
    st$snps$callrate[w.y] <- 1 - (st$snps$NAs[w.y] - st$snps$NAs.f[w.y])/nb.h
    n <- nrow(x) - st$snps$NAs
    pp <- st$snps$N1 / n
    st$snps$N0.f[!w.x & !w.y] <- NA
    st$snps$N1.f[!w.x & !w.y] <- NA
    st$snps$NAs.f[!w.x & !w.y] <- NA
    a <- st$snps$N1.f[w.x]
    b <- st$snps$N0.f[w.x] + st$snps$N1.f[w.x]
    pp[w.x] <- a/(a + b)
    st$snps$maf <- pmin(pp, 1 - pp)
    n.a <- sum(w.a)
    st$inds$callrate <- 1 - st$inds$NAs/n.a
    n.x <- sum(w.x)
    st$inds$callrate.x <- 1 - st$inds$NAs.x/n.x
    n.y <- sum(w.y)
    st$inds$callrate.y <- 1 - st$inds$NAs.y/n.y
    n.mt <- sum(w.mt)
    st$inds$callrate.mt <- 1 - st$inds$NAs.mt/n.mt
    x@snps[, names(st$snps)] <- st$snps
    x@ped[, names(st$inds)] <- st$inds
    if (verbose)
        cat("ped stats and snps stats have been set. \n")
    if (set.p) {
        x@p <- pp
        if (verbose)
            cat("'p' has been set. \n")
    }
    if (set.mu_sigma) {
        n <- nrow(x) - x@snps$NAs
        # mu <- (2 * x@snps$N2 + x@snps$N1)/n #à vérifier
        N <- nrow(x)
        # s <- sqrt((x@snps$N1 + 4 * x@snps$N2 + mu^2 * x@snps$NAs)/(N - 1) - N/(N - 1) * mu^2) #à vérifier
        # x@mu <- mu
        # x@sigma <- s
        if (verbose)
            cat("'mu' and 'sigma' have been set.\n")
    }
  drops <- c("N2","N2.x","N2.y","N2.mt","hz","hz.x","hz.y","hz.mt","N2.f")
	x <- list(x@ped[,!(names(x@ped) %in% drops)],x@snps[,!(names(x@snps) %in% drops)],x@bed,x@p,x@mu,x@sigma,x@standardize_p,x@standardize_mu_sigma)
	x <- new("bed.matrix", bed = x[[3]], snps = x[[2]], ped = x[[1]], p = x[[4]], mu = x[[5]], sigma = x[[6]], standardize_p = x[[7]], standardize_mu_sigma = x[[8]])
	x
}
}

library(gaston)
library(Mozza)
library(gaston.utils)

for(chr in 1:22) {
print(chr)
### Things you will need, vcf of your reference panel (eg. 1000 Genomes)
### Set of snps-array positions of interest (eg. a table of chr and position numbers)
### Information about groups of interest within the reference panel to trace (eg. the 1000G sample file)
### I have added lines for add an estimated genetic distance for each marker - not really necesary here but as you might
### want to simulat some mosaics with Mozza, this could be useful.
filename<-paste("/PUBLIC_DATA/ReferencePanels/1kG/beagle_vcf_ref/chr",chr,".1kg.phase3.v5a.vcf.gz",sep="")
namesChromR<-strsplit(system(paste("zcat ",filename," | grep -n -m 1 \"CHROM\"",sep=""),intern=T),":")[[1]][2]
namesChromR<-unlist(strsplit(namesChromR,"\t"));namesChromR[1]<-"CHROM";samples<-namesChromR[-c(1:9)]
KG_haplo_raw <- .Call("gg_read_vcf_chr_range_haplo", PACKAGE = "gaston.utils", filename, FALSE, -1L, 0L, 0L, samples)
KG <- read.vcf(file=filename)
KG.samples<- read.table("/PUBLIC_DATA/ReferencePanels/1kG/impute/Phase3/1000GP_Phase3/1000GP_Phase3.sample",header=T,as.is=T)
map <- read.table(paste("/PUBLIC_DATA/ReferencePanels/1kG/beagle_vcf_ref/plink.chr",chr,".GRCh37.map",sep=""),header=F,as.is=T)
ukbb_array<-read.table("/WORKING_DIRECTORY/share/Imputation_UKBB_array.txt",header=F,as.is=T)

### Here we are reading and organising the haplotypes as gaston ususally works with diploid data


m <- match(KG@snps$pos, map[,4]) 
distcM<-cbind(KG@snps$pos,map[m,3],rep(1,length(KG@snps$pos)))
map2 <- map[is.na(match(map[,4],distcM[,1])),c(4,3)] ; map2 <- cbind(map2,rep(2,nrow(map2)))
c("pos","dist","tag") -> colnames(distcM) -> colnames(map2)
map_tagged <- rbind(distcM,map2)
map_tagged <- map_tagged[order(map_tagged[,1], decreasing = FALSE),]
w2<-which(is.na(map_tagged[,2]))
d3<-approx(as.numeric(map_tagged[,1]),as.numeric(map_tagged[,2]),xout=map_tagged[w2,1],rule=1)
map_tagged[w2,2]<-d3$y
distcM <- map_tagged[which(map_tagged[,3] == 1),c(1,2)]
KG@snps$dist <- distcM[,2]

rm(list=c("map2","map_tagged","w2","d3","m"))

famid <- rep(samples, each = 2)
id <- paste0( famid, c(".1", ".2"))

ped <- data.frame(famid = famid, id = id, father = 0, mother = 0, sex = NA, pheno = NA, stringsAsFactors = FALSE)
m <- match(ped$famid, KG.samples$ID)  
ped$population       <- KG.samples$POP[m]
ped$super.population <- KG.samples$GROUP[m]

snp <- data.frame(chr = KG_haplo_raw$chr, id = ifelse( KG_haplo_raw$id ==".", paste(KG_haplo_raw$chr,KG_haplo_raw$pos,KG_haplo_raw$A1,KG_haplo_raw$A2,sep=":"),KG_haplo_raw$id), dist = distcM[,2], pos = KG_haplo_raw$pos , A1 = KG_haplo_raw$A1, A2 = KG_haplo_raw$A2,quality = KG_haplo_raw$quality, filter = factor(KG_haplo_raw$filter), stringsAsFactors = FALSE)
KG_haplo <- new("bed.matrix", bed = KG_haplo_raw$bed, snps = snp, ped = ped, p = NULL, mu = NULL, sigma = NULL, standardize_p = FALSE, standardize_mu_sigma = FALSE )

### Now that the haplotypes are storedm we remove some variants, eg. multiallelic

KG_haplo <- select.snps_haplotype(KG_haplo, (KG_haplo@snps$pos >= head(map[,4],1) & KG_haplo@snps$pos <= tail(map[,4],1)))
KG_haplo <- KG_haplo[,-which(duplicated(KG_haplo@snps$pos) | duplicated(KG_haplo@snps$pos,fromLast=TRUE))]
KG_haplo <- select.snps_haplotype(KG_haplo, (!grepl(",",KG_haplo@snps$A2)))
KG_haplo <- set.stats_haplotype(KG_haplo)

KG@snps$id <- ifelse( KG@snps$id ==".", paste(KG@snps$chr,KG@snps$pos,KG@snps$A1,KG@snps$A2,sep=":"),KG@snps$id)
KG <- select.snps(KG, (KG@snps$pos >= head(map[,4],1) & KG@snps$pos <= tail(map[,4],1)))
KG <- KG[,-which(duplicated(KG@snps$pos))]
KG <- select.snps(KG, (!grepl(",",KG@snps$A2)))
KG <- set.stats(KG)

### Choose some variants to tag

key_variants <- select.snps(KG,(paste(KG@snps$chr,KG@snps$pos,sep=":") %in% paste(ukbb_array[,1],ukbb_array[,2],sep=":")))
key_variants <- select.snps(key_variants, key_variants@snps$maf > 0.2)
key_variants <- LD.thin(key_variants, 0.02)
key_pos <- key_variants@snps$pos + 1
key_pos <- key_pos[!(key_pos %in% KG@snps$pos)]

### add barcodes, here is a version where one barcode is added for each of the 26 populations of 1000G - can easily make adjustments to tag only the 5 super populations

nb_tag<-26 
barcode <- rep(key_pos,each=nb_tag) + ((1:nb_tag)/(nb_tag+1))


## five super population version 
#nb_tag<-5 
#barcode <- rep(key_pos,each=nb_tag) + ((1:nb_tag)/(nb_tag+1))


snp2 <- data.frame(chr = chr, id = paste("mosaicbarcode",rep(1:length(key_pos),each=nb_tag),unique(KG_haplo@ped$population),sep=":"), dist = rep(KG_haplo@snps$dist[which(KG_haplo@snps$pos%in%(key_pos-1))],each=nb_tag), pos = barcode , A1 = "G", A2 = "T",quality = NA, filter = NA, stringsAsFactors = FALSE)

## five super population version 
#snp2 <- data.frame(chr = chr, id = paste("mosaicbarcode",rep(1:length(key_pos),each=nb_tag),unique(KG_haplo@ped$super.population),sep=":"), dist = rep(KG_haplo@snps$dist[which(KG_haplo@snps$pos%in%(key_pos-1))],each=nb_tag), pos = barcode , A1 = "G", A2 = "T",quality = NA, filter = NA, stringsAsFactors = FALSE)


bed2<-matrix(0,nrow(ped),nrow(snp2))

for (j in 1:nb_tag){
w<-seq(j,nb_tag*(length(key_pos)-1)+j,nb_tag)
pop<-unique(KG_haplo@ped$population)[j]
w2<-which(KG_haplo@ped$population==pop)
bed2[w2,w]<-1
}

## five super population version 
#for (j in 1:nb_tag){
#w<-seq(j,nb_tag*(length(key_pos)-1)+j,nb_tag)
#pop<-unique(KG_haplo@ped$super.population)[j]
#w2<-which(KG_haplo@ped$super.population==pop)
#bed2[w2,w]<-1
#}


synthetic_var_bedmat <- as.bed.matrix(x=bed2,fam=ped,bim=snp2)

KG_haplo_tagged <- cbind(KG_haplo,synthetic_var_bedmat)
m <- match(sort(KG_haplo_tagged@snps$pos),KG_haplo_tagged@snps$pos)
KG_haplo_tagged<-KG_haplo_tagged[,m]
KG_haplo_tagged@snps$pos <- as.integer(KG_haplo_tagged@snps$pos)
KG_haplo_tagged <- set.stats_haplotype(KG_haplo_tagged)

###Store the tagged reference panel 
write.bed.matrix(KG_haplo_tagged, paste("/WORKING_DIRECTORY/aherzig/tdekNOV22/haplo_tagged_chr",chr,sep=""))


}


