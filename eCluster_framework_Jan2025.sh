#!/bin/sh

#  eCluster_framework_Jan2025.sh
#  
#
#  Created by MolGenLab on 19.1.2025.
#  

########################################################################
####################### INFO for running: ##############################

### The enhancer coordinates can have distinct names. We might not have names at all (V1, V2, V3 or similar), or e.g. in the original dog DESeq2 file, they were:
# Chr_pl
# Start_pl
# End_pl

# I have changed these names to be:
# chr
# eStart
# eEnd

### The code checks, whether these names are present in the input file.
### If yes, nothing is done. If not columns 1:3 are named "chr", "eStart", "eEnd".
### Hence, the input files has to either:
###      i) have correct naming (anywhere in the columns),
###  or ii) have enhancer coordinates in cols 1:3.

### I'm planning to place that updated script to GitHub, so that we can state in the dog manuscript to provide a framework for calling eClusters & the individual enhnancers within them from PRO-seq data.

########################################################################



### FIRST PART, lines 50 - 182 ##
# You get eClusters in rows 165 - 166.
# We could provide the code only until this part.

### SECOND PART, below lines 166: ##
# Adds a threshold (here eRPK > 20) for the called clusters.
# Analyses the heat-induced change in eClusters. It requires the file to have eRNA values from the control and treatment conditions.
# These will work if the naming matches. We should harmonize these between our sample sets.








########################################################################
########### ENHANCER CLUSTERs #######

path_eClu="/Users/molgenlab/DH82_PROseq/eClusters_framework/"

#eData = read.table(paste(path_eClu,"DH82_CanFam_HS30min_vs_C_enhancers.txt",sep=""), header=T) ## DESeq2 data file provided by Samu.
eData = read.table(paste(path_eClu,"DH82_CanFam_HS30min_vs_C_enhancers_DESeq2.txt",sep=""), header=T) ## DESeq2 data file provided by Samu.



################################################
######## checking for correct column names ########

if (!("chr" %in% names(eData))) {
  # If not, rename the second column to "nameX"
  names(eData)[1] <- "chr"
}

if (!("eStart" %in% names(eData))) {
  # If not, rename the second column to "nameX"
  names(eData)[2] <- "eStart"
}

if (!("eEnd" %in% names(eData))) {
  # If not, rename the second column to "nameX"
  names(eData)[3] <- "eEnd"
}

################################################


eData = eData[order(eData$chr, eData$eStart),] ### important to order the enhancers based on coordinates.

write.table(eData, file=paste(path_eClu,"eData_noHeader.txt",sep=""), col.names=F, row.names=F, sep="\t", quote=F)

######################################################
###### Define the parameters for eClusters

initial_window <- 12500
extension_window <- 2000
min_enhancers <- 5
###################################################



# Function to find enhancer clusters
find_clusters <- function(cluData, init_win, ext_win, min_enh) {
  clusters <- list()
  start_pos <- 1
  n <- nrow(cluData)


while (start_pos <= n) {
    current_cluster <- cluData[start_pos, , drop = FALSE]
    cluster_end <- cluData$eEnd[start_pos] + init_win
    i <- start_pos + 1
    while (i <= n) {
      if (cluData$chr[i] != cluData$chr[start_pos]) break
      if (cluData$eStart[i] <= cluster_end) {
        current_cluster <- rbind(current_cluster, cluData[i, , drop = FALSE])
        cluster_end <- max(cluster_end, cluData$eEnd[i] + ext_win)
      } else {
        break
      }
      i <- i + 1
    }
    if (nrow(current_cluster) >= min_enh) {
      clusters <- append(clusters, list(current_cluster))
    }
    start_pos <- start_pos + 1
  }
  return(clusters)
}

enhancer_clusters <- find_clusters(eData, initial_window, extension_window, min_enhancers)


preClust = as.data.frame = enhancer_clusters[[1]]
preClust$clusterName = "pre1"
preClust$countEnhancers = as.numeric(0)
preClust$clusterStart = as.numeric(0)
preClust$clusterEnd = as.numeric(0)
preClust$clusterLength = preClust[nrow(preClust),"eEnd"] - preClust[1,"eStart"]
preClust$clusterCoords=""


preClust = preClust[1,]


for (cluster in enhancer_clusters) {
cluster$clusterName = ""
cluster$countEnhancers=nrow(cluster)
cluster$clusterStart = as.numeric(cluster[1,"eStart"])
cluster$clusterEnd = as.numeric(cluster[nrow(cluster),"eEnd"])
cluster$clusterLength = cluster[nrow(cluster),"eEnd"] - cluster[1,"eStart"]
cluster$clusterCoords=paste(cluster$chr,":",cluster[1,"eStart"],"-",cluster[nrow(cluster),"eEnd"],sep="")
preClust = rbind(preClust,cluster)
}


preClust = preClust[2:nrow(preClust),]


dim(preClust)
#[1] 10000    56
length(unique(preClust$clusterCoords))
#[1] 1789



preClust_u = preClust[!duplicated(preClust$clusterCoords),]
preClust_u = preClust_u[order(preClust_u$chr, preClust_u$clusterStart),]

write.table(preClust_u[,c("chr", "clusterStart", "clusterEnd", "clusterCoords")], file=paste(path_eClu,"preliminary_eClusters.txt", sep=""), col.names=F, row.names=F, quote=F, sep="\t")

#######################################################
## In these preliminary clusters, one enhancer can be part of multiple clusters. In other words, multiple clusters are called over a region with more than the required 5 enhancers.
## To unify these overlapping clusters, the next step is to conduct bedools merge. I use the default -d 0 . This means that clusters that are directly adjacent (0 nt gap) or overlap with at least 1 nt are merged into a single eCluster.

## Following the merge, the eCluster coordinates are intersected with enhancer coordinates.
## This assigns enhancers to each of the eCluster.

## This is the only place where I hop from R to shell. Perhaps in the future refined code, the merging could be done in R. The list there already contains the individual enhancers.

#################
### in shell: ###
################# ################# ################# ################# ################# ################# ################# ################# ################# #################

#cd to your working directory

bedtools merge -i preliminary_eClusters.txt > mergedClusters.bed       ### Reports eClusters. (The code could end here => eCluster coordinates)
bedtools intersect -wa -loj -a mergedClusters.bed -b eData_noHeader.txt > eClusters_withEnhancers.bed                   ### Assigns individual enhancers to each cluster (eClusters + their enhancers).

################# ################# ################# ################# ################# ################# ################# ################# ################# #################




#################
### back in R: ##
#################

########## Define here if an additional limit, eRNA minimum, is added for the 5 enahcners that consitute a cluster. I've used eRPK 20 for both pl and mn strand in the dog manuscript.
eRPK_limit = 20
##########


clust = read.table("~/DH82_PROseq/eClusters_framework/eClusters_withEnhancers.bed")

names(clust)=c("Chr", "clusterStart", "clusterEnd", names(eData))
clust$clusterCoords = paste(clust$chr,":", clust$clusterStart, "-", clust$clusterEnd, sep="")


dim(clust)
length(unique(clust$clusterCoords))

clusterList =unique(clust$clusterCoords)



preDF = clust[1,]

preDF$eCount = as.integer(0)
preDF$eCount_eRPKlimit = as.integer(0)
preDF$eCluSUM_C= as.numeric(0)
preDF$eCluMEAN_C= as.numeric(0)
preDF$eCluStDev_C= as.numeric(0)
preDF$eCluSUM_HS= as.numeric(0)
preDF$eCluMEAN_HS= as.numeric(0)
preDF$eCluStDev_HS= as.numeric(0)

preDF_ = preDF


for(i in 1:length(clusterList)) {
mClu = clusterList[i]
a = subset(clust, clusterCoords==mClu)
a_ = subset(a, gbRPK_C_pl >eRPK_limit & gbRPK_C_mn >eRPK_limit) ### eRPK_limit is taken in here.  I'm asking each enhancer in the cluster to have at least 20 eRPK on both pl and mn strands. Can be adjusted above.

a$eCount=nrow(a)
a$eCount_eRPKlimit = nrow(a_)
a$eCluSUM_C=sum(a$gbRPK_C_pl)+sum(a$gbRPK_C_mn)
a$eCluMEAN_C=(sum(a$gbRPK_C_pl)+sum(a$gbRPK_C_mn))/nrow(a)
a$eCluStDev_C=sd(c(a$gbRPK_C_pl, a$gbRPK_C_mn))

a$eCluSUM_HS=sum(a$gbRPK_HS30min_pl)+sum(a$gbRPK_HS30min_mn)
a$eCluMEAN_HS=(sum(a$gbRPK_HS30min_pl)+sum(a$gbRPK_HS30min_mn))/nrow(a)
a$eCluStDev_HS=sd(c(a$gbRPK_HS30min_pl, a$gbRPK_HS30min_mn))

preDF = rbind(preDF, a)

}

mClu = preDF[2:(nrow(preDF)),]

mClu$eClusterLength = mClu$clusterEnd - mClu$clusterStart +1


mClusters = subset(mClu, eCount_eRPKlimit>=5)


### In essence, I remove enhancers that have lower than eRPK 20 on pl or mn strand.
mClusters = subset(mClusters, gbRPK_C_pl>eRPK_limit & gbRPK_C_mn>eRPK_limit)  ### eRPK_limit is taken in here.


length(unique(mClusters$clusterCoords))
#[1] 404 ### nmber of eClusters


dim(mClusters)
#[1] 2482   61      ### 2482 enhcaners in the 404 clusters.




##### RE-DEFINING THE CLUSTER COORDINATES
### Not sure if this is actually required, seems that the coordinates were mostly the same before and after.
### Anyways, adds a unique cluster name as well.

refinedClusterList =unique(mClusters$clusterCoords)

preDF2 = mClusters[1,]
preDF2$eClusterStart=0
preDF2$eClusterEnd=0
preDF$eClusterLength=0
preDF2$eClusterCoordinates=""
preDF2$eCluster=""


for(i in 1:length(refinedClusterList)) {
eC = refinedClusterList[i]
cluster = subset(mClusters, clusterCoords==eC)

cluster$eClusterStart = as.numeric(cluster[1,"eStart"])
cluster$eClusterEnd = as.numeric(cluster[nrow(cluster),"eEnd"])
cluster$eClusterLength = cluster[nrow(cluster),"eEnd"] - cluster[1,"eStart"]
cluster$eClusterCoordinates=paste(cluster$chr,":",cluster[1,"eStart"],"-",cluster[nrow(cluster),"eEnd"],sep="")
cluster$eCluster = paste("DH82_canFam6_eCluster",i,sep="")

preDF2 = rbind(preDF2,cluster)

}


eCluster = preDF2[2:(nrow(preDF2)),]




########################################################################################################################################################
#### Here, I'm selecting just some of the columns. Some of the names are kinda specifc for the dog samples. We should harmonize the naming at some point.
########################################################################################################################################################


eCluster = eCluster[,c("chr", "eClusterStart", "eClusterEnd",  "eClusterCoordinates", "eCluster", "eCount_eRPKlimit", "eStart", "eEnd", "dTRE_name", "ppRPK_C_pl", "ppRPK_C_mn", "gbRPK_C_pl", "gbRPK_C_mn", "ppRPK_HS30min_pl", "ppRPK_HS30min_mn", "gbRPK_HS30min_pl", "gbRPK_HS30min_mn", "Reg_pl_mn_HS30min_to_C", "eRNA", "eRNA_C", "eRNA_HS", "elog2FC",  "eCount", "eCluSUM_C", "eCluMEAN_C", "eCluStDev_C", "eCluSUM_HS", "eCluMEAN_HS", "eCluStDev_HS", "eClusterLength")]

names(eCluster) = c("chr", "eClusterStart", "eClusterEnd",  "eClusterCoordinates", "eCluster", "eCount", "eStart", "eEnd", "dTRE_name", "ppRPK_C_pl", "ppRPK_C_mn", "gbRPK_C_pl", "gbRPK_C_mn", "ppRPK_HS30min_pl", "ppRPK_HS30min_mn", "gbRPK_HS30min_pl", "gbRPK_HS30min_mn", "enhancerReg_HS30_to_C", "eRNA", "eRNA_C", "eRNA_HS", "elog2FC",  "eCount_noRPKLimit", "eCluSUM_C", "eCluMEAN_C", "eCluStDev_C", "eCluSUM_HS", "eCluMEAN_HS", "eCluStDev_HS", "eClusterLength")


write.table(eCluster, file=paste(path_eClu,"DH82_eClusters_wClusterFrameWork.txt", sep=""), col.names=T, row.names=F, quote=F, sep="\t")


##############################################################################################################
# In the data frame, each row is an enhancer that belogs to a cluster.
# Which cluster it belongs to is indicated in eCluster and eClusterCoordinates columns.
# Some of the columns show info for the individual enhancers, others for the cluster the enhancer belongs to.
##############################################################################################################

dim(eCluster)
# [1] 2482   30     # number of enhancers that are part of an eCluster

length(unique(eCluster$eCluster))
# [1] 404           # number of eClusters

table(eCluster$eCount)
#   5    6    7    8    9   10   11   12   13   14   15  ## number of enhancers in a cluster
# 1010  552  322  224  108  140  44   12   26   14   30  ## total number of enhancers


table(eCluster$eCount)/as.numeric(names(table(eCluster$eCount)))
#   5   6   7   8   9  10  11  12  13  14  15           ## number of enhancers in a cluster
# 202  92  46  28  12  14   4   1   2   1   2           ## number of eClusters




##########################
### Generating the xy graph in Fig. 7D (likely becomes 7B).
##############################################################################################################


eCluster_plot = eCluster

eCluster_plot$eCluMEAN_dRPK = eCluster_plot$eCluMEAN_HS - eCluster_plot$eCluMEAN_C
eCluster_plot$eCluMEAN_log2FC = log2(eCluster_plot$eCluMEAN_HS/eCluster_plot$eCluMEAN_C)

#eCluster_plot_ = eCluster_plot_[order(abs(eCluster_plot_$eCluMEAN_log2FC)),]
eCluster_plot_increased = subset(eCluster_plot, eCluMEAN_log2FC>0)
eCluster_plot_decreased = subset(eCluster_plot, eCluMEAN_log2FC<0)



plot(eCluster_plot_increased$elog2FC, eCluster_plot_increased$eCluMEAN_log2FC, ylim=c(-1.6,1.6), xlim=c(-4.5,4.5), pch=21, col=rgb(1,0,0,0.5), bg=rgb(1,0,0,0.2), cex=0.6)
points(eCluster_plot_decreased$elog2FC, eCluster_plot_decreased$eCluMEAN_log2FC, ylim=c(-1.6,1.6), xlim=c(-4.5,4.5), pch=21, col=rgb(0,0,1,0.5), bg=rgb(0,0,1,0.2), cex=0.6)

pdf(paste(path_eClu,"eCluster_individual_enhancers_change_uponHS.pdf", sep=""))
plot(eCluster_plot_increased$elog2FC, eCluster_plot_increased$eCluMEAN_log2FC, ylim=c(-1.6,1.6), xlim=c(-4.5,4.5), pch=21, col=rgb(1,0,0,0.5), bg=rgb(1,0,0,0.2), cex=0.7)
points(eCluster_plot_decreased$elog2FC, eCluster_plot_decreased$eCluMEAN_log2FC, ylim=c(-1.6,1.6), xlim=c(-4.5,4.5), pch=21, col=rgb(0,0,1,0.5), bg=rgb(0,0,1,0.2), cex=0.7)
dev.off()


cor.test(eCluster_plot$elog2FC, eCluster_plot$eCluMEAN_log2FC, method="pearson")
cor.test(eCluster_plot$elog2FC, eCluster_plot$eCluMEAN_log2FC, method="spearman")


##############################################################################################################




### I will make a new script for the sequence similarity analyses.
