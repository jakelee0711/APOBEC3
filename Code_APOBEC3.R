setwd("~/Documents/SV_Enhancers/00_Analysis")
library(dplyr)
library(stringr)
library(data.table)

# Subset the amplicon data based on the mechanisms -- including FBI and TBA
df <- amphmf[((amphmf$mech_s == "fold_back_inversion" & amphmf$mcj_s != "no" & amphmf$mech_e != "translocation") | (amphmf$mech_e == "fold_back_inversion" & amphmf$mcj_e != "no" & amphmf$mech_s != "translocation")) & (amphmf$end - amphmf$start > 200000),]

# Mutations at the vicinty of the major fold-back inversions
specorder <- c("A[C>A]A", "A[C>A]C", "A[C>A]G", "A[C>A]T", "C[C>A]A", "C[C>A]C", "C[C>A]G", "C[C>A]T", "G[C>A]A", "G[C>A]C", "G[C>A]G", "G[C>A]T", "T[C>A]A", "T[C>A]C", "T[C>A]G", "T[C>A]T", 
               "A[C>G]A", "A[C>G]C", "A[C>G]G", "A[C>G]T", "C[C>G]A", "C[C>G]C", "C[C>G]G", "C[C>G]T", "G[C>G]A", "G[C>G]C", "G[C>G]G", "G[C>G]T", "T[C>G]A", "T[C>G]C", "T[C>G]G", "T[C>G]T", 
               "A[C>T]A", "A[C>T]C", "A[C>T]G", "A[C>T]T", "C[C>T]A", "C[C>T]C", "C[C>T]G", "C[C>T]T", "G[C>T]A", "G[C>T]C", "G[C>T]G", "G[C>T]T", "T[C>T]A", "T[C>T]C", "T[C>T]G", "T[C>T]T",
               "A[T>A]A", "A[T>A]C", "A[T>A]G", "A[T>A]T", "C[T>A]A", "C[T>A]C", "C[T>A]G", "C[T>A]T", "G[T>A]A", "G[T>A]C", "G[T>A]G", "G[T>A]T", "T[T>A]A", "T[T>A]C", "T[T>A]G", "T[T>A]T",
               "A[T>C]A", "A[T>C]C", "A[T>C]G", "A[T>C]T", "C[T>C]A", "C[T>C]C", "C[T>C]G", "C[T>C]T", "G[T>C]A", "G[T>C]C", "G[T>C]G", "G[T>C]T", "T[T>C]A", "T[T>C]C", "T[T>C]G", "T[T>C]T",
               "A[T>G]A", "A[T>G]C", "A[T>G]G", "A[T>G]T", "C[T>G]A", "C[T>G]C", "C[T>G]G", "C[T>G]T", "G[T>G]A", "G[T>G]C", "G[T>G]G", "G[T>G]T", "T[T>G]A", "T[T>G]C", "T[T>G]G", "T[T>G]T")
spec_vici <- data.frame(matrix(vector(), 0, 98, dimnames=list(c(), c("gw", specorder, "n_tumor"))), stringsAsFactors=F)

par(mar=c(5.1, 10.1, 4.1, 2.1))
#genomic_window <- 10
#genomic_window <- 50
#genomic_window <- 100
#genomic_window <- 1000
#genomic_window <- 10000
#genomic_window <- 100000
for (genomic_window in c(10, 30, 50, 100, 1000, 10000, 100000)){
  print(paste0("genomic window size = ", genomic_window))
  df <- amphmf[((amphmf$mech_s == "fold_back_inversion" & amphmf$mcj_s != "no" & amphmf$mech_e != "translocation") | (amphmf$mech_e == "fold_back_inversion" & amphmf$mcj_e != "no" & amphmf$mech_s != "translocation")) & (amphmf$end - amphmf$start > 200000),]
  for (w in unique(df$study_id)){
    vecdf <- read.csv(paste0("../27_EBI/final_data_022222/timing_veclonal/", w, ".sage30.snv.timeR.VEClonal.txt"), header=T, as.is=T, sep="\t")
    vecdf$study_id <- w
    if (w == unique(df$study_id)[1]){
      vici <- data.frame(matrix(vector(), 0, ncol(vecdf), dimnames=list(c(), colnames(vecdf))), stringsAsFactors=F)
    }
    for (i in which(df$study_id == w)){
      if (df$mech_s[i] == "fold_back_inversion" & df$mcj_s[i] != "no"){
        vici <- rbind(vici, vecdf[vecdf$chrom == df$chromosome[i] & vecdf$pos >= df$start[i] & vecdf$pos <= df$start[i]+genomic_window,])
      }
      if (df$mech_e[i] == "fold_back_inversion" & df$mcj_e[i] != "no"){
        vici <- rbind(vici, vecdf[vecdf$chrom == df$chromosome[i] & vecdf$pos <= df$end[i] & vecdf$pos >= df$end[i]-genomic_window,])
      }
    }
  }
  vici$mutspectra <- factor(vici$mutspectra, levels = specorder)
  spec_vici[(nrow(spec_vici)+1),] <- NA
  spec_vici$gw[nrow(spec_vici)] <- genomic_window
  spec_vici[nrow(spec_vici),c(2:97)] <- as.vector(table(vici$mutspectra))
  spec_vici$n_tumor[nrow(spec_vici)] <- length(unique(vici$study_id))
  
  pdf(paste0("FIP_revision/Vicinity_mutspectra/BFB/mutspectra.amphmf.boundary.fold_back.vicinity.", genomic_window, "bp.pdf"), height=4, width=8)
  barplot(as.vector(t(spec_vici[spec_vici$gw == genomic_window,c(2:97)])), border=NA, las=1, col=sigcolor, space=0.8, main=paste0("# mutations = ", nrow(vici), " # tumors = ", length(unique(vici$study_id))), ylab="Number of mutations", xlab="Trinucleotide context")
  dev.off()
}
bfb_vici <- vici

# Mutations at the vicinty of the major translocations
specorder <- c("A[C>A]A", "A[C>A]C", "A[C>A]G", "A[C>A]T", "C[C>A]A", "C[C>A]C", "C[C>A]G", "C[C>A]T", "G[C>A]A", "G[C>A]C", "G[C>A]G", "G[C>A]T", "T[C>A]A", "T[C>A]C", "T[C>A]G", "T[C>A]T", 
               "A[C>G]A", "A[C>G]C", "A[C>G]G", "A[C>G]T", "C[C>G]A", "C[C>G]C", "C[C>G]G", "C[C>G]T", "G[C>G]A", "G[C>G]C", "G[C>G]G", "G[C>G]T", "T[C>G]A", "T[C>G]C", "T[C>G]G", "T[C>G]T", 
               "A[C>T]A", "A[C>T]C", "A[C>T]G", "A[C>T]T", "C[C>T]A", "C[C>T]C", "C[C>T]G", "C[C>T]T", "G[C>T]A", "G[C>T]C", "G[C>T]G", "G[C>T]T", "T[C>T]A", "T[C>T]C", "T[C>T]G", "T[C>T]T",
               "A[T>A]A", "A[T>A]C", "A[T>A]G", "A[T>A]T", "C[T>A]A", "C[T>A]C", "C[T>A]G", "C[T>A]T", "G[T>A]A", "G[T>A]C", "G[T>A]G", "G[T>A]T", "T[T>A]A", "T[T>A]C", "T[T>A]G", "T[T>A]T",
               "A[T>C]A", "A[T>C]C", "A[T>C]G", "A[T>C]T", "C[T>C]A", "C[T>C]C", "C[T>C]G", "C[T>C]T", "G[T>C]A", "G[T>C]C", "G[T>C]G", "G[T>C]T", "T[T>C]A", "T[T>C]C", "T[T>C]G", "T[T>C]T",
               "A[T>G]A", "A[T>G]C", "A[T>G]G", "A[T>G]T", "C[T>G]A", "C[T>G]C", "C[T>G]G", "C[T>G]T", "G[T>G]A", "G[T>G]C", "G[T>G]G", "G[T>G]T", "T[T>G]A", "T[T>G]C", "T[T>G]G", "T[T>G]T")
spec_vici <- data.frame(matrix(vector(), 0, 98, dimnames=list(c(), c("gw", specorder, "n_tumor"))), stringsAsFactors=F)

par(mar=c(5.1, 10.1, 4.1, 2.1))
#genomic_window <- 10
#genomic_window <- 50
#genomic_window <- 100
#genomic_window <- 1000
#genomic_window <- 10000
#genomic_window <- 100000
for (genomic_window in c(10, 30, 50, 100, 1000, 10000, 100000)){
  print(paste0("genomic window size = ", genomic_window))
  df <- amphmf[((amphmf$mech_s == "translocation" & amphmf$mcj_s != "no" & amphmf$mech_e != "fold_back_inversion") | (amphmf$mech_e == "translocation" & amphmf$mcj_e != "no" & amphmf$mech_s != "fold_back_inversion")) & (amphmf$end - amphmf$start > 200000),]
  for (w in unique(df$study_id)){
    vecdf <- read.csv(paste0("../27_EBI/final_data_022222/timing_veclonal/", w, ".sage30.snv.timeR.VEClonal.txt"), header=T, as.is=T, sep="\t")
    vecdf$study_id <- w
    if (w == unique(df$study_id)[1]){
      vici <- data.frame(matrix(vector(), 0, ncol(vecdf), dimnames=list(c(), colnames(vecdf))), stringsAsFactors=F)
    }
    for (i in which(df$study_id == w)){
      if (df$mech_s[i] == "translocation" & df$mcj_s[i] != "no"){
        vici <- rbind(vici, vecdf[vecdf$chrom == df$chromosome[i] & vecdf$pos >= df$start[i] & vecdf$pos <= df$start[i]+genomic_window,])
      }
      if (df$mech_e[i] == "translocation" & df$mcj_e[i] != "no"){
        vici <- rbind(vici, vecdf[vecdf$chrom == df$chromosome[i] & vecdf$pos <= df$end[i] & vecdf$pos >= df$end[i]-genomic_window,])
      }
    }
  }
  vici$mutspectra <- factor(vici$mutspectra, levels = specorder)
  spec_vici[(nrow(spec_vici)+1),] <- NA
  spec_vici$gw[nrow(spec_vici)] <- genomic_window
  spec_vici[nrow(spec_vici),c(2:97)] <- as.vector(table(vici$mutspectra))
  spec_vici$n_tumor[nrow(spec_vici)] <- length(unique(vici$study_id))
  
  pdf(paste0("FIP_revision/Vicinity_mutspectra/BFB/mutspectra.amphmf.boundary.translocation.vicinity.", genomic_window, "bp.pdf"), height=4, width=8)
  barplot(as.vector(t(spec_vici[spec_vici$gw == genomic_window,c(2:97)])), border=NA, las=1, col=sigcolor, space=0.8, main=paste0("# mutations = ", nrow(vici), " # tumors = ", length(unique(vici$study_id))), ylab="Number of mutations", xlab="Trinucleotide context")
  dev.off()
}
tba_vici <- vici

# General comparison between APOBEC expression level vs APOBEC mutational signatures
df <- sigdf
df$AICDA <- NA
df$APOBEC1 <- NA
df$APOBEC2 <- NA
df$APOBEC3A <- NA
df$APOBEC3B <- NA
df$APOBEC3C <- NA
df$APOBEC3D <- NA
df$APOBEC3F <- NA
df$APOBEC3G <- NA
for (k in c("AICDA", "APOBEC1", "APOBEC2", "APOBEC3A", "APOBEC3B", "APOBEC3C", "APOBEC3D", "APOBEC3F", "APOBEC3G")){
  for (i in 1:nrow(df)){
    if (df$study_id[i] %in% colnames(ebasis)){
      if (!is.na(ebasis[which(ebasis$Name == k), which(colnames(ebasis) == df$study_id[i])])){
        df[i,which(colnames(df) == k)] <- ebasis[which(ebasis$Name == k), which(colnames(ebasis) == df$study_id[i])]
      }
    }
  }
}

for (i in 1:nrow(df)){
  if (df$study_id[i] %in% colnames(ebasis) & !is.na(ebasis[which(ebasis$Name == "AICDA"), which(colnames(ebasis) == df$study_id[i])])){
    df$AICDA[i] <- ebasis[which(ebasis$Name == "AICDA"), which(colnames(ebasis) == df$study_id[i])]
  }
}
df <- df[!is.na(df$APOBEC3D),]

pdf("Collab_APOBEC/A3A_SBS2+SBS13.pdf", height=5, width=5)
plot(log(df$SBS2[df$SBS2 != 0 & df$SBS13 != 0]+df$SBS13[df$SBS2 != 0 & df$SBS13 != 0]) ~ df$APOBEC3A[df$SBS2 != 0 & df$SBS13 != 0], xlim=c(-10,10), xlab="APOBEC3A expression (normalized)", ylab="Burden of SBS2+SBS13 (log10-transformed)", las=1, pch=19, col=rgb(0,0,0,.2), main=paste0(nrow(df[!is.na(df$APOBEC3A) & df$SBS2 != 0 & df$SBS13 != 0,]), " breast cancers (r = ", round(sqrt(summary(lm(log(df$SBS2[df$SBS2 != 0 & df$SBS13 != 0]+df$SBS13[df$SBS2 != 0 & df$SBS13 != 0]) ~ df$APOBEC3A[df$SBS2 != 0 & df$SBS13 != 0]))$adj.r.squared), 3), ")"))
abline(lm(log(df$SBS2[df$SBS2 != 0 & df$SBS13 != 0]+df$SBS13[df$SBS2 != 0 & df$SBS13 != 0]) ~ df$APOBEC3A[df$SBS2 != 0 & df$SBS13 != 0]))
dev.off()

summary(lm(log(df$SBS2[df$SBS2 != 0 & df$SBS13 != 0]+df$SBS13[df$SBS2 != 0 & df$SBS13 != 0]) ~ df$APOBEC3A[df$SBS2 != 0 & df$SBS13 != 0]))
sqrt(summary(lm(log(df$SBS2[df$SBS2 != 0 & df$SBS13 != 0]+df$SBS13[df$SBS2 != 0 & df$SBS13 != 0]) ~ df$APOBEC3A[df$SBS2 != 0 & df$SBS13 != 0]))$adj.r.squared)

pdf("Collab_APOBEC/A3B_SBS2+SBS13.pdf", height=5, width=5)
plot(log(df$SBS2[df$SBS2 != 0 & df$SBS13 != 0]+df$SBS13[df$SBS2 != 0 & df$SBS13 != 0]) ~ df$APOBEC3B[df$SBS2 != 0 & df$SBS13 != 0], xlim=c(-10,10), xlab="APOBEC3B expression (normalized)", ylab="Burden of SBS2+SBS13 (log10-transformed)", las=1, pch=19, col=rgb(0,0,0,.2), main=paste0(nrow(df[!is.na(df$APOBEC3B) & df$SBS2 != 0 & df$SBS13 != 0,]), " breast cancers (r = ", round(sqrt(summary(lm(log(df$SBS2[df$SBS2 != 0 & df$SBS13 != 0]+df$SBS13[df$SBS2 != 0 & df$SBS13 != 0]) ~ df$APOBEC3B[df$SBS2 != 0 & df$SBS13 != 0]))$adj.r.squared), 3), ")"))
abline(lm(log(df$SBS2[df$SBS2 != 0 & df$SBS13 != 0]+df$SBS13[df$SBS2 != 0 & df$SBS13 != 0]) ~ df$APOBEC3B[df$SBS2 != 0 & df$SBS13 != 0]))
dev.off()

summary(lm(log(df$SBS2[df$SBS2 != 0 & df$SBS13 != 0]+df$SBS13[df$SBS2 != 0 & df$SBS13 != 0]) ~ df$APOBEC3B[df$SBS2 != 0 & df$SBS13 != 0]))

pdf("Collab_APOBEC/A3C_SBS2+SBS13.pdf", height=5, width=5)
plot(log(df$SBS2[df$SBS2 != 0 & df$SBS13 != 0]+df$SBS13[df$SBS2 != 0 & df$SBS13 != 0]) ~ df$APOBEC3C[df$SBS2 != 0 & df$SBS13 != 0], xlim=c(-10,10), xlab="APOBEC3C expression (normalized)", ylab="Burden of SBS2+SBS13 (log10-transformed)", las=1, pch=19, col=rgb(0,0,0,.2), main=paste0(nrow(df[!is.na(df$APOBEC3C) & df$SBS2 != 0 & df$SBS13 != 0,]), " breast cancers (r = ", round(sqrt(summary(lm(log(df$SBS2[df$SBS2 != 0 & df$SBS13 != 0]+df$SBS13[df$SBS2 != 0 & df$SBS13 != 0]) ~ df$APOBEC3C[df$SBS2 != 0 & df$SBS13 != 0]))$adj.r.squared), 3), ")"))
abline(lm(log(df$SBS2[df$SBS2 != 0 & df$SBS13 != 0]+df$SBS13[df$SBS2 != 0 & df$SBS13 != 0]) ~ df$APOBEC3C[df$SBS2 != 0 & df$SBS13 != 0]))
dev.off()

pdf("Collab_APOBEC/AID_SBS2+SBS13.pdf", height=5, width=5)
plot(log(df$SBS2[df$SBS2 != 0 & df$SBS13 != 0]+df$SBS13[df$SBS2 != 0 & df$SBS13 != 0]) ~ df$AICDA[df$SBS2 != 0 & df$SBS13 != 0], xlim=c(-10,10), xlab="AID expression (normalized)", ylab="Burden of SBS2+SBS13 (log10-transformed)", las=1, pch=19, col=rgb(0,0,0,.2), main=paste0(nrow(df[!is.na(df$AICDA) & df$SBS2 != 0 & df$SBS13 != 0,]), " breast cancers (r = ", round(sqrt(summary(lm(log(df$SBS2[df$SBS2 != 0 & df$SBS13 != 0]+df$SBS13[df$SBS2 != 0 & df$SBS13 != 0]) ~ df$AICDA[df$SBS2 != 0 & df$SBS13 != 0]))$adj.r.squared), 3), ")"))
abline(lm(log(df$SBS2[df$SBS2 != 0 & df$SBS13 != 0]+df$SBS13[df$SBS2 != 0 & df$SBS13 != 0]) ~ df$AICDA[df$SBS2 != 0 & df$SBS13 != 0]))
dev.off()





df$bfb_vici <- 0
df$bfb_tcw <- 0
df$tba_vici <- 0
df$tba_tcw <- 0
for (i in 1:nrow(df)){
  df$bfb_vici[i] <- nrow(bfb_vici[bfb_vici$study_id == df$study_id[i],])
  df$bfb_tcw[i] <- nrow(bfb_vici[bfb_vici$study_id == df$study_id[i] & bfb_vici$mutspectra %in% c("T[C>T]A", "T[C>T]T", "T[C>G]A", "T[C>G]T", "T[C>A]A", "T[C>A]T"),])
  df$tba_vici[i] <- nrow(tba_vici[tba_vici$study_id == df$study_id[i],])
  df$tba_tcw[i] <- nrow(tba_vici[tba_vici$study_id == df$study_id[i] & tba_vici$mutspectra %in% c("T[C>T]A", "T[C>T]T", "T[C>G]A", "T[C>G]T", "T[C>A]A", "T[C>A]T"),])
}

pdf("Collab_APOBEC/A3A_BFB_all_mutations.pdf", height=5, width=5)
plot(log(df$bfb_vici[df$bfb_vici != 0]) ~ df$APOBEC3A[df$bfb_vici != 0], xlim=c(-10,10), ylim=c(0,3), xlab="APOBEC3A expression (normalized)", ylab="Number of mutations at the vicinity of BFB breakpoint (log10-transformed)", las=1, pch=19, col=rgb(0,0,0,.2), main=paste0(nrow(df[df$bfb_vici != 0,]), " breast cancers"))
dev.off()

pdf("Collab_APOBEC/A3B_BFB_all_mutations.pdf", height=5, width=5)
plot(log(df$bfb_vici[df$bfb_vici != 0]) ~ df$APOBEC3B[df$bfb_vici != 0], xlim=c(-10,10), ylim=c(0,3), xlab="APOBEC3B expression (normalized)", ylab="Number of mutations at the vicinity of BFB breakpoint (log10-transformed)", las=1, pch=19, col=rgb(0,0,0,.2), main=paste0(nrow(df[df$bfb_vici != 0,]), " breast cancers"))
dev.off()

pdf("Collab_APOBEC/A3A_BFB_tcw_mutations.pdf", height=5, width=5)
plot(log(df$bfb_tcw[df$bfb_tcw != 0]) ~ df$APOBEC3A[df$bfb_tcw != 0], xlim=c(-10,10), ylim=c(0,3), xlab="APOBEC3A expression (normalized)", ylab="Number of APOBEC-context mutations at the vicinity of BFB breakpoint (log10-transformed)", las=1, pch=19, col=rgb(0,0,0,.2), main=paste0(nrow(df[df$bfb_tcw != 0,]), " breast cancers"))
dev.off()

pdf("Collab_APOBEC/A3B_BFB_tcw_mutations.pdf", height=5, width=5)
plot(log(df$bfb_tcw[df$bfb_tcw != 0]) ~ df$APOBEC3B[df$bfb_tcw != 0], xlim=c(-10,10), ylim=c(0,3), xlab="APOBEC3B expression (normalized)", ylab="Number of APOBEC-context mutations at the vicinity of BFB breakpoint (log10-transformed)", las=1, pch=19, col=rgb(0,0,0,.2), main=paste0(nrow(df[df$bfb_tcw != 0,]), " breast cancers"))
dev.off()


pdf("Collab_APOBEC/A3A_TBA_all_mutations.pdf", height=5, width=5)
plot(log(df$tba_vici[df$tba_vici != 0]) ~ df$APOBEC3A[df$tba_vici != 0], xlim=c(-10,10), ylim=c(0,3), xlab="APOBEC3A expression (normalized)", ylab="Number of mutations at the vicinity of TBA breakpoint (log10-transformed)", las=1, pch=19, col=rgb(0,0,0,.2), main=paste0(nrow(df[df$tba_vici != 0,]), " breast cancers"))
dev.off()

pdf("Collab_APOBEC/A3B_TBA_all_mutations.pdf", height=5, width=5)
plot(log(df$tba_vici[df$tba_vici != 0]) ~ df$APOBEC3B[df$tba_vici != 0], xlim=c(-10,10), ylim=c(0,3), xlab="APOBEC3B expression (normalized)", ylab="Number of mutations at the vicinity of TBA breakpoint (log10-transformed)", las=1, pch=19, col=rgb(0,0,0,.2), main=paste0(nrow(df[df$tba_vici != 0,]), " breast cancers"))
dev.off()

pdf("Collab_APOBEC/A3A_TBA_tcw_mutations.pdf", height=5, width=5)
plot(log(df$tba_tcw[df$tba_tcw != 0]) ~ df$APOBEC3A[df$tba_tcw != 0], xlim=c(-10,10), ylim=c(0,3), xlab="APOBEC3A expression (normalized)", ylab="Number of APOBEC-context mutations at the vicinity of TBA breakpoint (log10-transformed)", las=1, pch=19, col=rgb(0,0,0,.2), main=paste0(nrow(df[df$tba_tcw != 0,]), " breast cancers"))
dev.off()

pdf("Collab_APOBEC/A3B_TBA_tcw_mutations.pdf", height=5, width=5)
plot(log(df$tba_tcw[df$tba_tcw != 0]) ~ df$APOBEC3B[df$tba_tcw != 0], xlim=c(-10,10), ylim=c(0,3), xlab="APOBEC3B expression (normalized)", ylab="Number of APOBEC-context mutations at the vicinity of TBA breakpoint (log10-transformed)", las=1, pch=19, col=rgb(0,0,0,.2), main=paste0(nrow(df[df$tba_tcw != 0,]), " breast cancers"))
dev.off()

# Extracting 5p and 3p bases
opposite.strand <- function(inputbase){
  if (inputbase == "A"){
    paste0("T")
  } else if (inputbase == "T"){
    paste0("A")
  } else if (inputbase == "C"){
    paste0("G")
  } else if (inputbase == "G"){
    paste0("C")
  } else if (inputbase == "N"){
    paste0("N")
  }
}

ext.TN <- function(chromname, poscoordinate, refbase, altbase){
  if (toupper(str_sub(fread(paste0("../../References/hg19_bychr/chr", chromname, ".fa"), skip = floor(poscoordinate/50) + 1, nrows = 1, header = F), start = poscoordinate - floor(poscoordinate/50)*50 - 0, end = poscoordinate - floor(poscoordinate/50)*50 + 0)) == "N"){
    paste0("not_available")
  } else {
    if (refbase != toupper(str_sub(fread(paste0("../../References/hg19_bychr/chr", chromname, ".fa"), skip = floor(poscoordinate/50) + 1, nrows = 1, header = F), start = poscoordinate - floor(poscoordinate/50)*50 - 0, end = poscoordinate - floor(poscoordinate/50)*50 + 0))){
      paste0("not_available")
    } else {
      if (refbase %in% c("C", "T")){
        tnc <- toupper(str_sub(fread(paste0("../../References/hg19_bychr/chr", chromname, ".fa"), skip = floor(poscoordinate/50) + 1, nrows = 1, header = F), start = poscoordinate - floor(poscoordinate/50)*50 - 1, end = poscoordinate - floor(poscoordinate/50)*50 + 1))
        paste0(str_sub(tnc, 1, 1), "[", refbase, ">", altbase, "]", str_sub(tnc, 3, 3))
      } else {
        tnc <- toupper(str_sub(fread(paste0("../../References/hg19_bychr/chr", chromname, ".fa"), skip = floor(poscoordinate/50) + 1, nrows = 1, header = F), start = poscoordinate - floor(poscoordinate/50)*50 - 1, end = poscoordinate - floor(poscoordinate/50)*50 + 1))
        paste0(opposite.strand(str_sub(tnc, 3, 3)), "[", opposite.strand(refbase), ">", opposite.strand(altbase), "]", opposite.strand(str_sub(tnc, 1, 1)))
      }
    }
  }
}

i <- 5
ext.TN(bfb_vici$chrom[i], bfb_vici$pos[i], bfb_vici$ref[i], bfb_vici$alt[i])

# Pan-cancer breakpoint association
tumortype <- c("Panc-AdenoCA", "Eso-AdenoCA")
tumortype <- c("Liver-HCC", "Prost-AdenoCA", "CNS-Medullo", "Ovary-AdenoCA", "Lymph-BNHL", "Skin-Melanoma")
tumortype <- c("Lymph-CLL", "CNS-PiloAstro", "Stomach-AdenoCA", "Head-SCC", "ColoRect-AdenoCA", "Lung-SCC", "Uterus-AdenoCA", "Lung-AdenoCA", "CNS-GBM", "Biliary-AdenoCA")
tumortype <- c("Breast-AdenoCA", "SoftTissue-Liposarc")
# Mutations at the vicinty of the major fold-back inversions

par(mar=c(5.1, 10.1, 4.1, 2.1))
#genomic_window <- 10
#genomic_window <- 50
#genomic_window <- 100
#genomic_window <- 1000
#genomic_window <- 10000
#genomic_window <- 100000

specorder <- c("A[C>A]A", "A[C>A]C", "A[C>A]G", "A[C>A]T", "C[C>A]A", "C[C>A]C", "C[C>A]G", "C[C>A]T", "G[C>A]A", "G[C>A]C", "G[C>A]G", "G[C>A]T", "T[C>A]A", "T[C>A]C", "T[C>A]G", "T[C>A]T", 
               "A[C>G]A", "A[C>G]C", "A[C>G]G", "A[C>G]T", "C[C>G]A", "C[C>G]C", "C[C>G]G", "C[C>G]T", "G[C>G]A", "G[C>G]C", "G[C>G]G", "G[C>G]T", "T[C>G]A", "T[C>G]C", "T[C>G]G", "T[C>G]T", 
               "A[C>T]A", "A[C>T]C", "A[C>T]G", "A[C>T]T", "C[C>T]A", "C[C>T]C", "C[C>T]G", "C[C>T]T", "G[C>T]A", "G[C>T]C", "G[C>T]G", "G[C>T]T", "T[C>T]A", "T[C>T]C", "T[C>T]G", "T[C>T]T",
               "A[T>A]A", "A[T>A]C", "A[T>A]G", "A[T>A]T", "C[T>A]A", "C[T>A]C", "C[T>A]G", "C[T>A]T", "G[T>A]A", "G[T>A]C", "G[T>A]G", "G[T>A]T", "T[T>A]A", "T[T>A]C", "T[T>A]G", "T[T>A]T",
               "A[T>C]A", "A[T>C]C", "A[T>C]G", "A[T>C]T", "C[T>C]A", "C[T>C]C", "C[T>C]G", "C[T>C]T", "G[T>C]A", "G[T>C]C", "G[T>C]G", "G[T>C]T", "T[T>C]A", "T[T>C]C", "T[T>C]G", "T[T>C]T",
               "A[T>G]A", "A[T>G]C", "A[T>G]G", "A[T>G]T", "C[T>G]A", "C[T>G]C", "C[T>G]G", "C[T>G]T", "G[T>G]A", "G[T>G]C", "G[T>G]G", "G[T>G]T", "T[T>G]A", "T[T>G]C", "T[T>G]G", "T[T>G]T")

# BFB
for (j in tumortype){
  print(paste0("processing ", j))
  for (genomic_window in c(10, 30, 50, 100, 1000, 10000, 100000)){
    print(paste0("genomic window size = ", genomic_window))
    df <- ampcawg[((ampcawg$mech_s == "fold_back_inversion" & ampcawg$mcj_s != "no" & ampcawg$mech_e != "translocation") | (ampcawg$mech_e == "fold_back_inversion" & ampcawg$mcj_e != "no" & ampcawg$mech_s != "translocation")) & (ampcawg$end - ampcawg$start > 200000) & ampcawg$histology_abbreviation == j,]
    if (nrow(df) != 0){
      for (w in unique(df$icgc_donor_id)){
        snvdf <- read.csv(paste0("../../ICGC_PCAWG_calls/final_consensus_12aug_passonly/snv_mnv/01.1_Parsed_txt/r01_", sumdf$tumor[sumdf$icgc_donor_id == w], ".consensus.snv_mnv.parsed.txt"), header=T, as.is=T, sep="\t")
        snvdf$icgc_donor_id <- w
        if (w == unique(df$icgc_donor_id)[1]){
          vici <- data.frame(matrix(vector(), 0, ncol(snvdf), dimnames=list(c(), colnames(snvdf))), stringsAsFactors=F)
        }
        for (i in which(df$icgc_donor_id == w)){
          if (df$mech_s[i] == "fold_back_inversion" & df$mcj_s[i] != "no"){
            vici <- rbind(vici, snvdf[snvdf$X.CHR == df$chromosome[i] & snvdf$POS >= df$start[i] & snvdf$POS <= df$start[i]+genomic_window,])
          }
          if (df$mech_e[i] == "fold_back_inversion" & df$mcj_e[i] != "no"){
            vici <- rbind(vici, snvdf[snvdf$X.CHR == df$chromosome[i] & snvdf$POS <= df$end[i] & snvdf$POS >= df$end[i]-genomic_window,])
          }
        }
      }
      if (nrow(vici) != 0){
        vici$mutspectra <- NA
        for (i in 1:nrow(vici)){
          vici$mutspectra[i] <- ext.TN(vici$X.CHR[i], vici$POS[i], vici$REF[i], vici$VAR[i])
        }
        vici$mutspectra <- factor(vici$mutspectra, levels = specorder)
        spec_vici <- data.frame(matrix(vector(), 0, 98, dimnames=list(c(), c("gw", specorder, "n_tumor"))), stringsAsFactors=F)
        spec_vici[(nrow(spec_vici)+1),] <- NA
        spec_vici$gw[nrow(spec_vici)] <- genomic_window
        spec_vici[nrow(spec_vici),c(2:97)] <- as.vector(table(vici$mutspectra))
        spec_vici$n_tumor[nrow(spec_vici)] <- length(unique(vici$icgc_donor_id))
        
        pdf(paste0("Collab_APOBEC/Mutspectra/PCAWG/mutspectra.ampcawg.", j, ".boundary.fold_back.vicinity.", genomic_window, "bp.pdf"), height=4, width=8)
        barplot(as.vector(t(spec_vici[spec_vici$gw == genomic_window,c(2:97)])), border=NA, las=1, col=sigcolor, space=0.8, main=paste0("# mutations = ", nrow(vici), " # tumors = ", length(unique(vici$icgc_donor_id))), ylab="Number of mutations", xlab="Trinucleotide context")
        dev.off()
        rm(spec_vici)
        rm(vici)
      }
    }
  }
}


# TB amplification
tumortype <- c("Panc-AdenoCA", "Eso-AdenoCA")
tumortype <- c("Liver-HCC", "Prost-AdenoCA", "CNS-Medullo", "Ovary-AdenoCA", "Lymph-BNHL", "Skin-Melanoma")
tumortype <- c("Lymph-CLL", "CNS-PiloAstro", "Stomach-AdenoCA", "Head-SCC", "ColoRect-AdenoCA", "Lung-SCC", "Uterus-AdenoCA", "Lung-AdenoCA", "CNS-GBM", "Biliary-AdenoCA")
tumortype <- c("Breast-AdenoCA", "SoftTissue-Liposarc")
# Mutations at the vicinty of the major fold-back inversions

par(mar=c(5.1, 10.1, 4.1, 2.1))
#genomic_window <- 10
#genomic_window <- 50
#genomic_window <- 100
#genomic_window <- 1000
#genomic_window <- 10000
#genomic_window <- 100000

specorder <- c("A[C>A]A", "A[C>A]C", "A[C>A]G", "A[C>A]T", "C[C>A]A", "C[C>A]C", "C[C>A]G", "C[C>A]T", "G[C>A]A", "G[C>A]C", "G[C>A]G", "G[C>A]T", "T[C>A]A", "T[C>A]C", "T[C>A]G", "T[C>A]T", 
               "A[C>G]A", "A[C>G]C", "A[C>G]G", "A[C>G]T", "C[C>G]A", "C[C>G]C", "C[C>G]G", "C[C>G]T", "G[C>G]A", "G[C>G]C", "G[C>G]G", "G[C>G]T", "T[C>G]A", "T[C>G]C", "T[C>G]G", "T[C>G]T", 
               "A[C>T]A", "A[C>T]C", "A[C>T]G", "A[C>T]T", "C[C>T]A", "C[C>T]C", "C[C>T]G", "C[C>T]T", "G[C>T]A", "G[C>T]C", "G[C>T]G", "G[C>T]T", "T[C>T]A", "T[C>T]C", "T[C>T]G", "T[C>T]T",
               "A[T>A]A", "A[T>A]C", "A[T>A]G", "A[T>A]T", "C[T>A]A", "C[T>A]C", "C[T>A]G", "C[T>A]T", "G[T>A]A", "G[T>A]C", "G[T>A]G", "G[T>A]T", "T[T>A]A", "T[T>A]C", "T[T>A]G", "T[T>A]T",
               "A[T>C]A", "A[T>C]C", "A[T>C]G", "A[T>C]T", "C[T>C]A", "C[T>C]C", "C[T>C]G", "C[T>C]T", "G[T>C]A", "G[T>C]C", "G[T>C]G", "G[T>C]T", "T[T>C]A", "T[T>C]C", "T[T>C]G", "T[T>C]T",
               "A[T>G]A", "A[T>G]C", "A[T>G]G", "A[T>G]T", "C[T>G]A", "C[T>G]C", "C[T>G]G", "C[T>G]T", "G[T>G]A", "G[T>G]C", "G[T>G]G", "G[T>G]T", "T[T>G]A", "T[T>G]C", "T[T>G]G", "T[T>G]T")

for (j in tumortype){
  print(paste0("processing ", j))
  for (genomic_window in c(10, 30, 50, 100, 1000, 10000, 100000)){
    print(paste0("genomic window size = ", genomic_window))
    df <- ampcawg[((ampcawg$mech_s == "translocation" & ampcawg$mcj_s != "no" & ampcawg$mech_e != "fold_back_inversion") | (ampcawg$mech_e == "translocation" & ampcawg$mcj_e != "no" & ampcawg$mech_s != "fold_back_inversion")) & (ampcawg$end - ampcawg$start > 200000) & ampcawg$histology_abbreviation == j,]
    if (nrow(df) != 0){
      for (w in unique(df$icgc_donor_id)){
        snvdf <- read.csv(paste0("../../ICGC_PCAWG_calls/final_consensus_12aug_passonly/snv_mnv/01.1_Parsed_txt/r01_", sumdf$tumor[sumdf$icgc_donor_id == w], ".consensus.snv_mnv.parsed.txt"), header=T, as.is=T, sep="\t")
        snvdf$icgc_donor_id <- w
        if (w == unique(df$icgc_donor_id)[1]){
          vici <- data.frame(matrix(vector(), 0, ncol(snvdf), dimnames=list(c(), colnames(snvdf))), stringsAsFactors=F)
        }
        for (i in which(df$icgc_donor_id == w)){
          if (df$mech_s[i] == "translocation" & df$mcj_s[i] != "no"){
            vici <- rbind(vici, snvdf[snvdf$X.CHR == df$chromosome[i] & snvdf$POS >= df$start[i] & snvdf$POS <= df$start[i]+genomic_window,])
          }
          if (df$mech_e[i] == "translocation" & df$mcj_e[i] != "no"){
            vici <- rbind(vici, snvdf[snvdf$X.CHR == df$chromosome[i] & snvdf$POS <= df$end[i] & snvdf$POS >= df$end[i]-genomic_window,])
          }
        }
      }
      if (nrow(vici) != 0){
        vici$mutspectra <- NA
        for (i in 1:nrow(vici)){
          vici$mutspectra[i] <- ext.TN(vici$X.CHR[i], vici$POS[i], vici$REF[i], vici$VAR[i])
        }
        vici$mutspectra <- factor(vici$mutspectra, levels = specorder)
        spec_vici <- data.frame(matrix(vector(), 0, 98, dimnames=list(c(), c("gw", specorder, "n_tumor"))), stringsAsFactors=F)
        spec_vici[(nrow(spec_vici)+1),] <- NA
        spec_vici$gw[nrow(spec_vici)] <- genomic_window
        spec_vici[nrow(spec_vici),c(2:97)] <- as.vector(table(vici$mutspectra))
        spec_vici$n_tumor[nrow(spec_vici)] <- length(unique(vici$icgc_donor_id))
        
        pdf(paste0("Collab_APOBEC/Mutspectra/PCAWG/mutspectra.ampcawg.", j, ".boundary.translocation.vicinity.", genomic_window, "bp.pdf"), height=4, width=8)
        barplot(as.vector(t(spec_vici[spec_vici$gw == genomic_window,c(2:97)])), border=NA, las=1, col=sigcolor, space=0.8, main=paste0("# mutations = ", nrow(vici), " # tumors = ", length(unique(vici$icgc_donor_id))), ylab="Number of mutations", xlab="Trinucleotide context")
        dev.off()
        rm(spec_vici)
        rm(vici)
      }
    }
  }
}


# Pan-cancer expression
load("../../ICGC_PCAWG_calls/RNA_Expression/GeneExpression_matrices.RData") ## gives us gexp.normal and gexp.tumor
tail(gexp.tumor[gexp.tumor$geneName == "APOBEC3B","DO1663"])
tail(gexp.tumor[gexp.tumor$geneName == "APOBEC3B","DO3037"])
length(tail(gexp.tumor[gexp.tumor$geneName == "APOBEC3B","DO3140"]))
length(tail(gexp.tumor[gexp.tumor$geneName == "APOBEC3B","DO0"]))
rm(gexp.normal)
rm(gexp.tumor)

pcsig <- read.csv("../../ICGC_PCAWG_calls/Mutational_signatures/SigProfiler_Sanger/PCAWG_sigProfiler_SBS_signatures_in_samples.csv", header=T)
pcsig$icgc_donor_id <- NA
for (i in 1:nrow(pcsig)){
  if (nrow(sumdf[sumdf$icgc_specimen_id == pcsig$Sample.Names[i],]) != 0){
    pcsig$icgc_donor_id[i] <- sumdf$icgc_donor_id[sumdf$icgc_specimen_id == pcsig$Sample.Names[i]]
  }
}
pcsig <- pcsig[!is.na(pcsig$icgc_donor_id),]
pcsig$APOBEC3A <- NA
pcsig$APOBEC3B <- NA
pcsig$APOBEC3C <- NA
pcsig$AICDA <- NA
for (i in 1:nrow(pcsig)){
  if (length(tail(gexp.tumor[gexp.tumor$geneName == "APOBEC3A", pcsig$icgc_donor_id[i]])) == 1){
    pcsig$APOBEC3A[i] <- tail(gexp.tumor[gexp.tumor$geneName == "APOBEC3A", pcsig$icgc_donor_id[i]])
  }
  if (length(tail(gexp.tumor[gexp.tumor$geneName == "APOBEC3B", pcsig$icgc_donor_id[i]])) == 1){
    pcsig$APOBEC3B[i] <- tail(gexp.tumor[gexp.tumor$geneName == "APOBEC3B", pcsig$icgc_donor_id[i]])
  }
  if (length(tail(gexp.tumor[gexp.tumor$geneName == "APOBEC3C", pcsig$icgc_donor_id[i]])) == 1){
    pcsig$APOBEC3C[i] <- tail(gexp.tumor[gexp.tumor$geneName == "APOBEC3C", pcsig$icgc_donor_id[i]])
  }
  if (length(tail(gexp.tumor[gexp.tumor$geneName == "AICDA", pcsig$icgc_donor_id[i]])) == 1){
    pcsig$AICDA[i] <- tail(gexp.tumor[gexp.tumor$geneName == "AICDA", pcsig$icgc_donor_id[i]])
  }
}

for (i in unique(pcsig$Cancer.Types)){
  if (nrow(pcsig[pcsig$Cancer.Types == i & !is.na(pcsig$APOBEC3A),]) >50){
    df <- pcsig[pcsig$Cancer.Types == i & !is.na(pcsig$APOBEC3A),]
    if (nrow(df[df$SBS2 >1 | df$SBS13 >1,]) != 0){
      pdf(paste0("Collab_APOBEC/Expression_APOBEC/PCAWG/", i, "_A3A_SBS2+SBS13.pdf"), height=5, width=5)
      plot(log(df$SBS2[df$SBS2 != 0 | df$SBS13 != 0]+df$SBS13[df$SBS2 != 0 | df$SBS13 != 0]) ~ df$APOBEC3A[df$SBS2 != 0 | df$SBS13 != 0], xlab="APOBEC3A expression", ylab="Burden of SBS2+SBS13 (log10-transformed)", las=1, pch=19, col=rgb(0,0,0,.2), main=paste0(nrow(df[!is.na(df$APOBEC3A) & df$SBS2 != 0 & df$SBS13 != 0,]), i, " (r = ", round(sqrt(summary(lm(log(df$SBS2[df$SBS2 != 0 & df$SBS13 != 0]+df$SBS13[df$SBS2 != 0 & df$SBS13 != 0]) ~ df$APOBEC3A[df$SBS2 != 0 & df$SBS13 != 0]))$adj.r.squared), 3), ")"))
      abline(lm(log(df$SBS2[df$SBS2 != 0 | df$SBS13 != 0]+df$SBS13[df$SBS2 != 0 | df$SBS13 != 0]) ~ df$APOBEC3A[df$SBS2 != 0 | df$SBS13 != 0]))
      dev.off()
    }
    if (nrow(df[df$SBS2 >1 | df$SBS13 >1,]) != 0){
      pdf(paste0("Collab_APOBEC/Expression_APOBEC/PCAWG/", i, "_A3B_SBS2+SBS13.pdf"), height=5, width=5)
      plot(log(df$SBS2[df$SBS2 != 0 | df$SBS13 != 0]+df$SBS13[df$SBS2 != 0 | df$SBS13 != 0]) ~ df$APOBEC3B[df$SBS2 != 0 | df$SBS13 != 0], xlab="APOBEC3B expression", ylab="Burden of SBS2+SBS13 (log10-transformed)", las=1, pch=19, col=rgb(0,0,0,.2), main=paste0(nrow(df[!is.na(df$APOBEC3B) & df$SBS2 != 0 & df$SBS13 != 0,]), i, " (r = ", round(sqrt(summary(lm(log(df$SBS2[df$SBS2 != 0 & df$SBS13 != 0]+df$SBS13[df$SBS2 != 0 & df$SBS13 != 0]) ~ df$APOBEC3B[df$SBS2 != 0 & df$SBS13 != 0]))$adj.r.squared), 3), ")"))
      abline(lm(log(df$SBS2[df$SBS2 != 0 | df$SBS13 != 0]+df$SBS13[df$SBS2 != 0 | df$SBS13 != 0]) ~ df$APOBEC3B[df$SBS2 != 0 | df$SBS13 != 0]))
      dev.off()
    }
    if (nrow(df[df$SBS2 >1 | df$SBS13 >1,]) != 0){
      pdf(paste0("Collab_APOBEC/Expression_APOBEC/PCAWG/", i, "_A3C_SBS2+SBS13.pdf"), height=5, width=5)
      plot(log(df$SBS2[df$SBS2 != 0 | df$SBS13 != 0]+df$SBS13[df$SBS2 != 0 | df$SBS13 != 0]) ~ df$APOBEC3C[df$SBS2 != 0 | df$SBS13 != 0], xlab="APOBEC3C expression", ylab="Burden of SBS2+SBS13 (log10-transformed)", las=1, pch=19, col=rgb(0,0,0,.2), main=paste0(nrow(df[!is.na(df$APOBEC3C) & df$SBS2 != 0 & df$SBS13 != 0,]), i, " (r = ", round(sqrt(summary(lm(log(df$SBS2[df$SBS2 != 0 & df$SBS13 != 0]+df$SBS13[df$SBS2 != 0 & df$SBS13 != 0]) ~ df$APOBEC3C[df$SBS2 != 0 & df$SBS13 != 0]))$adj.r.squared), 3), ")"))
      abline(lm(log(df$SBS2[df$SBS2 != 0 | df$SBS13 != 0]+df$SBS13[df$SBS2 != 0 | df$SBS13 != 0]) ~ df$APOBEC3C[df$SBS2 != 0 | df$SBS13 != 0]))
      dev.off()
    }
    if (nrow(df[df$SBS2 >1 | df$SBS13 >1,]) != 0){
      pdf(paste0("Collab_APOBEC/Expression_APOBEC/PCAWG/", i, "_AID_SBS2+SBS13.pdf"), height=5, width=5)
      plot(log(df$SBS2[df$SBS2 != 0 | df$SBS13 != 0]+df$SBS13[df$SBS2 != 0 | df$SBS13 != 0]) ~ df$AICDA[df$SBS2 != 0 | df$SBS13 != 0], xlab="AICDA expression", ylab="Burden of SBS2+SBS13 (log10-transformed)", las=1, pch=19, col=rgb(0,0,0,.2), main=paste0(nrow(df[!is.na(df$AICDA) & df$SBS2 != 0 & df$SBS13 != 0,]), i, " (r = ", round(sqrt(summary(lm(log(df$SBS2[df$SBS2 != 0 & df$SBS13 != 0]+df$SBS13[df$SBS2 != 0 & df$SBS13 != 0]) ~ df$AICDA[df$SBS2 != 0 & df$SBS13 != 0]))$adj.r.squared), 3), ")"))
      abline(lm(log(df$SBS2[df$SBS2 != 0 | df$SBS13 != 0]+df$SBS13[df$SBS2 != 0 | df$SBS13 != 0]) ~ df$AICDA[df$SBS2 != 0 | df$SBS13 != 0]))
      dev.off()
    }
  }
}

pdf("Collab_APOBEC/A3A_SBS2+SBS13.pdf", height=5, width=5)
plot(log(df$SBS2[df$SBS2 != 0 & df$SBS13 != 0]+df$SBS13[df$SBS2 != 0 & df$SBS13 != 0]) ~ df$APOBEC3A[df$SBS2 != 0 & df$SBS13 != 0], xlim=c(-10,10), xlab="APOBEC3A expression (normalized)", ylab="Burden of SBS2+SBS13 (log10-transformed)", las=1, pch=19, col=rgb(0,0,0,.2), main=paste0(nrow(df[!is.na(df$APOBEC3A) & df$SBS2 != 0 & df$SBS13 != 0,]), " breast cancers (r = ", round(sqrt(summary(lm(log(df$SBS2[df$SBS2 != 0 & df$SBS13 != 0]+df$SBS13[df$SBS2 != 0 & df$SBS13 != 0]) ~ df$APOBEC3A[df$SBS2 != 0 & df$SBS13 != 0]))$adj.r.squared), 3), ")"))
abline(lm(log(df$SBS2[df$SBS2 != 0 & df$SBS13 != 0]+df$SBS13[df$SBS2 != 0 & df$SBS13 != 0]) ~ df$APOBEC3A[df$SBS2 != 0 & df$SBS13 != 0]))
dev.off()

summary(lm(log(df$SBS2[df$SBS2 != 0 & df$SBS13 != 0]+df$SBS13[df$SBS2 != 0 & df$SBS13 != 0]) ~ df$APOBEC3A[df$SBS2 != 0 & df$SBS13 != 0]))
sqrt(summary(lm(log(df$SBS2[df$SBS2 != 0 & df$SBS13 != 0]+df$SBS13[df$SBS2 != 0 & df$SBS13 != 0]) ~ df$APOBEC3A[df$SBS2 != 0 & df$SBS13 != 0]))$adj.r.squared)






