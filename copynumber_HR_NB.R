
#' Packages and functions
#' ======================
library(copynumber)
library(dplyr)
library(biomaRt)
library(survival)
format_heatmap <- function(profiles, samples) {
  # need: sampleID  chrom  arm  start.pos  end.pos  n.probes  mean
  # not with %in% because should be in order of samples, rev because starts at bottom
  profiles <- bind_rows(profiles[rev(samples)])
  names(profiles) <- c("sampleID", "chrom", "start.pos", "end.pos", "mean", "annotation")
  profiles <- profiles[,1:5]
  profiles$arm <- c(rep("p",nrow(profiles)))
  profiles$n.probes <- c(rep(1,nrow(profiles)))
  profiles <- profiles[,c(1,2,6,3,4,7,5)] # would be better with select(..,.., arm, everything()), see fischer_data_TBX2_ampl.R
  profiles <- data.frame(profiles)
}
KM <- function(samples, group_var, S="OS") {
  if (S=="OS") {
    KM_curve <- do.call(survival::survfit, list(formula = Surv(samples$OStime, samples$OS) ~ group_var))
    KM_diff <- survdiff(Surv(samples$OStime, samples$OS) ~ group_var, rho=0)
  }
  if (S=="EFS") {
    KM_curve <- do.call(survival::survfit, list(formula = Surv(samples$EFStime, samples$EFS) ~ group_var))
    KM_diff <- survdiff(Surv(samples$EFStime, samples$EFS) ~ group_var, rho=0)
  }
  pKM <- 1 - pchisq(KM_diff$chisq, length(KM_diff$n) - 1)
  return(list(KM_curve=KM_curve, pKM=pKM, summ_tab=summary(KM_curve)$table, S=S))
}
plot_KM <- function(KM, red_first=T, p2=F, title="", chsize=0.85, psize=0.85, pposx=0.01, pposy=0.05, years=NULL) {
  # KM is list of KM_curve, pKM (2 kinds) and summary table to get number of samples/events
  KM$KM_curve$time <- KM$KM_curve$time/365.24
  print(names(KM$KM_curve$strata)) #to check red and blue on graph
  col <- c("firebrick2", "deepskyblue4", "chartreuse3", "orange", "brown", "darksalmon", "lightcyan3", "burlywood4", "chartreuse4")
  if(!red_first) {
    col <- col[c(2,1,3,4,5,6,7,8,9)] #change order if nec
  }
  par(cex=chsize, mgp=c(1.7,0.8,0))
  plot(KM$KM_curve, mark.time=TRUE, lwd=2.5, col=col, main=title, axes=FALSE, xlim=years, xlab="Years after diagnosis", ylab="")
  axis(1, las=1, pos=0, at=seq(0,20,2.5), labels=c(0,2.5,5,7.5,10,12.5,15,17.5,20))
  axis(2, las=1, pos=0)
  mtext(KM$S, side=2, line=2.6, cex=chsize)
  legend("topright", legend=paste(gsub("group_var=", "", names(KM$KM_curve$strata)), "  ", KM$summ_tab[,"records"], " (", KM$summ_tab[,"events"], ")"), lwd=2, col=col) 
  mylabel <- bquote(p == .(format(KM$pKM, digits=3)))
  par(cex=psize)
  text(pposx*max(KM$KM_curve$time), pposy, mylabel, pos=4)
  if (p2) {
    mylabel2 <- bquote(alt.p == .(format(KM$pKM2, digits=3)))
    par(cex=psize)
    text(pposx*max(KM$KM_curve$time), pposy+0.07, mylabel2, pos=4)
  }
  par(cex=1, mgp=c(3,1,0))
}
select_focal <- function(loc, ab="gain", prof_byplatf=profiles_byplatf, size=5e6, cutoff_plot=2, plot=T, save=FALSE) {
  chr <- gsub("chr","",loc$chromosome_name)
  start <- loc$start_position
  stop <- loc$end_position
  GOI <- loc$hgnc_symbol
  # prof is list of profiles within 1 platform, g/l is cutoff for 1 particular platform
  # gain
  profiles_gain_byplatf <- Map(function(prof, g) {
    Filter(function(y) {nrow(y)>0}, lapply(prof, function(x){x[x$chromosome%in%chr & x$min<=stop & x$max>=start &
                                                                 (x$max - x$min) <= size & x$logratio>=g,]}))
  }, prof_byplatf, cutoff_gain_platf)
  # loss
  profiles_loss_byplatf <- Map(function(prof, l) {
    Filter(function(y) {nrow(y)>0}, lapply(prof, function(x){x[x$chromosome%in%chr & x$min<=stop & x$max>=start &
                                                                 (x$max - x$min) <= size & x$logratio<=l,]}))
  }, prof_byplatf, cutoff_loss_platf)
  if (ab=="gain") profiles_GOI <- unlist(profiles_gain_byplatf, recursive=FALSE)
  if (ab=="loss") profiles_GOI <- unlist(profiles_loss_byplatf, recursive=FALSE)
  if (ab=="both") profiles_GOI <- c(unlist(profiles_gain_byplatf, recursive=FALSE), unlist(profiles_loss_byplatf, recursive=FALSE))
  if (length(profiles_GOI)==0) {print(paste("No aberrations of the selected size could be found for", GOI)); return(NULL)}
  names(profiles_GOI) <- gsub("^.+\\.", "", names(profiles_GOI))
  print(length(profiles_GOI))
  if (plot) {
    heatmap_colors <- c("deepskyblue4", "white", "firebrick")
    profiles_GOI_hm <- format_heatmap(profiles, names(profiles_GOI))
    if (save) pdf(file=paste0("figures/MYCNall/val_COG/heatmap_", GOI, ".pdf"), w=27, h=4.9)
    plotHeatmap(segments=profiles_GOI_hm, upper.lim=cutoff_plot, plot.ideo=TRUE, chrom=chr, sep.samples=0.1, cyto.text=TRUE,
                colors=heatmap_colors, ideo.frac=0.15
                #, layout=c(1,2)
    )
    abline(v=mean(c(start, stop))/1e6)
    text(mean(c(start, stop))/1e6, length(profiles_GOI)*1.01, labels=GOI, srt=90, adj=0, cex=0.8, xpd=TRUE)
  }
  return(profiles_GOI)
}
survival_focal <- function(profiles_GOI) {
  presence_GOI <- ifelse(samples$Name %in% names(profiles_GOI), "gain", "no gain")
  KM_GOI <- KM(samples, presence_GOI)
  plot_KM(KM_GOI)
  KM_GOI_EFS <- KM(samples, presence_GOI, S="EFS")
  plot_KM(KM_GOI_EFS)
}
cutoff_loss_platf <- list(Affymetrix=-0.25, Agilent=-0.3, Illumina=-0.25, NimbleGen=-0.3)
cutoff_gain_platf <- list(Affymetrix=0.15, Agilent=0.2, Illumina=0.15, NimbleGen=0.2)


#' Read data
#' =========
samples <- read.delim("samples_dd.txt")
profiles_df <- read.delim("profiles_dd.txt")
head(samples)
nrow(samples)
tail(profiles_df)

profiles <- split(profiles_df, f=profiles_df$Name)
names(profiles)
profiles_byplatf <- split(profiles, f=samples$Platform)
names(profiles_byplatf)


#' General view of the dataset
#' ===========================
#' Using frequency plot, heatmap and aberration plot from package copynumber
hm_col <- c("deepskyblue4", "white", "firebrick")
profiles_hm <- format_heatmap(profiles, names(profiles))
par(ask=FALSE)
plotHeatmap(segments=profiles_hm, upper.lim=0.8, plot.ideo=TRUE, sep.samples=0.1, cyto.text=TRUE,
            colors=hm_col, ideo.frac=0.1)
plotHeatmap(segments=profiles_hm, upper.lim=0.8, plot.ideo=TRUE, chrom=11, sep.samples=0.1, cyto.text=TRUE,
            colors=hm_col, ideo.frac=0.1)
plotFreq(segments=profiles_hm, thres.gain=0.2, thres.loss=-0.25, plot.ideo=TRUE, cyto.text=TRUE,
         col.gain="firebrick", col.loss="deepskyblue4", ideo.frac=0.1)
plotFreq(segments=profiles_hm, thres.gain=0.2, thres.loss=-0.25, plot.ideo=TRUE, chrom=11, cyto.text=TRUE,
         col.gain="firebrick", col.loss="deepskyblue4", ideo.frac=0.1)


#' Aberrations for gene of interest
#' ================================
#' By default, select_focal() looks for focal gains (max 5 Mb) overlapping with a certain region/gene, 
#' but the input parameters can be changed to look for larger aberrations ("size=..." in bp)
#' or for losses (ab="loss") or both gains and losses (ab="both")
GOI <- "ODC1"
ensembl_mart <- useMart("ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice")
ensembl_data <- useDataset("hsapiens_gene_ensembl", mart=ensembl_mart)
GOI_loc <- getBM(attributes=c("hgnc_symbol", "chromosome_name", "start_position", "end_position"), filters=c("hgnc_symbol"),
                 values=list(GOI), mart=ensembl_data)
GOI_loc
par(ask=FALSE)
profiles_GOI <- select_focal(GOI_loc, size=5e6, ab="both", cutoff_plot=0.8)
par(mfrow=c(1,2), mar=c(4,4,1,1))
survival_focal(profiles_GOI)


#' Get IDs for R2 platform (https://hgserver1.amc.nl/cgi-bin/r2/main.cgi)
#' ======================================================================
#' Create a custom sample list in R2 format, to use in the R2 genome browser
samples_MNA <- samples[samples$MYCN %in% 1,]
ID_R2_MNA <- paste0("KDP_1709_", formatC(samples_MNA$Name, width=3, flag="0"))
paste0(ID_R2_MNA, collapse=",")

