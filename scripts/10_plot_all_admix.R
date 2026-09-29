#!/usr/bin/env Rscript

# usage
#TIP_LABEL_COLUMN=phyloID \
#SAMPLE_LABEL_CEX=0.45 \
#PANEL_HEADER_CEX=0.68 \
#Rscript scripts/10_plot_all_admix.R \
#  results/admixture_10rep \
#  data/Dryas_sampledata.csv \
#  figures/09_dryas_tree_admixture_sample_order.txt \
#  figures/10_dryas_all_runs \
#  10 25
  
  
# Plot every ADMIXTURE run in active-tree order plus a separate CV figure.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4 || length(args) > 6) stop(paste("Usage: 10_plot_tree_and_cv.R",
  "ADMIXTURE_DIR METADATA.csv ORDER_FILE OUTPUT_PREFIX [K_MIN] [K_MAX]"))
admixture_dir <- normalizePath(args[1]); metadata_file <- args[2]
order_file <- args[3]; output_prefix <- sub("\\.pdf$", "", args[4], ignore.case=TRUE)
k_min <- if (length(args)>=5) as.integer(args[5]) else -Inf
k_max <- if (length(args)>=6) as.integer(args[6]) else Inf
label_column <- Sys.getenv("TIP_LABEL_COLUMN", "phyloID")
label_cex <- as.numeric(Sys.getenv("SAMPLE_LABEL_CEX", "0.40"))
header_cex <- as.numeric(Sys.getenv("PANEL_HEADER_CEX", "0.68"))
dir.create(dirname(output_prefix), recursive=TRUE, showWarnings=FALSE)

cv <- read.delim(file.path(admixture_dir,"cross_validation.tsv"), stringsAsFactors=FALSE)
needed <- c("K","replicate","CV_error")
if (length(setdiff(needed,names(cv)))) stop("CV table needs K, replicate, and CV_error")
cv <- cv[cv$K>=k_min & cv$K<=k_max, needed, drop=FALSE]
cv <- cv[order(cv$K,cv$replicate),,drop=FALSE]
if (!nrow(cv) || any(!is.finite(cv$CV_error))) stop("No valid CV records")
if (anyDuplicated(cv[c("K","replicate")])) stop("Duplicated K/replicate records")
q_order <- readLines(file.path(admixture_dir,"sample_order.txt")); plot_order <- readLines(order_file)
if (anyDuplicated(q_order)||anyDuplicated(plot_order)) stop("Duplicate IDs in an order file")
if (!setequal(q_order,plot_order)) stop("ORDER_FILE and sample_order.txt have different samples")
metadata <- read.csv(metadata_file,check.names=FALSE,fileEncoding="UTF-8-BOM",stringsAsFactors=FALSE)
if (!(label_column %in% names(metadata))) label_column <- "sampleID"
m <- metadata[match(plot_order,metadata$sampleID),,drop=FALSE]
if (any(is.na(m$sampleID))) stop("ORDER_FILE samples absent from metadata")
labels <- as.character(m[[label_column]]); bad <- is.na(labels)|labels==""; labels[bad] <- plot_order[bad]
groups <- sort(unique(as.character(m$spID))); gp <- setNames(hcl.colors(length(groups),"Dark 3"),groups)
label_colors <- unname(gp[as.character(m$spID)])

safe_cor <- function(x,y) { sx<-sd(x,na.rm=TRUE); sy<-sd(y,na.rm=TRUE)
  if(!is.finite(sx)||!is.finite(sy)||sx==0||sy==0) return(-Inf)
  z<-suppressWarnings(cor(x,y,use="pairwise.complete.obs")); if(is.finite(z)) z else -Inf }
align_q <- function(a,b) {
  s<-outer(seq_len(ncol(a)),seq_len(ncol(b)),Vectorize(function(i,j)safe_cor(a[,i],b[,j])))
  ui<-rep(FALSE,nrow(s)); uj<-rep(FALSE,ncol(s)); pairs<-matrix(integer(),0,2)
  for(z in seq_len(min(nrow(s),ncol(s)))) { x<-s; x[ui,]<--Inf; x[,uj]<--Inf
    ij<-which(x==max(x),arr.ind=TRUE)[1,]; pairs<-rbind(pairs,ij); ui[ij[1]]<-TRUE; uj[ij[2]]<-TRUE }
  b[,c(pairs[order(pairs[,1]),2],which(!uj)),drop=FALSE] }
find_q <- function(k,r) { x<-c(file.path(admixture_dir,paste0("rep",r),paste0("K",k,".rep",r,".Q")),
  file.path(admixture_dir,paste0("K",k,".rep",r,".Q"))); x<-x[file.exists(x)]
  if(!length(x)) stop("Q file not found for K=",k," R=",r); x[1] }
qs <- vector("list",nrow(cv))
for(i in seq_len(nrow(cv))) { q<-as.matrix(read.table(find_q(cv$K[i],cv$replicate[i]),header=FALSE))
  if(nrow(q)!=length(q_order)) stop("Q row count differs from sample_order.txt")
  rownames(q)<-q_order; qs[[i]]<-q[match(plot_order,q_order),,drop=FALSE] }
for(i in seq_along(qs)) { if(i==1) qs[[i]]<-qs[[i]][,order(apply(qs[[i]],2,which.max)),drop=FALSE]
  else qs[[i]]<-align_q(qs[[i-1]],qs[[i]]) }
cols <- if(requireNamespace("viridisLite",quietly=TRUE)) viridisLite::turbo(max(cv$K)) else hcl.colors(max(cv$K),"Turbo")

admix_pdf<-paste0(output_prefix,"_all_admixture_runs.pdf"); n<-length(plot_order); p<-nrow(cv); yy<-seq_len(n)
pdf(admix_pdf,width=max(18,7.2+.92*p),height=32,useDingbats=FALSE)
layout(matrix(seq_len(p+1),1),widths=c(7.2,rep(.92,p)))
par(mar=c(.6,.2,2.4,.1),xaxs="i",yaxs="i"); plot.new(); plot.window(c(0,1),c(.5,n+.5))
text(.99,yy,labels,adj=c(1,.5),cex=label_cex,col=label_colors); title("Samples in active-tree order",cex.main=.85,line=.25)
for(i in seq_len(p)) { q<-qs[[i]]; par(mar=c(.6,.02,2.4,.02),xaxs="i",yaxs="i")
  plot.new(); plot.window(c(0,1),c(.5,n+.5)); rect(0,.5,1,n+.5,col="#eeeeee",border=NA)
  for(j in seq_len(n)) { ends<-cumsum(q[j,]); starts<-c(0,head(ends,-1)); rect(starts,yy[j]-.49,ends,yy[j]+.49,col=cols[seq_along(ends)],border=NA) }
  title(sprintf("K=%d\nR%d\n%.5f",cv$K[i],cv$replicate[i],cv$CV_error[i]),cex.main=header_cex,line=.05); box(lwd=.7) }
dev.off()

cv_pdf<-paste0(output_prefix,"_cv.pdf"); reps<-sort(unique(cv$replicate)); repcols<-setNames(viridisLite::turbo(length(reps)),reps)
pdf(cv_pdf,width=9.5,height=5.8,useDingbats=FALSE); par(mar=c(4.2,4.5,2.8,7),xpd=NA)
plot(NA,xlim=range(cv$K),ylim=range(cv$CV_error),xlab="K",ylab="10-fold CV error",xaxt="n",main="ADMIXTURE cross-validation")
axis(1,at=sort(unique(cv$K))); for(r in reps) { z<-cv[cv$replicate==r,]; z<-z[order(z$K),]
  lines(z$K,z$CV_error,type="b",pch=16,cex=.5,col=repcols[as.character(r)]) }
legend("topright",inset=c(-.27,0),legend=paste0("R",reps),col=repcols,lty=1,pch=16,bty="n",ncol=2,cex=.75,title="Replicate")
dev.off()

# Summarize variation in CV error across replicate runs at each K. The outer
# light-gray ribbon is the observed minimum-to-maximum range; the darker ribbon
# is mean +/- one standard deviation; the black line is the replicate mean.
summary_k <- sort(unique(cv$K))
cv_summary <- do.call(rbind, lapply(summary_k, function(k) {
  values <- cv$CV_error[cv$K == k]
  data.frame(
    K = k,
    mean = mean(values),
    sd = if (length(values) > 1) sd(values) else 0,
    minimum = min(values),
    maximum = max(values),
    n = length(values)
  )
}))

summary_pdf <- paste0(output_prefix, "_cv_summary.pdf")
pdf(summary_pdf, width = 8.5, height = 5.8, useDingbats = FALSE)
par(mar = c(4.2, 4.5, 2.8, 1.2), xpd = FALSE)
plot(
  NA,
  xlim = range(cv_summary$K),
  ylim = range(c(cv_summary$minimum, cv_summary$maximum)),
  xlab = "K", ylab = "10-fold CV error", xaxt = "n",
  main = "ADMIXTURE cross-validation summary"
)
axis(1, at = cv_summary$K)
abline(h = pretty(range(cv_summary$minimum, cv_summary$maximum)),
       col = "#eeeeee", lwd = 0.7)

# Draw the full range first so the SD ribbon remains visible on top.
polygon(
  c(cv_summary$K, rev(cv_summary$K)),
  c(cv_summary$minimum, rev(cv_summary$maximum)),
  col = "#e5e5e5", border = NA
)
polygon(
  c(cv_summary$K, rev(cv_summary$K)),
  c(cv_summary$mean - cv_summary$sd,
    rev(cv_summary$mean + cv_summary$sd)),
  col = "#a6a6a6", border = NA
)
lines(cv_summary$K, cv_summary$mean, type = "b", pch = 16,
      cex = 0.72, lwd = 1.8, col = "black")
legend(
  "topright",
  legend = c("Mean", "Mean +/- 1 SD", "Replicate range"),
  col = c("black", "#a6a6a6", "#e5e5e5"),
  lty = c(1, NA, NA), pch = c(16, 15, 15),
  pt.cex = c(0.8, 1.4, 1.4), lwd = c(1.8, NA, NA),
  bty = "n", cex = 0.8
)
dev.off()

summary_table <- paste0(output_prefix, "_cv_summary.tsv")
write.table(cv_summary, summary_table, sep = "\t", row.names = FALSE,
            quote = FALSE)

message("ADMIXTURE figure: ", admix_pdf)
message("Replicate CV figure: ", cv_pdf)
message("Summary CV figure: ", summary_pdf)
message("Summary CV table: ", summary_table)
