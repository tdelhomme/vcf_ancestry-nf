args <- commandArgs(TRUE)
parseArgs <- function(x) {
  res = strsplit(sub("^--", "", x), "=")
  if(length(unlist(res))==1) res[[1]][2]=""
  return(res)
}

argsL <- as.list(do.call("cbind", parseArgs(args))[c(F,T)])
names(argsL) <- as.list(do.call("cbind", parseArgs(args))[c(T,F)])
args <- argsL;rm(argsL)

if(! is.null(args$help)) {
  cat("
      Mandatory arguments:
      --eigenvec_file            - Eigen vector file from plink PCA
      --PED                      - PED file for reference

      Optional arguments:

      --help \n\n")
  q(save="no")
}

if(is.null(args$eigenvec_file)) {stop("Option --eigenvec_file should be provided")} else{eigenvec_file=args$eigenvec_file}
if(is.null(args$PED)) {stop("Option --PED should be provided")} else{pedfile=args$PED}

eigenvec <- read.table(eigenvec_file, header = FALSE, skip=0, sep = ' ')
rownames(eigenvec) <- eigenvec[,2]
eigenvec <- eigenvec[,3:ncol(eigenvec)]
colnames(eigenvec) <- paste('Principal Component ', c(1:20), sep = '')

cohort_ids = which( !grepl("^HG0", rownames(eigenvec)) & !grepl("^NA", rownames(eigenvec)))

KGids = setdiff(1:length(rownames(eigenvec)), cohort_ids)
eigenvec = eigenvec[c(KGids, cohort_ids),] # correct order, just in case

PED <- read.table(pedfile, header = TRUE, skip = 0, sep = '\t')[,c("Individual.ID", "Population")]
PED = PED[which(PED$Individual.ID %in% rownames(eigenvec)[1:length(KGids)]), ]
PED1 = rbind(PED, data.frame(Individual.ID = rownames(eigenvec)[(length(KGids)+1):nrow(eigenvec)], Population = "unknown"))
PED2 <- PED1[match(rownames(eigenvec), PED1$Individual.ID),]
print(paste(date(), " INFO: does the individuals are concordant between new PED file and eigen vector from PCA? "))
all(PED2$Individual.ID == rownames(eigenvec)) == TRUE
#[1] TRUE

# set colours
require('RColorBrewer')

PED = PED2
# from: http://www.internationalgenome.org/category/population/
PED$Population <- factor(PED$Population, levels=c(
  "ACB","ASW","ESN","GWD","LWK","MSL","YRI",
  "CLM","MXL","PEL","PUR",
  "CDX","CHB","CHS","JPT","KHV",
  "CEU","FIN","GBR","IBS","TSI",
  "BEB","GIH","ITU","PJL","STU", "unknown"))

col <- colorRampPalette(c(
  "yellow","yellow","yellow","yellow","yellow","yellow","yellow",
  "forestgreen","forestgreen","forestgreen","forestgreen",
  "grey","grey","grey","grey","grey",
  "royalblue","royalblue","royalblue","royalblue","royalblue",
  "black","black","black","black","black", "pink"))(length(unique(PED$Population)))[factor(PED$Population)]

project.pca <- eigenvec
summary(project.pca)


pdf("plots.pdf", 15, 5)
par(mfrow = c(1,3))
plot(project.pca[,1], project.pca[,2], col = col, pch = 20, main = 'A', adj = 0.5, cex = 2,
     xlab = 'First component', ylab = 'Second component', font = 2, font.lab = 2)
legend('topright', bty = 'n', title = '', c('AFR', 'AMR', 'EAS', 'EUR', 'SAS', "Input cohort"),
       fill = c('yellow', 'forestgreen', 'grey', 'royalblue', 'black', 'pink'))
plot(project.pca[,1], project.pca[,3], col = col, pch=20, cex = 2, main="B", adj=0.5, xlab="First component", ylab="Third component", font=2, font.lab=2)
plot(project.pca[,2], project.pca[,3], col = col, pch=20, cex = 2, main="B", adj=0.5, xlab="Second component", ylab="Third component", font=2, font.lab=2)
dev.off()

project.pca$ID = rownames(project.pca)
write.table(project.pca[,c(1:3, match("ID", colnames(project.pca)))], file="table_3PCs.txt", quote = F, row.names = F, col.names = T, sep = "\t")

                        