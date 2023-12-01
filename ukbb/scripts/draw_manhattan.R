packages <- c("data.table", "dplyr", "qqman", "R.utils")

is.installed <- function(pkg) {
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) {
    install.packages(new.pkg, quiet=T)
  }
  sapply(pkg, require, character.only = TRUE)
}
is.installed(packages)

args = commandArgs(trailingOnly=TRUE)
if (length(args)!=2) {
  stop("Two arguments (infile, outfile) must be supplied (input file).n", call.=FALSE)
}

infile <- args[1]
outbase <- args[2]
#infile <- "assoc.resht_both.regenie.gz"
#outbase <- "resht_both"

df <- fread(infile)

png(paste0("manhattan.",outbase,".png"), width=1000, height=400)
manhattan(df, chr="CHROM", bp="GENPOS", p="P", snp="ID", 
              annotatePval=7.5, annotateTop=T) 
dev.off()

png(paste0("qqplot.",outbase,".png"), width=480, height=480)
qq(df$P) 
dev.off()

