#!/usr/bin/env Rscript
suppressMessages(library("dplyr"))
suppressMessages(library("reshape2"))

args   <- commandArgs(trailingOnly=TRUE)

#fn="control.varCount.CBUB.txt"
#fn="SJE2A037_D.varCount.CBUB.txt.tally.txt"
fn=args[1]
#nuc_var="tier1_nuclear.var.txt"
nuc_var=args[2]
#mt_var="mtDNA_variants.txt"
mt_var=args[3]



print(fn)

### load UMI count
data <- read.table(fn, sep="\t", header=TRUE)
data.nuc = data[which(data$chr != 'MT'),]
data.mt = data[which(data$chr == 'MT'),]

### get UMI summary at sample level
data.nuc.summary <- data.nuc %>% group_by(ID, chr, position, ref, alt) %>%
	summarize(total_CellN=n(), total_UMIn=sum(nUMI), total_UMIsum=(sum(ref_nUMI)+sum(alt_nUMI)), total_ReadN=(sum(ref_nReads) + sum(alt_nReads)),
		  total_refUMIn=sum(ref_nUMI), total_altUMIn=sum(alt_nUMI), total_altUMIratio=total_altUMIn/total_UMIn,
                  total_refReadN=sum(ref_nReads), total_altReadN=sum(alt_nReads), total_altReadRatio=total_altReadN/total_ReadN,
		  total_RR_cellN=length(CB[alt_nUMI_ratio==0]), total_RA_cellN=length(CB[alt_nUMI_ratio>0 & alt_nUMI_ratio<1]), total_AA_cellN=length(CB[alt_nUMI_ratio==1]))

data.mt.summary <- data.mt %>% group_by(ID, chr, position, ref, alt) %>%
        summarize(total_CellN=n(), total_UMIn=sum(nUMI), total_UMIsum=(sum(ref_nUMI)+sum(alt_nUMI)), total_ReadN=(sum(ref_nReads) + sum(alt_nReads)),
		  total_refUMIn=sum(ref_nUMI), total_altUMIn=sum(alt_nUMI), total_altUMIratio=total_altUMIn/total_UMIn,
		  total_refReadN=sum(ref_nReads), total_altReadN=sum(alt_nReads), total_altReadRatio=total_altReadN/total_ReadN,
                  total_RR_cellN=length(CB[alt_nUMI_ratio==0]), total_RA_cellN=length(CB[alt_nUMI_ratio>0 & alt_nUMI_ratio<1]), total_AA_cellN=length(CB[alt_nUMI_ratio==1]))


### load variants files
var.nuc =read.table(nuc_var, sep="\t", header=TRUE)
var.mt =read.table(mt_var, sep="\t", header=TRUE)

var.nuc.f <- var.nuc[,c('Sample','Gene','mRNA.Accession','Variant', 'Chr', 'Position','Class','Mutant.In.Tumor','Total.In.Tumor')]
var.nuc.f$MAF_Tumor <- round(var.nuc.f$Mutant.In.Tumor/var.nuc.f$Total.In.Tumor,2)
var.mt.f <- var.mt[, c('Sample','GroupID','Predicted.Impact','Location','mtDNA_pos','MAF_Tumor','Func..Impact')]
var.mt.f$Sample <- paste0(var.mt.f$Sample, "_D")

### annotate 
data.nuc.ann <- merge(data.nuc, var.nuc.f, by.x=c('ID', 'chr', 'position'), by.y=c('Sample','Chr','Position'))
data.mt.ann <- merge(data.mt, var.mt.f, by.x=c('ID', 'position'), by.y=c('Sample','mtDNA_pos'))

data.nuc.summary.ann <- merge(data.nuc.summary, var.nuc.f, by.x=c('ID', 'chr', 'position'), by.y=c('Sample','Chr','Position'))
data.mt.summary.ann <- merge(data.mt.summary, var.mt.f, by.x=c('ID', 'position'), by.y=c('Sample','mtDNA_pos'))

### cell with somatic mutations
data.nuc.ann.somatic <- data.nuc.ann[which(data.nuc.ann$alt_nUMI > 0),]
data.nuc.ann.somatic.cells <- unique(as.character(data.nuc.ann.somatic$CB))

data.mt.ann.somatic <- data.mt.ann[which(data.mt.ann$alt_nUMI > 0),]
data.mt.ann.somatic.cells <- unique(as.character(data.mt.ann.somatic$CB))

data.nuc.ann$Cell.with.somaticMut <- ifelse(data.nuc.ann$CB %in% data.nuc.ann.somatic.cells, 'y', '')
data.mt.ann$Cell.with.somaticMut <- ifelse(data.mt.ann$CB %in% data.nuc.ann.somatic.cells, 'y', '') ### use this for identifying somatic cell by nuclear or mitochondria
 


### output
output.nuc=paste0(fn, ".nuc.txt")
output.mt=paste0(fn, ".mt.txt")
output.nuc.summary=paste0(fn, ".nuc.summary.txt")
output.mt.summary=paste0(fn, ".mt.summary.txt")

write.table(data.nuc.ann, file=output.nuc, quote=FALSE, sep="\t", row.names=FALSE)
write.table(data.mt.ann, file=output.mt, quote=FALSE, sep="\t", row.names=FALSE)
write.table(data.nuc.summary.ann, file=output.nuc.summary, quote=FALSE, sep="\t", row.names=FALSE)
write.table(data.mt.summary.ann, file=output.mt.summary, quote=FALSE, sep="\t", row.names=FALSE)




