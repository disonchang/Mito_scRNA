#!/usr/bin/env Rscript
suppressMessages(library("Seurat"))
suppressMessages(library("dplyr"))
suppressMessages(library("reshape2"))
suppressMessages(library("ggrepel"))
suppressMessages(library("grid"))
suppressMessages(library("gridExtra"))
suppressMessages(library("cowplot"))

source("~/software/bin/R/paper_theme.R")

args=commandArgs(trailing=TRUE)

usage="seuratDataExc.R <seurat.rds> <MT tally> <NUC tally>"

rds=args[1]
mtfn=args[2]
nucfn=args[3]

#rds="../clustering2/SJE2A037_D.rds"
#mtfn="SJE2A015_D.varCount.CBUB.txt.tally.txt.mt.txt.pospoisson.0.2.pvalues.txt"
#mtfn="SJE2A037_D.varCount.CBUB.tally2.txt.mt.txt"
#nucfn="SJE2A037_D.varCount.CBUB.tally2.txt.nuc.txt"


d <- readRDS(rds)
print("tsne loaded")
tsne <- data.frame(d@dr$tsne@cell.embeddings)
cluster.id <- data.frame(d@ident)

data <- merge(tsne, cluster.id, by="row.names")
data$sample <- gsub("(SJ.*_[A-Z])_.*","\\1", data$Row.names)
data$BC <- gsub("SJ.*_[A-Z]_(.*)","\\1", data$Row.names)

output=paste(mtfn,"tSNE_cell_info", "txt", sep=".")
write.table(data, output, sep="\t", quote=FALSE, row.names=FALSE)

mt=read.table(mtfn, header=TRUE, sep="\t", stringsAsFactors=FALSE)
nuc=read.table(nucfn, header=TRUE, sep="\t", stringsAsFactors=FALSE)




print("table loaded")
### 
mt.somatic.all <- mt[which(grepl("TUMOR-specific|lowVAF", mt$GroupID) & mt$alt_nUMI > 0),]
mt.somatic.all2 <- mt[which(grepl("TUMOR-specific|lowVAF", mt$GroupID) > 0),] ### all cells, including alt_nUMI==0
mt.somatic.pos <- unlist(unique(mt.somatic.all$position))
print(dim(mt.somatic.all))


colnames(data) <- gsub("d\\.ident","Cluster",colnames(data))

p <- ggplot(data, aes(x=tSNE_1, y=tSNE_2, color=Cluster))+
	geom_point(alpha=1, size=2.5)
output=paste(mtfn,"tSNE","backbone","pdf", sep=".")
ggsave(output, plot=p, width=10, height=7)
print ("ann plot done")



for (site in 1:length(mt.somatic.pos)){
	print(mt.somatic.pos[site])
	mt.somatic <- mt.somatic.all[which(mt.somatic.all$position==mt.somatic.pos[site]),]
	mt.somatic2 <- mt.somatic.all2[which(mt.somatic.all2$position==mt.somatic.pos[site]),]
	mt.somatic$CB <- gsub("-1$","",mt.somatic$CB)
  	mt.somatic2$CB <- gsub("-1$","",mt.somatic2$CB)
	xmax=max(mt.somatic2$nUMI)
	binw=5
        if ( xmax < 60 ){ binw=2.5 } #else if (xmax < 100) { binw=10}
	bin_num=ceiling(xmax/binw)
	labWGS=round(max(mt.somatic2$MAF_Tumor),4)
	labX=max(mt.somatic2$nUMI)-10
	labY=max(mt.somatic2$MAF_Tumor)
	p0 <- ggplot(mt.somatic2, aes(x=nUMI, y=alt_nUMI_ratio))+
		geom_point()+
		geom_hline(yintercept=labY, color="red")+
		geom_text(label = labWGS, colour = "red", y = labY, x = labX, hjust = -.1, vjust=-.2, check_overlap = TRUE, size=4)+
		scale_x_continuous(expand=c(0,binw/2), breaks=seq(0,binw*bin_num,binw), limit=c(0,xmax+binw/2))
		#scale_y_continuous(expand=c(-0.05,0), limit=c(-0.1,1.1), breaks=seq(0,1,0.25))
	#pd <- mt.somatic2 %>% group_by(nUMI_rng=cut(nUMI, breaks=seq(0,binw*bin_num,by=5), alt_nUMI_ratio_rng=cut(alt_nUMI_ratio, breaks=seq(0,1,by=0.1)))) %>%
	#	summarise(n= n()) #%>%
		#arrange(as.numeric(nUMI_rng), as.numeric(alt_nUMI_ratio_rng))
	umi_step=1
	ratio_step=0.01
	mt.somatic2$nUMI_bin <- cut(mt.somatic2$nUMI, breaks=c(seq(0,binw*bin_num,by=umi_step), Inf), labels=seq(0,binw*bin_num,by=umi_step)+umi_step/2, include.lowest = TRUE)
	mt.somatic2$alt_nUMI_ratio_bin <- cut(mt.somatic2$alt_nUMI_ratio, breaks=c(seq(0,1,by=ratio_step), Inf), labels=seq(0,1,by=ratio_step)+ratio_step/2, include.lowest = TRUE)
	mt.somatic2$nUMI_bin <- as.numeric(as.character(mt.somatic2$nUMI_bin))
	mt.somatic2$alt_nUMI_ratio_bin <- as.numeric(as.character(mt.somatic2$alt_nUMI_ratio_bin))
	pd <- mt.somatic2 %>% group_by(nUMI_bin,alt_nUMI_ratio_bin) %>% summarise(cell_n= n()) %>% arrange(as.numeric(nUMI_bin), as.numeric(alt_nUMI_ratio_bin))


	p1 <- ggplot(pd, aes(x=nUMI_bin, y=alt_nUMI_ratio_bin, size=cell_n))+
                geom_point(color="dimgrey")+
                #geom_hline(yintercept=labY, color="red")+
                #geom_text(label = labWGS, colour = "red", y = labY, x = labX, hjust = -.1, vjust=-.2)+
                scale_x_continuous(expand=c(0,binw/2), breaks=seq(0,binw*bin_num,binw), limit=c(0,xmax+binw/2))+
		theme(legend.position="none")+
		geom_hline(yintercept=labY, color="red")+
		geom_text(label = labWGS, colour = "red", y = labY, x = labX, hjust = -.1, vjust=-.2, check_overlap = TRUE, size=4)


        ydensity0 <- ggplot(mt.somatic2, aes(alt_nUMI_ratio)) + 
		geom_histogram(aes(y=..count..), position="identity", binwidth=0.05, breaks=seq(0,1,0.05), closed="left") +
  		#geom_density(alpha=.5) +
  		theme(legend.position = "none")
	plotdata <- as.data.frame(ggplot_build(ydensity0)$data)


	#ydensity0 <- ydensity0 + geom_text(data = plotdata, aes(x=x, y= count, label = count), hjust=-0.15, size=4)
		#geom_vline(xintercept=unique(mt.somatic$MAF_Tumor, color="red", fill="red")) +
    	ymax1=max(plotdata$count)+800
	#ymax2=max(plotdata$count)+80000 ### for log10
	ymax2=200
	ydensity1 <- ydensity0 + geom_text(data = plotdata, aes(x=x, y= count, label = count), hjust=-0.15, size=4)
	ydensity1 <- ydensity1 +scale_x_continuous(expand=c(0,0), breaks = seq(0, 1, 0.25), limit=c(-0.05,1.05))
	ydensity1 <- ydensity1 + ylim(0, ymax1) 
	ydensity1 <- ydensity1 + coord_flip()+ theme(axis.title.y=element_blank())
	
	#ydensity2 <- ydensity + scale_y_continuous(trans='log10', limits=c(1,ymax2))
	
	ydensity2 <- ydensity0 + geom_text(data = subset(plotdata, count < 100), aes(x=x, y= count, label = count), hjust=-0.15, size=4)
	ydensity2 <- ydensity2 + scale_x_continuous(breaks = seq(0, 1, 0.25), limit=c(-0.05,1.05))
	if(nrow(subset(plotdata, count >= 100 & count <= ymax2)) > 0){
		ydensity2 <- ydensity2 + geom_text(data = subset(plotdata, count >= 100 & count <= ymax2), aes(x=x, y= count, label = count), hjust=1, color='white',size=4)
	}
	if(nrow(subset(plotdata, count > ymax2))>0){
		ydensity2 <- ydensity2 + geom_text(data = subset(plotdata, count > ymax2), aes(x=x, y= ymax2, label = count), hjust=1, color='white',size=4)
	}
	ydensity2 <- ydensity2 + coord_flip(ylim = c(0,ymax2), xlim = c(0, 1))+ theme(axis.title.y=element_blank()) #+ coord_cartesian(xlim=c(0, 200))
	
	xdensity <- ggplot(mt.somatic2, aes(nUMI)) +
                geom_histogram(aes(y=..count..), position="identity", binwidth=binw, breaks=seq(0,binw*bin_num, binw), closed = "left") +
                #geom_density(alpha=.5) +
                theme(legend.position = "none")
        plotdata <- as.data.frame(ggplot_build(xdensity)$data)
        xdensity <- xdensity + geom_text(data = plotdata, aes(x=x, y= count, label = count), vjust=-1, size=3.5)
                #geom_vline(xintercept=unique(mt.somatic$MAF_Tumor, color="red", fill="red")) +
        ymax=max(plotdata$count)+300
        xdensity <- xdensity +scale_x_continuous(expand=c(0,0), breaks = seq(0, binw*bin_num, binw), limit=c(-binw/2,xmax+binw))
        xdensity <- xdensity + ylim(0, ymax)
        xdensity <- xdensity + theme(axis.title.x=element_blank())

        
	#output=paste(mtfn,"MTsomatic","density",mt.somatic.pos[site], "pdf", sep=".")
	##pdf(output, width=10, height=7)
	##grid.arrange(p0, ydensity, nrow=1, layout_matrix=rbind(c(1,1,1,2)), top = textGrob(paste("position:", mt.somatic.pos[site], "cellN:", length(unique(mt.somatic2$CB)), sep=" "),gp=gpar(fontsize=14)))
	#comb <- plot_grid(xdensity,NULL,p0, ydensity1, align = 'hv', rel_heights = c(1, 3), rel_widths = c(3, 1))
	#title <- ggdraw() + draw_label(paste("position:", mt.somatic.pos[site], "cellN:", length(unique(mt.somatic2$CB)), sep=" "))
	#comb2 <- plot_grid (title, comb, ncol=1, rel_heights=c(0.1, 1))
	##comb
	##dev.off()
	#save_plot(output, comb2, base_height = 10, base_width=10)


	output=paste(mtfn,"MTsomatic","density2",mt.somatic.pos[site], "pdf", sep=".")
        #pdf(output, width=10, height=7)
        #grid.arrange(p0, ydensity, nrow=1, layout_matrix=rbind(c(1,1,1,2)), top = textGrob(paste("position:", mt.somatic.pos[site], "cellN:", length(unique(mt.somatic2$CB)), sep=" "),gp=gpar(fontsize=14)))
        
	comb <- plot_grid(xdensity,NULL,p1, ydensity2, align = 'hv', rel_heights = c(1, 3), rel_widths = c(3, 1))
	#comb <- plot_grid(p1,xdensity, align = 'v', rel_heights = c(1, 1))
        title <- ggdraw() + draw_label(paste("position:", mt.somatic.pos[site], "cellN:", length(unique(mt.somatic2$CB)), sep=" "))
        comb2 <- plot_grid (title, comb, ncol=1, rel_heights=c(0.1, 1))
        #comb
        #dev.off()
        save_plot(output, comb2, base_height = 10, base_width=10)


	mt.somatic.summary <- mt.somatic2 %>% group_by(CB) %>% ### for getting all cells using mt.somatic2, for getting cell with alt_nUMI > 0 use mt.somatic
		summarize(pos=paste(position, collapse=';'),
                  max_nUMI=max(nUMI),
                  min_nUMI=min(nUMI),
                  max_alt_nUMI=max(alt_nUMI),
                  min_alt_nUMI=min(alt_nUMI),
                  max_alt_nUMI_ratio=max(alt_nUMI_ratio),
                  min_alt_nUMI_ratio=min(alt_nUMI_ratio),
                  max_ref_nUMI=max(ref_nUMI),
                  min_ref_nUMI=min(ref_nUMI),
		  #pvalues=min(pvalues),
                  #fdr=min(FDR),
                  max_alt_nUMI_pos=paste(position[alt_nUMI==max_alt_nUMI], collapse=','),
                  min_alt_nUMI_pos=paste(position[alt_nUMI==min_alt_nUMI], collapse=','),
                  max_alt_nUMI_posn=length(position[alt_nUMI==max_alt_nUMI]),
                  min_alt_nUMI_posn=length(position[alt_nUMI==min_alt_nUMI])
		)
	mt.somatic.summary$BC <- mt.somatic.summary$CB; mt.somatic.summary$CB <- NULL

	### call cell as MT or NUC only
	determine_somatic_mode="ALL" ### UMI1/ALL
	mt.somaticCell=data.frame()
	nuc.somaticCell=data.frame()
	if (determine_somatic_mode=="UMI1"){
		### [option 1] cell with at least one read with somatic variation
		mt.somaticCell=data.frame(MT_cell=gsub("-1$","", unique(mt[which(grepl("TUMOR-specific", mt$GroupID) & mt$alt_nUMI > 0),][['CB']]))) ### for getting cell with alt_nUMI > 0
		nuc.somaticCell=data.frame(NUC_cell=gsub("-1$","", unique(nuc[which(nuc$alt_nUMI > 0),][['CB']]))) ### for getting cell with alt_nUMI > 0
	} else if (determine_somatic_mode=="ALL"){
		### all mt and nuc cells
		mt.somaticCell=data.frame(MT_cell=gsub("-1$","", unique(mt[which(grepl("TUMOR-specific", mt$GroupID)),][['CB']])))
        	nuc.somaticCell=data.frame(NUC_cell=gsub("-1$","", unique(nuc$CB)))
	}

	rownames(mt.somaticCell)=mt.somaticCell$MT_cell
	rownames(nuc.somaticCell)=nuc.somaticCell$NUC_cell
	all.somaticCell=merge(mt.somaticCell, nuc.somaticCell, by="row.names", all=TRUE)

	attach(all.somaticCell)
	all.somaticCell$SomaticCell <- ifelse(!is.na(MT_cell) & !is.na(NUC_cell), 'MT_NUC', ifelse(!is.na(MT_cell),'MT','NUC')) ### Determine based on MT and NUC status simultaneously. MT_UNC indicated that the cell was present in both MT and NUC count table at the same time
	all.somaticCell$MT_SomaticCell <- ifelse(!is.na(MT_cell), 'MT', NA) ### Determine based on MT_cell status only
	detach(all.somaticCell)

	colnames(all.somaticCell) <- gsub("Row.names","BC",colnames(all.somaticCell))
	#colnames(data) <- gsub("d\\.ident","Cluster",colnames(data))

	### merge tSNE and count data
	data.ann <- merge(data, all.somaticCell, by="BC", all.x=TRUE)
	data.ann <- merge(data.ann, mt.somatic.summary, by="BC", all.x=TRUE) 
	output=paste(mtfn,"tSNE","ann",mt.somatic.pos[site], determine_somatic_mode, "txt", sep=".")
	write.table(data.ann, output, row.names=FALSE, quote=FALSE, sep="\t")
	print ("ann summary done")
	#print(head(data.ann))

	### cluster summary
	### NOTE: the explanation for the plot here is different between UMI1 and ALL mode. 
	### For UMI1 mode, it implies the cell with matched barcodes in tSNE and at least one altUMI 
	### For ALL mode, it implies the cell with matched barcodes in tSNE 
	cluster.summary <- data.ann %>% group_by(sample, Cluster) %>%
		summarize(BCn=n(), 
		  uniqueBCn=length(unique(BC)), 
		  #uniqueBCn_wSomatic=length(unique(BC[!is.na(SomaticCell)])), 
                  CellN.woSomatic=length(unique(BC[is.na(SomaticCell)])),
                  #MT_CellN=length(BC[SomaticCell=='MT' & ! is.na(SomaticCell)]), 
		  CellN.MT=length(unique(BC[SomaticCell=='MT' & ! is.na(SomaticCell)])),
                  #NUC_CellN=length(BC[SomaticCell=='NUC' & ! is.na(SomaticCell)]), 
		  CellN.NUC=length(unique(BC[SomaticCell=='NUC' & ! is.na(SomaticCell)])),
                  #MT.NUC_CellN=length(BC[SomaticCell=='NUC' & ! is.na(SomaticCell)]), 
		  CellN.MT_NUC=length(unique(BC[SomaticCell=='MT_NUC' & ! is.na(SomaticCell)]))
		)

	output=paste(mtfn,"tSNE","ann",mt.somatic.pos[site], determine_somatic_mode, "summary", "txt", sep=".")
	write.table(cluster.summary, output, row.names=FALSE, quote=FALSE, sep="\t")
	print ("cluster summary done")


	### plot
	p <- ggplot(data.ann, aes(x=tSNE_1, y=tSNE_2, color=Cluster))+
		geom_point(data=subset(data.ann, !is.na(SomaticCell) & SomaticCell !="MT_NUC"), alpha=1, aes(shape=SomaticCell))+
		geom_point(data=subset(data.ann, !is.na(SomaticCell) & SomaticCell =="MT_NUC"), alpha=1, aes(shape=SomaticCell), size=2.5)+
		geom_point(data=subset(data.ann, is.na(SomaticCell)), alpha=0.3, shape=21)
	if(determine_somatic_mode=='UMI1'){
		p <- p + guides(shape=guide_legend(title="Cell with altUMIn > 1"))
	} else if (determine_somatic_mode=='ALL'){
		p <- p + guides(shape=guide_legend(title="Cell present at the site\n(including 0 altUMIn)"))
	}
	output=paste(mtfn,"tSNE","ann", mt.somatic.pos[site],determine_somatic_mode,"pdf", sep=".")
	ggsave(output, plot=p, width=10, height=7)
	print ("ann plot done")


	p <- ggplot(data.ann, aes(x=tSNE_1, y=tSNE_2, color=Cluster))+
		#geom_text(data=subset(data.ann, MT_SomaticCell %in% c('MT') & max_alt_nUMI_ratio > 0.5), aes(label=max_alt_nUMI))+
                geom_point(data=subset(data.ann, MT_SomaticCell %in% c('MT') & max_nUMI >= 3 & max_alt_nUMI_ratio > 0.4), alpha=0.8, shape=16, aes(size=max_alt_nUMI))+ ### solid circle
		scale_fill_manual(values=getPalette(8, n=8, p="Dark2"))+
                geom_point(data=subset(data.ann, MT_SomaticCell %in% c('MT') & max_nUMI > 0 & (max_nUMI < 3 | max_alt_nUMI_ratio <= 0.4)), alpha=0.8, shape=21, aes(size=max_alt_nUMI))+ ### empty circle
		geom_point(data=subset(data.ann, MT_SomaticCell %in% c('MT') & max_nUMI == 0), alpha=0.8, fill="grey", shape=21, size=1)+ ### empty circle
                geom_point(data=subset(data.ann, ! MT_SomaticCell %in% c('MT')), alpha=0.8, shape=22, fill="grey", size=1) ### empty square

	if(determine_somatic_mode=='UMI1'){
		p <- p + guides(shape=guide_legend(title="Cell with altUMIn > 1"))
	} else if(determine_somatic_mode=='ALL'){
		p <- p + guides(shape=guide_legend(title="Cell present at the site\n(including 0 altUMIn)"))
	}
	output=paste(mtfn,"tSNE","MTonly","ann", mt.somatic.pos[site], determine_somatic_mode, "pdf", sep=".")
	ggsave(output, plot=p, width=10, height=7)
	print ("ann MTonly plot done")




	#print(subset(data.ann, fdr <= 0.05))
	#p <- ggplot(data.ann, aes(x=tSNE_1, y=tSNE_2, color=Cluster))+
	#        #geom_point(data=subset(data.ann, MT_SomaticCell %in% c('MT') & fdr <= 0.05), alpha=0.8, shape=22, aes(fill=max_alt_nUMI_pos), size=4)+
	#	#geom_text_repel(data=subset(data.ann, MT_SomaticCell %in% c('MT') & fdr <= 0.05), aes(label=paste0(max_alt_nUMI, "/", max_nUMI, "(", round(fdr,2), ")")), force=10)+	
	#        geom_point(data=subset(data.ann, MT_SomaticCell %in% c('MT') & pvalues <= 0.05), alpha=0.8, shape=22, aes(fill=max_alt_nUMI_pos), size=4)+
	#	geom_text_repel(data=subset(data.ann, MT_SomaticCell %in% c('MT') & pvalues <= 0.05), aes(label=paste0(max_alt_nUMI, "/", max_nUMI, "(", round(fdr,2), ")")), force=10)+   
	#	scale_fill_manual(values=getPalette(8, n=8, p="Dark2"))+
	#       #geom_point(data=subset(data.ann, MT_SomaticCell %in% c('MT') & (is.na(fdr) | fdr > 0.05)), alpha=0.8, shape=16)+
	#        geom_point(data=subset(data.ann, MT_SomaticCell %in% c('MT') & (is.na(pvalues) | pvalues > 0.05)), alpha=0.8, shape=16)+
	#	geom_point(data=subset(data.ann, ! MT_SomaticCell %in% c('MT')), alpha=0.2, shape=21)+
	#        guides(shape=guide_legend(title="Cell with MT somatic\nmutations"))
	#output=paste(mtfn,"tSNE","MTonly","ann","pval", "pdf", sep=".")
	#ggsave(output, plot=p, width=10, height=7)


	### split position plots

	### 1. Update Tisney plots. 
	# Create separate plots for each mtDNA/RNA position of interest (ie somatic variants). 
	# All cells that pass qc should be color-coded based on group (as in the current plots). 
	# Cells with no umi’s covering the mtDNA position of interest should be empty circles (as in the current plots). 
	# The size of the filled circles (ie containing at least 1 mutant UMI) should reflect the total number of umi’s for the mtDNA position of interest in that cell (ref+alternate alleles). 
	# Would it be possible for these “filled circles” to actually be pie charts reflecting the mutant vafs (with alternate allele frequency in black)?
#positions=unique(data.ann$position[!is.na(data.ann$pos)])


	#for (i in 1:length(positions)){
        #var_cutoff=0.3
	umi_cutoff=1
	data.ann$alt_nUMI_ratio_bin <- ifelse(data.ann$max_alt_nUMI_ratio > 0.75, 75, ifelse(data.ann$max_alt_nUMI_ratio > 0.5, 50, ifelse(data.ann$max_alt_nUMI_ratio > 0.25, 25, 1)))
	print(head(data.ann))
	data0=subset(data.ann, is.na(data.ann$pos)) ### tSNE cell not present in count data, 0 coverage
	data0.1=subset(data.ann, ! is.na(data.ann$pos))
	#data1=subset(data0.1, pos==mt.somatic.pos[site] & max_alt_nUMI_ratio >= vaf_cutoff)
	#data2=subset(data0.1, pos!=mt.somatic.pos[site] | max_alt_nUMI_ratio < vaf_cutoff)
	data1=subset(data0.1, pos==mt.somatic.pos[site] & max_alt_nUMI >= umi_cutoff)
	data2=subset(data0.1, pos!=mt.somatic.pos[site] | max_alt_nUMI < umi_cutoff)
	print(paste0("total cell#: ",nrow(data.ann)))
	print(paste0("Non-present cells: ", nrow(data0), " Present cells: ", nrow(data0.1)))
	print(paste0("Match criteria cell#: ", nrow(data1)))
	print(paste0("Not match criteria cell#: ", nrow(data2)))

	scatterPie <- data1[, c("BC", "tSNE_1", "tSNE_2", "Cluster", "max_nUMI", "max_ref_nUMI","max_alt_nUMI")] ### NOTE: the max used here as it's unclear whether any cell has multiple position mapped ornot curretly
	colnames(scatterPie) <- c("BC", "tSNE_1", "tSNE_2", "Cluster", "UMIn", "ref_UMIn","alt_UMIn")
	scatterPie$tSNEtotal <- scatterPie$tSNE_1 + scatterPie$tSNE_2
	id_cols=c("tSNE_1", "tSNE_2", "tSNEtotal", "BC", "Cluster")
	scatterPie.w <- melt(scatterPie, id.vars=id_cols, measure.vars=c("ref_UMIn","alt_UMIn"), variable.name='Var', value.name='nUMI')

	### use plotrix
	suppressMessages(require(plotrix))
	Vmin=min(scatterPie[,c("tSNE_1", "tSNE_2")]) - 0.5
	Vmax=max(scatterPie[,c("tSNE_1", "tSNE_2")]) + 0.5
	output=paste(mtfn,"tSNE","MTonly","ann","pos",mt.somatic.pos[site], determine_somatic_mode,"pdf", sep=".")
	frac_cols=getPalette(2, p="Set3", n=9)
	clr_cols=getPalette(length(unique(data.ann$Cluster)), p="Set1", n=9)
	names(clr_cols) <- unique(data.ann$Cluster)
	pdf(output, width=10, height=7) 
		plot(Vmin:Vmax, xlim=c(Vmin, Vmax), ylim=c(Vmin, Vmax), type="n",main=paste("MT somatic mutation position:", mt.somatic.pos[site], sep=" "),xlab="tSNE_1",ylab="tSNE_2",bty='l')
		for (j in 1:nrow(data2)){
                        x=data2[j,"tSNE_1"]
                        y=data2[j,"tSNE_2"]
                        #r=log(data2[j,"max_nUMI"]+1,10)
                        r=0.15
			#r=data2$scale
                        cols=alpha(clr_cols[data2[j,"Cluster"]],1)
                        draw.circle(x, y, radius=r, col="white", border=cols)
                        #text(x, y-0.1, labels=scatterPie[j, "alt_UMIn"])
                }

		for (j in 1:nrow(data0)){
                        x=data0[j,"tSNE_1"]
                        y=data0[j,"tSNE_2"]
                        r=0.15
                        cols=alpha(clr_cols[data0[j,"Cluster"]],0.5)
                        draw.circle(x, y, radius=r, col="white", border="dimgrey")
                        #text(x, y-0.1, labels=scatterPie[j, "alt_UMIn"])
                }
		for (j in 1:nrow(data1)){
                        x=data1[j,"tSNE_1"]
                        y=data1[j,"tSNE_2"]
                        frac=unlist(data1[j, c("max_alt_nUMI","max_ref_nUMI")])
                        r=log(data1[j,"max_nUMI"]+1,10)/3
                        #r=data1$scale
                        cols=alpha(c(clr_cols[data1[j,"Cluster"]], "grey"),0.5)
                        floating.pie(x, y, frac, radius=r, col=cols)
                        #text(x, y-0.1, labels=scatterPie[j, "alt_UMIn"])
                }

	dev.off()


	clr_cols=ggplotColours(n = length(unique(data.ann$Cluster)))
	#clr_cols=getPalette(length(unique(data.ann$Cluster)), p="Set1", n=9)
	names(clr_cols) <- unique(data.ann$Cluster)
	p2 <- ggplot(data.ann, aes(x=tSNE_1, y=tSNE_2, color=Cluster, size=max_nUMI, fill=Cluster))+
        	geom_point(data=data1, shape=21, aes(alpha=max_alt_nUMI_ratio, fill=Cluster))+
		scale_size_continuous(name="UMI#", breaks=round(seq(min(data.ann$max_nUMI, na.rm=T), max(data.ann$max_nUMI, na.rm=T), length.out=10),0), range=c(1,5))+
		scale_alpha_continuous(name="altUMI ratio", breaks=seq(0,1,by=0.1), range=c(0.2,0.9))+
		scale_fill_manual(values=clr_cols)+
		scale_color_manual(values=clr_cols)+
		geom_point(data=data2, alpha=0.4, shape=21, size=1, fill=NA)+
		geom_point(data=data0, alpha=0.4, shape=21, size=1, fill=NA, color="dimgrey")+
		ggtitle(paste("MT somatic mutation position:", mt.somatic.pos[site], sep=" "))+
		guides(fill=FALSE, size=guide_legend(ncol=2), alpha=guide_legend(ncol=2))
	output=paste(mtfn,"tSNE","MTonly","ann","pos",mt.somatic.pos[site],determine_somatic_mode, "2", "pdf", sep=".")
	ggsave(output, plot=p2, width=10, height=7)

        p2 <- ggplot(data.ann, aes(x=tSNE_1, y=tSNE_2, color=Cluster, size=max_nUMI, fill=Cluster))+
                geom_point(data=data1, shape=21, aes(alpha=max_alt_nUMI_ratio, fill=Cluster))+
                scale_size_continuous(name="UMI#", breaks=round(seq(min(data.ann$max_nUMI, na.rm=T), max(data.ann$max_nUMI, na.rm=T), length.out=10),0), range=c(1,5))+
                scale_alpha_continuous(name="altUMI ratio", breaks=seq(0,1,by=0.1), range=c(0.2,0.9))+
                scale_fill_manual(values=clr_cols)+
                scale_color_manual(values=clr_cols)+
                geom_point(data=data2, alpha=0.4, shape=21, fill=NA)+
                geom_point(data=data0, alpha=0.4, shape=21, size=1, fill=NA, color="dimgrey")+
                ggtitle(paste("MT somatic mutation position:", mt.somatic.pos[site], sep=" "))+
                guides(fill=FALSE, size=guide_legend(ncol=2), alpha=guide_legend(ncol=2))
        output=paste(mtfn,"tSNE","MTonly","ann","pos",mt.somatic.pos[site],determine_somatic_mode, "3", "pdf", sep=".")
        ggsave(output, plot=p2, width=10, height=7)

	p2 <- ggplot(data.ann, aes(x=tSNE_1, y=tSNE_2, size=max_nUMI, fill=max_alt_nUMI_ratio))+
                #scale_color_manual(values=clr_cols)+
                geom_point(data=data2, alpha=0.5, shape=21, fill=NA, color="dimgrey")+ ### UMI in tSNE cells and alt count < cutoff
                geom_point(data=data0, alpha=0.8, shape=21, size=1, fill="black", color="black")+ ### tSNE cell not in count data 
		geom_point(data=data1, shape=21, color="dimgrey")+ ### UMI in tSNE cells and alt count >= cutoff
                scale_size_continuous(name="UMI# (total coverage)", breaks=round(seq(min(data.ann$max_nUMI, na.rm=T), max(data.ann$max_nUMI, na.rm=T), length.out=10),0), range=c(1,5))+
                scale_fill_gradient(name="VAF", low="blue", high="red")+
                ggtitle(paste("MT somatic mutation position:", mt.somatic.pos[site], sep=" "))+
                guides(size=guide_legend(ncol=2), alpha=guide_legend(ncol=2))
        output=paste(mtfn,"tSNE","MTonly","ann","pos",mt.somatic.pos[site],determine_somatic_mode, "4", "pdf", sep=".")
        ggsave(output, plot=p2, width=10, height=7)



	### plot2
	cluster.summary.p <- melt(cluster.summary, id.vars=c('sample','Cluster'), measure.vars=c('CellN.woSomatic','CellN.MT','CellN.NUC','CellN.MT_NUC'), variable.name='cell.type', value.name='count')
	cluster.summary.p <- merge(cluster.summary.p, cluster.summary[,c('sample','Cluster','uniqueBCn')], by=c('sample', 'Cluster'), all.x=TRUE)
	cluster.summary.p$proportion <- round(cluster.summary.p$count/cluster.summary.p$uniqueBCn, 4)
	cluster.summary.p$label <- paste0(round(cluster.summary.p$proportion*100,1),"%(", cluster.summary.p$count,")")
	p <- ggplot(cluster.summary.p, aes(x=factor(Cluster), y=proportion, fill=cell.type, label=label))+
		scale_fill_manual(values=getPalette(4, n=4, p="Dark2"))+	
        	geom_bar(stat='identity', position='fill')+
		geom_text(data=subset(cluster.summary.p, count !=0), size = 2, position = position_stack(vjust = 0.5), color="white")+
		theme_bw() + paper_theme+
		xlab("Cluster")+ylab("Proportion")

	output=paste(mtfn,"tSNE","ann",mt.somatic.pos[site],determine_somatic_mode, "summary","pdf", sep=".")
	ggsave(output, plot=p, width=10, height=7)
}
