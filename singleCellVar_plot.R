#!/usr/bin/env Rscript
suppressMessages(library("ggplot2"))
suppressMessages(library("ggrepel"))
suppressMessages(library("dplyr"))
suppressMessages(library("grid"))
suppressMessages(library("gridExtra"))
source("~/software/bin/R/paper_theme.R")

args=commandArgs(trailing=TRUE)

#fn="SJE2A015_D.varCount.CBUB.txt.tally.txt.mt.summary.txt"
fn=args[1]

d = read.table(fn, sep="\t", header=TRUE)
cut=0.03

title=paste(unique(d$ID,sep=','))

stats <- data.frame(table(d$GroupID))
print(stats)
colnames(stats) <- c('VarType','Frequency')

p1 <- ggplot(d, aes(x=MAF_Tumor, y=total_altUMIratio, color=GroupID))+
	geom_point(data=subset(d, (total_altUMIratio < (1-cut) & total_altUMIratio > cut) | (MAF_Tumor < (1-cut) & MAF_Tumor > cut)), shape=16, size=2)+
	geom_point(data=subset(d, (total_altUMIratio >= (1-cut) | total_altUMIratio <= cut) & (MAF_Tumor >= (1-cut) | MAF_Tumor <= cut)), shape=4)+
        theme_bw()+
	paper_theme+
	xlab("MAF_WGS")+
	ggtitle(title)+
	guides(color=FALSE)+
	geom_text_repel(data=subset(d, (total_altUMIratio < (1-cut) & total_altUMIratio > cut) | (MAF_Tumor < (1-cut) & MAF_Tumor > cut)),
			aes(label =paste(round(total_altUMIratio,2),round(MAF_Tumor,2), sep='/')), segment.color = 'grey50', force=100, show.legend=FALSE,max.iter=20)


p2 <- ggplot(d, aes(x=MAF_Tumor, y=total_altReadRatio, color=GroupID))+
	geom_point(data=subset(d, (total_altReadRatio < (1-cut) & total_altReadRatio > cut) | (MAF_Tumor < (1-cut) & MAF_Tumor > cut)), shape=16, size=2)+
        geom_point(data=subset(d, (total_altReadRatio >= (1-cut) | total_altReadRatio <= cut) & (MAF_Tumor >= (1-cut) | MAF_Tumor <= cut)), shape=4)+
        theme_bw()+
        paper_theme+
        xlab("MAF_WGS")+
        ggtitle(title)+
        geom_text_repel(data=subset(d, (total_altReadRatio < (1-cut) & total_altReadRatio > cut) | (MAF_Tumor < (1-cut) & MAF_Tumor > cut)), 
                        aes(label =paste(round(total_altReadRatio,2),round(MAF_Tumor,2), sep='/')), segment.color = 'grey50', force=100, show.legend=FALSE,max.iter=20)


stats$Type <- 'Type'
stats$Prop <- stats$Frequency/sum(stats$Frequency)
p3 <- ggplot(stats, aes(x=Type, y=Prop, fill=VarType, label=Frequency))+
	geom_bar(stat='identity', position='fill')+
	geom_text(size = 3, position = position_stack(vjust = 0.5))+
	theme_bw()+
	coord_flip()+
	theme(axis.title=element_blank(),
              axis.text=element_blank(),
              axis.ticks=element_blank(),
	      panel.border = element_blank(),
              legend.position="bottom")


lay <- rbind(c(rep(1,5),rep(2,5)),
	     c(rep(1,5),rep(2,5)),
	     c(rep(1,5),rep(2,5)),
	     c(rep(1,5),rep(2,5)),
	     c(rep(1,5),rep(2,5)),
	     c(rep(1,5),rep(2,5)),
             c(rep(1,5),rep(3,5))
		)



pdf(paste(fn, 'umi', 'pdf', sep='.'), width=12, height=6)
grid.arrange(grobs=list(p1,p2,p3), layout_matrix = lay)
dev.off()

#out <- arrangeGrob(p1, p2, layout_matrix = lay)

#ggsave(paste(fn, 'umi', 'pdf', sep='.'), plot=p2, width=9, height=7)


