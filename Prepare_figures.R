library(RiboseQC)
library(ggrepel)
library(ggridges)
library(cowplot)

trans <- function(x) {reso=-log(x, 10);reso[is.na(reso)]=0;return(reso)}
inv <- function(x) 10^(-x)
#setwd("/bfd/lcalviel/data/riboseq/new_ddx3/results/To_github/bef_submit/")
load("data/processed_data.RData")
dfa<-list_datasets$DE_siDDX3
dfa<-dfa[dfa$keep,]

dfa$padj2<-dfa$padj_RiboRNA
dfa$padj2[is.na(dfa$padj2)]<-1
dfa$padj2[dfa$padj2<10e-50]<-10e-50
dfa$experiment<-"siDDX3X vs control"

axm<-max(abs(c(dfa$log2FoldChange_RNA,dfa$log2FoldChange_Ribo)),na.rm = T)

toview<-c("DDX3X","ODC1","DVL1")
dfa$symb4<-NA
dfa$symb4[dfa$symb%in%toview]<-dfa$symb[dfa$symb%in%toview]

colsi<-alpha(c("blue","dark red","cornflowerblue","firebrick1","dark gray"),c(.8,.8,.8,.8,.3))
names(colsi)<-c("TE_down","TE_up","Concordant_down","Concordant_up","mixed_ns")
dfa$nudgex<-0
dfa$nudgey<-0
dfa$nudgex[which(dfa$symb3=="ODC1")]<-.7
dfa$nudgex[which(dfa$symb3=="DVL1")]<-.7
dfa$nudgey[which(dfa$symb3=="ODC1")]<-(-1)
dfa$nudgey[which(dfa$symb3=="DVL1")]<-(-.2)

a<-ggplot(dfa,aes(x=log2FoldChange_RNA,y=log2FoldChange_Ribo,color=tx_type,size=padj2,label=symb4,shape=gene_biotype)) + geom_point()
a<-a + theme_bw() +
    ylab(paste("Ribo-seq log2FC",sep = "")) +
    xlab(paste("RNA-seq log2FC",sep = "")) +
    theme(axis.title.x = element_text(size=22),axis.text.x  = element_text(angle=45, vjust=0.5, size=18)) +
    theme(axis.title.y = element_text(size=24),axis.text.y  = element_text(angle=45, vjust=0.5, size=18))  +
    theme(strip.background = element_rect(color="black", fill="white", size=1.5, linetype="solid")) +
    theme(strip.text.x = element_text(size = 15,face = "bold"))+
    scale_color_manual(values = colsi,"Tx_class") +
    scale_size_continuous("Adj. p-value (interaction)",trans = scales::trans_new(name = "test",transform = trans,inverse = inv, 
                                                                                 domain = c(1e-100, Inf)),breaks = c(10e-11,10e-5,10e-3,1))

a<-a + geom_hline(yintercept = 0,size=.5) + geom_vline(xintercept = 0,size=.5)
ariborna<-a + ylim(-axm,axm) + xlim(-axm,axm)
ariborna<-ariborna + facet_wrap(.~experiment,nrow = 1) 
ariborna<-ariborna + geom_text_repel(size=5,force=32,show.legend = F,
                                     nudge_x = dfa$nudgex[which(nchar(dfa$symb4)>2)],
                                     nudge_y = dfa$nudgey[which(nchar(dfa$symb4)>2)])

pdf(file = "figures/DE_plot.pdf",width = 9,height = 5)
ariborna
dev.off()


colsi<-alpha(c("blue","dark red","cornflowerblue","firebrick1","dark gray"),c(.8,.8,.8,.8,.3))
names(colsi)<-c("TE_down","TE_up","Concordant_down","Concordant_up","mixed_ns")

grd<-ggplot(dfa,aes(log2FC_5utrcds,group=tx_type,y=tx_type,fill=tx_type))  +theme_ridges()
grd<-grd + scale_fill_manual(values=colsi,"Tx_class")
grd<-grd +xlim(-4,4) + xlab("log2FC_5'UTR_skew")       
grd<-grd + stat_density_ridges(quantile_lines = TRUE, quantiles = 2) + geom_vline(xintercept = 0,lty=2,color="black")
grd<-grd + theme(legend.position="none")

pdf(file = "figures/DE_5riboskew.pdf",width = 7,height = 4)
grd
dev.off()


datasi<-dfa[,c("gene_id","symb","gene_biotype","log2FoldChange_RiboRNA","padj_RiboRNA")]
dfadeg<-list_datasets$DE_degron
dfadeg<-dfadeg[dfadeg$gene_id%in%datasi$gene_id,]
rownames(dfadeg)<-dfadeg$gene_id
rownames(datasi)<-datasi$gene_id
datasi<-datasi[rownames(dfadeg),]
dfadeg$padj_RiboRNA_si<-datasi$padj_RiboRNA
dfadeg$log2FoldChange_RiboRNA_si<-datasi$log2FoldChange_RiboRNA
dfadeg$pmin<-apply(dfadeg[,c("padj_RiboRNA","padj_RiboRNA_si")],1,FUN = min)
dfadeg$pmax<-apply(dfadeg[,c("padj_RiboRNA","padj_RiboRNA_si")],1,FUN = max)


df2<-dfadeg[,c("log2FoldChange_RiboRNA","log2FoldChange_RiboRNA_si","padj_RiboRNA","padj_RiboRNA_si","symb","gene_biotype")]

dfcp<-df2
nimok<-colnames(df2)[1:2]
trans <- function(x) -log(x, 10)
inv <- function(x) 10^(-x)
nmsok<-c("x_lfc","y_lfc","padj_RiboRNA","padj_RiboRNA_si","gene_name","gene_biotype")
colnames(df2)<-nmsok
df2$maxpv<-apply(df2[,3:4,drop=F],1,function(x){
    if(sum(is.na(x))==length(x)){return(NA)}
    else{return(mean(x,na.rm = T))}
})
df2<-df2[complete.cases(df2[,1:2]),]
df2<-df2[!is.infinite(df2[,1]),]
df2<-df2[!is.infinite(df2[,2]),]

df2$confidence<-"low"

df2$confidence<-factor(df2$confidence,levels=c("low","medium","high"))
df1<-df2
ok<-which(df2$maxpv<0.01)
df1$confidence[ok]<-"high"
if(length(ok)>2){
    df22<-df2[ok,]
    df22$confidence<-"high"
    df2<-rbind(df2,df22)
}
ok<-which(df2$maxpv<0.15)
df1$confidence[which(df1$maxpv<0.15)]<-"medium"

if(length(ok)>2){
    df22<-df2[ok,]
    df22$confidence<-"medium"
    df2<-rbind(df2,df22)
}

corrpe<-c(unlist(by(df2,df2$confidence,FUN = function(x){cor.test(x[,1],x[,2])$estimate})))
corrsp<-suppressWarnings(c(unlist(by(df2,df2$confidence,FUN = function(x){cor.test(x[,1],x[,2],method = "s")$estimate}))))
df3<-data.frame(corrs=c(corrpe,corrsp),stringsAsFactors = F)
df3$type<-rep(c("Pearson","Spearman"),c(length(corrpe),length(corrsp)))
df3$confidence<-factor(c(names(corrpe),names(corrsp)),levels = c("low","medium","high"))
df1$maxpv[df1$maxpv<1e-30]<-1e-30
df1$maxpv[is.na(df1$maxpv)]<-1

axm<-max(abs(c(df1$x_lfc,df1$y_lfc)),na.rm = T)
dfplot<-df1[complete.cases(df1[,1:2]),]
dfplot$symbok<-NA
dfplot$symbok[dfplot$gene_name=="ODC1"]<-"ODC1"
dfplot2<-dfplot
dfplot2$gene_biotype<-NULL
AAA<-ggplot(data = dfplot, aes(x=x_lfc,y=y_lfc,size=maxpv,alpha=maxpv,label=symbok,shape=gene_biotype)) + geom_point() + theme_bw() +
    xlab("TE log2FC\nDegron") +
    ylab("TE log2FC\nsiDDX3X") +
    theme(axis.title.x = element_text(size=22),axis.text.x  = element_text(angle=45, vjust=0.5, size=18)) +
    theme(axis.title.y = element_text(size=24),axis.text.y  = element_text(angle=45, vjust=0.5, size=18)) +
    scale_size_continuous("Adj. p-value (mean)",trans = scales::trans_new(name = "test",transform = trans,inverse = inv, 
                                                                          domain = c(1e-50, Inf)),breaks = c(10e-11,10e-5,10e-3,10e-2)) +
    scale_alpha_continuous("Adj. p-value (mean)",trans = scales::trans_new(name = "test",transform = trans,inverse = inv, 
                                                                           domain = c(1e-10, Inf)),breaks = c(10e-10,10e-3,10e-3,10e-1)) +
    geom_text_repel(size=5,force=32,show.legend = F)

AAA<-AAA + geom_hline(yintercept = 0,size=.5) + geom_vline(xintercept = 0,size=.5) 

suppressWarnings(AAA<-AAA  + geom_smooth(inherit.aes = F,data = dfplot2,mapping = aes(x=x_lfc,y=y_lfc),method = "lm",show.legend = F,se = F,formula = "y ~ x"))

AAA <- AAA  + ylim(-axm,axm) 
pdf(file = "figures/si_degron_plot.pdf",width = 9,height = 5)
AAA
dev.off()


dfa<-list_datasets$DE_mutwt
dfa<-dfa[dfa$keep,]

dfa$padj2<-dfa$padj_RiboRNA
dfa$padj2[is.na(dfa$padj2)]<-1
dfa$padj2[dfa$padj2<10e-50]<-10e-50
# add gene lengths to make rpkm column or like expression values
dfa$experiment<-"R326H vs wt"


axm<-max(abs(c(dfa$log2FoldChange_RNA,dfa$log2FoldChange_Ribo)),na.rm = T)

toview<-c("DDX3X","ODC1","DVL1")
dfa$symb4<-NA
dfa$symb4[dfa$symb%in%toview]<-dfa$symb[dfa$symb%in%toview]

colsi<-alpha(c("blue","dark red","cornflowerblue","firebrick1","dark gray"),c(.8,.8,.8,.8,.3))
names(colsi)<-c("TE_down","TE_up","Concordant_down","Concordant_up","mixed_ns")
dfa$nudgex<-0
dfa$nudgey<-0
dfa$nudgex[which(dfa$symb3=="ODC1")]<-.7
dfa$nudgex[which(dfa$symb3=="DVL1")]<-.7
dfa$nudgey[which(dfa$symb3=="ODC1")]<-(-1)
dfa$nudgey[which(dfa$symb3=="DVL1")]<-(-.2)

dfa$symbok<-NA
dfa$symbok[dfa$gene_name=="ODC1"]<-"ODC1"

a<-ggplot(dfa,aes(x=log2FoldChange_RNA,y=log2FoldChange_Ribo,color=tx_type,size=padj2,label=symbok,shape=gene_biotype)) + geom_point()
a<-a + theme_bw() +
    ylab(paste("Ribo-seq log2FC",sep = "")) +
    xlab(paste("RNA-seq log2FC",sep = "")) +
    theme(axis.title.x = element_text(size=22),axis.text.x  = element_text(angle=45, vjust=0.5, size=18)) +
    theme(axis.title.y = element_text(size=24),axis.text.y  = element_text(angle=45, vjust=0.5, size=18))  +
    theme(strip.background = element_rect(color="black", fill="white", size=1.5, linetype="solid")) +
    theme(strip.text.x = element_text(size = 15,face = "bold"))+
    scale_color_manual(values = colsi,"Tx_class") +
    scale_size_continuous("Adj. p-value (interaction)",trans = scales::trans_new(name = "test",transform = trans,inverse = inv, 
                                                                                 domain = c(1e-100, Inf)),breaks = c(10e-21,10e-11,10e-5,10e-3,1))

a<-a + geom_hline(yintercept = 0,size=.5) + geom_vline(xintercept = 0,size=.5)
ariborna<-a + ylim(-axm,axm) + xlim(-axm,axm)
ariborna<-ariborna + facet_wrap(.~experiment,nrow = 1) 
ariborna<-ariborna + geom_text_repel(size=5,force=32,show.legend = F)

pdf(file = "figures/DE_plot_mut.pdf",width = 9,height = 5)
ariborna
dev.off()

dfas<-dfa[order(dfa$log2FoldChange_RiboRNA,decreasing = T),]
dfas$posit<-1:dim(dfas)[1]

a<-ggplot(dfas,aes(x=posit,y=log2FoldChange_RiboRNA,alpha=padj2,size=padj2,label=symbok)) + geom_point(stroke=0)
a<-a + theme_classic() +
    ylab(paste("TE log2FC",sep = "")) +
    xlab("") +
    theme(axis.title.x = element_text(size=22),axis.text.x  = element_text(angle=45, vjust=0.5, size=18)) +
    theme(axis.title.y = element_text(size=24),axis.text.y  = element_text(angle=45, vjust=0.5, size=18))  +
    theme(strip.background = element_rect(color="black", fill="white", size=1.5, linetype="solid")) +
    theme(strip.text.x = element_text(size = 15,face = "bold"))+
    scale_color_manual(values = "dark grey") +
    scale_size_continuous("Adj. p-value",trans = scales::trans_new(name = "test",transform = trans,inverse = inv, 
                                                                   domain = c(1e-70, 1)),breaks = c(10e-4,10e-3,10e-2,10e-1)) +
    scale_alpha_continuous("Adj. p-value",trans = scales::trans_new(name = "test",transform = trans,inverse = inv, 
                                                                    domain = c(1e-7, 1)),breaks = c(10e-4,10e-3,10e-2,10e-1)) 

a<-a + facet_wrap(.~experiment,nrow = 1) + geom_hline(yintercept = 0,size=.5) + geom_hline(yintercept = c(-1,1),size=.2,lty=2) 
a<-a + geom_text_repel(size=5,force=122,show.legend = F)

pdf(file = "figures/DE_plot_mut_rank.pdf",width = 9,height = 5)
a
dev.off()






dfa<-list_datasets$DE_degron
dfa<-dfa[dfa$keep,]

dfa$padj2<-dfa$padj_RiboRNA
dfa$padj2[is.na(dfa$padj2)]<-1
dfa$padj2[dfa$padj2<10e-50]<-10e-50
# add gene lengths to make rpkm column or like expression values
dfa$experiment<-"IAA vs DMSO"


axm<-max(abs(c(dfa$log2FoldChange_RNA,dfa$log2FoldChange_Ribo)),na.rm = T)

toview<-c("DDX3X","ODC1","DVL1")
dfa$symb4<-NA
dfa$symb4[dfa$symb%in%toview]<-dfa$symb[dfa$symb%in%toview]

colsi<-alpha(c("blue","dark red","cornflowerblue","firebrick1","dark gray"),c(.8,.8,.8,.8,.3))
names(colsi)<-c("TE_down","TE_up","Concordant_down","Concordant_up","mixed_ns")
dfa$nudgex<-0
dfa$nudgey<-0
dfa$nudgex[which(dfa$symb3=="ODC1")]<-.7
dfa$nudgex[which(dfa$symb3=="DVL1")]<-.7
dfa$nudgey[which(dfa$symb3=="ODC1")]<-(-1)
dfa$nudgey[which(dfa$symb3=="DVL1")]<-(-.2)

a<-ggplot(dfa,aes(x=log2FoldChange_RNA,y=log2FoldChange_Ribo,color=tx_type,size=padj2,label=symb4,shape=gene_biotype)) + geom_point()
a<-a + theme_bw() +
    ylab(paste("Ribo-seq log2FC",sep = "")) +
    xlab(paste("RNA-seq log2FC",sep = "")) +
    theme(axis.title.x = element_text(size=22),axis.text.x  = element_text(angle=45, vjust=0.5, size=18)) +
    theme(axis.title.y = element_text(size=24),axis.text.y  = element_text(angle=45, vjust=0.5, size=18))  +
    theme(strip.background = element_rect(color="black", fill="white", size=1.5, linetype="solid")) +
    theme(strip.text.x = element_text(size = 15,face = "bold"))+
    scale_color_manual(values = colsi,"Tx_class") +
    scale_size_continuous("Adj. p-value (interaction)",trans = scales::trans_new(name = "test",transform = trans,inverse = inv, 
                                                                                 domain = c(1e-100, Inf)),breaks = c(10e-21,10e-11,10e-5,10e-3,1))

a<-a + geom_hline(yintercept = 0,size=.5) + geom_vline(xintercept = 0,size=.5)
ariborna<-a + ylim(-axm,axm) + xlim(-axm,axm)
ariborna<-ariborna + facet_wrap(.~experiment,nrow = 1) 
ariborna<-ariborna + geom_text_repel(size=5,force=20,show.legend = F,
                                     nudge_x = dfa$nudgex[which(nchar(dfa$symb4)>2)],
                                     nudge_y = dfa$nudgey[which(nchar(dfa$symb4)>2)])

pdf(file = "figures/DE_plot_degron.pdf",width = 9,height = 5)
ariborna
dev.off()


#####RF


rf_bad<-list_datasets$RF_results_badpreds
datasi<-list_datasets$DE_siDDX3
res_corrs<-lapply(rf_bad$siDDX3,function(x){
    xx<-t(sapply(x,function(y){as.numeric(y$res_tests[1:4])}))
    colnames(xx)<-names(x[[1]]$res_tests[1:4])
    xx
})
ress<-data.frame(do.call(res_corrs,what = rbind))
ress$feat<-rep(names(res_corrs),each=10)
ress<-melt(ress)

ress<-data.frame(do.call(lapply(rf_bad$siDDX3$delta_TE,function(x){
    x$res_tests[[5]]
}),what = rbind))

ress<-data.frame(scale(ress))
ress$padj<-rf_bad$df_list$siDDX3$RiboRNA_padj[match(rownames(ress),rf_bad$df_list$siDDX3$gene_id)]
ress$padj2<-ress$padj
ress$padj2[is.na(ress$padj2)]<-1
ress$padj2[ress$padj2<10e-50]<-10e-50

ress$TPM<-rf_bad$df_list$siDDX3$TPM_RNA[match(rownames(ress),rf_bad$df_list$siDDX3$gene_id)]

rf_teplot<-ggplot(ress,aes(y=predicted,x=original,size=padj2,alpha=padj2)) +  geom_point() +theme_bw()+
    ylab(paste("TE log2FC\n(predicted)",sep = "")) +
    xlab(paste("TE log2FC\n(measured)",sep = "")) +
    theme(axis.title.x = element_text(size=22),axis.text.x  = element_text(angle=45, vjust=0.5, size=18)) +
    theme(axis.title.y = element_text(size=24),axis.text.y  = element_text(angle=45, vjust=0.5, size=18))  +
    theme(strip.background = element_rect(color="black", fill="white", size=1.5, linetype="solid")) +
    theme(strip.text.x = element_text(size = 15,face = "bold"))+
    scale_size_continuous("Adj. p-value",trans = scales::trans_new(name = "test",transform = trans,inverse = inv, 
                                                                   domain = c(1e-100, Inf)),breaks = c(10e-21,10e-11,10e-5,10e-3,1)) +
    scale_alpha_continuous("Adj. p-value",trans = scales::trans_new(name = "test",transform = trans,inverse = inv, 
                                                                    domain = c(1e-100, Inf)),breaks =  c(10e-21,10e-11,10e-5,10e-3)) 

pdf(file = "figures/si_RFbad_plot.pdf",width = 9,height = 5)
rf_teplot + stat_cor(label.x = -8, label.y = 5,show.legend = F)
dev.off()
rf_few<-list_datasets$RF_results_fewpreds
datasi<-list_datasets$DE_siDDX3
res_corrs<-lapply(rf_few$siDDX3,function(x){
    xx<-t(sapply(x,function(y){as.numeric(y$res_tests[1:4])}))
    colnames(xx)<-names(x[[1]]$res_tests[1:4])
    xx
})
ress<-data.frame(do.call(res_corrs,what = rbind))
ress$feat<-rep(names(res_corrs),each=10)
ress<-melt(ress)

ress<-data.frame(do.call(lapply(rf_few$siDDX3$delta_TE,function(x){
    x$res_tests[[5]]
}),what = rbind))

ress<-data.frame(scale(ress))
ress$padj<-rf_few$df_list$siDDX3$RiboRNA_padj[match(rownames(ress),rf_few$df_list$siDDX3$gene_id)]
ress$padj2<-ress$padj
ress$padj2[is.na(ress$padj2)]<-1
ress$padj2[ress$padj2<10e-50]<-10e-50

ress$TPM<-rf_few$df_list$siDDX3$TPM_RNA[match(rownames(ress),rf_few$df_list$siDDX3$gene_id)]

rf_teplot<-ggplot(ress,aes(y=predicted,x=original,size=padj2,alpha=padj2)) +  geom_point() +theme_bw()+
    ylab(paste("TE log2FC\n(predicted)",sep = "")) +
    xlab(paste("TE log2FC\n(measured)",sep = "")) +
    theme(axis.title.x = element_text(size=22),axis.text.x  = element_text(angle=45, vjust=0.5, size=18)) +
    theme(axis.title.y = element_text(size=24),axis.text.y  = element_text(angle=45, vjust=0.5, size=18))  +
    theme(strip.background = element_rect(color="black", fill="white", size=1.5, linetype="solid")) +
    theme(strip.text.x = element_text(size = 15,face = "bold"))+
    scale_size_continuous("Adj. p-value",trans = scales::trans_new(name = "test",transform = trans,inverse = inv, 
                                                                   domain = c(1e-100, Inf)),breaks = c(10e-21,10e-11,10e-5,10e-3,1)) +
    scale_alpha_continuous("Adj. p-value",trans = scales::trans_new(name = "test",transform = trans,inverse = inv, 
                                                                    domain = c(1e-100, Inf)),breaks =  c(10e-21,10e-11,10e-5,10e-3)) 

pdf(file = "figures/si_RFfew_plot.pdf",width = 9,height = 5)
rf_teplot + stat_cor(label.x = -8, label.y = 5,show.legend = F)
dev.off()


rfone<-list_datasets$RF_results
datasi<-list_datasets$DE_siDDX3
res_corrs<-lapply(rfone$siDDX3,function(x){
    xx<-t(sapply(x,function(y){as.numeric(y$res_tests[1:4])}))
    colnames(xx)<-names(x[[1]]$res_tests[1:4])
    xx
})
ress<-data.frame(do.call(res_corrs,what = rbind))
ress$feat<-rep(names(res_corrs),each=10)
ress<-melt(ress)

rf_corplot<-ggplot(ress,aes(y=value,x=feat,group=feat,color=feat)) + 
    geom_jitter(alpha=.06) + 
    stat_summary(position=position_dodge(.6),size=.7) +
    ylab("Correlation") +
    xlab("") +
    #scale_color_manual(values = c("blue","dark red","cornflowerblue","firebrick1","dark gray"),"Tx_class") + 
    theme_classic() +
    theme(axis.title.x = element_text(size=22),axis.text.x  = element_text(angle=45, vjust=0.5, size=13)) +
    theme(axis.title.y = element_text(size=18),axis.text.y  = element_text(angle=0, vjust=0.5, size=15))  +
    theme(strip.text.x = element_text(size=10, face="bold"),strip.text.y = element_text(size=24),strip.background = element_rect(colour="black", fill="white"))+
    facet_wrap(.~variable)


ress<-data.frame(do.call(lapply(rfone$siDDX3$delta_TE,function(x){
    x$res_tests[[5]]
}),what = rbind))

ress<-data.frame(scale(ress))
ress$padj<-rfone$df_list$siDDX3$RiboRNA_padj[match(rownames(ress),rfone$df_list$siDDX3$gene_id)]
ress$padj2<-ress$padj
ress$padj2[is.na(ress$padj2)]<-1
ress$padj2[ress$padj2<10e-50]<-10e-50

ress$TPM<-rfone$df_list$siDDX3$TPM_RNA[match(rownames(ress),rfone$df_list$siDDX3$gene_id)]

rf_teplot<-ggplot(ress,aes(y=predicted,x=original,size=padj2,alpha=padj2)) +  geom_point() +theme_bw()+
    ylab(paste("TE log2FC\n(predicted)",sep = "")) +
    xlab(paste("TE log2FC\n(measured)",sep = "")) +
    theme(axis.title.x = element_text(size=22),axis.text.x  = element_text(angle=45, vjust=0.5, size=18)) +
    theme(axis.title.y = element_text(size=24),axis.text.y  = element_text(angle=45, vjust=0.5, size=18))  +
    theme(strip.background = element_rect(color="black", fill="white", size=1.5, linetype="solid")) +
    theme(strip.text.x = element_text(size = 15,face = "bold"))+
    scale_size_continuous("Adj. p-value",trans = scales::trans_new(name = "test",transform = trans,inverse = inv, 
                                                                   domain = c(1e-100, Inf)),breaks = c(10e-21,10e-11,10e-5,10e-3,1)) +
    scale_alpha_continuous("Adj. p-value",trans = scales::trans_new(name = "test",transform = trans,inverse = inv, 
                                                                    domain = c(1e-100, Inf)),breaks =  c(10e-21,10e-11,10e-5,10e-3)) 

pdf(file = "figures/si_RF_plot.pdf",width = 9,height = 5)
rf_teplot + stat_cor(label.x = -8, label.y = 5,show.legend = F)
dev.off()

ress<-melt(lapply(rfone$siDDX3$delta_TE,function(x){
    x$importance
}))

ress<-ress[ress$Var2=="%IncMSE",]
ress$tocol<-F
ress$tocol[ress$Var1%in%c("base_TE","GCpct_CDS","GCpct_fiveUT","scores_cert")]<-T

ress$Var1<-gsub(ress$Var1,pattern = "scores_pqs" ,replacement = "G4_propensity")
ress$Var1<-gsub(ress$Var1,pattern = "scores_top" ,replacement = "TOP_mRNA")
ress$Var1<-gsub(ress$Var1,pattern = "scores_" ,replacement = "motifscores_")
ress$Var1<-gsub(ress$Var1,pattern = "posdens_base" ,replacement = "Ribo_density")
ress$Var1<-gsub(ress$Var1,pattern = "len" ,replacement = "_length")
ress$Var1<-gsub(ress$Var1,pattern = "base_TE" ,replacement = "baseline_TE")
ress$Var1<-gsub(ress$Var1,pattern = "base_intrex" ,replacement = "intron_exonratio")

ress$Var2<-NULL
ress$Var1<-factor(ress$Var1,levels=unique(ress$Var1))

vars<-unique(ress$Var1)


aaa<-ifelse(!vars%in%c("baseline_TE","GCpct_CDS","GCpct_fiveUT","motifscores_cert"), "black", "blue")
bbb<-ifelse(!vars%in%c("baseline_TE","GCpct_CDS","GCpct_fiveUT","motifscores_cert"), 13, 21)

rf_impor<-ggplot(ress,aes(y=Var1,x=value,color=tocol)) +  stat_summary() + geom_jitter(alpha=.5) +theme_bw()+
    ylab("") +
    xlab("Feature Importance") +
    scale_color_manual(values = c("black","blue"),"") +
    theme(axis.title.x = element_text(size=22),axis.text.x  = element_text(angle=45, vjust=0.5, size=18)) +
    theme(axis.title.y = element_text(size=24),axis.text.y  = element_text(angle=0, vjust=0.5, size=bbb,colour=aaa))  +
    theme(strip.background = element_rect(color="black", fill="white", size=1.5, linetype="solid")) +
    theme(strip.text.x = element_text(size = 15,face = "bold")) + theme(legend.position="none")

pdf(file = "figures/si_RF_import.pdf",width = 8,height = 12)
rf_impor + theme( # remove the vertical grid lines
    panel.grid.major.x = element_blank() ,
    # explicitly set the horizontal lines (or they will disappear too)
    panel.grid.major.y = element_line( size=.1, color="black" ) 
)
dev.off()    


ress<-rfone$df_list$siDDX3
ress[,colnames(ress)[grep(colnames(ress),pattern = "GCpc")]]<-ress[,colnames(ress)[grep(colnames(ress),pattern = "GCpc")]]*100

ress<-ress[,colnames(ress)%in%c("base_TE","GCpct_CDS","GCpct_fiveUT","scores_cert","gene_id")]
ress$tx_type<-"mixed_ns"
ress$tx_type<-datasi$tx_type[match(ress$gene_id,datasi$gene_id)]
ress<-ress[!is.na(ress$tx_type),]
ress<-melt(ress)
ress<-ress[complete.cases(ress$value),]
ress<-ress[!is.infinite(ress$value),]
ress$Var1<-ress$variable
ress$Var1<-gsub(ress$Var1,pattern = "scores_pqs" ,replacement = "G4_propensity")
ress$Var1<-gsub(ress$Var1,pattern = "scores_top" ,replacement = "TOP_mRNA")
ress$Var1<-gsub(ress$Var1,pattern = "scores_" ,replacement = "motifscores_")
ress$Var1<-gsub(ress$Var1,pattern = "posdens_base" ,replacement = "Ribo_density")
ress$Var1<-gsub(ress$Var1,pattern = "len" ,replacement = "_length")
ress$Var1<-gsub(ress$Var1,pattern = "base_TE" ,replacement = "baseline_TE")
ress$Var1<-gsub(ress$Var1,pattern = "base_intrex" ,replacement = "intron_exonratio")

ress$Var1<-factor(ress$Var1,levels=unique(ress$Var1))
colsi<-alpha(c("blue","dark red","cornflowerblue","firebrick1","dark gray"),c(.8,.8,.8,.8,.3))
names(colsi)<-c("TE_down","TE_up","Concordant_down","Concordant_up","mixed_ns")

rf_dens<-ggplot(ress,aes(y=tx_type,x=value,fill=tx_type)) +  stat_density_ridges(quantile_lines = TRUE, quantiles = 2) +
    theme_bw() +
    ylab("") +
    xlab("") +
    scale_fill_manual(values = colsi,"Tx_type",guide = guide_legend(reverse = TRUE)) +
    theme(axis.title.x = element_text(size=22),axis.text.x  = element_text(angle=45, vjust=0.5, size=18)) +
    theme(axis.title.y = element_blank(),axis.text.y  = element_blank())  +
    theme(strip.background = element_rect(color="black", fill="white", size=1.5, linetype="solid")) +
    theme(strip.text.x = element_text(size = 15,face = "bold")) 

pdf(file = "figures/si_RF_feats.pdf",width = 13,height = 6)
rf_dens + facet_wrap(.~Var1,scales = "free")
dev.off() 


#rRNA


minpos_rdna<-3000
maxpos_rdna<-13000

list_clips_rRNA<-list_datasets$rRNA_parclip_cov
list_clips_rRNA_frac<-list_datasets$rRNA_parclip_tc
list_covs_iclip<-list_datasets$rRNA_iclip_cov
signall<-do.call(lapply(list_clips_rRNA,function(x){melt(lapply(x,function(y){as.vector(y[minpos_rdna:maxpos_rdna])}))}),what = rbind)
signall$position<-rep(minpos_rdna:maxpos_rdna,3)
signall$dataset<-sapply(strsplit(rownames(signall),"[.]"),"[[",1)
signall$mutation<-signall$L1

totpos<-aggregate(signall$value,list(signall$position),sum)
maxxa<-max(signall$value)

colrna<-c("grey","orange","navyblue")
pos_18s<-c(3657,5527)
pos_5s<-c(6623,6779)
pos_28s<-c(7935,12969)
parrnapl_reps<-ggplot(signall,aes(x=position,y=value,fill=mutation)) + geom_bar(stat="identity",width=1) +
    facet_wrap(dataset~.,nrow = 2,scales = "free") +
    ylab("Read coverage") +
    xlab("") +
    annotate("rect",xmin = pos_18s[1],xmax = pos_18s[2],ymin = 0,ymax=Inf,color=alpha("black",.5),alpha=0,lty=2) +
    annotate("rect",xmin = pos_5s[1],xmax = pos_5s[2],ymin = 0,ymax=Inf,color=alpha("black",.5),alpha=0,lty=2) +
    annotate("rect",xmin = pos_28s[1],xmax = pos_28s[2],ymin = 0,ymax=Inf,color=alpha("black",.5),alpha=0,lty=2) +
    scale_fill_manual(values = colrna,"Mutations") + 
    theme_classic() +
    theme(axis.title.x = element_text(size=18),axis.text.x  = element_text(angle=0, vjust=0.5, size=18)) +
    theme(axis.title.y = element_text(size=15),axis.text.y  = element_text(angle=0, vjust=0.5, size=15))  +
    theme(strip.text.x = element_text(size=10, face="bold"),strip.text.y = element_text(size=24),strip.background = element_rect(colour="black", fill="white"))  +
    scale_x_continuous(breaks=c(pos_18s[1],(pos_18s[1]+pos_18s[2])/2,pos_18s[2],pos_5s[1],(pos_5s[1]+pos_5s[2])/2,pos_5s[2],pos_28s[1],(pos_28s[1]+pos_28s[2])/2,pos_28s[2]), labels=c("","18S","","","5.8S","","","28S",""))


signall<-do.call(lapply(list_clips_rRNA_frac,function(x){melt(lapply(x,function(y){as.vector(y[minpos_rdna:maxpos_rdna])}))}),what = rbind)
signall$position<-rep(minpos_rdna:maxpos_rdna,3)
signall$dataset<-sapply(strsplit(rownames(signall),"[.]"),"[[",1)
signall$mutation<-signall$L1
signall$totpos<-unname(setNames(totpos$x,totpos$Group.1)[as.character(signall$position)])
signall_cp<-signall
#signall$totpos[match(signall$position,totpos[,1])]<-totpos[match(signall$position,totpos[,1]),2]
maxxa<-max(signall$value)

colrna<-c("orange","navyblue")
pos_18s<-c(3657,5527)
pos_5s<-c(6623,6779)
pos_28s<-c(7935,12969)

signall<-signall[which(signall$totpos>2000),]
#signall<-signall[which(signall$mutation!="None"),]
signall<-signall[which(signall$mutation=="T>C"),];colrna<-c("dark grey","navyblue")
signall$dataset[signall$dataset=="parclip_1"]<-"Replicate_1"
signall$dataset[signall$dataset=="parclip_2"]<-"Replicate_2"
signall$value<-signall$value*100
parrnapl_reps_fr<-ggplot(signall,aes(x=position,y=value)) + geom_bar(stat="identity",width=1) +
    facet_wrap(dataset~.,nrow = 2,scales = "free") +
    ylab("% of PAR-CLIP coverage\nwith T>C conversions ") +
    xlab("") +
    annotate("rect",xmin = pos_18s[1],xmax = pos_18s[2],ymin = 0,ymax=Inf,color=alpha("black",.5),alpha=0,lty=2) +
    annotate("rect",xmin = pos_5s[1],xmax = pos_5s[2],ymin = 0,ymax=Inf,color=alpha("black",.5),alpha=0,lty=2) +
    annotate("rect",xmin = pos_28s[1],xmax = pos_28s[2],ymin = 0,ymax=Inf,color=alpha("black",.5),alpha=0,lty=2) +
    #scale_fill_manual(values = colrna,"Mutations") + 
    theme_classic() +
    theme(axis.title.x = element_text(size=18),axis.text.x  = element_text(angle=0, vjust=0.5, size=18)) +
    theme(axis.title.y = element_text(size=15),axis.text.y  = element_text(angle=0, vjust=0.5, size=15))  +
    theme(strip.text.x = element_text(size=10, face="bold"),strip.text.y = element_text(size=24),strip.background = element_rect(colour="black", fill="white"))  +
    scale_x_continuous(breaks=c(pos_18s[1],(pos_18s[1]+pos_18s[2])/2,pos_18s[2],pos_5s[1],(pos_5s[1]+pos_5s[2])/2,pos_5s[2],pos_28s[1],(pos_28s[1]+pos_28s[2])/2,pos_28s[2]), labels=c("","18S","","","5.8S","","","28S",""))

sign_one<-aggregate(signall$value,list(signall$position),mean)
colnames(sign_one)<-c("position","value")
parrnapl_fr<-ggplot(sign_one,aes(x=position,y=value)) + geom_bar(stat="identity",width=1) +
    #facet_wrap(dataset~.,nrow = 2,scales = "free") +
    ylab("% PAR-CLIP coverage\nT>C reads") +
    xlab("") +
    annotate("rect",xmin = pos_18s[1],xmax = pos_18s[2],ymin = 0,ymax=Inf,color=alpha("black",.5),alpha=0,lty=2) +
    annotate("rect",xmin = pos_5s[1],xmax = pos_5s[2],ymin = 0,ymax=Inf,color=alpha("black",.5),alpha=0,lty=2) +
    annotate("rect",xmin = pos_28s[1],xmax = pos_28s[2],ymin = 0,ymax=Inf,color=alpha("black",.5),alpha=0,lty=2) +
    #scale_fill_manual(values = colrna,"Mutations") + 
    theme_classic() +
    theme(axis.title.x = element_text(size=18),axis.text.x  = element_text(angle=0, vjust=0.5, size=18)) +
    theme(axis.title.y = element_text(size=15),axis.text.y  = element_text(angle=0, vjust=0.5, size=15))  +
    theme(strip.text.x = element_text(size=10, face="bold"),strip.text.y = element_text(size=24),strip.background = element_rect(colour="black", fill="white"))  +
    scale_x_continuous(breaks=c(pos_18s[1],(pos_18s[1]+pos_18s[2])/2,pos_18s[2],pos_5s[1],(pos_5s[1]+pos_5s[2])/2,pos_5s[2],pos_28s[1],(pos_28s[1]+pos_28s[2])/2,pos_28s[2]), labels=c("","18S","","","5.8S","","","28S",""))



aa<-cbind(signall[signall$mutation=="T>C",]$position,signall[signall$mutation=="T>C",]$value)



covicl<-list_covs_iclip$DDX3_WT_iCLIP_rep1$by_rl$all_reads[[1]]+list_covs_iclip$DDX3_WT_iCLIP_rep2$by_rl$all_reads[[1]][]
signall<-melt(covicl[minpos_rdna:maxpos_rdna])
signall$position<-rep(minpos_rdna:maxpos_rdna)
signall<-signall[signall$position>=min(parrnapl_reps_fr$data$position),]
icl<-ggplot(signall,aes(x=position,y=value)) + geom_bar(stat="identity",width=1) +
    ylab("iCLIP coverage") +
    xlab("") +
    annotate("rect",xmin = pos_18s[1],xmax = pos_18s[2],ymin = 0,ymax=Inf,color=alpha("black",.5),alpha=0,lty=2) +
    annotate("rect",xmin = pos_5s[1],xmax = pos_5s[2],ymin = 0,ymax=Inf,color=alpha("black",.5),alpha=0,lty=2) +
    annotate("rect",xmin = pos_28s[1],xmax = pos_28s[2],ymin = 0,ymax=Inf,color=alpha("black",.5),alpha=0,lty=2) +
    theme_classic() +
    theme(axis.title.x = element_text(size=18),axis.text.x  = element_text(angle=0, vjust=0.5, size=18)) +
    theme(axis.title.y = element_text(size=15),axis.text.y  = element_text(angle=0, vjust=0.5, size=15))  +
    theme(strip.text.x = element_text(size=10, face="bold"),strip.text.y = element_text(size=24),strip.background = element_rect(colour="black", fill="white"))  +
    scale_x_continuous(breaks=c(pos_18s[1],(pos_18s[1]+pos_18s[2])/2,pos_18s[2],pos_5s[1],(pos_5s[1]+pos_5s[2])/2,pos_5s[2],pos_28s[1],(pos_28s[1]+pos_28s[2])/2,pos_28s[2]), labels=c("","18S","","","5.8S","","","28S",""))




pdf(file = "figures/rDNA_clips.pdf",width = 8,height = 4)
plot_grid(NULL,parrnapl_fr,icl,ncol = 1,align = "v",rel_heights =  c(.3,1,1))
dev.off()


pdf(file = "figures/suppl_rDNA_par_clips.pdf",width = 8,height = 4)
parrnapl_reps_fr
dev.off()

signall_sup<-signall_cp
signall_sup$position2<-"other"

signall_sup$position2[signall_sup$position%in%4173:4220]<-"h16"
signall_sup<-signall_sup[which(signall_sup$mutation=="T>C"),];colrna<-c("dark grey","navyblue")
signall_sup$value<-signall_sup$value*100
parrnapl_fr2<-ggplot(signall_sup,aes(x=totpos,y=value,color=position2,shape=dataset)) + geom_point(size=2.5) +
    #facet_wrap(dataset~.,nrow = 2,scales = "free") +
    ylab("% PAR-CLIP coverage\nT>C reads") +
    xlab("Read coverage") +
    theme_classic() +
    scale_color_manual(values = alpha(rev(c("grey","black")),.8),"") +
    theme(axis.title.x = element_text(size=18),axis.text.x  = element_text(angle=0, vjust=0.5, size=18)) +
    theme(axis.title.y = element_text(size=15),axis.text.y  = element_text(angle=0, vjust=0.5, size=15))  +
    theme(strip.text.x = element_text(size=10, face="bold"),strip.text.y = element_text(size=24),strip.background = element_rect(colour="black", fill="white"))  
pdf(file = "figures/suppl_rDNA_coveragepar.pdf",width = 8,height = 4)
parrnapl_fr2+xlim(0,10000)
dev.off()
#PARCLIP mRNA


peak_gr_ok<-list_datasets$peaks_gen
ddd<-data.frame(mcols(peak_gr_ok))
ddd$regions<-factor(ddd$region,levels = names(sort(table(ddd$region),decreasing = T)))
colssi<-c("red","orange","gold","cornflowerblue","blue","dark blue")
barpeaks<-ggplot(ddd,aes(x=regions,fill=regions)) + geom_bar() +
    scale_fill_manual(values = colssi) +
    theme_classic() +
    ylab("Number of peaks") +
    xlab("") +
    theme(axis.title.x = element_text(size=22),axis.text.x  = element_text(angle=45, vjust=0.5, size=15)) +
    theme(axis.title.y = element_text(size=18),axis.text.y  = element_text(angle=0, vjust=0.5, size=15))  +
    theme(strip.text.x = element_text(size=10, face="bold"),strip.text.y = element_text(size=24),strip.background = element_rect(colour="black", fill="white")) + 
    theme(legend.position="none")

pdf(file = "figures/peak_stats_ddx3.pdf",width = 8,height = 4)
print(barpeaks)
dev.off()


maptx<-list_datasets$peaks_tx

bin5<-15
bin_c<-30
bin_u<-20


aggg<-aggregate(maptx$score,list(as.vector(seqnames(maptx))),sum)
agg_sc<-aggg[,2]
names(agg_sc)<-aggg[,1]
maptx$txreads_sum<-agg_sc[as.vector(seqnames(maptx))]
maptx$norm_cov_sum<-maptx$score/maptx$txreads_sum

maptx$tx_bin_ok<-maptx$tx_bin

maptx$tx_bin_ok[maptx$tx_position=="CDS"]<-maptx$tx_bin[maptx$tx_position=="CDS"]+bin5
maptx$tx_bin_ok[maptx$tx_position=="3_UTR"]<-maptx$tx_bin[maptx$tx_position=="3_UTR"]+(bin5+bin_c)

maptx2<-maptx
maptx2<-maptx2[maptx2$txreads_sum>10]
agga2<-aggregate(list(maptx2$score,maptx2$norm_cov_sum,maptx2$conv_spec),by=list(maptx2$tx_bin_ok),sum)
colnames(agga2)<-c("tx_bin","score","norm_cov_sum","conv_spec")


aggplot<-ggplot(agga2,aes(x=tx_bin,y=norm_cov_sum)) + 
    geom_line() +  
    ylab("Aggregate peak coverage") +
    xlab("") +
    geom_vline(xintercept =bin5+1,col="black",lty=3)+
    #scale_color_manual(values = alpha(c("black","cornflowerblue","orange","blue","dark red"),.7),"Tx_class") + 
    geom_vline(xintercept =bin5+bin_c+1,col="black",lty=3)+ 
    theme_classic() +
    theme(axis.title.x = element_text(size=22),axis.text.x  = element_text(angle=0, vjust=0.5, size=13)) +
    theme(axis.title.y = element_text(size=18),axis.text.y  = element_text(angle=0, vjust=0.5, size=15))  +
    theme(strip.text.x = element_text(size=10, face="bold"),strip.text.y = element_text(size=24),strip.background = element_rect(colour="black", fill="white")) +
    scale_x_continuous(breaks=c(0,bin5+1,bin5+bin_c+1,bin5+bin_c+bin_u+1), labels=c("TSS","start codon","stop codon","TES"))
pdf(file = "figures/peak_aggregate.pdf",width = 8,height = 4)
aggplot
dev.off()


maptx2<-maptx
maptx2<-maptx2[maptx2$txreads_sum>10]

agga2<-data.frame(mcols(maptx2),stringsAsFactors = F)
conv_specplot1<-ggplot(agga2,aes(x=tx_bin_ok,y=conv_spec)) + 
    geom_smooth(method="loess",span=.2,se = F) +  
    stat_summary() + 
    ylab("T>C conversion specificity") +
    xlab("") +
    geom_vline(xintercept =bin5+1,col="black",lty=3)+
    #scale_color_manual(values = alpha(c("black","cornflowerblue","orange","blue","dark red"),.7),"Tx_class") + 
    geom_vline(xintercept =bin5+bin_c+1,col="black",lty=3)+ 
    theme_classic() +
    theme(axis.title.x = element_text(size=22),axis.text.x  = element_text(angle=0, vjust=0.5, size=13)) +
    theme(axis.title.y = element_text(size=18),axis.text.y  = element_text(angle=0, vjust=0.5, size=15))  +
    theme(strip.text.x = element_text(size=10, face="bold"),strip.text.y = element_text(size=24),strip.background = element_rect(colour="black", fill="white")) +
    scale_x_continuous(breaks=c(0,bin5+1,bin5+bin_c+1,bin5+bin_c+bin_u+1), labels=c("TSS","start codon","stop codon","TES"))

pdf(file = "figures/peak_meanconv.pdf",width = 8,height = 4)
conv_specplot1
dev.off()


maptx2<-maptx
maptx2<-maptx2[maptx2$txreads_sum>10]
agga2<-data.frame(mcols(maptx2),stringsAsFactors = F)
agga2$tx_position[agga2$dist_start%in%(-25:25)]<-"Start_codon"
agga2$tx_position<-factor(agga2$tx_position,levels=c("5_UTR","Start_codon","CDS","3_UTR"))


btps<-unique(agga2$tx_type)
comps<-unname(split(cbind(rep("mixed_ns",length(btps)-1),as.character(btps[-1])),f = 1:(length(btps)-1)))
comps<-comps[c(1,3,4,2)]
conv_specplot3<-ggplot(agga2,aes(y=conv_spec,x=tx_type,group=tx_type,color=tx_type)) + 
    geom_jitter(alpha=.06) + 
    stat_summary(position=position_dodge(.6),size=.7) +
    ylab("T>C conversion specificity") +
    xlab("") +
    geom_vline(xintercept =bin5+1,col="black",lty=3)+
    scale_color_manual(values = c("blue","dark red","cornflowerblue","firebrick1","dark gray"),"Tx_class") + 
    geom_vline(xintercept =bin5+bin_c+1,col="black",lty=3)+ 
    theme_classic() +
    theme(axis.title.x = element_text(size=22),axis.text.x  = element_text(angle=45, vjust=0.5, size=13)) +
    theme(axis.title.y = element_text(size=18),axis.text.y  = element_text(angle=0, vjust=0.5, size=15))  +
    theme(strip.text.x = element_text(size=10, face="bold"),strip.text.y = element_text(size=24),strip.background = element_rect(colour="black", fill="white")) 
conv_specplot3<-conv_specplot3 + stat_compare_means(label="p.signif",method = "wilcox.test",hide.ns = F,comparisons = comps,label.y = .43 + seq(0,.1,length.out = 4)) + facet_wrap(.~tx_position,ncol = 4) 
pdf(file = "figures/peak_diffconv.pdf",width = 8,height = 4)
conv_specplot3+ylim(0,0.55) + theme(axis.title.x=element_blank(),
                                    axis.text.x=element_blank(),
                                    axis.ticks.x=element_blank())
dev.off()


res_clips_meta_ok<-list_datasets$peaks_parclips
prots=c("DDX3X","EIF3B","FMR1","MOV10")
dat<-melt(data.frame(res_clips_meta_ok))
dat$position<-rep(1:(bin5+bin_c+bin_u),dim(res_clips_meta_ok)[2])
dat2<-dat[dat$variable%in%prots,]
dat2$variable<-factor(dat2$variable,prots)
covpl<-ggplot(dat2,aes(x=position,y=value,group=variable,color=variable)) +
    geom_vline(xintercept =bin5+1,col="black",lty=3)+
    scale_color_manual(values = alpha(c("red","orange","forestgreen","blue"),.7),"RBP") + 
    geom_line(size=1.1) +
    geom_vline(xintercept =bin5+bin_c+1,col="black",lty=3)+ 
    theme_classic() +
    ylab("Aggregate\nPAR-CLIP peak score") +
    xlab("") +
    theme(axis.title.x = element_text(size=22),axis.text.x  = element_text(angle=0, vjust=0.5, size=13)) +
    theme(axis.title.y = element_text(size=18),axis.text.y  = element_text(angle=0, vjust=0.5, size=15))  +
    theme(strip.text.x = element_text(size=10, face="bold"),strip.text.y = element_text(size=24),strip.background = element_rect(colour="black", fill="white")) +
    scale_x_continuous(breaks=c(0,bin5+1,bin5+bin_c+1,bin5+bin_c+bin_u+1), labels=c("TSS","start codon","stop codon","TES"))

pdf(file = "figures/peak_aggregate_PARCLIPs.pdf",width = 8,height = 4)
print(covpl)
dev.off()


res_clips_meta_ok<-list_datasets$peaks_eclips
prots=c("DDX3X.K562","GEMIN5.K562","RPS3.K562","DDX6.K562")
dat<-melt(data.frame(res_clips_meta_ok))
dat$position<-rep(1:(bin5+bin_c+bin_u),dim(res_clips_meta_ok)[2])
dat2<-dat[dat$variable%in%prots,]
dat2$variable<-factor(dat2$variable,prots)
covpl<-ggplot(dat2,aes(x=position,y=value,group=variable,color=variable)) +
    geom_vline(xintercept =bin5+1,col="black",lty=3)+
    scale_color_manual(values = alpha(c("red","orange","forestgreen","blue"),.7),"RBP") + 
    geom_line(size=1.1) +
    geom_vline(xintercept =bin5+bin_c+1,col="black",lty=3)+ 
    theme_classic() +
    ylab("Aggregate\nPAR-CLIP peak score") +
    xlab("") +
    theme(axis.title.x = element_text(size=22),axis.text.x  = element_text(angle=0, vjust=0.5, size=13)) +
    theme(axis.title.y = element_text(size=18),axis.text.y  = element_text(angle=0, vjust=0.5, size=15))  +
    theme(strip.text.x = element_text(size=10, face="bold"),strip.text.y = element_text(size=24),strip.background = element_rect(colour="black", fill="white")) +
    scale_x_continuous(breaks=c(0,bin5+1,bin5+bin_c+1,bin5+bin_c+bin_u+1), labels=c("TSS","start codon","stop codon","TES"))
pdf(file = "figures/peak_aggregate_ECLIPs.pdf",width = 8,height = 4)
print(covpl)
dev.off()


seqse<-list_datasets$peaks_seq
certt<-list_datasets$peaks_cert
pqs<-list_datasets$peaks_pqs
de_tab<-list_datasets$DE_siDDX3

aa<-list_datasets$peaks_shuf

nt_freqs_g<-t(sapply(seqse,function(x){letterFrequencyInSlidingView(x,view.width = 1,letters = "G")}))
nt_freqs_c<-t(sapply(seqse,function(x){letterFrequencyInSlidingView(x,view.width = 1,letters = "C")}))


agga2<-melt(nt_freqs_g)
colnames(agga2)<-c("tx","position","value")
agga2$type="G_freq"
agga3<-melt(nt_freqs_c)
colnames(agga3)<-c("tx","position","value")
agga3$type="C_freq"
agga2<-rbind(agga2,agga3)
agga2$position<-agga2$position-50


dgplot1<-ggplot(agga2,aes(x=position,y=value,color=type,group=type)) + 
    geom_smooth(method="loess",span=.2,se = F) +  
    stat_summary() + 
    ylab("nt frequency") +
    xlab("Distance from peak summit") +
    
    scale_color_manual(values = alpha(c("blue","orange"),.7),"nt") + 
    
    theme_classic() +
    theme(axis.title.x = element_text(size=22),axis.text.x  = element_text(angle=0, vjust=0.5, size=13)) +
    theme(axis.title.y = element_text(size=18),axis.text.y  = element_text(angle=0, vjust=0.5, size=15))  +
    theme(strip.text.x = element_text(size=10, face="bold"),strip.text.y = element_text(size=24),strip.background = element_rect(colour="black", fill="white")) 

pdf(file = "figures/peak_gc_ddx3.pdf",width = 8,height = 4)
print(dgplot1)
dev.off()

aa2<-list_datasets$peaks_shuf

dg<-pqs
mt<-as(dg,"matrix")

dg2<-certt
mt2<-as(dg2,"matrix")

agga2<-melt(mt)
colnames(agga2)<-c("tx","position","value")
agga2$value<-scale(agga2$value)
agga2$type="G4 propensity"
agga3<-melt(mt2)
colnames(agga3)<-c("tx","position","value")
agga3$value<-scale(agga3$value)
agga3$type="CERT motif scores"
agga2<-rbind(agga2,agga3)
agga2$position<-agga2$position-50
agga2$control<-"positive"
agga2$control[agga2$tx>length(which(nchar(names(aa2))>2))]<-"shuffled"

agga2$type[agga2$type=="CERT motif score" & agga2$control=="shuffled"]="CERT motif score\n(shuffled)"
agga2$type[agga2$type=="G4 propensity" & agga2$control=="shuffled"]="G4 propensity\n(shuffled)"
agga2$motiff<-"CERT motif score"
agga2$motiff[grep(agga2$type,pattern = "G4")]<-"G4 propensity"

dgplot2_shuf<-ggplot(agga2,aes(x=position,y=value,color=control,group=control)) + 
    geom_smooth(method="loess",span=.2,se = F) +
    stat_summary() + 
    ylab("z-score") +
    xlab("Distance from peak summit") +
    
    #scale_color_manual(values = alpha(c("dark grey","grey","dark grey","grey"),c(.8,.5,.8,.5)),"") +
    scale_color_manual(values = alpha(c("dark blue","grey"),c(.9,.8)),"") + 
    
    theme_classic() +
    theme(axis.title.x = element_text(size=22),axis.text.x  = element_text(angle=0, vjust=0.5, size=13)) +
    theme(axis.title.y = element_text(size=18),axis.text.y  = element_text(angle=0, vjust=0.5, size=15))  +
    theme(strip.text.x = element_text(size=10, face="bold"),strip.text.y = element_text(size=24),strip.background = element_rect(colour="black", fill="white")) 

pdf(file = "figures/peak_motif_shuf.pdf",width = 13,height = 4)

dgplot2_shuf + facet_wrap(.~motiff,nrow = 1,scales="free")
dev.off()



agga2<-melt(as(list_datasets$peaks_deltaG,"matrix"))
colnames(agga2)<-c("tx","position","value")
agga2$position<-agga2$position-50
dgplot3<-ggplot(agga2,aes(x=position,y=value)) + 
    geom_smooth(method="loess",span=.2,se = F) +  
    stat_summary() + 
    ylab(paste(expression(Delta),"G",sep = "")) +
    xlab("Distance from peak summit") +
    
    #scale_color_manual(values = alpha(c("blue","orange"),.7),"nt") + 
    
    theme_classic() +
    theme(axis.title.x = element_text(size=22),axis.text.x  = element_text(angle=0, vjust=0.5, size=13)) +
    theme(axis.title.y = element_text(size=18),axis.text.y  = element_text(angle=0, vjust=0.5, size=15))  +
    theme(strip.text.x = element_text(size=10, face="bold"),strip.text.y = element_text(size=24),strip.background = element_rect(colour="black", fill="white")) 

pdf(file = "figures/peak_deltaG_ddx3.pdf",width = 8,height = 4)
print(dgplot3)
dev.off()


dfa<-list_datasets$DE_siDDX3
dfa<-dfa[dfa$keep,]
datasi<-dfa[,c("gene_id","symb","gene_biotype","log2FoldChange_RiboRNA","padj_RiboRNA")]
dfadeg<-list_datasets$DE_degron
dfadeg<-dfadeg[dfadeg$gene_id%in%datasi$gene_id,]
rownames(dfadeg)<-dfadeg$gene_id
rownames(datasi)<-datasi$gene_id
datasi<-datasi[rownames(dfadeg),]
dfadeg$padj_RiboRNA_si<-datasi$padj_RiboRNA
dfadeg$log2FoldChange_RiboRNA_si<-datasi$log2FoldChange_RiboRNA
dfadeg$pmin<-apply(dfadeg[,c("padj_RiboRNA","padj_RiboRNA_si")],1,FUN = min)
dfadeg$pmax<-apply(dfadeg[,c("padj_RiboRNA","padj_RiboRNA_si")],1,FUN = max)


df2<-dfadeg[,c("log2FoldChange_RiboRNA","log2FoldChange_RiboRNA_si","padj_RiboRNA","padj_RiboRNA_si","symb","gene_biotype")]

dfcp<-df2
nimok<-colnames(df2)[1:2]
trans <- function(x) -log(x, 10)
inv <- function(x) 10^(-x)
nmsok<-c("x_lfc","y_lfc","padj_RiboRNA","padj_RiboRNA_si","gene_name","gene_biotype")
colnames(df2)<-nmsok
df2$maxpv<-apply(df2[,3:4,drop=F],1,function(x){
    if(sum(is.na(x))==length(x)){return(NA)}
    else{return(mean(x,na.rm = T))}
})
df2<-df2[complete.cases(df2[,1:2]),]
df2<-df2[!is.infinite(df2[,1]),]
df2<-df2[!is.infinite(df2[,2]),]

df1<-df2

df2$confidence<-"low"

ok<-which(df2$maxpv<.5)
df2$confidence[which(df2$maxpv<.5)]<-"medium"

if(length(ok)>2){
    df22<-df2[ok,]
    df22$confidence<-"medium"
    df2<-rbind(df2,df22)
}

df2$confidence<-factor(df2$confidence,levels=c("low","medium","high"))

ok<-which(df2$maxpv<.2)
df2$confidence[ok]<-"high"
if(length(ok)>2){
    df22<-df2[ok,]
    df22$confidence<-"high"
    df2<-rbind(df2,df22)
}


corrpe<-c(unlist(by(df2,df2$confidence,FUN = function(x){cor.test(x[,1],x[,2])$estimate})))
corrsp<-suppressWarnings(c(unlist(by(df2,df2$confidence,FUN = function(x){cor.test(x[,1],x[,2],method = "s")$estimate}))))
df3<-data.frame(corrs=c(corrpe,corrsp),stringsAsFactors = F)
df3$type<-rep(c("Pearson","Spearman"),c(length(corrpe),length(corrsp)))
df3$confidence<-factor(c(names(corrpe),names(corrsp)),levels = c("low","medium","high"))
df1$maxpv[df1$maxpv<1e-30]<-1e-30
df1$maxpv[is.na(df1$maxpv)]<-1

axm<-max(abs(c(df1$x_lfc,df1$y_lfc)),na.rm = T)
dfplot<-unique(df1[complete.cases(df1[,1:2]),])
dfplot$symbok<-NA
dfplot$symbok[dfplot$gene_name=="ODC1"]<-"ODC1"
dfplot2<-df2
dfplot2$gene_biotype<-NULL

dfplot$confidence="low"
dfplot$confidence[dfplot$maxpv<.5]="medium"
dfplot$confidence[dfplot$maxpv<.2]="high"
dfplot$confidence<-factor(dfplot$confidence,levels=c("low","medium","high"))

dfplot$pvalue_cutoff=">.5"
dfplot$pvalue_cutoff[dfplot$maxpv<.5]="<.5"
dfplot$pvalue_cutoff[dfplot$maxpv<.2]="<.2"
dfplot$pvalue_cutoff<-factor(dfplot$pvalue_cutoff,levels=c(">.5","<.5","<.2"))


AAA<-ggplot(data = dfplot, aes(x=x_lfc,y=y_lfc,size=maxpv,alpha=maxpv,label=symbok)) + geom_point() + theme_bw() +
    xlab("TE log2FC\nDegron") +
    ylab("TE log2FC\nsiDDX3X") +
    theme(axis.title.x = element_text(size=22),axis.text.x  = element_text(angle=45, vjust=0.5, size=18)) +
    theme(axis.title.y = element_text(size=24),axis.text.y  = element_text(angle=45, vjust=0.5, size=18)) +
    scale_size_continuous("Adj. p-value (mean)",trans = scales::trans_new(name = "test",transform = trans,inverse = inv, 
                                                                          domain = c(1e-50, Inf)),breaks = c(10e-11,10e-5,10e-3,10e-2)) +
    scale_alpha_continuous("Adj. p-value (mean)",trans = scales::trans_new(name = "test",transform = trans,inverse = inv, 
                                                                           domain = c(1e-10, Inf)),breaks = c(10e-10,10e-3,10e-3,10e-1)) +
    geom_text_repel(size=5,force=32,show.legend = F)

AAA<-AAA + geom_hline(yintercept = 0,size=.5) + geom_vline(xintercept = 0,size=.5) 

suppressWarnings(AAA<-AAA  + geom_smooth(inherit.aes = F,data = dfplot,mapping = aes(x=x_lfc,y=y_lfc,group=pvalue_cutoff,color=pvalue_cutoff),method = "lm",show.legend = F,se = T,formula = "y ~ x"))

AAA <- AAA  + ylim(-axm,axm) 
AAA<-AAA + stat_cor(method = "p",inherit.aes = F,data = dfplot,aes(x=x_lfc,y=y_lfc,group=pvalue_cutoff,color=pvalue_cutoff),label.x = c(-6,-6,-6), label.y = c(1,2,3),show.legend = T)
AAA<-AAA + scale_color_manual(values = rev(c("dark blue","steelblue","grey12")),"pvalue_cutoff")
pdf(file = "figures/si_degron_plot2.pdf",width = 9,height = 5)
AAA
dev.off()



#SUPPL RF


# 
# rfone<-get(load("/bfd/lcalviel/data/riboseq/new_ddx3/results/siRNA_annotated/siRNA_annotated_results_DEDEX_res_RFreduced_onedelta"))
# #datasi<-list_datasets$DE_siDDX3
# res_corrs<-lapply(rfone$siDDX3,function(x){
#     xx<-t(sapply(x,function(y){as.numeric(y$res_tests[1:4])}))
#     colnames(xx)<-names(x[[1]]$res_tests[1:4])
#     xx
# })
# ress<-data.frame(do.call(res_corrs,what = rbind))
# ress$feat<-rep(names(res_corrs),each=10)
# ress<-melt(ress)
# 
# rf_corplot<-ggplot(ress,aes(y=value,x=feat,group=feat,color=feat)) + 
#     geom_jitter(alpha=.06) + 
#     stat_summary(position=position_dodge(.6),size=.7) +
#     ylab("Correlation") +
#     xlab("") +
#     #scale_color_manual(values = c("blue","dark red","cornflowerblue","firebrick1","dark gray"),"Tx_class") + 
#     theme_classic() +
#     theme(axis.title.x = element_text(size=22),axis.text.x  = element_text(angle=45, vjust=0.5, size=13)) +
#     theme(axis.title.y = element_text(size=18),axis.text.y  = element_text(angle=0, vjust=0.5, size=15))  +
#     theme(strip.text.x = element_text(size=10, face="bold"),strip.text.y = element_text(size=24),strip.background = element_rect(colour="black", fill="white"))+
#     facet_wrap(.~variable)
# 
# 
# ress<-data.frame(do.call(lapply(rfone$siDDX3$delta_TE,function(x){
#     x$res_tests[[5]]
# }),what = rbind))
# 
# ress<-data.frame(scale(ress))
# ress$padj<-rfone$df_list$siDDX3$RiboRNA_padj[match(rownames(ress),rfone$df_list$siDDX3$gene_id)]
# ress$padj2<-ress$padj
# ress$padj2[is.na(ress$padj2)]<-1
# ress$padj2[ress$padj2<10e-50]<-10e-50
# 
# ress$TPM<-rfone$df_list$siDDX3$TPM_RNA[match(rownames(ress),rfone$df_list$siDDX3$gene_id)]
# 
# rf_teplot<-ggplot(ress,aes(y=predicted,x=original,size=padj2,alpha=padj2)) +  geom_point() +theme_bw()+
#     ylab(paste("TE log2FC\n(predicted)",sep = "")) +
#     xlab(paste("TE log2FC\n(measured)",sep = "")) +
#     theme(axis.title.x = element_text(size=22),axis.text.x  = element_text(angle=45, vjust=0.5, size=18)) +
#     theme(axis.title.y = element_text(size=24),axis.text.y  = element_text(angle=45, vjust=0.5, size=18))  +
#     theme(strip.background = element_rect(color="black", fill="white", size=1.5, linetype="solid")) +
#     theme(strip.text.x = element_text(size = 15,face = "bold"))+
#     scale_size_continuous("Adj. p-value",trans = scales::trans_new(name = "test",transform = trans,inverse = inv, 
#                                                                    domain = c(1e-100, Inf)),breaks = c(10e-21,10e-11,10e-5,10e-3,1)) +
#     scale_alpha_continuous("Adj. p-value",trans = scales::trans_new(name = "test",transform = trans,inverse = inv, 
#                                                                     domain = c(1e-100, Inf)),breaks =  c(10e-21,10e-11,10e-5,10e-3)) 
# 
# pdf(file = "figures/si_RFred_plot.pdf",width = 9,height = 5)
# rf_teplot + stat_cor(label.x = -8, label.y = 5,show.legend = F)
# dev.off()
# 
# ress<-melt(lapply(rfone$siDDX3$delta_TE,function(x){
#     x$importance
# }))
# 
# ress<-ress[ress$Var2=="%IncMSE",]
# ress$tocol<-F
# ress$tocol[ress$Var1%in%c("base_TE","GCpct_CDS","GCpct_fiveUT","scores_cert")]<-T
# 
# ress$Var1<-gsub(ress$Var1,pattern = "scores_pqs" ,replacement = "G4_propensity")
# ress$Var1<-gsub(ress$Var1,pattern = "scores_top" ,replacement = "TOP_mRNA")
# ress$Var1<-gsub(ress$Var1,pattern = "scores_" ,replacement = "motifscores_")
# ress$Var1<-gsub(ress$Var1,pattern = "posdens_base" ,replacement = "Ribo_density")
# ress$Var1<-gsub(ress$Var1,pattern = "len" ,replacement = "_length")
# ress$Var1<-gsub(ress$Var1,pattern = "base_TE" ,replacement = "baseline_TE")
# ress$Var1<-gsub(ress$Var1,pattern = "base_intrex" ,replacement = "intron_exonratio")
# 
# ress$Var2<-NULL
# ress$Var1<-factor(ress$Var1,levels=unique(ress$Var1))
# 
# vars<-unique(ress$Var1)
# 
# 
# aaa<-ifelse(!vars%in%c("baseline_TE","GCpct_CDS","GCpct_fiveUT","motifscores_cert"), "black", "blue")
# bbb<-ifelse(!vars%in%c("baseline_TE","GCpct_CDS","GCpct_fiveUT","motifscores_cert"), 13, 21)
# 
# rf_impor<-ggplot(ress,aes(y=Var1,x=value,color=tocol)) +  stat_summary() + geom_jitter(alpha=.5) +theme_bw()+
#     ylab("") +
#     xlab("Feature Importance") +
#     scale_color_manual(values = c("black","blue"),"") +
#     theme(axis.title.x = element_text(size=22),axis.text.x  = element_text(angle=45, vjust=0.5, size=18)) +
#     theme(axis.title.y = element_text(size=24),axis.text.y  = element_text(angle=0, vjust=0.5, size=bbb,colour=aaa))  +
#     theme(strip.background = element_rect(color="black", fill="white", size=1.5, linetype="solid")) +
#     theme(strip.text.x = element_text(size = 15,face = "bold")) + theme(legend.position="none")
# 
# pdf(file = "figures/si_RFred_import.pdf",width = 8,height = 12)
# rf_impor + theme( # remove the vertical grid lines
#     panel.grid.major.x = element_blank() ,
#     # explicitly set the horizontal lines (or they will disappear too)
#     panel.grid.major.y = element_line( size=.1, color="black" ) 
# )
# dev.off()    

