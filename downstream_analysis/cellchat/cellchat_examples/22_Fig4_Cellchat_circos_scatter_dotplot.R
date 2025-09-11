library(Seurat)
library(ggplot2)
library(dplyr)
library(Polychrome) #color palettes
library(cowplot) #stitch together plots
library(ggpubr) #rotate axis options
library(pals) #more color pallets
library(eyedroppeR)
library(ggalt)
library(CellChat)
library(patchwork)
library(circlize)
library(CCPlotR)


##MERFISH neighborhoods scatter plot circos plot unfiltered

object.list.merfish <- readRDS("20240105_combined_NN_IAN_MERFISH_trimmed_object.list.rds")

##scatter plot
for (i in 1:length(object.list.merfish)) {
  num.link <- rowSums(object.list.merfish[[i]]@net$count) + colSums(object.list.merfish[[i]]@net$count)-diag(object.list.merfish[[i]]@net$count)
  weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
  gg <- list()
  p1 <-  netAnalysis_signalingRole_scatter(object.list.merfish[[i]], title = names(object.list.merfish)[i], weight.MinMax = weight.MinMax, color.use = celltype_colors, show.legend = F)
  pdf(paste0("20250105_",names(object.list.merfish)[i],"_scatter_sender_receiver_MERFISH.pdf"), height =3.75, width = 3.75)
  print(p1)
  dev.off()
}

##Circos Plots
for (i in 1:length(object.list.merfish)) {
  strwidth <- function(x) {0.5}
  pathways <- object.list.merfish[[i]]@netP$pathways
  p1 <- netVisual_chord_cell(object.list.merfish[[i]], signaling = pathways, thresh = 0.0001, remove.isolate = T, title.name = NULL, show.legend = F, color.use = celltype_colors, lab.cex = 1.2)
  pdf(paste("Neighborhood.chordmap_MERFISH_",names(object.list.merfish)[i],".pdf"), height = 6, width = 6)
  print(p1)
  dev.off()  
}



##SC neighborhoods scatter plot circos plot unfiltered

object.list.sc <- readRDS("20240105_combined_NN_IAN_sc_trimmed_object.list.rds")

##scatter plot
for (i in 1:length(object.list.sc)) {
  num.link <- rowSums(object.list.sc[[i]]@net$count) + colSums(object.list.sc[[i]]@net$count)-diag(object.list.sc[[i]]@net$count)
  weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
  gg <- list()
  p1 <-  netAnalysis_signalingRole_scatter(object.list.sc[[i]], title = names(object.list.sc)[i], weight.MinMax = weight.MinMax, color.use = celltype_colors, show.legend = F)
  pdf(paste0("20250105_",names(object.list.sc)[i],"_scatter_sender_receiver_SC.pdf"), height =3.75, width = 3.75)
  print(p1)
  dev.off()
}

##Circos Plots
for (i in 1:length(object.list.sc)) {
  strwidth <- function(x) {0.5}
  pathways <- object.list.sc[[i]]@netP$pathways
  p1 <- netVisual_chord_cell(object.list.sc[[i]], signaling = pathways, thresh = 0.0001, remove.isolate = T, title.name = NULL, show.legend = F, color.use = celltype_colors, lab.cex = 1.2)
  pdf(paste("Neighborhood.chordmap_SC_",names(object.list.sc)[i],".pdf"), height = 6, width = 6)
  print(p1)
  dev.off()  
}


##list FIB and immune populations in neighborhoods
###########################################
object.list.merfish <- readRDS("20240105_combined_NN_IAN_MERFISH_trimmed_object.list.rds")

NN1.celltypes <- unique(object.list.merfish$NN1@idents)
NN1.celltypes
NN1.fibs <- c("FIB_Reticular","FIB_Papillary")
NN1.immune <- c("DC","TC_Treg","MAC","LC","TC_Th17","TC_Cd4","TC_gdt","ILC2","Basophil","TC_Cd8")
NN2.celltypes <- unique(object.list.merfish$NN2@idents)
NN2.celltypes
NN2.fibs <- c("FIB_Reticular","FIB_ProInf")
NN2.immune <- c("TC_Cd4")
NN3.celltypes <- unique(object.list.merfish$NN3@idents)
NN3.celltypes
NN3.fibs <- c("FIB_Subcutis","FIB_Reticular","FIB_Papillary")
NN3.immune <- c("MAC","Mast","DC","ILC2","TC_Cd4","Basophil")
NN4.celltypes <- unique(object.list.merfish$NN4@idents)
NN4.celltypes
NN4.fibs <- c("FIB_Subfascia","FIB_Subcutis")
NN4.immune <- c("DC","MAC","Neutrophil","NK","TC_Cd4")
NN5.celltypes <- unique(object.list.merfish$NN5@idents)
NN5.celltypes
NN5.fibs <- c("FIB_ProInf","FIB_Subcutis","FIB_Reticular")
NN5.immune <- c("MAC","Neutrophil","TC_gdt","TC_Th17","NK","Basophil")

IAN1.celltypes <- unique(object.list.merfish$IAN1@idents)
IAN1.celltypes
IAN1.fibs <- c()
IAN1.immune <- c("LC")
IAN2.celltypes <- unique(object.list.merfish$IAN2@idents)
IAN2.celltypes
IAN2.fibs <- c("FIB_Reticular","FIB_Papillary","FIB_ProInf")
IAN2.immune <- c("MAC","DC","Basophil","LC","Mast","TC_Cyc","ILC2","TC_gdt","TC_Cd4")
IAN3.celltypes <- unique(object.list.merfish$IAN3@idents)
IAN3.celltypes
IAN3.fibs <- c("FIB_Subfascia","FIB_ProInf","FIB_Subcutis","FIB_Reticular","FIB_Papillary")
IAN3.immune <- c("MAC","TC_Th17","DC","Basophil","TC_Cd4","Mast","TC_Treg","Neutrophil","NK","ILC2","TC_gdt","TC_Cd8")
IAN4.celltypes <- unique(object.list.merfish$IAN4@idents)
IAN4.celltypes
IAN4.fibs <- c("FIB_Reticular","FIB_ProInf")
IAN4.immune <- c("Neutrophil","MAC","DC","Mast","Basophil","NK","TC_Cd4","TC_gdt")
IAN5.celltypes <- unique(object.list.merfish$IAN5@idents)
IAN5.celltypes
IAN5.fibs <- c()
IAN5.immune <- c("Neutrophil","TC_gdt")
IAN6.celltypes <- unique(object.list.merfish$IAN6@idents)
IAN6.celltypes
IAN6.fibs <- c("FIB_Subcutis","FIB_ProInf")
IAN6.immune <- c()


###########################################

##SC FIB and immune filtered scatter and circos plots
###########################################
object.list.sc <- readRDS("20240105_combined_NN_IAN_SC_trimmed_object.list.rds")
object.list.sc.filtered <- list()

object.list.sc.filtered$NN1.WT <- subsetCellChat(object.list.sc$NN1.WT, idents.use = c(NN1.immune, NN1.fibs))
object.list.sc.filtered$NN2.WT <- subsetCellChat(object.list.sc$NN1.WT, idents.use = c(NN2.immune, NN2.fibs))
object.list.sc.filtered$NN3.WT <- subsetCellChat(object.list.sc$NN1.WT, idents.use = c(NN3.immune, NN3.fibs))
object.list.sc.filtered$NN4.WT <- subsetCellChat(object.list.sc$NN1.WT, idents.use = c(NN4.immune, NN4.fibs))
object.list.sc.filtered$NN5.WT <- subsetCellChat(object.list.sc$NN1.WT, idents.use = c(NN5.immune, NN5.fibs))

object.list.sc.filtered$IAN1.MC903 <- subsetCellChat(object.list.sc$IAN1.MC903, idents.use = c(IAN1.immune, IAN1.fibs))
object.list.sc.filtered$IAN1.OXA <- subsetCellChat(object.list.sc$IAN1.OXA, idents.use = c(IAN1.immune, IAN1.fibs))
object.list.sc.filtered$IAN2.MC903 <- subsetCellChat(object.list.sc$IAN2.MC903, idents.use = c(IAN2.immune, IAN2.fibs))
object.list.sc.filtered$IAN2.OXA <- subsetCellChat(object.list.sc$IAN2.OXA, idents.use = c(IAN2.immune, IAN2.fibs))
object.list.sc.filtered$IAN3.MC903 <- subsetCellChat(object.list.sc$IAN3.MC903, idents.use = c(IAN3.immune, IAN3.fibs))
object.list.sc.filtered$IAN3.OXA <- subsetCellChat(object.list.sc$IAN3.OXA, idents.use = c(IAN3.immune, IAN3.fibs))
object.list.sc.filtered$IAN4.MC903 <- subsetCellChat(object.list.sc$IAN4.MC903, idents.use = c(IAN4.immune, IAN4.fibs))
object.list.sc.filtered$IAN4.OXA <- subsetCellChat(object.list.sc$IAN4.OXA, idents.use = c(IAN4.immune, IAN4.fibs))
object.list.sc.filtered$IAN5.MC903 <- subsetCellChat(object.list.sc$IAN5.MC903, idents.use = c(IAN5.immune, IAN5.fibs))
object.list.sc.filtered$IAN5.OXA <- subsetCellChat(object.list.sc$IAN5.OXA, idents.use = c(IAN5.immune, IAN5.fibs))
object.list.sc.filtered$IAN6.MC903 <- subsetCellChat(object.list.sc$IAN6.MC903, idents.use = c(IAN6.immune, IAN6.fibs))
object.list.sc.filtered$IAN6.OXA <- subsetCellChat(object.list.sc$IAN6.OXA, idents.use = c(IAN6.immune, IAN6.fibs))


saveRDS(object.list.sc.filtered, "20240105_combined_NN_IAN_sc_FIB_Immune_only_object.list.rds")


##scatter plot
for (i in 1:length(object.list.sc.filtered)) {
  num.link <- rowSums(object.list.sc.filtered[[i]]@net$count) + colSums(object.list.sc.filtered[[i]]@net$count)-diag(object.list.sc.filtered[[i]]@net$count)
  weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
  gg <- list()
  p1 <-  netAnalysis_signalingRole_scatter(object.list.sc.filtered[[i]], title = names(object.list.sc.filtered)[i], weight.MinMax = weight.MinMax, color.use = celltype_colors, show.legend = F)
  pdf(paste0("20250105_filtered",names(object.list.sc.filtered)[i],"_scatter_sender_receiver_SC.pdf"), height =3.75, width = 3.75)
  print(p1)
  dev.off()
}

##Circos Plots
for (i in 1:length(object.list.sc.filtered)) {
  strwidth <- function(x) {0.5}
  pathways <- object.list.sc.filtered[[i]]@netP$pathways
  p1 <- netVisual_chord_cell(object.list.sc.filtered[[i]], signaling = pathways, thresh = 0.0001, remove.isolate = T, title.name = NULL, show.legend = F, color.use = celltype_colors, lab.cex = 1.2)
  pdf(paste("Neighborhood.chordmap_SC_filtered_",names(object.list.sc.filtered)[i],".pdf"), height = 6, width = 6)
  print(p1)
  dev.off()  
}

###########################################

##SC FIB sender circos 
###########################################
##NN1
pathways <- object.list.sc.filtered$NN1.WT@netP$pathways
p1 <- netVisual_chord_cell(object.list.sc.filtered$NN1.WT, signaling = pathways, thresh = 0.0001, remove.isolate = T, title.name = NULL, show.legend = F, color.use = celltype_colors, lab.cex = 1.2, sources.use = NN1.fibs, targets.use = NN1.immune)
pdf(paste("Neighborhood.chordmap_SC_filtered_fibs_outgoing_NN1.WT.pdf"), height = 6, width = 6)
print(p1)
dev.off()  


##NN3
pathways <- object.list.sc.filtered$NN3.WT@netP$pathways
p1 <- netVisual_chord_cell(object.list.sc.filtered$NN3.WT, signaling = pathways, thresh = 0.0001, remove.isolate = T, title.name = NULL, show.legend = F, color.use = celltype_colors, lab.cex = 1.2, sources.use = NN3.fibs, targets.use = NN3.immune)
pdf(paste("Neighborhood.chordmap_SC_filtered_fibs_outgoing_NN3.WT.pdf"), height = 6, width = 6)
print(p1)
dev.off()  


##IAN2
pathways <- object.list.sc.filtered$IAN2.MC903@netP$pathways
p1 <- netVisual_chord_cell(object.list.sc.filtered$IAN2.MC903, signaling = pathways, thresh = 0.0001, remove.isolate = T, title.name = NULL, show.legend = F, color.use = celltype_colors, lab.cex = 1.2, sources.use = IAN2.fibs, targets.use = IAN2.immune)
pdf(paste("Neighborhood.chordmap_SC_filtered_fibs_outgoing_IAN2.MC903.pdf"), height = 6, width = 6)
print(p1)
dev.off()  

pathways <- object.list.sc.filtered$IAN2.OXA@netP$pathways
p1 <- netVisual_chord_cell(object.list.sc.filtered$IAN2.OXA, signaling = pathways, thresh = 0.0001, remove.isolate = T, title.name = NULL, show.legend = F, color.use = celltype_colors, lab.cex = 1.2, sources.use = IAN2.fibs, targets.use = IAN2.immune)
pdf(paste("Neighborhood.chordmap_SC_filtered_fibs_outgoing_IAN2.OXA.pdf"), height = 6, width = 6)
print(p1)
dev.off()  

##IAN3
pathways <- object.list.sc.filtered$IAN3.MC903@netP$pathways
p1 <- netVisual_chord_cell(object.list.sc.filtered$IAN3.MC903, signaling = pathways, thresh = 0.0001, remove.isolate = T, title.name = NULL, show.legend = F, color.use = celltype_colors, lab.cex = 1.2, sources.use = IAN3.fibs, targets.use = IAN3.immune)
pdf(paste("Neighborhood.chordmap_SC_filtered_fibs_outgoing_IAN3.MC903.pdf"), height = 6, width = 6)
print(p1)
dev.off()  

pathways <- object.list.sc.filtered$IAN3.OXA@netP$pathways
p1 <- netVisual_chord_cell(object.list.sc.filtered$IAN3.OXA, signaling = pathways, thresh = 0.0001, remove.isolate = T, title.name = NULL, show.legend = F, color.use = celltype_colors, lab.cex = 1.2, sources.use = IAN3.fibs, targets.use = IAN3.immune)
pdf(paste("Neighborhood.chordmap_SC_filtered_fibs_outgoing_IAN3.OXA.pdf"), height = 6, width = 6)
print(p1)
dev.off()  

##IAN4
pathways <- object.list.sc.filtered$IAN4.MC903@netP$pathways
p1 <- netVisual_chord_cell(object.list.sc.filtered$IAN4.MC903, signaling = pathways, thresh = 0.0001, remove.isolate = T, title.name = NULL, show.legend = F, color.use = celltype_colors, lab.cex = 1.2, sources.use = IAN4.fibs, targets.use = IAN4.immune)
pdf(paste("Neighborhood.chordmap_SC_filtered_fibs_outgoing_IAN4.MC903.pdf"), height = 6, width = 6)
print(p1)
dev.off()  

pathways <- object.list.sc.filtered$IAN4.OXA@netP$pathways
p1 <- netVisual_chord_cell(object.list.sc.filtered$IAN4.OXA, signaling = pathways, thresh = 0.0001, remove.isolate = T, title.name = NULL, show.legend = F, color.use = celltype_colors, lab.cex = 1.2, sources.use = IAN4.fibs, targets.use = IAN4.immune)
pdf(paste("Neighborhood.chordmap_SC_filtered_fibs_outgoing_IAN4.OXA.pdf"), height = 6, width = 6)
print(p1)
dev.off()  

###########################################

##SC FIB reciever circos
###########################################
##NN1
pathways <- object.list.sc.filtered$NN1.WT@netP$pathways
p1 <- netVisual_chord_cell(object.list.sc.filtered$NN1.WT, signaling = pathways, thresh = 0.0001, remove.isolate = T, title.name = NULL, show.legend = F, color.use = celltype_colors, lab.cex = 1.2, sources.use = NN1.immune, targets.use = NN1.fibs)
pdf(paste0("Neighborhood.chordmap_SC_filtered_fibs_recieving_NN1.WT.pdf"), height = 6, width = 6)
print(p1)
dev.off()  


##NN3
pathways <- object.list.sc.filtered$NN3.WT@netP$pathways
p1 <- netVisual_chord_cell(object.list.sc.filtered$NN3.WT, signaling = pathways, thresh = 0.0001, remove.isolate = T, title.name = NULL, show.legend = F, color.use = celltype_colors, lab.cex = 1.2, sources.use = NN3.immune, targets.use = NN3.fibs)
pdf(paste("Neighborhood.chordmap_SC_filtered_fibs_recieving_NN3.WT.pdf"), height = 6, width = 6)
print(p1)
dev.off()  


##IAN2
pathways <- object.list.sc.filtered$IAN2.MC903@netP$pathways
p1 <- netVisual_chord_cell(object.list.sc.filtered$IAN2.MC903, signaling = pathways, thresh = 0.0001, remove.isolate = T, title.name = NULL, show.legend = F, color.use = celltype_colors, lab.cex = 1.2, sources.use = IAN2.immune, targets.use = IAN2.fibs)
pdf(paste("Neighborhood.chordmap_SC_filtered_fibs_recieving_IAN2.MC903.pdf"), height = 6, width = 6)
print(p1)
dev.off()  

pathways <- object.list.sc.filtered$IAN2.OXA@netP$pathways
p1 <- netVisual_chord_cell(object.list.sc.filtered$IAN2.OXA, signaling = pathways, thresh = 0.0001, remove.isolate = T, title.name = NULL, show.legend = F, color.use = celltype_colors, lab.cex = 1.2, sources.use = IAN2.immune, targets.use = IAN2.fibs)
pdf(paste("Neighborhood.chordmap_SC_filtered_fibs_recieving_IAN2.OXA.pdf"), height = 6, width = 6)
print(p1)
dev.off()  

##IAN3
pathways <- object.list.sc.filtered$IAN3.MC903@netP$pathways
p1 <- netVisual_chord_cell(object.list.sc.filtered$IAN3.MC903, signaling = pathways, thresh = 0.0001, remove.isolate = T, title.name = NULL, show.legend = F, color.use = celltype_colors, lab.cex = 1.2, sources.use = IAN3.immune, targets.use = IAN3.fibs)
pdf(paste("Neighborhood.chordmap_SC_filtered_fibs_recieving_IAN3.MC903.pdf"), height = 6, width = 6)
print(p1)
dev.off()  

pathways <- object.list.sc.filtered$IAN3.OXA@netP$pathways
p1 <- netVisual_chord_cell(object.list.sc.filtered$IAN3.OXA, signaling = pathways, thresh = 0.0001, remove.isolate = T, title.name = NULL, show.legend = F, color.use = celltype_colors, lab.cex = 1.2, sources.use = IAN3.immune, targets.use = IAN3.fibs)
pdf(paste("Neighborhood.chordmap_SC_filtered_fibs_recieving_IAN3.OXA.pdf"), height = 6, width = 6)
print(p1)
dev.off()  

##IAN4
pathways <- object.list.sc.filtered$IAN4.MC903@netP$pathways
p1 <- netVisual_chord_cell(object.list.sc.filtered$IAN4.MC903, signaling = pathways, thresh = 0.0001, remove.isolate = T, title.name = NULL, show.legend = F, color.use = celltype_colors, lab.cex = 1.2, sources.use = IAN4.immune, targets.use = IAN4.fibs)
pdf(paste("Neighborhood.chordmap_SC_filtered_fibs_recieving_IAN4.MC903.pdf"), height = 6, width = 6)
print(p1)
dev.off()  

pathways <- object.list.sc.filtered$IAN4.OXA@netP$pathways
p1 <- netVisual_chord_cell(object.list.sc.filtered$IAN4.OXA, signaling = pathways, thresh = 0.0001, remove.isolate = T, title.name = NULL, show.legend = F, color.use = celltype_colors, lab.cex = 1.2, sources.use = IAN4.immune, targets.use = IAN4.fibs)
pdf(paste("Neighborhood.chordmap_SC_filtered_fibs_recieving_IAN4.OXA.pdf"), height = 6, width = 6)
print(p1)
dev.off()  

###########################################

##Basophil PIFB signaling
###########################################
WT.cellchat <- readRDS("20241111_mouse_merfish_merge_WT_OXA_MC903_timecourse/Cellchat/20250102_WT_v3_sc_ji_gs_cellchat_nodoublet.rds")

MC903.cellchat <- readRDS("20241111_mouse_merfish_merge_WT_OXA_MC903_timecourse/Cellchat/20250102_MC903_v3_sc_ji_gs_cellchat_nodoublet.rds")

OXA.cellchat <- readRDS("20241111_mouse_merfish_merge_WT_OXA_MC903_timecourse/Cellchat/20250102_OXA_v3_sc_ji_gs_cellchat_nodoublet.rds")

object.list <- list(WT = WT.cellchat, MC903 = MC903.cellchat, OXA= OXA.cellchat)  


old_celltype_colors <- c( "Basophil" ='#EE1289',
                          "FIB_ProInf"="red3")

for (i in 1:length(object.list)) {
  p1 <- netVisual_chord_gene(object.list[[i]], sources.use = "Basophil", targets.use = "FIB_ProInf", lab.cex = 1.5 ,legend.pos.y = 30, thresh = 0.00000001, color.use = old_celltype_colors, small.gap = 3, big.gap = 15)
  pdf(paste("Basophil_PIFB_",names(object.list)[i],".pdf"), height = 8, width = 8)
  print(p1)
  dev.off()  
}
###########################################

##PIFB Immune signaling
###########################################
for (i in 1:length(object.list)) {
  p1 <- netVisual_chord_gene(object.list[[i]], sources.use = "FIB_ProInf", targets.use = IAN3.immune, lab.cex = 1.5 ,legend.pos.y = 30, thresh = 0.00000001, color.use = celltype_colors, signaling = c("CXCL","CCL"), show.legend = F)
  pdf(paste("PIFB_immune_CXCL_CCL",names(object.list)[i],".pdf"), height = 8, width = 8)
  print(p1)
  dev.off()  
}

p1 <- netVisual_chord_cell(object.list.sc.filtered$IAN3.MC903, signaling = "CCL", thresh = 0.0001, remove.isolate = T, title.name = NULL, show.legend = F, color.use = celltype_colors, lab.cex = 1.2, sources.use = "FIB_ProInf", targets.use = IAN3.immune)
pdf(paste("PIFB_CCL_SC_IAN3.MC903.pdf"), height = 6, width = 6)
print(p1)
dev.off()  

p1 <- netVisual_chord_cell(object.list.sc.filtered$IAN3.MC903, signaling = "CXCL", thresh = 0.0001, remove.isolate = T, title.name = NULL, show.legend = F, color.use = celltype_colors, lab.cex = 1.2, sources.use = "FIB_ProInf", targets.use = IAN3.immune)
pdf(paste("PIFB_CXCL_SC_IAN3.MC903.pdf"), height = 6, width = 6)
print(p1)
dev.off()  

p1 <- netVisual_chord_cell(object.list.sc.filtered$IAN3.OXA, signaling = "CCL", thresh = 0.0001, remove.isolate = T, title.name = NULL, show.legend = F, color.use = celltype_colors, lab.cex = 1.2, sources.use = "FIB_ProInf", targets.use = IAN3.immune)
pdf(paste("PIFB_CCL_SC_IAN3.OXA.pdf"), height = 6, width = 6)
print(p1)
dev.off()  

p1 <- netVisual_chord_cell(object.list.sc.filtered$IAN3.OXA, signaling = "CXCL", thresh = 0.0001, remove.isolate = T, title.name = NULL, show.legend = F, color.use = celltype_colors, lab.cex = 1.2, sources.use = "FIB_ProInf", targets.use = IAN3.immune)
pdf(paste("PIFB_CXCL_SC_IAN3.OXA.pdf"), height = 6, width = 6)
print(p1)
dev.off()  

pathways <- object.list.sc.filtered$IAN3.OXA@netP$pathways
p1 <- netVisual_chord_cell(object.list.sc.filtered$IAN3.OXA, signaling = pathways, thresh = 0.0001, remove.isolate = T, title.name = NULL, show.legend = F, color.use = celltype_colors, lab.cex = 1.2, sources.use = IAN3.fibs, targets.use = IAN3.immune)
pdf(paste("Neighborhood.chordmap_SC_filtered_fibs_outgoing_IAN3.OXA.pdf"), height = 6, width = 6)
print(p1)
dev.off()  

###########################################

##Basophil to PIFB DotPlot
###########################################
WT.cellchat <- readRDS("20241111_mouse_merfish_merge_WT_OXA_MC903_timecourse/Cellchat/20250102_WT_v3_sc_ji_gs_cellchat_nodoublet.rds")

MC903.cellchat <- readRDS("20241111_mouse_merfish_merge_WT_OXA_MC903_timecourse/Cellchat/20250102_MC903_v3_sc_ji_gs_cellchat_nodoublet.rds")

OXA.cellchat <- readRDS("20241111_mouse_merfish_merge_WT_OXA_MC903_timecourse/Cellchat/20250102_OXA_v3_sc_ji_gs_cellchat_nodoublet.rds")

object.list <- list(WT = WT.cellchat, MC903 = MC903.cellchat, OXA= OXA.cellchat)  

cellchat <- mergeCellChat(object.list, add.names = names(object.list), cell.prefix = F)

cytokines <- subsetCommunication(cellchat, sources.use = c("Basophil","DC","Mast","MAC"), targets.use = c("FIB_ProInf"))

pairLR.use.MC903 = cytokines$MC903[, "interaction_name", drop = F]
pairLR.use.OXA = cytokines$OXA[, "interaction_name", drop = F]
pairLR.use <- unique(rbind(pairLR.use.MC903, pairLR.use.OXA))
pairLR.use.subset <-  subset(pairLR.use, interaction_name %in% c("TGFB1_TGFBR1_TGFBR2","TGFB1_ACVR1B_TGFBR2","TGFB1_ACVR1_TGFBR1",
                                                                 "IL4_IL4R","IL4_IL4R_IL13RA1","IL4_IL4R_IL13RA2","IL4_IL4R_IL2RG","IL13_IL4R_IL13RA1","IL13_IL4R_IL13RA2","IL13_IL13RA1","IL13_IL13RA2",
                                                                 "OSM_LIFR_IL6ST","OSM_OSMR_IL6ST","IL6_IL6R_IL6ST", "IL1A_IL1R1_IL1RAP","IL1B_IL1R1_IL1RAP","TNF_TNFRSF1A","TNF_TNFRSF1B","IFNG_IFNGR1_IFNGR2"))
pairLR.use.subset$interaction_name <- sort(pairLR.use.subset$interaction_name)
p1 <- netVisual_bubble(cellchat, sources.use = rev(c("Basophil","DC","Mast","MAC")), targets.use = "FIB_ProInf", remove.isolate = F, pairLR.use = pairLR.use.subset, comparison = c(2,3), sort.by.source = F)+ coord_flip() + rotate_x_text(angle = 45, vjust = 1)
pdf(paste("Basophil_PIFBSignaling.pdf"), height = 2.5, width = 6.25)
print(p1)
dev.off() 
###########################################

##PIFB to immune DotPlot
###########################################
WT.cellchat <- readRDS("20241111_mouse_merfish_merge_WT_OXA_MC903_timecourse/Cellchat/20250102_WT_v3_sc_ji_gs_cellchat_nodoublet.rds")

MC903.cellchat <- readRDS("20241111_mouse_merfish_merge_WT_OXA_MC903_timecourse/Cellchat/20250102_MC903_v3_sc_ji_gs_cellchat_nodoublet.rds")

OXA.cellchat <- readRDS("20241111_mouse_merfish_merge_WT_OXA_MC903_timecourse/Cellchat/20250102_OXA_v3_sc_ji_gs_cellchat_nodoublet.rds")

object.list <- list(WT = WT.cellchat, MC903 = MC903.cellchat, OXA= OXA.cellchat)  

cellchat <- mergeCellChat(object.list, add.names = names(object.list), cell.prefix = F)

chemokines <- subsetCommunication(cellchat, signaling = c("CXCL","CCL"), sources.use = c("FIB_ProInf"), targets.use = c("Basophil","DC","ILC2","MAC"))

pairLR.use.MC903 = chemokines$MC903[, "interaction_name", drop = F]
pairLR.use.OXA = chemokines$OXA[, "interaction_name", drop = F]
pairLR.use <- unique(rbind(pairLR.use.MC903, pairLR.use.OXA))
pairLR.use.subset <-  subset(pairLR.use, !interaction_name %in% c("CCL6_CCR2", "CCL6_CCR3",
                                                                  "CCL9_CCR1",
                                                                  "CCL3_CCR1","CCL3_CCR5",
                                                                  "CCL4_CCR5",
                                                                  "CCL27A_CCR3","CCL28_CCR3","CCL27A_CCR10","CCL27A_CCR2",
                                                                  "CCL24_CCR3",
                                                                  "CCL25_CCR9",
                                                                  "CCL21A_CCR7",
                                                                  "CCL20_CCR6",
                                                                  "CXCL4_CXCR3",
                                                                  "CXCL13_CXCR5","CXCL13_CXCR3",
                                                                  "CXCL10_CXCR3",
                                                                  "CXCL9_CXCR3",
                                                                  "CXCL16_CXCR6"))
pairLR.use.subset$interaction_name <- sort(as.character(pairLR.use.subset$interaction_name))

p1 <- netVisual_bubble(cellchat, sources.use = c("FIB_ProInf"), targets.use = c("TC_Cd4","MAC", "DC","Basophil"), remove.isolate = F, comparison = c(2,3), pairLR.use = pairLR.use.subset, sort.by.source.priority = F, sort.by.target = F) + coord_flip() + rotate_x_text(angle = 45, vjust = 1)
pdf(paste("PIFB_Immune_Signaling_rotated.pdf"), height = 2, width = 5.25)
print(p1)
dev.off() 
###########################################


