library(Seurat)
library(FastIntegration)
library(harmony)
library(pbmcapply)
library(Matrix)
library(SeuratData)

InstallData("panc8")
data("panc8")
rna.list = SplitObject(panc8, split.by = "tech")
rna.list = rna.list[c("celseq", "celseq2", "fluidigmc1", "smartseq2")]

FastIntegration::BuildIntegrationFile(rna.list = rna.list, 
                                      tmp.dir = "./", nCores = 5)
FastIntegration::FastFindAnchors(tmp.dir = "./", nCores = 5)

features = SelectIntegrationFeatures(rna.list)
rna.bind = rna.list[[1]]@assays$RNA@data[features,]
for (i in 2:4) {
  rna.bind = cbind(rna.bind, rna.list[[i]]@assays$RNA@data[features,])
}
rna.bind = ScaleData(rna.bind)
pca = RunPCA(rna.bind)
meta = data.frame(sample = c(rna.list[[1]]$tech, rna.list[[2]]$tech, rna.list[[3]]$tech, rna.list[[4]]$tech))
harmony = HarmonyMatrix(pca@cell.embeddings, meta, "sample", do_pca=FALSE)

rna.integrated = FastIntegration::FastIntegration(tmp.dir = "./", 
                                                  input.pca = harmony[,1:30], slot = "data",
                                                  features.to.integrate = rownames(rna.list[[1]]))

rna.integrated = CreateSeuratObject(rna.integrated)
rna.integrated = ScaleData(rna.integrated, features = features)
rna.integrated = RunPCA(rna.integrated, features = features)
rna.integrated = RunUMAP(rna.integrated, dims = 1:30)

meta = rbind(rna.list[[1]]@meta.data, rna.list[[2]]@meta.data, rna.list[[3]]@meta.data, rna.list[[4]]@meta.data)
rna.integrated = AddMetaData(rna.integrated, meta)
rna.integrated@active.ident = as.factor(rna.integrated$celltype)
names(rna.integrated@active.ident) = Cells(rna.integrated)

DimPlot(rna.integrated, group.by = "tech", label = T)

markers.integrated = pbmclapply(
  unique(rna.integrated$celltype), function(i) {
    marker = FindMarkers(rna.integrated, ident.1 = i)
    marker$gene = rownames(marker)
    marker$cell = i
    return(marker)
  }, mc.cores = 13
)
markers.integrated = do.call(rbind, markers.integrated)
markers.integrated = markers.integrated[which(abs(markers.integrated$avg_log2FC) > 0.5),]

rna.data = merge(rna.list[[1]], rna.list[2:4], merge.data = TRUE)
rna.data@active.ident = as.factor(rna.data$celltype)
names(rna.data@active.ident) = Cells(rna.data)

markers.conserve = pbmclapply(
  unique(rna.integrated$celltype), function(i) {
    tryCatch(
      {
        marker = FindConservedMarkers(rna.data, ident.1 = i, grouping.var = "tech")
        marker$gene = rownames(marker)
        marker$cell = i
        return(marker)
      },
      error = function(e) {
        return(NULL)
      }
    )
  }, mc.cores = 13
)

markers.conserve = FindConservedMarkers(rna.data, ident.1 = "gamma", grouping.var = "tech")
markers.conserve = do.call(rbind, markers.conserve)

setdiff(unique(rna.data$celltype), unique(markers.conserve$cell))
# Four cell types (schwann, mast, macrophage, epsilon and quiescent_stellate) has too few cells to detect DEGs.
# Error Info: Cell group 1 has fewer than 3 cells 

unique(markers.conserve$cell)[1]
marker = markers.integrated[which(markers.integrated$cell == "gamma"),]
marker = marker[order(-(marker$avg_log2FC)),]
setdiff(marker$gene[1:10], markers.conserve$gene[which(markers.conserve$cell == "gamma")]) # GPC5-AS1


unique(markers.conserve$cell)[2]
marker = markers.integrated[which(markers.integrated$cell == "alpha"),]
marker = marker[order(-(marker$avg_log2FC)),]
setdiff(marker$gene[1:10], markers.conserve$gene[which(markers.conserve$cell == "alpha")]) # same


unique(markers.conserve$cell)[3]
marker = markers.integrated[which(markers.integrated$cell == "delta"),]
marker = marker[order(-(marker$avg_log2FC)),]
setdiff(marker$gene[1:10], markers.conserve$gene[which(markers.conserve$cell == "delta")]) # same


unique(markers.conserve$cell)[4]
marker = markers.integrated[which(markers.integrated$cell == "beta"),]
marker = marker[order(-(marker$avg_log2FC)),]
setdiff(marker$gene[1:10], markers.conserve$gene[which(markers.conserve$cell == "beta")]) # FAM159B


unique(markers.conserve$cell)[5]
marker = markers.integrated[which(markers.integrated$cell == "acinar"),]
marker = marker[order(-(marker$avg_log2FC)),]
setdiff(marker$gene[1:10], markers.conserve$gene[which(markers.conserve$cell == "acinar")]) # same

unique(markers.conserve$cell)[6]
marker = markers.integrated[which(markers.integrated$cell == "ductal"),]
marker = marker[order(-(marker$avg_log2FC)),]
setdiff(marker$gene[1:10], markers.conserve$gene[which(markers.conserve$cell == "ductal")]) # same

unique(markers.conserve$cell)[7]
marker = markers.integrated[which(markers.integrated$cell == "activated_stellate"),]
marker = marker[order(-(marker$avg_log2FC)),]
setdiff(marker$gene[1:10], markers.conserve$gene[which(markers.conserve$cell == "activated_stellate")]) # same

unique(markers.conserve$cell)[8]
marker = markers.integrated[which(markers.integrated$cell == "endothelial"),]
marker = marker[order(-(marker$avg_log2FC)),]
setdiff(marker$gene[1:10], markers.conserve$gene[which(markers.conserve$cell == "endothelial")]) # PECAM1
