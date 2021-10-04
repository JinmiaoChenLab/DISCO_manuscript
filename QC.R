## This script is QC and cell selection for global atlas, QC of other subatlas is
## similar to it

######################## QC ########################
pbmclapply(
  1:nrow(all.sample), function(i) {
    p = all.sample[i,1]
    s = all.sample[i,2]
    
    rna = readRDS(paste0("/data_dir/", p, "/", s, ".rds"))
    rna[["percent.mt"]] = PercentageFeatureSet(rna, pattern = "^MT-")
    
    mt.cut = qnorm(0.95, mean = mean(rna$percent.mt), sd = sd(rna$percent.mt))
    nFeature_RNA.cut = qnorm(0.95, mean = mean(rna$nFeature_RNA), sd = sd(rna$nFeature_RNA))
    nFeature_RNA.cutmin = qnorm(0.05, mean = mean(rna$nFeature_RNA), sd = sd(rna$nFeature_RNA))
    
    if(length(which(rna$nFeature_RNA > max(500, nFeature_RNA.cutmin) & rna$nFeature_RNA < nFeature_RNA.cut & rna$percent.mt < min(mt.cut, 15))) > 500) {
      rna = subset(rna, subset = nFeature_RNA > max(500, nFeature_RNA.cutmin) & nFeature_RNA < nFeature_RNA.cut & percent.mt < min(mt.cut, 15))
      rna = NormalizeData(rna)
      d = computeDoubletDensity(rna@assays$RNA@data)
      db.cut = qnorm(0.95, mean = mean(d), sd = sd(d))
      rna = subset(rna,  cells = Cells(rna)[which(d < db.cut)])
      rna@assays$RNA@data = matrix(0)
      dir.create(paste0("/clean_rds/", p))
      saveRDS(rna, file = paste0("/clean_rds/", p, "/", s, ".rds"), compress = F)
    }
  }, mc.cores = 80
)

data.list = pbmclapply(
  1:nrow(res),
  function(i) {
    f = paste0("/clean_rds/", res$project_id[i], "/", res$sample_id[i], ".rds") #read data
    if (file.exists(f)) {
      d = readRDS(f)
      d@assays$RNA@data = d@assays$RNA@counts
      d$sample = res$sample_id[i]
      if (ncol(d) < 700) {
        return(NULL)
      }
      d = NormalizeData(d)
      d = FindVariableFeatures(d)
      if (ncol(d) > 1500) {
        d = ScaleData(d, features = VariableFeatures(d))
        d = RunPCA(d, features = VariableFeatures(d))
        d = FindNeighbors(d, nn.method = "rann")
        d = FindClusters(d, resolution = 0.5)
        d = subset(d, downsample = ceiling(1500/max(as.numeric(d@active.ident))))
        d@assays$RNA@scale.data = matrix(0)
      }
      return(d)
    }else {
      return(NULL)
    }
  }, mc.cores = 100
)
