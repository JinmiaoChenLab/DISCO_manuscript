args.input=commandArgs(T)
library(Seurat)
library(data.table)
library(Matrix)
library(uwot)
library(RMySQL)
library(RJSONIO)

job.id = args.input[1]

con <- dbConnect(
  MySQL(), 
  host="", 
  dbname="", 
  user="", 
  password=""
)

sql = paste0('update mapper_job set job_status="', 'load reference...', '" where job_id="', job.id, '"')
dbSendQuery(con, sql)

tryCatch(
  {
    ref = dbSendQuery(con, paste0('SELECT ref from mapper_job where job_id="', job.id, '"'))
    ref = dbFetch(ref, n = 1)
    ref = ref[1,1]
    
    rna.ref = readRDS(paste0("/home/eu0316/hscrm/mapping/ref/", ref, "_ref.rds"))
    umap.mod = uwot::load_uwot(paste0("/home/eu0316/hscrm/mapping/ref/", ref, "_umap.rds"))
    rna.ref$cluster = as.character(rna.ref$cluster)
    rna.ref@assays$RNA@counts = rna.ref@assays$RNA@data
    rna.ref@assays$RNA@counts@x = exp(rna.ref@assays$RNA@data@x) - 1
  }, 
  error=function(e) {
    sql = paste0('update mapper_job set job_status="', 'fail to load ref', '" where job_id="',
                 job.id, '"')
    dbSendQuery(con, sql)
    quit(save = "no")
  }
)


sql = paste0('update mapper_job set job_status="', 'check input...', '" where job_id="',
             job.id, '"')
dbSendQuery(con, sql)

tryCatch(
  {
    rna.query = readRDS(paste0("~/hscrm/mapping/jobs/", job.id,"/input.rds"))
    rna.query@assays$RNA@counts = rna.query@assays$RNA@data
    rna.query = UpdateSeuratObject(rna.query)
    if (length(setdiff(rownames(rna.ref), rownames(rna.query))) > 0) {
      sql = paste0('update mapper_job set warning_msg = "Some feature genes are lacked in input data. This will inflence the accuray of UMAP" where job_id="', job.id, '"')
      dbSendQuery(con, sql)
    }
    overlap.gene = intersect(rownames(rna.query), rownames(rna.ref@reductions$pca@feature.loadings))
    rna.query = ScaleData(rna.query, features = overlap.gene)
    pca = t(rna.query@assays$RNA@scale.data[overlap.gene,]) %*% rna.ref@reductions$pca@feature.loadings[overlap.gene,]
    rna.query@reductions$pca = CreateDimReducObject(pca)
  }, 
  error=function(e) {
    sql = paste0('update mapper_job set job_status="', 'fail to load ref', '",error_msg="',
                 e,'" where job_id="',job.id, '"')
    dbSendQuery(con, sql)
    quit(save = "no")
  }
)

sql = paste0('update mapper_job set job_status="', 'find anchors...', '" where job_id="',
             job.id, '"')
dbSendQuery(con, sql)

tryCatch(
  {
    anchors = FindTransferAnchors(reference = rna.ref, query = rna.query, reduction = "cca",
                                  features = rownames(rna.ref), 
                                  dims = 1:50, k.filter = 100)
  }, 
  error=function(e) {
    sql = paste0('update mapper_job set job_status="', 'unsupported data', '",error_msg="',
                 e,'" where job_id="',job.id, '"')
    dbSendQuery(con, sql)
    quit(save = "no")
  }
)

sql = paste0('update mapper_job set job_status="', 'Map data...', '" where job_id="',
             job.id, '"')
dbSendQuery(con, sql)

tryCatch(
  {
    
    anchors = data.frame(anchors@anchors)
    nn.cells2 = Cells(rna.query)
    anchors.cells2 = unique(x = nn.cells2[anchors[, "cell2"]])
    data.use = rna.query@reductions$pca@cell.embeddings[nn.cells2, ]
    knn_2_2 <- Seurat:::NNHelper(
      data = data.use[anchors.cells2, ],
      query = data.use,
      k = 50,
      method = "annoy",
      n.trees = 50,
      eps = 0
    )
    
    distances <- Distances(object = knn_2_2)
    distances <- 1 - (distances / distances[, ncol(x = distances)])
    cell.index <- Indices(object = knn_2_2)
    
    weights <- Seurat:::FindWeightsC(
      cells2 = 0:(length(x = nn.cells2) - 1),
      distances = as.matrix(x = distances),
      anchor_cells2 = anchors.cells2,
      integration_matrix_rownames = nn.cells2[anchors$cell2],
      cell_index = cell.index,
      anchor_score = anchors[, "score"],
      min_dist = 0,
      sd = 1,
      display_progress = TRUE
    )
    
    integration.matrix <- rna.query@reductions$pca@cell.embeddings[anchors$cell2, ] - 
      rna.ref@reductions$pca@cell.embeddings[anchors$cell1, ]
    pca = rna.query@reductions$pca@cell.embeddings - t(weights) %*% integration.matrix
    
    
    #### cell type
    knn_2_2 <- Seurat:::NNHelper(
      data = rna.ref@reductions$pca@cell.embeddings,
      query = pca,
      k = 10,
      method = "annoy",
      n.trees = 50,
      eps = 0
    )
    
    cell.index <- Indices(object = knn_2_2)
    
    predicted.cluster = as.character(apply(cell.index, 1, function(i) {
      return(names(which.max(table(rna.ref$cluster[i]))))
    }))
    #######
    
    
    umap.query = umap_transform(as.matrix(pca), umap.mod, verbose = TRUE)
    
    res = cbind(as.matrix(umap.query)[,c(1,2)], predicted.cluster)
    rownames(res) = Cells(rna.query)
    colnames(res) = c("disco_umap1", "disco_umap2", "disco_cluster")
    
    res = data.frame(res)
    res$disco_umap1 = round(as.numeric(res$disco_umap1), 3)
    res$disco_umap2 = round(as.numeric(res$disco_umap2), 3)
    
    cell.type = fromJSON(paste0("https://www.immunesinglecell.org/api/subatlas/getCellType?atlas=", ref)) 
    cell.type = lapply(cell.type, data.frame)
    cell.type = do.call(rbind, cell.type)
    res$disco_celltype = cell.type$standardizedName[as.numeric(res$disco_cluster)+1]
    res$disco_anchor = "False"
    res[anchors.cells2, "disco_anchor"] = "True"
    
    res = data.frame(cbind(res), rna.query@meta.data)
    res$cell = rownames(res)
    
    write.table(res, paste0("~/hscrm/mapping/jobs/", job.id,"/output.txt"), 
                sep = "\t", quote = F, row.names = F, col.names = T)
    
    sql = paste0('update mapper_job set mapped=', length(anchors.cells2), ',total=', 
                 length(nn.cells2),' where job_id="',
                 job.id, '"')
    dbSendQuery(con, sql)
  }, 
  error=function(e) {
    sql = paste0('update mapper_job set job_status="fail to map",error_msg="',
                 e,'" where job_id="',job.id, '"')
    dbSendQuery(con, sql)
    quit(save = "no")
  }
)

sql = paste0('update mapper_job set job_status="', 'finished', '" where job_id="',
             job.id, '"')
dbSendQuery(con, sql)




