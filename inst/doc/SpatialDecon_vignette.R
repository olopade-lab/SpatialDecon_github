## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(SpatialDecon)

## ----loaddata-----------------------------------------------------------------
data("mini_geomx_dataset")
norm = mini_geomx_dataset$normalized
raw = mini_geomx_dataset$raw
annot = mini_geomx_dataset$annot
dim(raw)
head(annot)
raw[seq_len(5), seq_len(5)]

# better segment names:
colnames(norm) = colnames(raw) = rownames(annot) = paste0(annot$ROI, annot$AOI.name)


## ----estimateBG---------------------------------------------------------------
# use the NegProbe to estimate per-observation background
per.observation.mean.neg = norm["NegProbe", ]
# and define a background matrix in which each column (observation) is the appropriate 
#  value of per-observation background:
bg = sweep(norm * 0, 2, per.observation.mean.neg, "+")
dim(bg)


## ----derivebg-----------------------------------------------------------------
bg2 = derive_GeoMx_background(norm = norm,
                             probepool = rep(1, nrow(norm)),
                             negnames = "NegProbe")

## ----showsafetme, fig.height=5, fig.width=10, fig.cap = "The safeTME cell profile matrix"----
signif(safeTME[seq_len(3), seq_len(3)], 2)

heatmap(sweep(safeTME, 1, apply(safeTME, 1, max), "/"),
        labRow = NA, margins = c(10, 5))


## ----downloadmatrix, fig.height=7, fig.width=10, fig.cap = "The Mouse Brain profile matrix"----
mousebrain = download_profile_matrix(matrixname = "Mouse_Brain")
dim(mousebrain)

heatmap(sweep(mousebrain, 1, apply(mousebrain, 1, max), "/"),
        labRow = NA, margins = c(12, 5), cexCol = 0.7)


## ----runiss-------------------------------------------------------------------
res = spatialdecon(norm = norm,
                   bg = bg,
                   X = safeTME,
                   align_genes = TRUE)
str(res)

## ----plotissres, fig.height = 5, fig.width = 8, fig.cap = "Cell abundance estimates"----
heatmap(res$beta, cexCol = 0.5, cexRow = 0.7, margins = c(10,7))

## ----showmatches--------------------------------------------------------------
str(safeTME.matches)

## ----runisstils---------------------------------------------------------------
# vector identifying pure tumor segments:
annot$istumor = (annot$AOI.name == "Tumor")

# run spatialdecon with all the bells and whistles:
restils = spatialdecon(norm = norm,                     # normalized data
                       raw = raw,                       # raw data, used to down-weight low-count observations 
                       bg = bg,                         # expected background counts for every data point in norm
                       X = safeTME,                     # safeTME matrix, used by default
                       cellmerges = safeTME.matches,   # safeTME.matches object, used by default
                       cell_counts = annot$nuclei,      # nuclei counts, used to estimate total cells
                       is_pure_tumor = annot$istumor,   # identities of the Tumor segments/observations
                       n_tumor_clusters = 5)            # how many distinct tumor profiles to append to safeTME

str(restils)

## ----shownewX, fig.height=5, fig.width=8, fig.cap = "safeTME merged with newly-derived tumor profiles"----
heatmap(sweep(restils$X, 1, apply(restils$X, 1, max), "/"),
         labRow = NA, margins = c(10, 5))


## ----compareresults, fig.height=6, fig.width=8, fig.cap = "Cell abundance estimates with and without modelling tumor profiles"----
par(mfrow = c(2, 3))
par(mar = c(5,7,2,1))
for (i in seq_len(6)) {
  cell = rownames(res$beta)[i]
  plot(res$beta[cell, ], restils$beta.granular[cell, ],
       xlab = paste0(cell, " score under basic setting"), 
       ylab = paste0(cell, " score when tumor\ncells are modelled"), 
       pch = 16,
       col = c(rgb(0,0,1,0.5), rgb(1,0,0,0.5))[1 + annot$istumor],
       xlim = range(c(res$beta[cell, ], restils$beta.granular[cell, ])),
       ylim = range(c(res$beta[cell, ], restils$beta.granular[cell, ])))
  abline(0,1)
  if (i == 1) {
    legend("topleft", pch = 16, col = c(rgb(0,0,1,0.5), rgb(1,0,0,0.5)),
           legend = c("microenv.", "tumor"))
  }
}


## ----barplot, fig.width=9, fig.height=6, fig.cap="Barplots of TIL abundance"----
# For reference, show the TILs color data object used by the plotting functions when safeTME has been used:
cellcols

# show just the TME segments, since that's where the immune cells are:
layout(mat = (matrix(c(1, 2), 1)), widths = c(7, 3))
TIL_barplot(restils$cell.counts$cell.counts, draw_legend = TRUE, cex.names = 0.5)
# or the proportions of cells:
TIL_barplot(restils$prop_of_nontumor[, annot$AOI.name == "TME"], draw_legend = TRUE, cex.names = 0.75)


## ----florets, fig.width=8, fig.height=6, fig.cap = "TIL abundance plotted on PC space"----
# PCA of the normalized data:
pc = prcomp(t(log2(pmax(norm, 1))))$x[, c(1, 2)]

# run florets function:
par(mar = c(5,5,1,1))
layout(mat = (matrix(c(1, 2), 1)), widths = c(6, 2))
florets(x = pc[, 1], y = pc[, 2],
        b = restils$beta, cex = 2,
        xlab = "PC1", ylab = "PC2")
par(mar = c(0,0,0,0))
frame()
legend("center", fill = cellcols[rownames(restils$beta)], legend = rownames(restils$beta), cex = 0.7)

## ----collapse, fig.width=5, fig.height=5, fig.cap="Cell abundance estimates with related cell types collapsed"----
matching = list()
matching$myeloid = c( "macrophages", "monocytes", "mDCs")
matching$T.NK = c("CD4.T.cells","CD8.T.cells", "Treg", "NK")
matching$B = c("B")
matching$mast = c("mast")
matching$neutrophils = c("neutrophils")
matching$stroma = c("endothelial.cells", "fibroblasts")


collapsed = collapseCellTypes(fit = restils, 
                              matching = matching)

heatmap(collapsed$beta, cexRow = 0.85, cexCol = 0.75)

## ----appendtumor, fig.width = 10, fig.height= 5, fig.cap = "safeTME merged with newly-derived tumor profiles"----
pure.tumor.ids = annot$AOI.name == "Tumor"
safeTME.with.tumor = mergeTumorIntoX(norm = norm, 
                                     bg = bg, 
                                     pure_tumor_ids = pure.tumor.ids, 
                                     X = safeTME, 
                                     K = 3) 

heatmap(sweep(safeTME.with.tumor, 1, apply(safeTME.with.tumor, 1, max), "/"),
        labRow = NA, margins = c(10, 5))


## ----reverse, fig.height=6, fig.width=6, fig.cap="Residuals from reverseDecon"----
rdecon = reverseDecon(norm = norm,
                      beta = res$beta)
str(rdecon)

# look at the residuals:
heatmap(pmax(pmin(rdecon$resids, 2), -2))

## ----reverse2, fig.height=6, fig.width=6, fig.cap="Genes' dependency on cell mixing"----
# look at the two metrics of goodness-of-fit:
plot(rdecon$cors, rdecon$resid.sd, col = 0)
showgenes = c("CXCL14", "LYZ", "NKG7")
text(rdecon$cors[setdiff(names(rdecon$cors), showgenes)], rdecon$resid.sd[setdiff(names(rdecon$cors), showgenes)], 
     setdiff(names(rdecon$cors), showgenes), cex = 0.5)
text(rdecon$cors[showgenes], rdecon$resid.sd[showgenes], showgenes, cex = 0.75, col = 2)


