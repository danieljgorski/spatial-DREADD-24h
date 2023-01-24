# In this script we will annotate the molecular niches

# Load libraries
library(Seurat)
library(ggplot2)
library(ape)

# Set up output dirs
output_dirs <- c("results",
                 "results/annotation")
for (i in output_dirs) {
  if (!dir.exists(i)) {
    dir.create(i)
    print(paste0("made ", i, " directory"))
  } else {
    print(paste0(i, " directory already exists."))
  }
}

# Load object
load(file = "results/objects/obj_integrated.Rdata")

# Set Idents to the optimum resolution
Idents(obj) <- "integrated_snn_res.0.3"

# Calculate niche markers
niche_markers <- FindAllMarkers(obj,
                                assay = "Spatial",
                                only.pos = T)
write.csv(niche_markers,
          file = "results/annotation/niche_markers.csv",
          row.names = F)
top5 <- niche_markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)

# Dendrogram
obj <- BuildClusterTree(obj, dims = 1:25)
PlotClusterTreeDJG <- function(object, ...) {
  if (is.null(x = Tool(object = object, slot = "BuildClusterTree"))) {
    stop("Phylogenetic tree does not exist, build using BuildClusterTree")
  }
  data.tree <- Tool(object = object, slot = "BuildClusterTree")
  plot.phylo(x = data.tree, font = 2, direction = "rightwards", label.offset = 2, edge.width = 1)
}
pdf(file = "results/annotation/Dendrogram.pdf",
    useDingbats = F)
PlotClusterTreeDJG(obj)
dev.off()

###############################################################################
# Annotation notes
###############################################################################

# 0 = Ischemic zone (IZ): S100a8, S100a9 (neutrophils), Arg1, Ccl6, Ccl9, Ly2,
#   Pf4, classic inflammatory response

# 1 = Border zone (BZ-1): Ankrd1, Nppa, Nppb Similar to Calcagno et al. 2022 BZ-1
# 2 = Remote zone (RZ-1): mt-genes, Lpl, Myh6, surviving myocardium
# 3 = Border zone (BZ-2): Cryab, Des, Xirp2, Flnc, Similar to Calcagno et al. 2022 BZ-2,
#   very thing layer, similar to Alcagno et al. 2022 BZ-2
# 4 = Remote zone (RZ-2): Tcap, Tnni3, Mb, Tnnt2, (sarcomere proteins), surviving myocardium
# 5 = Remote zone RBC (RZ-3): Myh6, Mb, Lpl, Tcap, Tnni surviving myocardium, with
# hemoglobins from erythrocytes (Hbb-bt, Hbb-bs, Hba-a2, Hba-a1, Alas2)
# 6 = Fibrotic zone (FZ): Mfap4, Sfrp2, Cthrc1, Ltbp2, Ccn5, Lox, collagens, strong
# fibrotic signature. Previously we annotated a similar niche as an ischemic zone
# with a firbotic signature, but the dendrogram in this experiment shows a much
# weaker relationship to the ischemic zone (Niche 0).

###############################################################################

# Rename Idents to annotations
obj <- RenameIdents(obj,
                    "0" = "IZ",
                    "1" = "BZ-1",
                    "2" = "RZ-1",
                    "3" = "BZ-2",
                    "4" = "RZ-2",
                    "5" = "RZ-3",
                    "6" = "FZ")

# Store renamed idents as a new meta data column, set as Idents
obj@meta.data$basic_annotation <- Idents(obj)

# Refactor annotation levels
source("scripts/dimplotlevels.R")
obj@meta.data$basic_annotation <- factor(obj@meta.data$basic_annotation,
                                         levels = dimplotlevels)
DimPlot(obj,
        group.by = "basic_annotation",
        label = T,
        repel = T)

# Set Idents as re-factored basic_annotation identities
Idents(obj) <- "basic_annotation"

# Save object with basic annotations
save(obj, file = "results/objects/obj_annotated.Rdata")
