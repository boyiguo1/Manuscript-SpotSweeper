# Load library ----


# Load data ----
data(DLPFC_artifact)
spe <- DLPFC_artifact

# Caclualte stat across orders----

## Define orders ----


## Iterate over the function ----


spe <- findArtifacts(spe,
    mito_percent = "expr_chrM_ratio",
    mito_sum = "expr_chrM",
    n_order = 2, # 5 recommended, using 2 for time
    shape = "hexagonal",
    name = "artifact"
)

# Make Plots -----

# Session Info
sessioninfo::session_info()