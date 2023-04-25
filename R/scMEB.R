#' scMEB: A fast and clustering-independent method for detecting differentially
#' expressed genes in single-cell RNA-seq data
#' @description Using the Minimum Enclosing Ball (MEB) method to discriminate
#' differential expression genes (DEGs) without prior cell clustering results.
#' @usage scMEB(sce, stable_idx, filtered = FALSE,
#' gamma = seq(1e-04,0.001,1e-05), nu = 0.01, reject_rate = 0.1)
#' @param sce A SingleCellExperiment class scRNA-seq data.
#' @param stable_idx A vector shows the name of stable expressed gene in sce.
#' @param filtered A logical value to show if the data have been filtered.
#' @param gamma A parameter needed for all kernels except linear.
#' @param nu A parameter needed for one-classification.
#' @param reject_rate A value used in controling the scale of ball, default is
#' 0.1.
#' @return list(.) A list of results, "model" represents the model of scMEB,
#' which could be used to discriminate a new gene, "dat_pca" represents the
#' first 50 PCs of each genes, "gamma" represents the selected gamma parameters
#' in model scMEB, "train_error" represents the corresponding train_error when
#' the value of gamma changed, "dist" shows the distance between the points and
#' the sphere in feature space.
#' @examples
#' ## Simulation data for scRNA-seq data are generated from splatter package.
#' library(SingleCellExperiment)
#' data(sim_scRNA_data)
#' data(stable_gene)
#' sim_scRNA <- scMEB(sce=sim_scRNA_data, stable_idx=stable_gene,
#' filtered = FALSE, gamma = seq(1e-04,0.001,1e-05), nu = 0.01,
#' reject_rate = 0.1)
#' @export
#' @importFrom e1071 svm
#' @importFrom wrswoR sample_int_expj
#' @importFrom stats median na.omit
#' @importFrom scater calculatePCA
#' @importFrom edgeR cpm
#' @rawNamespace import(SingleCellExperiment, except = cpm)



scMEB <- function(sce, stable_idx, filtered = FALSE,
                  gamma = seq(1e-04,0.001,1e-05),
                  nu = 0.01, reject_rate = 0.1){

    counts <- counts(sce)
    dat_cpm <- edgeR::cpm(counts)
    cpm_threshold <- apply(dat_cpm, 1, stats::median)
    dat_detection <- dat_cpm > cpm_threshold

    if (filtered == FALSE){
        gene_ind_filter <- rowSums(dat_detection) >= 10
        cell_ind_filter <- colSums(dat_detection) >= 100

        sce <- SingleCellExperiment(list(counts = counts[gene_ind_filter,
                                                         cell_ind_filter]),
                                    colData = colData(sce)[cell_ind_filter,],
                                    rowData = rowData(sce)[gene_ind_filter,])
        dat_detection <- dat_detection[gene_ind_filter,cell_ind_filter]
        stable_idx <- intersect(stable_idx, rownames(sce)[gene_ind_filter])
    }

    if(ncol(sce)>1000){
        general.detection <- apply(dat_detection, 2, sum)
        sp <- wrswoR::sample_int_expj(ncol(sce), size=1000,
                                      prob=general.detection)
    }else{
        sp <- seq_len(ncol(sce))
    }

    dat_sample <- log2(counts(sce)+1)[,sp]
    dat.pca_scMEB <- scater::calculatePCA(t(dat_sample),ncomponents=50)

    trainData <- dat.pca_scMEB[stable_idx,]
    train_error <- numeric(length(gamma))
    for (k in seq_len(length(gamma))) {
        model <- e1071::svm(trainData, y = NULL, scale = FALSE,
                            type = "one-classification",
                            kernel = "radial", gamma = gamma[k],
                            nu = nu, tolerance = 0.001,
                            shrinking = TRUE, cross = 5, probability = FALSE,
                            fitted = TRUE,
                            na.action = stats::na.omit)
        train_error[k] <- 1 - model$tot.accuracy/100
    }
    gamma_num_new <- which.min(abs(train_error - reject_rate))

    model_new <- e1071::svm(trainData, y = NULL, scale = FALSE,
                            type = "one-classification",
                            kernel = "radial", gamma = gamma[gamma_num_new],
                            nu = nu, tolerance = 0.001,
                            shrinking = TRUE, cross = 5, probability = FALSE,
                            fitted = TRUE,
                            na.action = stats::na.omit)

    check <- numeric(nrow(dat.pca_scMEB))
    for (k in seq_len(nrow(dat.pca_scMEB))) {
        check[k] <-
            .decision_function(dat.pca_scMEB[k,],
                               model = model_new,
                               gamma = gamma[gamma_num_new]) - model_new$rho}

    list(model=model_new, dat_pca=dat.pca_scMEB, gamma = gamma[gamma_num_new],
         train_error = train_error, dist = check)
}



.kernel <- function(x,sv,gamma){
    k <- numeric(nrow(sv))
    for(i in seq_len(nrow(sv))){
        k[i] <- exp(-gamma*sum((sv[i,]-x)^2))
    }
    return(k)
}



.decision_function <- function(x,model,gamma){
    index <- model$coefs
    sv <- model$SV
    dec_val <- sum(index*.kernel(x,sv=sv,gamma))
    return(dec_val)
}







