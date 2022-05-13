#' Plot the gene-specific variances
#'
#' This function returns the ggplots of the distribution of the gene-specific variances within the spot clusters discovered by SpaRTaCo.
#' In addition, it returns a data frame with the most variable genes in each of the spot clusters requested.
#'
#' @import ggplot2
#' @export
#'
#' @param x an object of class `spartaco.genes`;
#' @param r the spot clusters to be displayed (by default, they are all displayed);
#' @param g the number of highly variable genes in each spot cluster to be highligthed.
#' @param return.plots if `FALSE`, it returns only the data frame containing the `g` most variable genes into the spot clusters given by `r`.
#' @param plot.labels a list containing the plot parameters of the most variable genes labels (see **Details**).
#'
#' @return If `return.plots == T`, it returns the requested plots as `ggplot` objects. In addition, if `g > 0`, it returns a list of the top `g` highly variable genes in each spot cluster.
#'
#' @details When `g > 0`, the posterior distribution of top `g` genes with the largest expectations are displayed in red. In addition, a label with the gene names is placed close to each distribution.
#' The plot parameters of the labels are passed to `plot.labels`, that is a list of four elements:
#'  - `angle` gives the rotation of the labels;
#'  - `size` gives the size of the labels;
#'  - `hjust` gives the horizontal adjustment;
#'  - `vjust` gives the vertical adjustment.
#'
#' You can modify just some of the four parameters, without need to redefine the entire list.
#'
#' @references
#'
#'Sottosanti, A. and Risso, D. (2021+) Co-clustering of Spatially Resolved Transcriptomic Data [(preprint)](https://arxiv.org/abs/2110.04872)
#'
#'
#' @examples
#' library(spartaco)
#' # Not run:
#' # x <- matrix(runif(n*p), n, p)
#' # coordinates <- matrix(runif(2*p), p, 2)
#' # output <- spartaco(x = x, coordinates = coordinates, K = K, R = R)
#'
#' sigma2 <- GeneVariances(output)
#'
#' # To plot the first 10 highly variable genes within each spot cluster:
#' plot(sigma2, g = 10)
#'
#' # To plot the first 10 highly variable genes within the spot clusters 1 and 3 (let us assume `p>=3`):
#' plot(sigma2, g = 10, r = c(1, 3))


plot.spartaco.genes <- function(x, r = 1:ncol(x$Expectation), g = 5, return.plots = T, plot.labels = list(angle = 45, size = 4, hjust = 0, vjust = 0)){
    if(is.null(plot.labels$angle)) plot.labels$angle <-  45
    if(is.null(plot.labels$size)) plot.labels$size <- 4
    if(is.null(plot.labels$hjust)) hjust <- 0
    if(is.null(plot.labels$vjust)) vjust <- 0
    gene.names <- factor(row.names(x$Expectation), levels = row.names(x$Expectation))
    K <- length(unique(x$Cs))
    R <- ncol(x$Expectation)
    r.values <- r
    if(is.null(g)) g <- 0
    if(g > 0){
        first.g.genes <- data.frame(matrix(NA, g, length(r.values)))
        names(first.g.genes) <- paste("r=", r.values, sep="")
    }
    Plots <- list()
    theme.settings <- theme(plot.title = element_text(hjust = 0.5),
                            axis.text=element_text(size=18),
                            axis.title=element_text(size=18),
                            legend.text = element_text(size = 15),
                            legend.title = element_text(size = 15),
                            legend.position = "bottom",
                            title = element_text(size=18),
                            axis.ticks.x = element_blank(),
                            axis.text.x = element_text(
                                angle = plot.labels$angle,
                                vjust = 1, hjust=1,
                                size = plot.labels$size))
    j <- 1
    for(r in r.values){
            gene.names.first.g <- as.character(rep("", length(gene.names)))
            if(g > 0){
                ranked <- order(x$Expectation[,r], decreasing = T )[1:g]
                gene.names.first.g[ranked] <- first.g.genes[,j] <- as.character(gene.names[ranked])
            }
            df <- data.frame(
                Genes = gene.names,
                y = x$Expectation[,r],
                Cs = as.factor(x$Cs))
                df$v <- x$HPD.left[,r]
                df$w <- x$HPD.right[,r]

            if(return.plots){
                Plots[[j]] <- local({
                    p <- ggplot(df, aes(Genes, y, col = Cs))+
                        geom_pointrange(data = df, aes(ymin = v, ymax = w))+
                        scale_color_grey(start=0.4, end=0.7) + theme_classic()+
                        labs(x = "gene", y = expression(sigma^2~"|"~data), col = "cluster")+
                        ggtitle(paste("r = ",r,sep=""))+
                            scale_x_discrete(labels = rep("",length(gene.names)))+
                        theme.settings
                    if(g > 0){
                        gene.names.first.g <- gene.names.first.g
                        ranked <- ranked
                        p <- p+
                        #geom_pointrange(data = df[ranked,], aes(ymin = v, ymax = w), col = "darkred")+
                        geom_point(data = df[ranked,], aes(x = Genes, y = y), col = "darkred")+
                        geom_text(data = df[ranked,], mapping = aes(x = Genes, y = w, label = gene.names.first.g[ranked]),
                                  angle = plot.labels$angle,
                                  vjust = plot.labels$vjust,
                                  hjust = plot.labels$hjust,
                                  size = plot.labels$size, col = "black")}
                    p})
                #plot(Plots[[j]])
                } else Plots <- NULL
            j <- j + 1
    }
    if(g > 0) output <- list(first.g.genes = first.g.genes, Plots = Plots) else
        output <- Plots
    return(output)
}
