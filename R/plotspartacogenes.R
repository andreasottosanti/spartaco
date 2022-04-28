#' Plot the gene-specific variances given by SpaRTaCo
#'
#' This function returns the ggplots of the distribution of the gene variances estimated by SpaRTaCo.
#'
#' @import ggplot2
#' @export
#'
#' @param x a `spartaco.genes` object;
#' @param use.gene.clusters if `TRUE`, it returns the plots for every row cluster.
#' @param g if not `NULL`, it is the number of genes in decreasing order of the expected value whose names are displayed.
#' @return The requested plots are displayed. In addition, if assigned to an object, it will return the `ggplot` object.
#'
#'
plot.spartaco.genes <- function(x, r = 1:ncol(x$Expectation), display.plots = T, g = 5, angle = 45, gene.label.size = 4, hjust = 0, vjust = 0){
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
                            legend.text = element_text(size = 18),
                            legend.title = element_text(size = 22),
                            legend.position = "bottom",
                            title = element_text(size=18),
                            axis.ticks.x = element_blank(),
                            axis.text.x = element_text(
                                angle = angle, vjust = 1, hjust=1, size = gene.label.size))
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
            if(is.null(x$HPD.left)){
                    df$v <- x$Expectation[,r] - 1.96*sqrt(x$Variance[,r])
                    df$w <- x$Expectation[,r] + 1.96*sqrt(x$Variance[,r])} else {
                        df$v <- x$HPD.left[,r]
                        df$w <- x$HPD.right[,r]
                    }
            if(display.plots){
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
                                  angle = angle, vjust = vjust, hjust = hjust, size = gene.label.size, col = "black")}
                    p})
                plot(Plots[[j]])
                } else Plots <- NULL
            j <- j + 1
    }
    if(g > 0) output <- list(first.g.genes = first.g.genes, Plots = Plots) else
        output <- Plots
    return(output)
}
