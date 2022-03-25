#' Plot SpaRTaCo
#'
#' This function returns the ggplots of the mean and the spatial signal-to-noise ratios from a SpaRTaCo model.
#'
#' @import ggplot2
#' @export
#'
#' @param x a `spartaco` object;
#' @param type the type of plot you want to show.
#' `1` displays the co-cluster mean levels,
#' `2` displays the co-cluster spatial signal-to-noise ratios,
#' `3` displays the map of the spots colored with respect to the estimated column clusters,
#' and `4` displays the sample mean of some or of all the gene clusters.
#' @param k if `type = 4`, it displays the sample mean of each gene cluster is shown.
#' @param r if `type = 4`, it displays the sample mean of each gene cluster is shown.
#' @param title the plot title. If `NULL`, it does not add any title.
#' @param manual.palette is a vector of colors used when `type = "spots"`.
#' @return The requested plot is displayed. In addition, if assigned to an object, it will return the `ggplot` object.
#'
plot.spartaco <- function(x, type = 1, k = NULL, r = NULL, title = NULL, manual.palette = NULL, ...){
    if(class(x) != "spartaco") stop("the input file is not a spartaco object")
    if(length(type) > 1) type <- 1
    K <- nrow(x$mu)
    R <- ncol(x$mu)
    k.lab <- 1:K
    r.lab <- 1:R
    gr <- expand.grid(k.lab,r.lab)
    gr$Mu <- as.vector(x$mu)
    gr$Ratio <- as.vector(x$tau/x$xi)
    gr <- cbind(gr, expand.grid(as.vector(table(x$Cs))/nrow(x$x),as.vector(table(x$Ds))/ncol(x$x)))
    names(gr)[-c(3,4)] <- c("Y","X","height","width")

    prop.x <- as.vector(table(x$Ds))/ncol(x$x)
    prop.y <- as.vector(table(x$Cs))/nrow(x$x)
    xlim.left <- c(0, cumsum(prop.x)[-R])
    xlim.right <- cumsum(prop.x)
    ylim.left <- c(0, cumsum(prop.y)[-K])
    ylim.right <- cumsum(prop.y)


    # ---plot mu
    if(type == 1){
        Plots <- ggplot(gr, aes(xmin = as.vector(sapply(1:R, function(i) rep(xlim.left[i],K))),
                                     xmax = as.vector(sapply(1:R, function(i) rep(xlim.right[i],K))),
                                     ymin = rep(ylim.left,R),
                                     ymax = rep(ylim.right,R),
                                     fill = Mu)
        )+geom_rect()+theme_bw()+
            scale_fill_distiller(palette = "RdPu")+
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  axis.text=element_text(size=18),
                  axis.title=element_text(size=18),
                  legend.text = element_text(size = 18),
                  legend.title = element_text(size = 22),
                  #legend.position = "bottom",
                  plot.title = element_text(hjust = 0.5),
                  title = element_text(size=18),
                  plot.margin=grid::unit(c(3,2,3,2), "mm"))+
            labs(fill=expression(hat(mu)[kr]))+
            scale_x_continuous(breaks=(xlim.left+xlim.right)/2,
                               labels=paste("r =",1:R)
            )+
            scale_y_continuous(breaks=(ylim.left+ylim.right)/2,
                               labels=paste("k =",1:K)
            )
        if(!is.null(title)) Plots <- Plots + ggtitle(label = title)
        plot(Plots)
        invisible(Plots)
    }
    # ---plot tau/xi

    if(type == 2){
        Plots <- ggplot(gr, aes(xmin = as.vector(sapply(1:R, function(i) rep(xlim.left[i],K))),
                                     xmax = as.vector(sapply(1:R, function(i) rep(xlim.right[i],K))),
                                     ymin = rep(ylim.left,R),
                                     ymax = rep(ylim.right,R),
                                     fill = Ratio)
        )+geom_rect()+theme_bw()+
            viridis::scale_fill_viridis(discrete=FALSE)+
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  axis.text=element_text(size=18),
                  axis.title=element_text(size=18),
                  legend.text = element_text(size = 18),
                  legend.title = element_text(size = 22),
                  #legend.position = "bottom",
                  plot.title = element_text(hjust = 0.5),
                  title = element_text(size=18),
                  plot.margin=grid::unit(c(3,2,3,2), "mm"))+
            labs(fill=expression(hat(tau)[kr]/hat(xi)[kr]))+
            scale_x_continuous(breaks=(xlim.left+xlim.right)/2,
                               labels=paste("r =",1:R)
            )+
            scale_y_continuous(breaks=(ylim.left+ylim.right)/2,
                               labels=paste("k =",1:K)
            )
        if(!is.null(title)) Plots <- Plots + ggtitle(label = title)
        plot(Plots)
        invisible(Plots)
    }

    if(type == 3){
        # ---plot column clusters
        manual.palette <- c("red","yellow","lightblue","green","blue","purple","salmon","black","grey")
        Coord <- data.frame(x = x$coordinates[,2], y = -x$coordinates[,1], z = as.factor(x$Ds))
        Plots <- ggplot(Coord, aes(x, y, color = z))+geom_point(size = 3)+theme_bw()+
            labs(col = "")+
            scale_color_manual(values = manual.palette)+#scale_color_brewer(palette="Set1")
            labs(col = expression(D[r]))+
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  axis.text.x = element_blank(),#element_text(size=18),
                  axis.text.y = element_blank(),
                  axis.ticks.x = element_blank(),
                  axis.ticks.y = element_blank(),
                  axis.title=element_blank(),#element_text(size=18),
                  legend.text = element_text(size = 18),
                  legend.title = element_text(size = 22),
                  plot.title = element_text(hjust = 0.5),
                  #legend.position = "bottom",
                  #legend.spacing.x = unit(0.3, 'cm'),
                  title = element_text(size=18),
                  plot.margin=grid::unit(c(3,2,3,2), "mm"))+
            geom_point(shape = 1,size = 3,colour = "black")
        if(!is.null(title)) Plots <- Plots + ggtitle(label = title)
        plot(Plots)
        invisible(Plots)
    }

    if(type == 4){
        # ---plot sample means
        Coord <- data.frame(x = x$coordinates[,2], y = -x$coordinates[,1], z = as.factor(x$Ds))
        if(is.null(k)) k <- 1:K
        if(is.null(r)) r <- 1:R
        if(length(k) == 1){
            x.bar <- colMeans(x$x[x$Cs == k, which(x$Ds %in% r)])
            Plots <- ggplot(Coord[which(x$Ds %in% r),], aes(x, y, color = x.bar))+
                geom_point(size = 3)+theme_bw()+
                labs(col = "")+
                theme(panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      axis.text.x = element_blank(),#element_text(size=18),
                      axis.text.y = element_blank(),
                      axis.ticks.x = element_blank(),
                      axis.ticks.y = element_blank(),
                      axis.title=element_blank(),#element_text(size=18),
                      legend.text = element_text(size = 18),
                      legend.title = element_text(size = 22),
                      plot.title = element_text(hjust = 0.5),
                      #legend.position = "bottom",
                      #legend.spacing.x = unit(0.3, 'cm'),
                      title = element_text(size=18),
                      plot.margin=grid::unit(c(3,2,3,2), "mm"))+
                geom_point(shape = 1,size = 3,colour = "black")+
                ggtitle(label = paste("k =",k))
            plot(Plots)
            invisible(Plots)
        } else {
            Plots <- list()
            for(k.ind in 1:length(k)){
                x.bar <- colMeans(x$x[x$Cs == k[k.ind], which(x$Ds %in% r)])
                Plots[[k.ind]] <- ggplot(Coord[which(x$Ds %in% r),], aes(x, y, color = x.bar))+
                    geom_point(size = 3)+theme_bw()+
                    labs(col = "")+
                    theme(panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank(),
                          axis.text.x = element_blank(),#element_text(size=18),
                          axis.text.y = element_blank(),
                          axis.ticks.x = element_blank(),
                          axis.ticks.y = element_blank(),
                          axis.title=element_blank(),#element_text(size=18),
                          legend.text = element_text(size = 18),
                          legend.title = element_text(size = 22),
                          plot.title = element_text(hjust = 0.5),
                          title = element_text(size=18),
                          plot.margin=grid::unit(c(3,2,3,2), "mm"))+
                    geom_point(shape = 1,size = 3,colour = "black")+
                    ggtitle(label = paste("k =",k[k.ind]))
                plot(Plots[[k.ind]])
            }
            invisible(Plots)
        }
    }
}
