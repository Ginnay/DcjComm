#' Plot a chord diagram showing the amount of cell-cell communication
#'
#' @description Create an interactive bipartite chord diagram showing number of observation among two categories of groups (var1-var2).
#'
#' @param data A dataframe containing the results of signaling analysis.
#' @param var1 A character string representing the first category group. It must be a subset of \code{colnames(data)}.
#' @param var2 A character string representing the second category group. It must be a subset of \code{colnames(data)}.
#' @param group_fontsize Numeric font size in pixels for the group name labels.
#' @param margin Numeric margin in pixels between the outer diagram radius and the edge of the display.
#' @param var1_label A character string to be used for category var1 label
#' @param var2_label A character string to be used for category var2 label
#'
#' @return The function returns an image in the plot window of Rstudio (Viewer)
#'
#' @export

scSeqComm_chorddiag_cardinality <- function(data, var1 = "cluster_R", var2 = "cluster_L", group_fontsize = 12, margin = 90,
                                            var1_label = "Receptor-expressing cell types", var2_label ="Ligand-expressing cell types")
{
  #select variables from data.frame
  data <- data[,c(var1,var2)]

  #create a square matrix containing number of observations of (var1 - var2) couple: var1 on rows (right side) and var2 on columns (left side)
  df <- as.matrix(igraph::as_adjacency_matrix(tidygraph::as_tbl_graph(data)))

  #number of groups names
  n <- ncol(df)

  #set color palette for each group of var2 (left side)
  if (n>9)
  {
    palette1 <- RColorBrewer::brewer.pal(11, "Spectral")
    palette1 <- colorRampPalette(palette1)(n)
  }
  else{
    palette1 <- RColorBrewer::brewer.pal(n, "Spectral")
  }

  #set color palette for each group of var1 (right side)
  palette2 <- rep("#838383",n)

  #chordiagram
  c <- chorddiag::chorddiag(df, type="bipartite",groupColors = c(palette2,palette1),
                 showTicks = F, groupnameFontsize = group_fontsize, groupnamePadding = 5, margin = margin,
                 categoryNames = c(var1_label, var2_label),
                 categorynamePadding = 110, categorynameFontsize = 20)
  return(c)
}

#' Plot an heatmap showing the amount of cell-cell communication
#'
#' @description Creates an heatmap showing the number of observations among two categories of groups (x-y)
#'
#' @param data A dataframe containing results of signaling analysis.
#' @param x A character string describing the variable to be plotted on x axis. It must be a subset of \code{colnames(data)}.
#' @param y A character string describing the variable to be plotted on y axis. It must be a subset of \code{colnames(data)}.
#' @param cell_types A character vector containing the cell types to be plotted in the heatmap.
#' @param limit_fill A numeric vector of length two providing limits of the scale. Use NULL to refer to the existing minimum or maximum in \code{data}.
#' @param xlab A character string to be used as title of x axis.
#' @param ylab A character string to be used as title of y axis.
#' @param title A character string to be used as title of plot.
#'
#' @return The function returns an image in the plot window of Rstudio
#'
#' @export

scSeqComm_heatmap_cardinality <- function(data, y = "cluster_L", x = "cluster_R", limit_fill = NULL,
                                          cell_types = NULL,
                                          ylab = "Ligand-expressing cell types", xlab = "Receptor-expressing cell types", title)
{
  windowsFonts(myFont = windowsFont("Times New Roman"))
  #select variables from data
  data <- data[,c(x,y)]

  #count the number of observations of x-y pair
  df <- data %>% count(data[[x]],data[[y]])

  if (is.null(cell_types)) cell_types = unique(c(data[[x]],data[[y]]))
  res <- gtools::permutations(n=length(cell_types),r=2,v=cell_types,repeats.allowed=T)
  colnames(df) <- c(x,y,"cardinality")
  colnames(res) <- c(x,y)
  df <- left_join(as.data.frame(res[,1:2]), as.data.frame(df))
  df[is.na(df)]=0

  #str(df)

  if (is.null(limit_fill)) limit_fill=c(0, max(df[,"cardinality"]))

  g <- ggplot(df,aes_string(x = x, y = y, fill = "cardinality"))+
    geom_tile(colour="white",size=0.5)+
    geom_text(aes(label = cardinality),colour="black") +
    labs(x= xlab,y=ylab, title=title)+ #labels
    scale_fill_continuous(type = "viridis",alpha=0.7, limits = limit_fill, na.value = "white")+
    theme_grey(base_size=10)+
    theme(legend.position="right",legend.direction="vertical",
          legend.title=element_text(colour="black"),
          legend.margin=margin(grid::unit(0.8,"cm")),
          legend.text=element_text(colour="black",size=9,family = 'myFont'),
          legend.key.height=grid::unit(0.5,"cm"),
          legend.key.width=grid::unit(0.5,"cm"),
          axis.text.x=element_text(size=9,colour="black",angle=45, hjust=1, vjust = 1,family = 'myFont'),
          axis.text.y=element_text(size=9,vjust=0.2,colour="black",family = 'myFont'),
          axis.title= element_text(size=12, vjust=0.5,family = 'myFont'),
          axis.ticks=element_line(size=0.4),
          plot.background=element_blank(),
          panel.border=element_blank(),
          plot.margin=margin(0.7,0.4,0.1,0.2,"cm"),
          plot.title=element_text(colour="black",hjust=0,size=10,face="bold"),
          axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
          axis.title.x = element_text(margin = margin(t = 15, r = 0, b = 0, l = 0)))

  return(g)
}


#' Summary plot (heatmap)
#'
#' @description The function creates a heatmap based on specified parameters.
#'
#' @param data A data frame containing results of signaling analysis.
#' @param x A character string describing the variable to be mapped on x axis. It must be a subset of \code{colnames(data)}.
#' @param y A character string describing the variable to be mapped on y axis. It must be a subset of \code{colnames(data)}.
#' @param fill A character string describing the variable to be mapped in color scale. It must be a subset of \code{colnames(data)}.
#' @param size a character string describing the variable to be mapped in point size scale. It must be a subset of \code{colnames(data)}.
#' @param title A character string to be used as title of plot.
#' @param xlab A character string to be used as title of x axis.
#' @param ylab A character string to be used as title of y axis.
#' @param fill_lab A character string to be used as label for fill scale.
#' @param size_lab A character string to be used as label for size scale.
#' @param limit_fill A numeric vector of length two describing the fill scale limits.
#' @param limit_size A numeric vector of length two describing the size scale limits.
#' @param breaks_size A numeric vector of positions of breaks.
#' @param facet_grid_x 	A variable defining faceting group on columns dimension.
#' @param facet_grid_y A variable defining faceting group on row dimension.
#' @param annotation_GO A data.frame or a list of data.frame with results of enrichment analysis.
#' @param cutoff A numeric value defining the significance value for GO terms.
#' @param topGO A numeric value defining the maximum number of terms to be visualized.
#' @param GO_ncol A numeric value defining the number of columns on which show the GO terms table.
#'
#' @return The function returns the image in the plot window of Rstudio
#'
#' @export
scSeqComm_heatmap <- function(data, x, y, fill = "S_intra", size = "S_inter", title, xlab, ylab="",
                              fill_lab ="S_intra", size_lab = "S_inter",
                              limit_fill=c(0,1),limit_size= c(0,1),breaks_size=c(0.2,0.5,0.8),
                              facet_grid_x=NULL, facet_grid_y=NULL,
                              annotation_GO = NULL, cutoff = 0.05, topGO = 5, GO_ncol = 1)
{
  #if user does not provide limits values for fill and/or size, they are computed by the function
  if (is.null(limit_fill))
  {
    limit_fill <- c(min(data[[fill]]),max(data[[fill]]))
  }
  if (is.null(limit_size) & !is.null(size))
  {
    limit_size <- c(min(data[[size]]),max(data[[size]]))
  }

  if(!is.null(annotation_GO))
  {
    # order LR pairs by receptor (for visualization purposes)
    data$LR_pair <- factor(data$LR_pair, levels = unique(data$LR_pair[order(data$receptor)]))
  }

  #if heatmap is a 3 or 4 dimensional plot
  if(!is.null(size))
  {
    g <- ggplot(data, aes_string(x = x, y = y, col = fill, size = size))+ geom_point() +
      labs(x = xlab,y = ylab, col = fill_lab, size = size_lab, title = title)+ #labels
      scale_size_binned(range = c(0.01,5), breaks = breaks_size, trans ="exp", limits = limit_size)+ add2ggplot::theme_white()+ #size scale
      scico::scale_color_scico(palette = "lajolla", limit = limit_fill, na.value = "grey")  #color scale
  }else{
    g <- ggplot(data, aes_string(x = x, y = y, fill = fill))+
      labs(x = xlab,y = ylab, fill = fill_lab, title = title)+ #labels
      geom_tile(colour = "white",size = 0.2) +
      scico::scale_fill_scico(palette = "lajolla", limit = limit_fill, na.value = "grey") #color scale
  }

  #theme settings
  g <- g +
    #theme_grey(base_size=10)+
    add2ggplot::theme_white()+
    theme(legend.position="right",legend.direction="vertical",
          legend.title=element_text(colour="black"),
          legend.margin=margin(grid::unit(0.5,"cm")),
          legend.text=element_text(colour="black",size=7),
          legend.key.height=grid::unit(0.6,"cm"),
          legend.key.width=grid::unit(0.5,"cm"),
          strip.text.x = element_text(size=7, color="black",
                                      face="bold.italic"),
          strip.text.y = element_text(size=7, color="black",
                                      face="bold.italic",angle=0),
          axis.text.x=element_text(size=8,colour="black",angle=45, hjust=1, vjust = 1),
          axis.text.y=element_text(vjust=0.2,colour="black",size=8),
          axis.ticks=element_line(size=0.4),
          plot.background=element_blank(),
          panel.border=element_blank(),
          plot.margin=margin(0.7,0.4,0.1,0.2,"cm"),
          plot.title=element_text(colour="black",hjust=0,size=12,face="bold"),
          axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
          axis.title.x = element_text(margin = margin(t = 15, r = 0, b = 0, l = 0)))

  #based on faceting variables, the panel is divided in grid

  if (!is.null(facet_grid_x) & !is.null(facet_grid_y)){
    g <- g + facet_grid(data[[facet_grid_y]] ~ data[[facet_grid_x]],switch = "x", scales = "free", space = "free")
  }
  if (is.null(facet_grid_x) & !is.null(facet_grid_y)){
    g <- g + facet_grid(data[[facet_grid_y]] ~ .,scales = "free_y", space = "free_y")
  }
  if (!is.null(facet_grid_x) & is.null(facet_grid_y)){
    g <- g + facet_grid(~ data[[facet_grid_x]], switch = "x", scales = "free_x", space = "free_x")
  }

  #if user specify annotation results, GO terms are visualized

  if (!is.null(annotation_GO))
  {
    #if annotation results is list, more than one table are displayed
    if(class(annotation_GO)=="list")
    {
      #gtables <- sapply(names(annotation_GO)[1:3],function(y) create_gtable(annotation_GO[[y]],y, cutoff, topGO))
      gtables <- sapply(names(annotation_GO),function(y) create_gtable(annotation_GO[[y]],y, cutoff, topGO))
      #g1 <- do.call("gridExtra::grid.arrange", c(gtables, list(ncol=2)))
      #print(class(gtables))
      #g1 <- gridExtra::grid.arrange(gtables, ncol=2)
      g1 <- gridExtra::arrangeGrob(grobs = gtables, ncol=GO_ncol)

      #print(class(g1))
    }
    else{
      g1 <- create_gtable(annotation_GO," ",cutoff, topGO)
      #print(class(g1))
    }

    #arrange plots
    #print(class(g))
    g <- gridExtra::grid.arrange(g,g1, ncol=2)
  }

  return(g)
}


#function that create a gtable containing text grobs representing a character matrix.
create_gtable <- function(x, name, cutoff, topGO){

  #significant GO terms
  data_GO <- x[x$pval.adjusted < cutoff,c("GO.ID","Term")]
  if (nrow(data_GO) > topGO) data_GO <- data_GO[1:topGO,]

  #if there are none significant GO terms
  if(nrow(data_GO)==0) {return(grid::textGrob(paste(name,": none GO terms"), gp=grid::gpar(fontsize=14, fontface=4)))}
  #otherwise the gtable is built
  t <- grid::textGrob(paste(name), gp=grid::gpar(fontsize=14, fontface=4))#, hjust=0, x=0)
  g <- gridExtra::tableGrob(data_GO, theme = gridExtra::ttheme_default(base_size=8))
  g <- g %>%
    gtable::gtable_add_rows(heights = grid::grobHeight(t) + unit(2,"mm"), pos = 0L) %>%
    gtable::gtable_add_grob(t, t=1,b=1,l=2,r=ncol(g))

  return(g)
}





#' Plot intercellular and intracellular scores (heatmap)
#'
#' Function to plot the intercellular (mapped in point size scale) and intracellular (mapped in color scale) scores of the given input data with cluster pairs on x axis and ligand-receptor pairs on y axis.
#' Grey color is associated to NA values of intracellular score.
#'
#' @param data A data frame containing results of signaling analysis.
#' @param title A character string to be used as title of plot.
#' @param facet_grid_x 	A variable defining faceting group on columns dimension.
#' @param facet_grid_y A variable defining faceting group on row dimension.
#' @param annotation_GO A data.frame or a list of data.frame with results of enrichment analysis.
#' @param cutoff A numeric value defining the significance value for GO terms.
#' @param topGO A numeric value defining the maximum number of terms to be visualized.
#' @param GO_ncol A numeric value defining the number of columns on which show the GO terms table.
#'
#' @return Plot showing intracellular and intercellular scores
#'
#' @export
scSeqComm_plot_scores <- function(data, title,
                                  facet_grid_x="cluster_R", facet_grid_y="receptor",
                                  annotation_GO = NULL, cutoff = 0.05, topGO = 5, GO_ncol = 1){

  return(scSeqComm_heatmap(data = data, title = title,
                           y = "LR_pair", x = "interaction",
                           ylab = "Ligand-receptor pairs" , xlab = "L expressing cluster --> R expressing cluster",
                           fill = "S_intra", size = "S_inter",
                           limit_fill=c(0,1),limit_size= c(0,1),breaks_size=c(0.2,0.5,0.8),
                           facet_grid_x=facet_grid_x, facet_grid_y=facet_grid_y,
                           annotation_GO = annotation_GO, cutoff = cutoff, topGO = topGO, GO_ncol = GO_ncol))


}

#' Plot intercellular and intracellular scores with pathway information (heatmap)
#'
#' Function to plot the intercellular (mapped in point size scale) and intracellular scores (mapped in color scale) of the given input data with pathway information on y axis and cluster pairs on x axis.
#' Grey color is associated to NA values of intracellular score.
#'
#' @param data A data frame containing results of signaling analysis.
#' @param title A character string to be used as title of plot.
#' @param facet_grid_x 	A variable defining faceting group on columns dimension.
#' @param annotation_GO A data.frame or a list of data.frame with results of enrichment analysis.
#' @param cutoff A numeric value defining the significance value for GO terms.
#' @param topGO A numeric value defining the maximum number of terms to be visualized.
#' @param GO_ncol A numeric value defining the number of columns on which show the GO terms table.
#'
#' @return Plot showing intracellular and intercellular scores and pathway
#'
#' @export
scSeqComm_plot_scores_pathway <- function(data, title,
                                          facet_grid_x="cluster_R",
                                          annotation_GO = NULL, cutoff = 0.05, topGO = 5, GO_ncol = 1){

  return(scSeqComm_heatmap(data = data, title = title,
                           y = "pathway", x = "interaction",
                           ylab = "Pathway" , xlab = "L expressing cluster --> R expressing cluster",
                           fill = "S_intra", size = "S_inter",
                           limit_fill=c(0,1),limit_size= c(0,1),breaks_size=c(0.2,0.5,0.8),
                           facet_grid_x=facet_grid_x, facet_grid_y="LR_pair",
                           annotation_GO = annotation_GO, cutoff = cutoff, topGO = topGO, GO_ncol = GO_ncol))


}


#' Plot intercellular and intracellular scores of a given ligand-receptor pair across cell clusters (heatmap)
#'
#' Function to plot intercellular (mapped in point size scale) and intracellular (mapped in color scale) scores of a given ligand-receptor pair across cell clusters.
#' If input data contains more than one ligand-receptor pair, a plot will be generated for each ligand-receptor pair.
#' Grey color is associated to NA values of intracellular score.
#'
#' @param data A data frame containing results of signaling analysis.
#' @param title A character string to be used as title of plot.
#' @param selected_LR_pair String or array of strings containing ligand-receptor pairs to plot.
#'
#' @return A plot (or a list of plots)
#'
#' @export
scSeqComm_plot_LR_pair <- function(data, title, selected_LR_pair = NULL){

  if(!is.null(selected_LR_pair)){
    pairs <- selected_LR_pair
  }else{
    pairs <- unique(data$LR_pair)
  }

  if(length(pairs)>1){
    pairs_list <- list()
    for(pair in pairs){
      pairs_list[[pair]] <- scSeqComm_heatmap(data = scSeqComm_select(data, LR_pair = pair), title = pair,
                                              y = "cluster_L", x = "cluster_R",
                                              ylab = "Ligand-expressing cluster" , xlab = "Receptor-expressing cluster",
                                              fill = "S_intra", size = "S_inter",
                                              limit_fill=c(0,1),limit_size= c(0,1),breaks_size=c(0.2,0.5,0.8),
                                              facet_grid_x=NULL, facet_grid_y=NULL,
                                              annotation_GO = NULL, cutoff = 0.05, topGO = 5, GO_ncol = 1)
    }
    return(pairs_list)
  }else{
    if(!is.null(selected_LR_pair)){
      data <- scSeqComm_select(data, LR_pair = selected_LR_pair)
    }
    return(scSeqComm_heatmap(data = data, title = title,
                             y = "cluster_L", x = "cluster_R",
                             ylab = "Ligand-expressing cluster" , xlab = "Receptor-expressing cluster",
                             fill = "S_intra", size = "S_inter",
                             limit_fill=c(0,1),limit_size= c(0,1),breaks_size=c(0.2,0.5,0.8),
                             facet_grid_x=NULL, facet_grid_y=NULL,
                             annotation_GO = NULL, cutoff = 0.05, topGO = 5, GO_ncol = 1))
  }

}


#' Plot intercellular and intracellular scores of ligand-receptor pairs between a given cell clusters pair (heatmap)
#'
#' Function to plot intercellular (mapped in point size scale) and intracellular scores (mapped in color scale) of ligand-receptor pairs between a given cell clusters pair.
#' If input data contains more than one cell clusters pair, a plot will be generated for each cell clusters pair.
#' Grey color is associated to NA values of intracellular score.
#'
#' @param data A data frame containing results of signaling analysis.
#' @param title A character string to be used as title of plot.
#' @param selected_cluster_pair String or array of strings containg cell clusters pairs to plot.
#'
#' @return A plot (or a list of plots)
#'
#' @export
scSeqComm_plot_cluster_pair <- function(data, title, selected_cluster_pair = NULL){

  if(!is.null(selected_cluster_pair)){
    cl_pairs <- selected_cluster_pair
  }else{
    cl_pairs <- unique(data$interaction)
  }

  if(length(cl_pairs)>1){
    cl_pairs_list <- list()
    for(cl_pair in cl_pairs){
      cl_pairs_list[[cl_pair]] <- scSeqComm_heatmap(data = scSeqComm_select(data, interaction = cl_pair), title = cl_pair,
                                                 y = "ligand", x = "receptor",
                                                 ylab = "Ligand" , xlab = "Receptor",
                                                 fill = "S_intra", size = "S_inter",
                                                 limit_fill=c(0,1),limit_size= c(0,1),breaks_size=c(0.2,0.5,0.8),
                                                 facet_grid_x=NULL, facet_grid_y=NULL,
                                                 annotation_GO = NULL, cutoff = 0.05, topGO = 5, GO_ncol = 1)
    }
    return(cl_pairs_list)
  }else{
    if(!is.null(selected_cluster_pair)){
      data <- scSeqComm_select(data, interaction = selected_cluster_pair)
    }
    return(scSeqComm_heatmap(data = data, title = title,
                             y = "ligand", x = "receptor",
                             ylab = "Ligand" , xlab = "Receptor",
                             fill = "S_intra", size = "S_inter",
                             limit_fill=c(0,1),limit_size= c(0,1),breaks_size=c(0.2,0.5,0.8),
                             facet_grid_x=NULL, facet_grid_y=NULL,
                             annotation_GO = NULL, cutoff = 0.05, topGO = 5, GO_ncol = 1))
  }
}
