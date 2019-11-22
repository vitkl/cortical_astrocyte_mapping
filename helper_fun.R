# helper functions
library(ggplot2)
load_gene_loadings = function(dir = "./data/Holt_astrocytes/out/", prefix = "res", k = 3){
  
  gene_load = fread(paste0(dir, prefix, k, ".gene_score.txt"),
                    header = F, data.table = F)
  gene_load_names = fread(paste0(dir, "genes.txt"),
                          header = F)
  rownames(gene_load) = gene_load_names$V2
  # expand table to match rows in sce object
  gene_load = gene_load[rownames(normalised), ]
  rownames(gene_load) = rownames(normalised)
  colnames(gene_load) = paste0("k", k, "_", seq_len(ncol(gene_load)))
  
  gene_load
}

plot_factor_loadings = function(normalised, k = 3,
                                gene_list = c("Chrdl1", "Id1", "Il33", "Eogt", "Gfap", "Id3")) {
  
  factors_dt = as.data.table(as.data.frame(rowData(normalised)[gene_list,
                                                               paste0("k",k,"_", seq_len(k))]), keep.rownames = "gene")
  factors_dt = melt.data.table(factors_dt, "gene")
  
  ggplot(factors_dt, aes(gene, variable, color = value, fill = value)) +
    geom_tile() + scale_color_viridis_c() + scale_fill_viridis_c() +
    ggtitle(paste0(k," factors"))
}

# plot zonation

plot_gradient = function(zonation_mean, gene_list, gene_order = NULL,
                         layer_levels =  c("L1", "L2-3", "L4", "L5", "L6"),
                         lm_genes = NULL, #c("Chrdl1", "Id3", "Il33"),
                         n_matching_genes = 3, include_lm = F,
                         ref_zonation_mean = NULL,
                         log_offset = min(as.matrix(zonation_mean, rownames = "genes")[as.matrix(zonation_mean, rownames = "genes") > 0]),
                         limits = NULL, normalise = c(TRUE, FALSE, "log10")[1]){
  
  zonation_mean_dt = copy(zonation_mean[genes %in% gene_list])
  
  # or center of the peak
  #gene_order_dt = as.data.table(t(apply(-zonation_mean_mat, 1, frank)))
  #gene_order_dt$genes = rownames(zonation_mean_mat)
  #setorder(gene_order_dt, -V1, -V2, -V3)
  #gene_order_peak = gene_order_dt$genes
  
  # if reference profiles provided add them to all profiles
  if(!is.null(ref_zonation_mean)){
    ref_zonation_mean_dt = copy(ref_zonation_mean[genes %in% lm_genes])
    zonation_mean_dt = rbind(zonation_mean_dt, ref_zonation_mean_dt)
  }
  # convert to matrix
  zonation_mean_mat = as.matrix(zonation_mean_dt, rownames = "genes")
  
  # create data.table for plotting
  for_plot = melt.data.table(zonation_mean_dt, id.vars = "genes", variable.name = "layer", value.name = "expression")
  
  # create order of layers in a plot as factor
  for_plot[, layer := factor(layer, levels = layer_levels[length(layer_levels):1])]
  for_plot[, layer_num := as.integer(layer)]
  
  # add log offset
  for_plot[, expression := expression + log_offset, by = .(genes)]
  
  # normalise within each gene
  if(isTRUE(as.logical(normalise[1]))){
    for_plot[, expression := expression / sum(expression), by = .(genes)]
  } else if(normalise[1] == "log10") {
    for_plot[, expression := log10(expression), by = .(genes)]
  }
  
  if(!is.null(lm_genes)){
    # compute distances to landmark genes
    ref_dist_mat_1 = kl_dist(for_plot, symmetric = F, type = "kl",
                             ids1 = unique(for_plot$genes), ids2 = lm_genes,
                             id_col = "genes", freq_col = "expression") / 2
    ref_dist_mat_2 = kl_dist(for_plot, symmetric = F, type = "kl",
                             ids2 = unique(for_plot$genes), ids1 = lm_genes,
                             id_col = "genes", freq_col = "expression") / 2
    ref_dist_mat = ref_dist_mat_1 + t(ref_dist_mat_2)
    
    # find genes that match profiles of landmark genes
    matching_genes = sapply(lm_genes, function(lm_gene) {
      shortest = sort(ref_dist_mat[lm_gene, ])
      if(include_lm) start = 1 else start = 2
      names(shortest)[seq(start, n_matching_genes+1)]
    })
    
    # set order by similarity to markers
    gene_order = unique(as.character(matching_genes))
    
    # select genes matching those landmark profiles for plotting
    for_plot = for_plot[genes %in% as.character(matching_genes) & 
                          (!genes %in% as.character(lm_genes) | include_lm)]
  } else matching_genes = NULL
  
  if(is.null(gene_order)){
    
    # compute distances between profiles & do hierarchical clustering
    dist_mat = kl_dist(for_plot, symmetric = T, type = "kl",
                       id_col = "genes", freq_col = "expression")
    dist = as.dist(dist_mat)
    hcl = hclust(dist)
    
    # create gene order by hierarchical clustering
    gene_order_hcl = rownames(dist_mat)[hcl$order]
    
  } else {
    # use order provided
    gene_order_hcl = gene_order
  }
  
  gene_order_hcl = gene_order_hcl[!is.na(gene_order_hcl)]
  
  # create order of genes in a plot as factor
  for_plot[, genes := factor(genes, levels = gene_order_hcl[seq(1, length(gene_order_hcl))])]
  
  # set colors and theme
  col = colorRampPalette(RColorBrewer::brewer.pal(9, "YlOrRd"))(100)
  theme_change = theme(
    strip.background = element_blank(),
    strip.text = element_blank(),
    plot.background = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.border =element_blank(),
    #panel.border = element_rect(colour = "black", fill=NA, size=0.5),
    axis.line = element_blank(),
    panel.background = element_blank(),
    legend.key.size = unit(1, "lines"),
    #legend.text=element_text(size=6),
    #legend.title=element_text(size=6),
    plot.title = element_text(hjust = 0.5, size=12),
    axis.title = element_blank(),
    #axis.text.y = element_text(size=7),
    axis.ticks = element_blank(),
    axis.text.x = element_text(angle = 90,
                               #size=7,
                               hjust=1,
                               vjust=0.6),
    #plot.margin = margin(0,-1,0,-1,unit="pt"),
    legend.position = "right",
    legend.direction = "vertical"
  )
  
  plot = ggplot(for_plot, aes(genes, layer,
                              fill = expression, color = expression)) +
    geom_tile() 
  if(!is.null(limits)) {
    plot  = plot + scale_fill_gradientn(colours=col, limits = limits) +
      scale_color_gradientn(colours=col, limits = limits) 
  } else {
    plot  = plot + scale_fill_gradientn(colours=col) +
      scale_color_gradientn(colours=col) 
  }
  
  plot = plot + scale_y_discrete(expand=c(0,0)) +
    theme_change +
    guides(fill = guide_colorbar(title.position = "top")) +
    theme(panel.spacing = unit(1.2,"lines"))
  
  list(plot = plot,
       gene_order = gene_order_hcl, dist = dist, 
       matching_genes = matching_genes,
       theme_change = theme_change, col = col)
}

plot_gradient2 = function(spatial, gene_list, normalise = FALSE) {
  sel_spatial = spatial[genes %in% gene_list]
  sel_spatial[, normalisedDepth_bin := round(normalisedDepth, 1)]
  sel_spatial[, spotcounts_bin := mean(spotcounts), by = .(genes, normalisedDepth_bin)]
  for_plot = unique(sel_spatial[, .(genes, normalisedDepth_bin, spotcounts_bin)])
  for_plot[, spot_frequency_bin := spotcounts_bin / sum(spotcounts_bin), by = .(genes)]
  
  dist_mat = kl_dist(for_plot, symmetric = T, type = "kl")
  set.seed(1)
  hcluster = hclust(as.dist(dist_mat))
  hierarchy = as.dendrogram(hcluster)
  
  gplots::heatmap.2(dist_mat,
                    Rowv = (hierarchy),
                    Colv = (hierarchy),
                    trace = "none",
                    col = viridis::viridis_pal())
  
  for_plot[, genes := factor(genes, levels = hcluster$labels[hcluster$order])]
  
  if(isTRUE(normalise)){
    for_plot[, spotcounts_bin := spotcounts_bin / sum(spotcounts_bin),
             by = .(genes)]
  } else {
    for_plot[, spotcounts_bin := log2(spotcounts_bin)]
  }
  
  print(ggplot(for_plot, aes(x = genes, y = -normalisedDepth_bin,
                             fill = (spotcounts_bin), colour = (spotcounts_bin))) + 
          geom_raster() +
          scale_fill_viridis_c() + scale_colour_viridis_c()+
          theme(axis.text.x = element_text(angle = -45)))
  list(hierarchy = hierarchy, hcluster = hcluster)
}



kl_dist = function(dt, id_col = "genes",
                   ids1 = unique(dt[, get(id_col)]), ids2 = ids1,
                   freq_col = "spot_frequency_bin", log_fun = log2,
                   symmetric = TRUE,
                   type = c("kl", "cor", "Jensen_Shannon")){
  
  kl_ass = vapply(ids1, function(id1){
    vapply(ids2, function(id2){
      ind1 = dt[, get(id_col) %in% id1]
      ind2 = dt[, get(id_col) %in% id2]
      d1 = dt[ind1, get(freq_col)]
      d2 = dt[ind2, get(freq_col)]
      if(type == "kl") {
        
        LR = ifelse(d1 > 0, log_fun(d1) - log_fun(d2), 0)
        return(d1 %*% LR)
        
      } else if(type == "cor"){
        
        return(cor(d1, d2, method = "spearman"))
        
      } else if(type == " Jensen_Shannon") {
        # doesn't work
        m = 0.5 * (d1 + d2)
        JS = ifelse(d1 > 0 | m > 0 , 0.5 * (sum(d1 * log(d1 / m)) +
                                              sum(d2 * log(d2 / m))), 0)
        return(JS)
      }
      
    }, FUN.VALUE = numeric(1))
  }, FUN.VALUE = numeric(length(ids2)))
  
  if(symmetric & type == "kl"){
    if(length(ids1) != length(ids2)) stop("When symmetric kl length(id1) should be equal to length(id2)")
    kl_ass = (Matrix::tril(kl_ass) + Matrix::t(Matrix::triu(kl_ass))) / 2
    kl_ass = as.matrix(kl_ass + t(kl_ass))
  }
  rownames(kl_ass) = ids2
  colnames(kl_ass) = ids1
  kl_ass
}


plot_cell_prob = function(prob_mat, spread_of_profiles, sce, sce_slot = "counts",
                          min = 1.2, max = 1.3, n_cells = 10, 
                          limits_num_reads = c(6.5e5, 2.2e6), limits_prob = c(0, 1),
                          layer_levels = paste0("L", seq_along(ncol(prob_mat)))) {
  cells_ind = which(spread_of_profiles > min & spread_of_profiles < max)
  cells_ind = cells_ind[sample.int(length(cells_ind), n_cells)]
  cells = as.data.table(prob_mat[names(cells_ind), ], keep.rownames = "genes")
  
  c_p = plot_gradient(cells, gene_list = cells$genes, limits = limits_prob,
                      layer_levels = layer_levels) 
  
  sc_cov = data.table(num_reads = colSums(assay(sce[, c_p$gene_order], sce_slot)),
                      gene = factor(c_p$gene_order, levels = c_p$gene_order))
  
  num_reads = ggplot(sc_cov, aes(gene, "Marker gene\nreads",
                                 fill = num_reads, color = num_reads)) +
    geom_tile() +
    scale_fill_gradientn(colours=c_p$col, limits = limits_num_reads) +
    scale_color_gradientn(colours=c_p$col, limits = limits_num_reads) +
    scale_y_discrete(expand=c(0, 0)) +
    c_p$theme_change +
    guides(fill = guide_colorbar(title.position = "top")) +
    labs(fill = "Marker gene reads", color = "Marker gene reads") +
    theme(axis.text.x = element_blank(),
          legend.direction = "horizontal", legend.position = "top",
          legend.key.width = unit(30, "pt"))
  
  cowplot::plot_grid(num_reads,
                     c_p$plot + labs(fill = "probability", color = "probability") +
                       theme(legend.direction = "horizontal", legend.position = "bottom",
                             legend.key.width = unit(30, "pt")),
                     nrow = 2, 
                     rel_heights = c(1.2, 3), align = "v")
}

plot_cell_prob_and_mark = function(prob_mat, spread_of_profiles, sce, sce_slot = "counts",
                                   min = 1.2, max = 1.3, n_cells = 10, limits_prob = c(0, 1),
                                   layer_levels = paste0("L", seq_along(ncol(prob_mat))),
                                   markers = c("Chrdl1", "Gfap", "Il33"),
                                   normalise =  c(T, F, "log10"), flip_cells = F) {
  if(length(min) > 1){
    cells_ind = integer()
    for (i in 1:length(min)) {
      cells_ind_1 = which(spread_of_profiles > min[i] & spread_of_profiles < max[i])
      cells_ind_1 = cells_ind_1[sample.int(length(cells_ind_1), n_cells)]
      cells_ind = c(cells_ind, cells_ind_1)
    }
  } else {
    cells_ind = which(spread_of_profiles > min & spread_of_profiles < max)
    cells_ind = cells_ind[sample.int(length(cells_ind), n_cells)]
  }
  
  cells = as.data.table(prob_mat[names(cells_ind), ], keep.rownames = "genes")
  
  # generate nice order of cells
  # compute distances between profiles & do hierarchical clustering
  gene_order_hcl = character()
  for (j in (1:length(min) - 1)) {
    cells_ord = melt.data.table(cells[seq(j * n_cells + 1, j * n_cells + n_cells), ],
                                id.vars = "genes", variable.name = "layer", value.name = "expression")
    
    log_offset = min(as.matrix(cells, rownames = "genes")[as.matrix(cells, rownames = "genes") > 0])
    
    # add log offset
    cells_ord[, expression := expression + log_offset, by = .(genes)]
    # normalise within each gene
    cells_ord[, expression := expression / sum(expression), by = .(genes)]
    
    dist_mat = kl_dist(cells_ord, symmetric = T, type = "kl",
                       id_col = "genes", freq_col = "expression")
    dist = as.dist(dist_mat)
    hcl = hclust(dist)
    # create gene order by hierarchical clustering
    gene_order_hcl = c(gene_order_hcl, rownames(dist_mat)[hcl$order])
    
  }
  
  
  # plot probabilities
  c_p = plot_gradient(cells,gene_list = gene_order_hcl, gene_order = gene_order_hcl,
                      limits = limits_prob, layer_levels = layer_levels)
  
  # plot markers
  mark_mat = assay(sce[markers, c_p$gene_order], sce_slot)
  mark = as.data.table(as.matrix(t(mark_mat)), 
                       keep.rownames = "genes")
  c_mark = plot_gradient(mark, gene_list = gene_order_hcl, gene_order = gene_order_hcl,
                         #limits = c(min(mark_mat) - min(mark_mat) * 0.01, max(mark_mat) + max(mark_mat) * 0.01),
                         layer_levels = markers, normalise = normalise,
                         log_offset = 1)
  
  c_mark = c_mark$plot + labs(fill = "expression", color = "expression") + 
    theme(legend.direction = "horizontal", legend.position = "top",
          legend.key.width = unit(30, "pt"),
          axis.text.x = element_blank())
  
  if(flip_cells) {
    c_p = c_p$plot + coord_flip() + theme(axis.text.y = element_blank()) 
  } else {
    c_p = c_p$plot + theme(axis.text.x = element_blank()) 
  }
  
  c_p = c_p + labs(fill = "posterior probability", color = "posterior probability") + 
    theme(legend.direction = "horizontal", legend.position = "top",
          legend.key.width = unit(30, "pt"))
  
  
  cowplot::plot_grid(c_mark, c_p, nrow = 2, rel_heights = c(3, 2), align = "v")
}

# set colors and theme
ob5_col = colorRampPalette(RColorBrewer::brewer.pal(9, "YlOrRd"))(100)
ob5_theme_change = theme(
  strip.background = element_blank(),
  strip.text = element_blank(),
  plot.background = element_blank(),
  panel.grid.minor = element_blank(),
  panel.grid.major = element_blank(),
  panel.border =element_blank(),
  #panel.border = element_rect(colour = "black", fill=NA, size=0.5),
  axis.line = element_blank(),
  panel.background = element_blank(),
  legend.key.size = unit(1, "lines"),
  #legend.text=element_text(size=6),
  #legend.title=element_text(size=6),
  plot.title = element_text(hjust = 0.5, size=12),
  axis.title = element_blank(),
  #axis.text.y = element_text(size=7),
  axis.ticks = element_blank(),
  axis.text.x = element_text(angle = 90,
                             #size=7,
                             hjust=1,
                             vjust=0.6),
  #plot.margin = margin(0,-1,0,-1,unit="pt"),
  legend.position = "right",
  legend.direction = "vertical"
)