
tmpFile <- function(){ paste0(paste0(stringi::stri_rand_strings(30, 1, '[A-Za-z0-9]'), collapse = ''), '.tmp') }

ppNum <- function (n)  format(n, big.mark = ",", scientific = FALSE, trim = TRUE)

make_square <- function(p, fudge=1) {
  dims <- heatmap_dims(p)
  p + ggplot2::theme(aspect.ratio = (dims$nrows/dims$ncols)*fudge)
}

heatmap_dims <- function(p) {
  .x <- as.character(p$mapping$x)
  .y <- as.character(p$mapping$y)
  
  .x <- .x[-grep('~', .x)]
  .y <- .y[-grep('~', .y)]
  
  ncols <- length(unique(p$data[[.x]]))
  nrows <- length(unique(p$data[[.y]]))
  return(list(ncols=ncols, nrows=nrows))
}




createEpiGenomicHeatMapData <- function(d, Rscript_path = '/home/opt/R-3.4.0/bin/Rscript', maxSites = 25000, outputDir = NULL){
  if(is.null(outputDir)){
    stop('Output dir not defined.')
  }
  
  if(! file.exists('epiCellTypes')) stop('"epiCellTypes" file missing.')
  
  write(c('sampleName,GTSP,patient', paste0(unique(d$patient), ',', unique(d$patient), ',', 'xxx', collapse = '\n')), file = 'samples')
  write('seqnames,strand,position,sampleName,refGenome', file = 'sites')
  
  d <- subset(d, seqnames %in% c(paste0('chr', 1:23), 'chrX', 'chrY'))
  
  if(nrow(d) > maxSites){
    n <- ceiling(maxSites / n_distinct(d$patient))
    set.seed(1)
    d <- bind_rows(lapply(split(d, d$patient), function(x){
      if(nrow(x) > n){
        x <- sample_n(x, n)
      }
      x
    }))
  }
  
  d$seqnames <- factor(d$seqnames, levels = unique(d$seqnames))
  
  write.table(dplyr::select(d, seqnames, strand, start, patient, refGenome), 
              file = 'sites', col.names = FALSE, row.names = FALSE, sep = ',', append = TRUE, quote = FALSE)
  
  
  comm <- paste(Rscript_path, '~/ext/epigeneticHeatmapMaker/epi_heatmap_from_file.R samples ', 
                '-c ~/ext/genomicHeatMapMaker/INSPIIRED.yml ',
                '-t epiCellTypes',
                '-o', outputDir,
                '-f sites ',
                '-r', d$refGenome[1])
  
  system(comm)
  system(paste0('mv samples sites ', outputDir, '/'))
}


createGenomicHeatMapData <- function(d, Rscript_path = '/home/opt/R-3.4.0/bin/Rscript', maxSites = 25000, outputDir = NULL){
  if(is.null(outputDir)){
    stop('Output dir not defined.')
  }
  
  write(c('sampleName,GTSP,patient', paste0(unique(d$patient), ',', unique(d$patient), ',', 'xxx', collapse = '\n')), file = 'samples')
  write('seqnames,strand,position,sampleName,refGenome', file = 'sites')
  
  d <- subset(d, seqnames %in% c(paste0('chr', 1:23), 'chrX', 'chrY'))
  
  if(nrow(d) > maxSites){
    n <- ceiling(maxSites / n_distinct(d$patient))
    set.seed(1)
    d <- bind_rows(lapply(split(d, d$patient), function(x){
           if(nrow(x) > n){
             x <- sample_n(x, n)
           }
           x
         }))
  }
  
  d$seqnames <- factor(d$seqnames, levels = unique(d$seqnames))
  
  write.table(dplyr::select(d, seqnames, strand, start, patient, refGenome), 
              file = 'sites', col.names = FALSE, row.names = FALSE, sep = ',', append = TRUE, quote = FALSE)
  
  
  comm <- paste(Rscript_path, '~/ext/genomicHeatMapMaker/genomic_heatmap_from_file.R samples ', 
                '-c ~/ext/genomicHeatMapMaker/INSPIIRED.yml ',
                '-o', outputDir,
                '-f sites ',
                '-r', d$refGenome[1])
  
  system(comm)
  system(paste0('mv samples sites ', outputDir, '/'))
}


createGenomicHeatMapPlot <- function(outputDir, sampleOrder = NULL, featuresToDrop = NULL, disableSigMarks = FALSE){
  if(! dir.exists(outputDir)){
    stop('Output dir missing or not defined.')
  }

  genomicHeatmap <- within(
    list(), {
      heatmap_sample_info <- read.csv(file.path(outputDir, 'samples'))
      gen_heatmap <- readRDS(paste0(outputDir, '/roc.res.rds'))
      
      heatmap_scale <- seq(0.2, 0.8, 0.1)
      gen_heatmap_colors <- colorspace::diverge_hsv(
        length(heatmap_scale), h = c(240, 0), v = 1, power = 1)

      select_gen_features <- row.names(gen_heatmap$ROC)
      if(! is.null(featuresToDrop)) select_gen_features <- select_gen_features[! select_gen_features %in% featuresToDrop]
      
      gen_heatmap$ROC <- gen_heatmap$ROC[select_gen_features,]
      
      heatmap_sample_levels <- unique(d$patient) # c("H19", "J60", "Linus", "M06", "M50", "M66")
      if(! is.null(sampleOrder)) heatmap_sample_levels <- sampleOrder
      
      heatmap_figure_labels <- heatmap_sample_levels
      
      stat_cuts <- c(0, 0.001, 0.01, 0.05, 1)
      gen_comp_stats <- structure(cut(
        gen_heatmap$pvalues$op[select_gen_features, 1],
        stat_cuts,
        labels = c("***", " **", " * ", "   "),
        include.lowest = TRUE),
        names = select_gen_features)
      gen_row_names <- paste0(names(gen_comp_stats), " - ", gen_comp_stats)
      
      plot_data <- gen_heatmap$ROC %>%
        reshape2::melt() %>%
        mutate(
          feat = Var1,
          comp.sym = gen_comp_stats[Var1],
          Var1 = paste0(Var1, " - ", comp.sym),
          Var1 = factor(Var1, levels = gen_row_names),
          Var2 = factor(Var2, levels = heatmap_sample_levels),
          grp = " ",
          sig = as.vector(gen_heatmap$pvalues$np[select_gen_features,]),
          sym = cut(
            sig, stat_cuts, labels = c("***", " **", " * ", "   "),
            include.lowest = TRUE))
      
      levels(plot_data$Var2) <- heatmap_figure_labels
      
      
      plot_data$Var1 <- gsub('\\*', '', as.character(plot_data$Var1))
      plot_data$Var1 <- gsub('\\-', '', as.character(plot_data$Var1))
      levels_var1= gsub('\\*', '',gen_row_names)
      levels_var1= gsub('\\-', '',levels_var1)
      plot_data$Var1=factor(plot_data$Var1,levels=levels_var1)
      
      
      if(disableSigMarks){
        plot_data$sym  <- ''
      }
      
      gen_plot <- ggplot(plot_data, aes(x = Var2, y = Var1, fill = value)) +
        geom_tile(color = 'black') +
        geom_text(aes(label = sym), color = "black", size = 3, nudge_y = -0.15) +
        scale_x_discrete(position = "top") +
        scale_fill_gradient2(
          breaks = c(0.2, 0.4, 0.6, 0.8),
          low = gen_heatmap_colors[1],
          mid = gen_heatmap_colors[round(length(heatmap_scale)/2)],
          high = gen_heatmap_colors[length(heatmap_scale)],
          midpoint = 0.5) +
        guides(fill = guide_colorbar(
          title.position = "left", title.hjust = 0.5,
          direction = "horizontal")) +
        labs(x = NULL, y = NULL, fill = "ROC\nScore") +
        theme(
          legend.position = "bottom",
          axis.text.y = element_text( angle = 0, hjust = 0, vjust = 0.5, size = 12),
          axis.text.x.top = element_text(
            angle = 90, hjust = 0.5, vjust = 0.5, size = 12),
          strip.placement = "outside",
          axis.line = element_blank(),
          axis.ticks.y = element_blank(),
          aspect.ratio = nrow(gen_heatmap$ROC)/ncol(gen_heatmap$ROC))
    })
  
  genomicHeatmap
}



