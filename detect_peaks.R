library(GenomicRanges)
library(igraph)
library(qtl)
library(pracma)
# quartz( )
# lod_col <- match(row$trait, colnames((single_scan,))-2
# plot(single_scan, lodcolumn=lod_col, type = "l")


get_peak_marker_lod <- function(peaks, single_scan) {
  apply(peaks, 1, function(r) {
    marker <- r[["marker"]]
    trait <- r[["trait"]] 
    single_scan[marker, trait] })
}



get_peak_table <- function(perms = NULL, single_scan = NULL){
  perm_sum = summary(single_scan, perms=perms, alpha=0.05,
                     format = "tabByChr", pvalues=TRUE, ci.function="lodint",
                     expandtomarkers = TRUE) 

  lod_sum <- summary(perms, alpha = 0.05)
  
  lod_thresh <- data.frame(
    trait = dimnames(lod_sum)[[2]],
    thresh = as.vector(lod_sum))

  lod_adj_sum <- summary(perms , controlAcrossCol=TRUE, alpha = 0.05)
  
  lod_thresh$thresh_adj <- as.vector(lod_adj_sum)
  
  # perm_sum <- perm_sum[-which(sapply(perm_sum, is.null))]  #remove empty elements

  
  peaks <- lapply(perm_sum, `[`)  %>%
    dplyr::bind_rows() %>%
    mutate(id = unlist(lapply(perm_sum, rownames))) %>%
    tidyr::separate(id, c("trait","pseudom"), sep = " : ") %>%
    left_join(lod_thresh) %>%
    mutate(trait = factor(trait))
  
  peaks$marker <- find.marker(cross,chr = peaks$chr, pos = peaks$pos)
  peaks$marker_lod <-  get_peak_marker_lod(peaks,single_scan)
  peaks$ci_left  <- find.marker(cross,chr = peaks$chr, pos = peaks$ci.low) 
  peaks$ci_left  <-  gsub(".*_","", peaks$ci_left, perl =TRUE)   %>% as.numeric()
  peaks$ci_right <- find.marker(cross,chr = peaks$chr, pos = peaks$ci.high)
  peaks$ci_right <-  gsub(".*_","", peaks$ci_right, perl =TRUE) %>% as.numeric()
#  peaks$width <- (peaks$ci_right - peaks$ci_left)/1000000
  peaks$width <- (peaks$ci_right - peaks$ci_left)
  peaks
}


refine_peaks<- function(peaks = NULL, single_scan =NULL, perms =NULL, cross = NULL){
  out_peaks <- peaks[0,]
  for(idx in 1:nrow(peaks)) {
    
    row <- peaks[idx,]
    trait_name <- as.character(row$trait)
    data <- as.data.frame(single_scan)[single_scan$chr == row$chr,c("pos",trait_name)]
    
    data$pos
    data[trait_name]
    intensity <- 10
    akima_lod <- akimaInterp(data$pos,data[,trait_name],
                             seq(from = 1/intensity, to = max(data$pos), by = 1/intensity))
    
    akima_lod[1] <- data[1,trait_name]
    akima_lod[length(akima_lod)] <- data[nrow(data),2]
    
    row
    subpeaks <- pracma::findpeaks(akima_lod, minpeakdistance = 5*intensity, 
                                  minpeakheight = row$thresh,
                                  npeaks = 5 )
    
    if(length(subpeaks)==0){
      trio <- floor(intensity*(c(row$pos, row$ci.low, row$ci.high)))
      subpeaks <- matrix(c(row$lod, trio), nrow = 1)
    }
    
    # get the positions inside each subpeak with a drop of at least 1.5 lod
    
    in_drop <- t(apply(subpeaks,1,function(x) {
      d <- round(which((x[1] - akima_lod) > 1.5) - x[2],0)
      y <-  c( x[2] + max(d[d < 0]), x[2] + min(d[d > 0]))
      if(xor(y[1] == -Inf,y[2] == Inf)) {
        inf_idx <- which(y == Inf | y == -Inf)
        y[[inf_idx]] <- c(1,length(akima_lod))[inf_idx]} 
      else if (y[1] == -Inf & y[2] == Inf){
        inf_idx <- which( y == -Inf)
        y <- c(1,length(akima_lod))} 
      y
    }))
    row
    subpeaks
    in_drop 
    # Consider multiple subpeaks with a lod drop of 1.5
    if(nrow(subpeaks) > 1 ) {
      edge_list <- as.matrix(IRanges::findOverlaps(IRanges(in_drop[,1], in_drop[,2])))
      
      g <- graph_from_edgelist(edge_list ,directed  = FALSE)
      # Iam asuumig small number of subpeaks
      # if subpeaks >10 this algorithm may never end
      
      ivs_list <- ivs(g, min=2)
      if(length(ivs_list) == 0){
        ivs_list <- ivs(g, min=1)
      }
      # Consider multiple solutions with the same vertex size
      vsize <- sapply(ivs_list, function(x) length(x)) 
      
      vs <- ivs_list[vsize == max(vsize)]
      
      is_max_lod <- sapply(ivs_list, function(x) max(subpeaks[x,1]) == max(subpeaks[,1]))
      
      # Select the solution containning the peak with max lod
      
      max_peak <- vs[which(is_max_lod)]    
    } else {
      max_peak <-list()
      max_peak[1] <- 1 
    }
    length(max_peak[1]) 
    length(max_peak[[1]]) 
    
    # Consider  multiple solutions with the same vertex size and including max lod
    # get the solution with the greatest lod sum.
    # Hideous!!!!
    
    if (length(max_peak) > 1 ){
      sum_lod <- sapply(max_peak, function(x) sum(subpeaks[x,1]))
      new_peaks <- cbind(subpeaks, in_drop)[ max_peak[[which(sum_lod == max(sum_lod))]],]
    }else if(length(max_peak)== 1 & length(max_peak[[1]]) > 1 ){
      new_peaks <- cbind(subpeaks, in_drop)[ max_peak[[1]],]
    } else if((length(max_peak) == 1) & (length(max_peak[[1]]) ==  1 )){
      # Remember this is an array of multiple vertices
      new_peaks <- matrix(cbind(subpeaks, in_drop)[max_peak[[1]],], nrow=1)
    }
    
    #finding the markers just outside
    chr_map <-  cross$geno[[row$chr]]$map
    
    bound_m <- t(apply(new_peaks,1, function(x) {
      left_pos <- max(chr_map[ chr_map  <= x[5]/intensity]) 
      right_pos <- min(chr_map[ chr_map >= x[6]/intensity]) 
      find.marker(cross, chr =row$chr, pos= c(left_pos,right_pos))
    }))
    
    out_row <- data.frame(
      chr = row$chr,
      pos = new_peaks[,2]/intensity,
      ci.low  = find.markerpos(cross, bound_m[,1])[,2],
      ci.high = find.markerpos(cross, bound_m[,2])[,2],
      lod = new_peaks[,1],
      pval =  sum(perms[,trait_name] > new_peaks[,1]) / length(perms[,trait_name]),
      trait = trait_name,
      pseudom = find.pseudomarker(cross, chr = row$chr, 
                                  pos = new_peaks[,2]/intensity,
                                  where = "prob"),
      thresh = row$thresh,
      thresh_adj = row$thresh_adj,
      marker = find.marker(cross, chr = row$chr, pos = new_peaks[,2]/intensity)
    )
    
    np_marker <- as.character(out_row$marker)
    np_trait <- as.character(row$trait)
    
    out_row$marker_lod <- single_scan[np_marker,np_trait]
    out_row$ci_left <- gsub(".*_","", bound_m[,1], perl =TRUE) %>% as.integer()
    out_row$ci_right <- gsub(".*_","", bound_m[,2], perl =TRUE) %>% as.integer()
    out_row$width <- (out_row$ci_right - out_row$ci_left) /1000000
    out_peaks <- rbind(out_peaks,out_row)
    
    # quartz()
    # plot(akima_lod, main = row$trait, type="l",
    #      xlab = paste("Chromosome", row$chr),
    #      ylab = "Akima Interpolated LOD")
    # points(subpeaks[, 2], akima_lod[subpeaks[, 2]], pch =2)
    # points(subpeaks[, 3], akima_lod[subpeaks[, 3]], pch = 0)
    # points(subpeaks[, 4], akima_lod[subpeaks[, 4]], pch = 0)
    # points(floor(intensity*out_row$ci.low), akima_lod[floor(intensity*out_row$ci.low)], pch = 15)
    # points(floor(intensity*out_row$ci.high), akima_lod[floor(intensity*out_row$ci.high)], pch = 15)
  }
  
  out_peaks <- out_peaks[out_peaks$lod > out_peaks$thresh,]
  out_peaks[out_peaks$lod > out_peaks$thresh_adj,]
  nrow(out_peaks)
  out_peaks
  
}

