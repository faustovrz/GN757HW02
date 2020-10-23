source("detect_peaks.R")
library(tidyr)
library(dplyr)
library(purrr)
library(qtl)
library(ggplot2)

# Load trait data

pheno <-  read.table("NAM-Project/NAM_populations/population_Z021/RIL_traitMatrix.txt", sep ="\t", header =TRUE)

trait_names <- pheno %>% 
  dplyr::select(!tidyselect::contains("_"), -SampleID) %>% 
  colnames()

# Calculate mean over environments, i.e. field/year
ncol(mean_pheno)
mean_pheno <- lapply(trait_names, function(trait) {
    pheno %>%
      dplyr::select(SampleID, tidyselect::starts_with(trait)) %>%
      tidyr::pivot_longer(
        cols = starts_with(trait),
        names_to = "env",
        names_prefix = trait,
        values_to = trait,
        values_drop_na = FALSE
      ) %>%
      dplyr::group_by(SampleID) %>%
      dplyr::summarise(!!trait :=  mean(!! rlang::sym(trait), na.rm = TRUE)) 
    } 
  ) %>% 
  purrr::reduce(left_join) %>%
  dplyr::select_if(~sum(!is.na(.)) > 0) %>%        #  select non empty traits
  dplyr::select_if(~table(.)[1]/length(.) < 0.9)   #  select traits with at least 10% different plants

# Visualize trait correlation structure



# Add id column and write csv for rqtl

colnames(mean_pheno)[1] <- "RIL"
write.csv(mean_pheno, "cross_pheno.csv", row.names = FALSE, quote = FALSE)

# Read genotype data

geno <-  read.table("NAM-Project/NAM_populations/population_Z021/NAM_genos_imputed.txt", sep ="\t", header =FALSE)
geno[1:5,1:5]
# Add id column, modify table and write csv for rqtl

chr <- as.integer(geno[1,])
cM <- geno[2,]
ids<- geno[3,]

geno[1,] <-ids
geno[2,] <- chr
geno[3,] <- cM
geno[2:3,1] <- ""
geno[1:5,1:5]


write.table(geno, "cross_geno.csv", col.names = FALSE, 
            row.names = FALSE, quote = FALSE, sep =",")


# QTL  calculations ------------------------------------------------------------


cross <- read.cross(
    "csvs", ".", 
    "cross_geno.csv", 
    "cross_pheno.csv", 
    genotypes = c("0.0","2.0"),
    na.strings = c("NaN","NA")
  ) %>%
  convert2riself() %>%
  jittermap(amount = 0.1) %>%
  fill.geno() %>% # for simplicity I fill the genotypes.
  calc.genoprob() %>%
  sim.geno()

pull.map(cross)
rownames(cross$pheno) <- cross$pheno$RIL

cross$pheno <- cross$pheno[,-1]

# Check missing genotypes
quartz()
geno.image(cross)


# Single marker analysis ------------------------------------------------------

n_perm <- 1000
n_pheno <-  nphe(cross)

singleqtl <- scanone(cross, pheno.col=1:n_pheno, method="em")

singleqtl.perm  <- scanone(cross, method="hk",  pheno.col=1:n_pheno, n.perm=n_perm, n.cluster = 10)


threshold <-summary(singleqtl.perm, alpha=0.05)

by_chr <- summary(singleqtl, threshold=3, format="tabByChr") 
class(by_chr)
with_peaks <- !grepl("NULL",lapply(by_chr, FUN = class) %>% unlist())

# Detect Peaks/ select initial marker set for MIM  -----------------------------

one_peak <- get_peak_table(singleqtl.perm, singleqtl) %>%
  dplyr::mutate(
    marker = pseudom,
    marker_lod = lod,
    ci_left  = ci.low*1e6,
    ci_right = ci.high*1e6,
    width = (ci.high-ci.low)*1e6
  )

one_peak %>% dplyr::group_by(trait) %>%
  dplyr::summarise(n = length(trait)) %>%
  dplyr::arrange(-n)


new_peaks <-refine_peaks(peaks = one_peak, single_scan =singleqtl, perms = singleqtl.perm, cross =cross) %>%
  dplyr::mutate(
    marker = pseudom,
    marker_lod = lod,
    ci_left  = ci.low*1e6,
    ci_right = ci.high*1e6,
    width = (ci.high-ci.low)*1e6
  )

new_peaks  %>% dplyr::group_by(trait) %>%
  dplyr::summarise(n = length(trait)) %>%
  dplyr::arrange(-n)

# Single marker QTL visulaization #############################################

# Convert QTLs to genomic ranges
gr <- GRanges(seqnames= new_peaks$chr, 
                IRanges(start=new_peaks$ci_left, end = new_peaks$ci_right),
                chr = new_peaks$chr,
                trait = new_peaks$trait,
                lod = new_peaks$lod,
                peak_pos = new_peaks$pos,
                peak_marker = new_peaks$marker,
                peak_marker_lod = new_peaks$marker_lod)

qtl_count <- 
  as.data.frame(gr) %>%
  dplyr::left_join(
    as.data.frame(gr) %>%
    dplyr::group_by(trait) %>%
    dplyr::summarise(n = length(trait)) %>%
    dplyr::arrange(n) %>% 
    dplyr::mutate(trait_n = factor(paste0(trait,"_",n), levels = factor(paste0(trait,"_",n))))
  ) %>% 
  dplyr::arrange(trait_n,lod) %>% as.data.frame()
 

qtl_gr <- makeGRangesFromDataFrame (qtl_count, seqnames = "chr")
mcols(qtl_gr) <-    qtl_count[,7:13]

write.csv(qtl_count %>% dplyr::arrange(-n,-lod) , row.names=FALSE,
          file = paste0("single_marker_intervals.csv"))


grl <- GenomicRanges::split(qtl_gr,qtl_gr$trait_n)

pdf(file = paste0("single_marker_intervals.pdf"),
    height = 7, width = 14)

# Single marker track plot 

p <-  ggbio::autoplot(grl,  group.selfish = TRUE,
                      aes(fill = lod, col = seqnames)) + 
  ggtitle("NAM NC358xB73 QTLs (1.5 LOD interval, n=186) Single Marker Analysis" ) +
  scale_x_continuous(labels=function(x)x/1000000) +
  scale_fill_viridis_c() +
  scale_color_manual(values = rep("black",10), guide = 'none') +
  theme(axis.text.x = element_text(angle = 90))

print(p)
dev.off()

# Composite interval mapping ##########################################################

# Reload the data and calculate sim.geno without filling genotypes for MIM
cross <- read.cross(
  "csvs", ".", 
  "cross_geno.csv", 
  "cross_pheno.csv", 
  genotypes = c("0.0","2.0"),
  na.strings = c("NaN","NA")
) %>%
  convert2riself() %>%
  jittermap(amount = 0.1) %>%
  #fill.geno() %>% #  for MIM I can't use this option
  calc.genoprob() %>%
  sim.geno()
rownames(cross$pheno) <- cross$pheno$RIL
cross$pheno <- cross$pheno[,-1]


# Initialize tables for  MIM QTL visualization and summaries

mim_qtl <- data.frame()
mim_fit <- list()
mim_rqtl <- list()
tidx <- 1


with_qtl <-  gsub("_.*","", levels(qtl_count$trait_n) ,perl =TRUE)
with_qtl

for(t in with_qtl){
nt <- which(colnames(cross$pheno) == t)
print(t)
print( paste(tidx, "out of", length(with_qtl)))

# Getting initial QTLs from single marker analysis

chr <- qtl_count %>% 
  dplyr::filter(trait ==t) %>% 
  dplyr::arrange(trait_n,-lod) %>%
  dplyr::pull("chr")
pos <- qtl_count %>% 
  dplyr::filter(trait ==t) %>% 
  dplyr::arrange(trait_n,-lod) %>%
  dplyr::pull("peak_pos")

qtl <- makeqtl(sim.geno(cross), chr, pos) # for some reason cross object was not enough

# Building an initial additive model

add_str <- paste(qtl$altname, collapse = " + ")
add_formula <- paste("y ~ ",add_str) %>% as.formula()

# Fit the additive model
fitqtl(cross, pheno.col=nt, qtl=qtl, formula= add_formula )

# Calculate the two marker scan for estimating lilelihood penalties.

perm2 <- scantwo(cross,pheno.col=t,method="hk", n.perm =n_perm)

# Calculate likelihood  penalties
lpen <- calc.penalties(perm2)

# Finaly stepwise selection including forward and  backward selection

stepout1 <- stepwiseqtl(cross, pheno.col=nt,
                        additive.only=TRUE, 
                        qtl = qtl,
                        penalties=lpen,
                        verbose=TRUE)


if(attr(stepout1,"pLOD") == 0){ 
  print("All QTLs dropped in multiqtl model!")
  next
}

# Initial additive multilocus model fit
# I need it in order to add interations for epistasis analysis

add_str <- paste(stepout1$altname, collapse = " + ")
add_formula <- paste("y ~ ",add_str) %>% as.formula()

additive <- fitqtl(cross, pheno.col=nt, qtl=stepout1, formula= add_formula )

# Check that I retained peaks after the multilocus analysis

if(!is.null(additive$result.drop)){

by_var <- additive$result.drop %>% as.data.frame() %>% arrange(-`%var`) %>% row.names()

qtl_step <- makeqtl(sim.geno(cross), 
               stepout1$chr[match(by_var,stepout1$name)], 
               stepout1$pos[match(by_var,stepout1$name)]
               )

} else {
  qtl_step <- stepout1
}

# Fit again the full additive model with the selected markers

full <- fitqtl(cross, pheno.col=nt, qtl=qtl_step, 
               formula = paste("y ~ ", paste(qtl_step$altname, collapse = " + ")) )

# QTL interaction analysis (Epistatic models) ---------------------------------

if(length(qtl_step$name)>1){
# Add interactions to the model and assess them
  intx <- addint(cross, pheno.col=nt, qtl=qtl_step , formula = add_formula)
  attributes(intx)
  is_significant <- intx$`Pvalue(F)` < 0.05
  
# Retain significant interaction
  if (any(is_significant)){
    intx_df <- data.frame("intx" = attr(intx,"row.names")[is_significant]) %>%
      separate(intx, c("qtl1", "qtl2"), ":") 
    sign_int <- apply(intx_df,1, function(x){
      paste0(
        qtl_step$altname[which(x[1] == qtl_step$name)], "*",
        qtl_step$altname[which(x[2] == qtl_step$name)]
      )
    })
    sign_int_str <- paste(sign_int, collapse  = " + ")
    
# Reevaluate model with significant interactions added
   full <- fitqtl(
     cross, 
     pheno.col=nt,
     qtl = qtl_step, 
     formula = as.formula(paste0("y ~ ", add_str, " + ", sign_int_str )))
   
   is_significant <- full$result.drop[,"Pvalue(F)"] < 0.05
   print( is_significant )
 }
}


# Reevaluate model with or without interactions depending what I found in previous steps

full.fit <- fitqtl(cross, 
       pheno.col=nt,
       qtl=qtl_step, formula = attr(full,"formula"))


mim_fit[[t]] <- full.fit

# Get the MIM scan for plotting !
rqtl <- refineqtl(cross, pheno.col=nt, qtl=qtl_step)
mim_rqtl[[t]] <- rqtl

# Get the peaks for visualization  with lodint
# My custom function refine_peaks frm the file detect_peaks.R fails with the MIM scans

mim_peaks <- 
  lapply(1:rqtl$n.qtl,function(x) { 
  lodpeak <- lodint(rqtl, qtl.index = x, drop=1.5)
  lastidx <- nrow(lodpeak)
  data.frame(
          IRanges(start=lodpeak$pos[1]*1e6 , end = lodpeak$pos[lastidx]*1e6),
          chr = lodpeak$chr[1],
          trait = t,
          lod = lodpeak$lod[2],
          peak_pos = lodpeak$pos[2]*1e6,
          peak_marker = row.names(lodpeak)[2],
          peak_marker_lod = lodpeak$lod[2]) 
  }) %>% dplyr::bind_rows()
mim_qtl <- rbind(mim_peaks, mim_qtl)
tidx <- tidx + 1
}

# Check what qTLs where discarded with MIM in comparison with the single marker analysis
left_over <- with_qtl[!with_qtl %in% names(mim_fit)]
# "SecondaryBranchNumber"

mim_gr <- GRanges(mim_qtl) 
mim_gr$chr <- seqnames(mim_gr )

mim_fit[["SouthernLeafBlight"]]

# MIM QTL visualization -------------------------------------------------------

# Build GenomicRanges for QTLs
mim_qtl_count <- 
  as.data.frame(mim_gr) %>%
  dplyr::left_join(
    as.data.frame(mim_gr) %>%
      dplyr::group_by(trait) %>%
      dplyr::summarise(n = length(trait)) %>%
      dplyr::arrange(n) %>% 
      dplyr::mutate(trait_n = factor(paste0(trait,"_",n), levels = factor(paste0(trait,"_",n))))
  ) %>% 
  dplyr::arrange(trait_n,lod) %>%
  tibble() %>%
  print(n=100)

mim_qtl_gr <- makeGRangesFromDataFrame (mim_qtl_count, seqnames = "chr")
mcols(mim_qtl_gr) <-    mim_qtl_count[,7:13]


write.csv(as.data.frame(mim_qtl_gr) %>% dplyr::arrange(-n,-lod) , row.names=FALSE,
          file = paste0("mim_intervals.csv"))


mim_grl <- GenomicRanges::split(mim_qtl_gr, mim_qtl_gr$trait_n)

# Build  single marker scan list for plot comparing Single Marker and MIM

sm_list <- list()

for(t in with_qtl){
  nt <- which(colnames(cross$pheno) == t)
  print(t)
  chr <- qtl_count %>% 
    dplyr::filter(trait ==t) %>% 
    dplyr::arrange(trait_n,-lod) %>%
    dplyr::pull("chr")
  pos <- qtl_count %>% 
    dplyr::filter(trait ==t) %>% 
    dplyr::arrange(trait_n,-lod) %>%
    dplyr::pull("peak_pos")
  
  qtl <- makeqtl(sim.geno(cross), chr, pos)
  
  add_str <- paste(qtl$altname, collapse = " + ")
  add_formula <- paste("y ~ ",add_str) %>% as.formula()
  sm_list[[t]] <- fitqtl(cross, pheno.col=nt, qtl=qtl, formula= add_formula )
}


# MIM track plot vizualization

pdf(file = paste0("mim.pdf"),
    height = 7, width = 14)

p <-  ggbio::autoplot(mim_grl,  group.selfish = TRUE,
                      aes(fill = lod, col = chr )) + 
  ggtitle("NAM NC358xB73 QTLs (1.5 LOD interval, n=186) Multiplie Interval Mapping" ) +
  scale_x_continuous(labels=function(x)x/1000000) +
  scale_fill_viridis_c() +
  scale_color_manual(values = rep("black",10), guide = 'none') +
  theme(axis.text.x = element_text(angle = 90))

print(p)
dev.off()

# Single marker model summary. 
# I can  compare SM with MIM in variance explained (R2) and Model LOD

sm_sum <- lapply(names(sm_list), function(x){
  myfit <- sm_list[[x]]
  data.frame(
    trait = x,
    interactions = 0,
    model_lod =  myfit$lod,
    R2  = myfit$result.full["Model","%var"]/100
  )
}
) %>% dplyr::bind_rows()

# This is the single Marker summary for the stats figures
SingleMarker <- cbind( data.frame(method ="SingleMarker"),
       as.data.frame(gr) %>%
         dplyr::group_by(trait) %>%
         dplyr::summarise(QTL_count = length(trait), max_LOD = max(lod), mean_LOD = mean(lod), mean_width = mean(width)/1e6) %>%
         dplyr::arrange(-QTL_count) 
) %>%
  dplyr::left_join(sm_sum, by ="trait") 


# Interactions  summary. 
# So I can see wether the more QTLs in the single marker analysis produce more
# Interactions in the MIM ananlysis

intx_sum <- lapply(names(mim_fit), function(x){
  myfit <- mim_fit[[x]]
  data.frame(
    trait = x,
    interactions = as.character(attr(mim_fit[[x]],"formula")) %>% 
      stringr::str_count( pattern = ":"),
    model_lod =  myfit$lod,
    R2  = myfit$result.full["Model","%var"]/100
  )
}
) %>% dplyr::bind_rows()


# This is the MIM summary for the stats figures
MIM <- cbind( data.frame(method ="MIM"),
  as.data.frame(mim_gr) %>%
  dplyr::group_by(trait) %>%
  dplyr::summarise(QTL_count = length(trait), max_LOD = max(lod),  mean_LOD = mean(lod), mean_width = mean(width)/1e6) %>%
  dplyr::arrange(-QTL_count) 
) %>%
  dplyr::left_join(intx_sum, by ="trait") 


# Single Marker vs MIM comparison ---------------------------------------------
df1 <- SingleMarker %>% 
  dplyr::left_join(MIM, by ="trait") %>%
  dplyr::mutate(QTL_count_diff = QTL_count.y - QTL_count.x)

df2 <- rbind(SingleMarker,MIM) %>%
  dplyr::filter(!trait == "SecondaryBranchNumber")



# Stat comparison plot stats1.png paired mean diferences

p1 <- ggpubr::ggarrange( 
  
  ggpubr::ggpaired(df2, x = "method", y = "QTL_count",
                   color = "method", line.color = "gray", line.size = 0.4,
                   palette = "jco") +
    ggplot2::xlab("") +
    ggplot2::ylab("QTL count per trait") +
    ggpubr::stat_compare_means(paired = TRUE) +
    ggplot2::theme(legend.position = "none"),
  
 
  ggpubr::ggpaired(df2, x = "method", y = "model_lod",
                   color = "method", line.color = "gray", line.size = 0.4,
                   palette = "jco") +
    ggplot2::xlab("") +
    ggplot2::ylab("model LOD") +
    ggpubr::stat_compare_means(paired = TRUE) +
    ggplot2::theme(legend.position = "none"),
 
  ggpubr::ggpaired(df2, x = "method", y = "R2",
                   color = "method", line.color = "gray", line.size = 0.4,
                   palette = "jco") +
    ggplot2::xlab("") +
    ggplot2::ylab("model R2") +
    ggpubr::stat_compare_means(paired = TRUE) +
    ggplot2::theme(legend.position = "none"),
  
  ggpubr::ggpaired(df2, x = "method", y = "mean_width",
                   color = "method", line.color = "gray", line.size = 0.4,
                   palette = "jco") +
    ggplot2::xlab("") +
    ggplot2::ylab("mean QTL width cM") +
    ggpubr::stat_compare_means(paired = TRUE) +
    ggplot2::theme(legend.position = "none"),
  
  
  nrow = 2,ncol=2,  align = "v"
)

# Stat comparison plot stats2.png  
#how more initial QTL in singlemarker analysis produces more 
# additive QTLs and interaction QTLs in MIM analysis

p2 <- 
ggpubr::ggarrange(   
  ggpubr::ggscatter(df1, x = "QTL_count.x", y = "QTL_count_diff",
                  add = "reg.line",                                 # Add regression line
                  conf.int = TRUE,                                  # Add confidence interval
                  add.params = list(color = "blue",
                                    fill = "lightgray"),
                  
) + 
  ggplot2::xlab("Single Marker QTL count")+
  ggplot2::ylab("MIM QTL count - SM QTL count") +
  ggpubr::stat_cor(method = "pearson", label.x = 1, label.y = 10),  # Add correlation coefficient

ggpubr::ggscatter(df1, x = "QTL_count.x", y = "interactions.y",
                        add = "reg.line",                                 # Add regression line
                        conf.int = TRUE,                                  # Add confidence interval
                        add.params = list(color = "blue",
                                          fill = "lightgray"),
                        
) + 
  ggplot2::xlab("Single Marker QTL count")+
  ggplot2::ylab("Pairwise QTL interactions count") +
  ggpubr::stat_cor(method = "pearson", label.x = 1, label.y = 8),  # Add correlation coefficient
nrow = 2, align = "v"
)

quartz()
ggpubr::ggarrange(p1,p2, ncol =2)

# Example Analysis of a single trait: GDDDaystoTassel -------------------------

singleqtl$GDDDaystoTassel

singleqtl$GDDDaystoTassel

# Get mim LOD profile
mim_dtascan <- attr(mim_rqtl[["GDDDaystoTassel"]],"lodprofile")
quartz()

dta.col <- which(colnames(cross$pheno) == "GDDDaystoTassel")

# Compare LOD profiles of SM and MIM analysis 

#  Using basegraphics to plot each profile by its own
plot(singleqtl, lodcolumn=dta.col)
abline(h=threshold[dta.col], lty =2, col ="turquoise")

plotLodProfile(mim_rqtl[["GDDDaystoTassel"]])
pull.map(cross)[[x]] 

# Get the whole gene map for merging the scans
map <-  lapply(1:10, FUN = function(x){
  data.frame(
    marker = names(pull.map(cross)[[x]]),
    chr = factor(x),
    pos = pull.map(cross)[[x]] %>% as.numeric()
  )
}) %>% dplyr::bind_rows() %>%as.data.frame()

# Merging LOD profiles for comparing them with ggplot

toplot  <- append(
  list( SM = map %>%
    dplyr::left_join(
      as.data.frame(singleqtl["GDDDaystoTassel"]) %>% 
        tibble::rownames_to_column("marker") %>%
        dplyr::rename(lod = "GDDDaystoTassel") %>%
        dplyr::mutate(method = "Single Marker", peak = NA)
    )) ,
  lapply( names(attr(mim_rqtl[["GDDDaystoTassel"]],"lodprofile")),
          function(x){
            qtl_peak <- attr(mim_rqtl[["GDDDaystoTassel"]],"lodprofile")[[x]]
            qtl_peak %>% tibble::rownames_to_column("marker") %>%
              dplyr::mutate(method = "MIM",peak = x)
          }),
) %>% dplyr::bind_rows() %>%as.data.frame()

# The ANOVA table corresponding to the full MIM model  and the drop one out comparison
dta.fit <- mim_fit[["GDDDaystoTassel"]]

write.csv(dta.fit$result.full,file = "dta_fit_full.csv")
write.csv(dta.fit$result.drop,file = "dta_fit_drop.csv")

# Adding information about interactions to each model peak label (like Tuckey groups)
int_label <- data.frame(
  interaction = rownames(dta.fit$result.drop)[grepl(":",rownames(dta.fit$result.drop))]
) %>%
  tidyr::separate(interaction, c("qtl1","qtl2"),":", remove = FALSE) %>%
  tidyr::pivot_longer(qtl1:qtl2, values_to ="peak") %>%
  dplyr::select(-name) %>%
  dplyr::mutate(int_id = letters[as.integer(factor(interaction))]) %>%
  dplyr::select(-interaction) %>%
  dplyr::arrange(int_id) %>%
  # dplyr::mutate(label = paste0(qtl,"*",int_id)) %>%
  tidyr::pivot_wider( names_from = int_id ,values_from = int_id) %>%
  tidyr::unite("label",a:f, sep ="", na.rm = TRUE)%>%
  dplyr::mutate(label = paste0(peak,"*",label))


# Check chromosomes to plot
chr_toplot <- toplot  %>% 
  dplyr::group_by(chr,method,peak) %>%
  dplyr::filter(method == "MIM" & !all(is.na(lod))) %>%
  dplyr::pull(chr) %>% unique()

# Adding interaction labels to peaks 
peak_toplot <- toplot  %>% 
  dplyr::group_by(chr,method,peak) %>%
  dplyr::filter(method == "MIM" & !all(is.na(lod))) %>%
  dplyr::select(chr,pos,peak,lod) %>%
  dplyr::filter(lod == max(lod)) %>%
  dplyr::arrange(chr,pos) %>%
  dplyr::left_join(int_label) %>%
  dplyr::mutate(label = coalesce(label, peak))

# Adding hung 2012 previous detected QTLs info to the plot

attributes(mim_fit$GDDDaystoTassel)
hung2012 <- read.csv("hung2012.tab")


hung2012_gr <- GRanges(
  data.frame(
    IRanges(start=hung2012$pos*1e6 , end = hung2012$pos*1e6),
    pos = hung2012$pos*1e6,
    chr = hung2012$chr,
    trait = "GDDDaystoTassel",
    R2 = hung2012$R2/100
  )
)

dta.gr <- mim_gr[mim_gr$trait =="GDDDaystoTassel"]

data.frame(SNP = dta.gr$peak_marker) %>%
  dplyr::inner_join(
    read.table(file = "map.tab", sep = "\t", header= TRUE)
  )

olap <- findOverlaps(dta.gr,hung2012_gr)

cbind(
  as.data.frame(dta.gr)[queryHits(olap),],
  as.data.frame(hung2012_gr)[subjectHits(olap),]
)

#  Profile comparison of SM and MIM  using ggplot: dta_scan.png

quartz()
toplot  %>% 
  dplyr::filter(chr %in% chr_toplot ) %>%
  ggplot2::ggplot() + 
  ggplot2::ggtitle("GDDDaystoTassel") +
  ggplot2::xlab("Chromosome") +
  ggplot2::ylab("LOD") +
  ggplot2::ylim(c(NA,15)) +
  ggplot2::geom_hline(aes(yintercept = threshold[dta.col]), col = "turquoise", lty =2) +
  ggplot2::geom_line(aes(x = pos, y = lod, group = peak, col = method), ) +
  ggplot2::facet_wrap(.~ as.factor(chr) , 
                      #scales = "free_x",
                      strip.position = "bottom",
                      ncol = length(chr_toplot)) + 
  ggplot2::geom_text(
    data = peak_toplot, 
    aes(x =pos, y = lod,label = label, angle = 90), 
    size =3, nudge_y = 1.5) +
  ggplot2::geom_text(
    data = data.frame(chr = c(8,9), pos = c(75.4,44.1),label = c("FT","vgt1")),
    aes(x =pos, y = c(6.5,10),label = label), size =5) +
  ggplot2::geom_rug(aes(x = pos),sides="b", size =0.01) +
  ggpubr::theme_classic2(base_size = 15) +
  ggplot2::theme( strip.background = ggplot2::element_blank(),
                  strip.placement = "outside",
                  axis.text.x = ggplot2::element_text(size =8, angle = 90),
                  legend.position = "top")

# Comparison of LOD profiles among traits--------------------------------------
# Are QTLs of correlated traits  correlated as well?

#  Compare heatmaps of trait correlation and LOD profiles

  # lapply( names(attr(mim_rqtl[[t]],"lodprofile")),
  #         function(x){
  #           qtl_peak <- attr(mim_rqtl[[t]],"lodprofile")[[x]]
  #           qtl_peak %>% tibble::rownames_to_column("marker") %>%
  #             dplyr::mutate(method = "MIM",peak = x)
  #         })

# Make  LOD profile matrix

m <- lapply(names(mim_rqtl), function(t){
  colnew <- "lod"
  names(colnew) <- t
  
  dplyr::left_join( map,
  lapply( names(attr(mim_rqtl[[t]],"lodprofile")),
          function(x){
            qtl_peak <- attr(mim_rqtl[[t]],"lodprofile")[[x]]
            qtl_peak %>% tibble::rownames_to_column("marker") %>%
              dplyr::mutate(peak = x)
          }) %>%
    dplyr::bind_rows() %>%
    dplyr::group_by(marker, chr, pos) %>%
    dplyr::filter(lod >0) %>%     #  just select positive values for ploting
    dplyr::summarise( lod = mean(lod, na.rm = TRUE)) %>%
    dplyr::rename(!!!colnew) %>% 
    dplyr::ungroup() 
  ) %>%
    dplyr::select(tail(names(.), 1))
  }
) %>% bind_cols() %>% t() 

m[is.na(m)] <-0 # for plotting empty values ploted as 0s

empty_col <- colSums(m, na.rm = TRUE)  == 0 # for skipping over empty columns

# Make trait correlation matrix

traitcor <- cor(cross$pheno[,rownames(m)], use ="pairwise.complete.obs") %>%abs()%>% as.matrix()


library(circlize)
library(ComplexHeatmap)

# Get  MIM QTL count for adding to LOD profile labels
peak_count <- mim_qtl_count %>% 
  dplyr::ungroup() %>%
  dplyr::group_by(trait,n) %>%
  dplyr::summarise(trait = trait[1],n =n[1]) %>%
  tibble::column_to_rownames("trait") %>%
  dplyr::select("n")

peak_count<- setNames(as.integer(peak_count$n), rownames(peak_count))

# Trait correlation hetmap and trait clustering

hm1 <- ComplexHeatmap::Heatmap(
  traitcor, 
  show_column_dend = FALSE,
  show_column_names = FALSE,
  column_title = "Trait Cluster Chromosome",
  column_title_side = "bottom",
  row_title_rot = 0,
  col = circlize::colorRamp2(c(0, 1), c("red", "yellow")),
  heatmap_legend_param = list(title = "Cor",
                              direction = "horizontal"),
  column_gap = unit(0.1, "mm"),
  row_split =  10, column_split =  10,
  width = unit(8, "cm")) 
ComplexHeatmap::draw(hm1)

order <- ComplexHeatmap::row_order(hm1) %>%unlist()


# Anntotation of MIM QTL QTL count 
ra <- ComplexHeatmap::rowAnnotation( 
        n = anno_text(peak_count[rownames(cm)],
                      gp = grid::gpar(fontsize = 8))
      )

# MIM LOD profiles matrix
hm2 <- ComplexHeatmap::Heatmap(m[,!empty_col],
                              cluster_columns = FALSE,
                              col = circlize::colorRamp2(c(0, 5), c("black", "green")),
                              column_split = map$chr[!empty_col],
                              column_title_side = "bottom",
                              heatmap_legend_param = list(title = "MIM LOD",
                                                          direction = "horizontal"),
                              width = unit(12, "cm"),
                              right_annotation = ra,
                              row_names_gp  = grid::gpar(fontsize = 8))


# Plot  heatmap comparison of trait correlation an LOD profile correlation

ComplexHeatmap::draw(hm1 + hm2,  # merge_legend = TRUE,
                      heatmap_legend_side = "top",
                      ht_gap = unit(0, "mm")
                      # align_heatmap_legend = "heatmap_top"
                     ) 



# Compare  Trait correlations vs QTL correlations ------------
trait_adj <- traitcor
diag(trait_adj) <- 0

lodcor <- cor(m[,!empty_col] %>% t()) %>% abs()
lod_adj <- lodcor
diag(lod_adj) <- 0

lod_trait_cor <-  dplyr::inner_join(
igraph::graph_from_adjacency_matrix(
  trait_adj,mode ="undirected", weighted = TRUE
) %>%
  igraph::as_data_frame() %>%
  dplyr::rename(traitcor = weight),

# Pairwise comparisons as a two column table
igraph::graph_from_adjacency_matrix(
  lod_adj,mode ="undirected", weighted = TRUE
  ) %>%
  igraph::as_data_frame()  %>%
  dplyr::rename(lodcor = weight)
) %>% 
  tibble::rowid_to_column("id")

order <- row_order(hm1)

# Pairwise comparisons corresponding to trait clusters
trait_cluster <- lapply(1:length(order), FUN = function(x) {
  v <- rownames(m[,!empty_col])[order[[x]]]
  df <- combn(v,2) %>% t() %>% as.data.frame()
  colnames(df) =c("from","to")
  df$cluster = x
  df
  }
) %>%
 dplyr::bind_rows()

elist <- split(lod_trait_cor[,1:3],f = as.factor(lod_trait_cor$id))

# Sorting the comparison labels (to and from ) so the information on trait clusters can be added

sorted_edge <- lapply(elist, function(x){
    out <- as.character(x[2:3]) %>% sort
    data.frame(from = out[1], to  = out[2])
  }
  ) %>% dplyr::bind_rows()

lod_trait_cor$from <-  sorted_edge$from
lod_trait_cor$to <-  sorted_edge$to


# Cluster labels according to the most common class of traits
# These I made manually:

cluster.label <- c("1 Flowering", 
                   "2 Leaf Width",
                   "3 Ear",
                   "4 Height",
                   "5 NIR",
                   "6 Misc",
                   "7 Ear Length",
                   "8 N Leaf Blight",
                   "9 Leaf Angle",
                   "10 Tassel" )

# Add trait cluster info to pairwise comparisons

toplot <- lod_trait_cor %>%
  dplyr::left_join(trait_cluster  %>%
    left_join(
      data.frame (cluster = 1:10, 
                  cluster.label =factor( cluster.label, levels = cluster.label)) 
      )
  ) %>%
  dplyr::mutate(cluster = factor(cluster, levels = 1:10)) %>%
  dplyr::mutate(dot.size = ifelse(is.na(cluster), 1,5))


# Marginal dot plot

pt2 <- toplot %>%
  # dplyr::filter(!is.na(cluster)) %>%
  ggplot2::ggplot() +
  ggplot2::geom_jitter(aes(y =lodcor, 
                          x = forcats::fct_reorder(cluster, lodcor),
                          color = cluster),
                          #size = dot.size),
                          width = 0.1)+ 
  ggplot2::xlab("Trait Cluster") +
  ggplot2::ylab("Correlation between MIM LOD profiles") +
  ggplot2::scale_color_discrete( name = "Trait Cluster",
                                 labels = cluster.label[c(-2,-9)]) + 
  ggplot2::scale_size(guide=FALSE) +
  ggpubr::theme_pubr(base_size = 10, legend = "none")

# Main Scatterplot of trait correlation vs LOD correlation last.png

pt1 <- toplot %>%
  ggplot2::ggplot() +
  ggplot2::geom_point(aes(y =lodcor, x = traitcor), col = "grey75") +
  ggplot2::geom_point(data = toplot %>% # 
                               dplyr::filter(!is.na(cluster)),
                      aes(y =lodcor, x = traitcor, color = cluster, size = 2)) +
  ggplot2::xlab("Correlation between trait values") +
  ggplot2::ylab("Correlation between MIM LOD profiles") +
  ggplot2::scale_size(guide=FALSE) +
  ggplot2::scale_color_discrete( name = "Trait Cluster",
                                 labels = cluster.label[c(-2,-9)]) + 
  ggpubr::theme_pubr(base_size = 15, legend = "top")

# ggpubr::ggarrange(
#   pt1 + theme(plot.margin = unit(c(1,0,0,0), "lines")),
#   pt2 + theme(axis.text.y = element_blank(),
#           #axis.ticks.y = element_blank(),
#           axis.title.y = element_blank(),
#           plot.margin = unit(c(0,0,0,0), "lines") ),
#    widths = c(1.5, 1),
#   common.legend = TRUE, align = "hv", nrow = 1 )

pr <- cowplot::insert_yaxis_grob(pt1, pt2, grid::unit(.2, "null"), position = "right")

cowplot::ggdraw(prpr)

# Failed QTL Overlap analysis -------------------------------------
# Not much to analyse here 
# overlap count or percentage as calcultaed here does not tell
# anything about QTL correlation between traits

# Calculate QTL range overlap from granges
# does not work wth lists and the pct overlap may be
from <- mim_gr
names(from) <- gsub("_\\d+","",names(from),perl =TRUE)
to <- from

mim_olap <- findOverlaps(from, to)
overlaps <- pintersect(from[queryHits(mim_olap)], to[subjectHits(mim_olap)])
percentOverlap <- width(overlaps) / width(to[subjectHits(mim_olap)])

overlaps$pct_olap <- percentOverlap 
overlaps$from <- as.data.frame(from[queryHits(mim_olap)])[,"trait"]
overlaps$to <- as.data.frame(from[subjectHits(mim_olap)])[,"trait"]

# Scatter plot, nothing to see here
lod_trait_cor %>%
  dplyr::left_join(trait_cluster)  %>%
  dplyr::mutate(cluster = factor(cluster, levels =as.character(1:10))) %>%
  dplyr::inner_join(
  overlaps %>% as.data.frame() %>%
    dplyr::filter(from != to) %>%
    dplyr::select(from,to,pct_olap)
) %>% dplyr::filter(!is.na(cluster)) %>%
  ggplot2::ggplot() +
  ggplot2::geom_point(aes(x = traitcor, y = pct_olap, color = cluster ))

# Multiple trait ------------------------------------------------------
# mqm <- mqmscan(cross, autocofactors, pheno.col = 1:nphe,  n.cluster=4)
# 
# mqmplot.cofactors(cross.full, autocofactors, justdots=TRUE)
# plot(mqmgetmodel(mqm))
# 
# require(snow)
# 
# results <- mqmpermutation(cross, scanfunction=mqmscan, cofactors=autocofactors,
#                                       n.cluster=2, n.perm=n_perm, batchsize=25)
# plot(mqm)
