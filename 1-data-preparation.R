##-----------------------------------------------##
##  NOVEL MILLENNIUM - PREPARATION               ##
##  Ecological novelty in Queensland corals      ##
##  Emer Cunningham 2024                         ##
##-----------------------------------------------##

# load necessary packages
library(analogue)   #for join()
library(ciTools)    #for genus compositional change, add_ci()
library(FD)         #for trait data
library(picante)    #for trait data, comdist()
library(tidyverse)  #for various helpful functions

# load colour for cumulative novelty
cumul.col <- rgb(0.373, 0.651, 0.765)

# load novelty detection framework
identify.novel.gam <- function(site.sp.mat, alpha, metric, site,
                               plot = TRUE, plot.data = FALSE, 
                               reverse.mat = TRUE){
  
  require(vegan) # community disimilarity metrics
  require(mgcv) # additive modelling
  
  # check to see if site species matrix rows run oldest -> youngest
  # if so, rotate site-species matrix to get oldest sites first
  if(reverse.mat){
    if(as.numeric(as.character(rownames(site.sp.mat)[1])) <
       as.numeric(as.character(rownames(site.sp.mat)[dim(site.sp.mat)[1]]))){
      
      site.sp.mat <- site.sp.mat[dim(site.sp.mat)[1]:1, ]
      
    }
  }
  
  # calculate dissimilarity matrix, splitting those that can be
  # incorporated into vegdist from those that can't
  if(metric %in% c("chord", "SQchord")){
    
    require(analogue)
    
    site.dist <- as.matrix(distance(site.sp.mat,
                                    method=metric))
    
    if(metric=="chord"){site.dist <- site.dist / sqrt(2)}
    if(metric=="SQchord"){site.dist <- site.dist / 2}
    
  }
  
  if(metric == "hellinger"){
    
    hell.mat <- decostand(site.sp.mat, method="hellinger")
    site.dist <- as.matrix(dist(hell.mat))
    site.dist <- site.dist / sqrt(2)
    
  }
  
  if(!metric %in% c("chord", "SQchord", "hellinger")){
    
    site.dist <- as.matrix(vegdist(site.sp.mat,
                                   method=metric))
  }
  
  # option for trait-based dissimilarity (community MPD)
  # site.sp.mat = site.dist
  # ensure that both objects are in matrix format
  if(metric == "gower"){
    
    site.sp.mat <- as.matrix(site.sp.mat)
    site.dist <- as.matrix(site.sp.mat)
    
  }
  
  # rescale distance matrices that scale beyond 1
  if(metric %in% c("euclidean")){
    site.dist <- site.dist / max(site.dist)
  }
  
  # check to make sure there are some non-0 or non-1 values
  # so we can actually estimate novelty
  if(var(as.vector(site.dist)) == 0){
    return(NA)
  }
  
  # convert bin names to numbers to use as continuous variable
  bins <- as.numeric(rownames(site.sp.mat))
  
  # This sets the maximum knot number for the additive model. One of the biggest
  # problems with additive models is overfitting. Ideally we want penalized
  # splines, where knot number is weighted against variance explained. This will
  # provide a good measure of local mean trends in dissimilarity. But if we want
  # this method to run on smaller time-series, setting penalized splines doesn't work
  # very well. So I've set this to set maximum knots to half the number of bins, rounded
  # down, for time-series with between 4 and 10 bins. Less than 4 will not run using
  # splines, and is set to run as a linear regression below.
  set.k <- ifelse(dim(site.sp.mat)[1] > 10, #if >10
                  -1,
                  ifelse(dim(site.sp.mat)[1] > 4, #for those <10,
                         floor(dim(site.sp.mat)[1])/2, #if >4
                         NA)) #if <4
  
  # SEQUENTIAL DISIMILARITY ####
  
  # calculate differenced sequential dissimilarities between
  # time t and t-1
  seq.dist <- c(NA, diag(site.dist[-1,-dim(site.dist)[2]]))
  
  # transform to remove 0s and 1s for beta regression
  seq.dist.tr <- (seq.dist * (length(seq.dist)-1) + 0.5) / length(seq.dist)
  
  # model localised trend over time and extract residuals to get dissimilarity
  # compared to local mean.
  if(var(seq.dist, na.rm=TRUE)==0){return(NULL)}
  
  if(!is.na(set.k)){
    seq.gam <- gam(seq.dist.tr ~ s(bins, bs="cr", k= set.k),
                   family=betar(),
                   method="REML")
  } else{
    seq.gam <- gam(seq.dist.tr ~ bins, family=betar(), method="REML")
  }
  
  # MINIMUM DISIMILARITY ####
  
  # calculate minimum dissimilarity from time 1 to time t (time 1 being
  # earliest time point)
  min.dist <- sapply(1:dim(site.dist)[1], function(n){
    if(n==1){return(NA)}
    min(site.dist[n,1:n][-n])
  })
  
  # transform to remove 0s and 1s for beta regression
  min.dist.tr <- (min.dist * (length(min.dist)-1) + 0.5) / length(min.dist)
  
  if(var(min.dist, na.rm=TRUE)==0){return(NULL)}
  # model localised trend over time
  if(!is.na(set.k)){
    min.gam <- gam(min.dist.tr ~ s(bins, bs="cr", k= set.k),
                   family=betar(),
                   method="REML")
  } else{
    min.gam <- gam(min.dist.tr ~ bins, family=betar(), method="REML")
  }
  
  # COMPARING OBSERVED TO EXPECTED SEQUENTIAL DISIMILARITY ####
  
  # This process calculates the p-value of the observed disimilarity score
  # being part of the expected distribution at the point in the time-series.
  
  # convert mu & phi beta parameters to A & B to use in qbeta.
  # mu = the mean predicted dissimilarity at each time point based on the
  #      additive model.
  # phi = a dispersion parameter, which is set globally across the model.
  #       phi is unfortunately hidden deep in the gam model as a text string.
  
  seq.mu <- c(NA, seq.gam$fitted.values)
  phi <- as.numeric(substr(seq.gam$family$family,
                           regexpr("\\(", seq.gam$family$family)+1,
                           nchar(seq.gam$family$family)-1))
  
  # shape parameters used in qbeta.
  A = seq.mu * phi
  B = phi - A
  
  # predict 5% and 95% prediction intervals from beta distribution parameters.
  # We use 95% rather than 97.5% because this is a one-tailed probability test.
  # We are not concerned about communities that are MORE similar than predictions.
  # This is done for each bin along the time-series.
  seq.p <- do.call("rbind", lapply(1:length(A), function(n){
    data.frame(lwr=qbeta(0.05, shape1=A[n], shape2=B[n]),
               upr=qbeta(0.95, shape1=A[n], shape2=B[n]),
               seq.p = pbeta(seq.dist[n], shape1=A[n], shape2=B[n], lower.tail=FALSE))
  }))
  
  # EXPECTED OBSERVED TO MINIMUM DISIMILARITY ####
  
  # convert mu & phi beta parameters to A & B to use in qbeta
  min.mu <- c(NA, min.gam$fitted.values)
  phi <- as.numeric(substr(min.gam$family$family,
                           regexpr("\\(", min.gam$family$family)+1,
                           nchar(min.gam$family$family)-1))
  A = min.mu * phi
  B = (1 - min.mu) * phi
  
  # predict 5% and 95% prediction intervals from beta distribution parameters
  min.p <- do.call("rbind", lapply(1:length(A), function(n){
    data.frame(lwr=qbeta(0.05, shape1=A[n], shape2=B[n]),
               upr=qbeta(0.95, shape1=A[n], shape2=B[n]),
               min.p = pbeta(min.dist[n], shape1=A[n], shape2=B[n], lower.tail=FALSE))
  }))
  
  # CREATE RETURN DATA-FRAME ####
  return.data <- data.frame(site = site,
                            bins = rownames(site.sp.mat),
                            bin.lag = c(NA, abs(diff(as.numeric(rownames(site.sp.mat))))),
                            seq.dist = seq.dist,
                            seq.resid = c(NA, resid(seq.gam)),
                            raw.min.dist = min.dist,
                            min.resid = c(NA, resid(min.gam)),
                            seq.p = seq.p$seq.p,
                            min.p = min.p$min.p,
                            instant = seq.p$seq.p <= alpha,
                            cumul = min.p$min.p <= alpha,
                            novel = seq.p$seq.p <= alpha & min.p$min.p <= alpha,
                            seq.exp=seq.mu,
                            min.exp=min.mu)
  
  return.data$cat = "back"
  return.data[return.data$instant &
                !is.na(return.data$instant) &
                !is.na(return.data$cumul) &
                !return.data$cumul, "cat"] = "instant"
  return.data[!return.data$instant &
                !is.na(return.data$instant) &
                !is.na(return.data$cumul) &
                return.data$cumul, "cat"] = "cumul"
  return.data[return.data$instant &
                !is.na(return.data$cumul) &
                !is.na(return.data$instant) &
                return.data$cumul, "cat"] = "novel"
  
  return.data$cat.bef <- c(NA, return.data$cat[-nrow(return.data)])
  
  # PLOT TIME SERIES ####
  
  # this section of code generates the plot seen when running the function.
  if(plot){
    par(mfrow=c(2,1), mar=c(0,0,0,0), oma=c(3,6,0.5,9), las=1)
    
    plot(seq.dist ~ bins, type="n",
         ylim=c(max(seq.dist, na.rm=TRUE)+0.1,
                min(seq.dist, na.rm=TRUE)),
         axes=FALSE)
    
    polygon(x=c(bins, rev(bins)),
            y=c(seq.p[,1], rev(seq.p[,2])),
            col="grey80", border=NA)
    lines(plogis(predict(seq.gam)) ~ bins[-1], col="grey50", lwd=2)
    lines(seq.dist ~ bins)
    points(y=seq.dist[seq.p$seq.p< 0.05],
           x=bins[seq.p$seq.p < 0.05], pch=21, bg="red")
    lims <- par("usr")
    
    sapply(which(return.data$novel), function(x){
      segments(x0 = bins[x],
               x1 = bins[x],
               y0 = seq.dist[x] + (0.05 * (par("usr")[3] - par("usr")[4])),
               y1 = par("usr")[3], col="orange", lwd=1)
    })
    
    segments(x0=par("usr")[1], x1=par("usr")[1], y0=par("usr")[3], y1=par("usr")[4])
    segments(x0=par("usr")[1], x1=par("usr")[2], y0=par("usr")[4], y1=par("usr")[4])
    segments(x0=par("usr")[2], x1=par("usr")[2], y0=par("usr")[3], y1=par("usr")[4])
    axis(side=2)
    mtext(side=2, text = "Sequential\ndissimilarity", line=3.5, las=0)
    
    
    plot(min.dist ~ bins, type="n", axes=FALSE,
         ylim=c(min(min.dist, na.rm=TRUE),
                max(min.dist, na.rm=TRUE)+0.1))
    polygon(x=c(bins, rev(bins)),
            y=c(min.p[,1], rev(min.p[,2])),
            col="grey80", border=NA)
    lines(plogis(predict(min.gam)) ~ bins[-1], col="grey50", lwd=2)
    lines(min.dist ~ bins)
    
    points(y=min.dist[min.p$min.p< 0.05],
           x=bins[min.p$min.p < 0.05], pch=21, bg="skyblue")
    
    sapply(which(return.data$novel), function(x){
      segments(x0 = bins[x],
               x1 = bins[x],
               y0 = min.dist[x] + (0.05 * (par("usr")[4] - par("usr")[3])),
               y1 =par("usr")[4], col="orange", lwd=1)
    })
    
    par(xpd=NA)
    points(y=rep(par("usr")[4], sum(min.p$min.p< 0.05 & seq.p$seq.p < 0.05, na.rm=TRUE)),
           x=bins[min.p$min.p< 0.05 & seq.p$seq.p < 0.05 & !is.na(min.p$min.p)],
           pch=21, bg="orange")
    par(xpd=TRUE)
    
    segments(x0=par("usr")[1], x1=par("usr")[1], y0=par("usr")[3], y1=par("usr")[4])
    segments(x0=par("usr")[1], x1=par("usr")[2], y0=par("usr")[3], y1=par("usr")[3])
    segments(x0=par("usr")[2], x1=par("usr")[2], y0=par("usr")[3], y1=par("usr")[4])
    
    axis(side=1)
    mtext(side=1, "Time series age", line=2)
    
    axis(side=2)
    mtext(side=2, text = "Minimum\ndissimilarity", line=3.5, las=0)
    
    par(xpd=NA)
    legend(x=par("usr")[2], y= par("usr")[4],
           legend=c("", "", ""),
           pch=c(NA,NA,NA), lwd=c(NA,NA,1),
           pt.bg=c("red","skyblue","orange"), col=c("black","black","orange"),
           yjust=0.5, xjust=0, bty="n", seg.len=1, x.intersp=0.5)
    
    legend(x=par("usr")[2], y= par("usr")[4],
           legend=c("Instantaneous", "Cumulative", "Novel"),
           pch=c(21,21,21), lwd=c(NA,NA,NA),
           pt.bg=c("red","skyblue","orange"), col=c("black","black","black"),
           yjust=0.5, xjust=0, bty="n", seg.len=1, x.intersp=0.5)
    par(xpd=FALSE)
  }
  
  if(plot.data){return(list(return.data,
                            cbind(seq.p, c(NA, plogis(predict(seq.gam)))),
                            cbind(min.p, c(NA, plogis(predict(min.gam))))))
  } else {return(return.data)}
}

# load the taxonomic community classification matrices
# sections 1. and 2. below show code used to create this object;
# please see README.md for information on raw data access
sorted_matrices <- readRDS("output/sorted-matrices.RDS")

#

#### 1. wrangling the coral community data ####

# # import the long-form core data
# dated_community_data <- read.csv("data/dated-community-data.csv")
# 
# # reduce these data
# core_dat_all <- dated_community_data[,c("coreID", "taxon", "raw.abund", "region",
#                                         "reef", "site", "core", "age.mean")]
# 
# # update some genus names (given "worms" package, 
# # and insider knowledge i.e., advice from Dr Hannah Markham Summers)
# core_dat_all$taxon <- str_replace_all(core_dat_all$taxon, 
#                                       c("Carophyllia" = "Caryophyllia",
#                                         "Caulastrea" = "Caulastraea",
#                                         "Favia" = "Dipsastraea"))
# 
# # make a list containing data for all individual cores
# core_dat <- split(core_dat_all, f = core_dat_all$coreID)
# 
# # there are some duplicate entries of taxa at the same time point
# # but! they are not strict duplicates - their abundances differ
# # so find the mean raw abundance for these pseudo-duplicates
# # and remove coral abundance observations with missing ages
# 
# # fix double-entries for Halfway 3.2 and Pandora 1H
# core_dat <- lapply(core_dat, function(core){
#   
#   # create column with enough detail to extract pseudo-duplicates
#   core$dupes <- str_c(core$age.mean, core$taxon, sep = "_")
#   
#   # return their mean raw abundance
#   new_abunds <- sapply(split(core$raw.abund, f = core$dupes), mean)
#   
#   # remove the duplicates and match up the new mean raw abundances with core data
#   new_coredat <- core[!duplicated(core$dupes),]
#   new_coredat$raw.abund <- new_abunds[match(new_coredat$dupes, names(new_abunds))]
#   
#   # finally get rid of the pseudo-duplicates-identifier
#   return(new_coredat[,-9])
#   
# })
# 
# # remove ageless observations for Halfway 3.1, 3.2, North Keppel 1.3
# core_dat <- lapply(core_dat, function(core){
#   
#   core <- core[!is.na(core$age.mean),]
#   return(core)
#   
# })
# 
# ### create raw abundance matrices
# 
# # make a list of raw coral abundances (by mass) over time
# raw_matrices_all <- lapply(core_dat, function(core){
#   
#   # the "id_cols" argument in pivot_wider stopped working after R updated...
#   # so here's a temporary manual fix - give pivot_wider only the columns of interest
#   core <- core[,c("taxon", "raw.abund", "age.mean")]
#   
#   # create a matrix of raw abundance (mass) of genera over time
#   core_mass <- pivot_wider(core, #id_cols = c(taxon, raw.abund, age.mean),
#                            names_from = taxon, values_from = raw.abund, values_fill = 0)
#   
#   # order and name by time
#   core_mass <- core_mass[order(core_mass$age.mean, decreasing = TRUE),]
#   core_mass <- column_to_rownames(core_mass, var = "age.mean")
#   
# })
# 
# # we now have all 86 cores in the correct format
# # but some will have to be cut, to leave suitable cores for analysis
# # so let's remove unwanted time series and taxa based on selection criteria
# 
# # trim the >1500 ybp dates from Big Woody 1.8 (because younger is usable)
# bw18 = raw_matrices_all$Hervey_Bay.Big_Woody.1.8
# 
# bw18 <- bw18[which(as.numeric(row.names(bw18)) < 1500), ] #trim dates
# bw18 <- bw18[ ,-which(colSums(bw18) == 0)] #remove empty Pocillopora
# 
# raw_matrices_all$Hervey_Bay.Big_Woody.1.8 = bw18
# 
# # remove cores that stretch beyond 1500 ybp
# core_max <- sapply(raw_matrices_all, function(core){
#   
#   max(as.numeric(row.names(core)))
#   
# })
# raw_matrices <- raw_matrices_all[-which(core_max > 1500)]
# 
# # and remove Middle Island 3.1 - distal reef flat core (not reef slope)
# raw_matrices$Keppel_Islands.Middle_Island.3.1 <- NULL
# 
# # remove columns with unwanted or undefined genera
# genus.remove <- function(core, genus){
#   
#   strTF <- str_detect(names(core), pattern = genus)
#   
#   if(length(unique(strTF)) == 2){
#     rm.column <- which(strTF == "TRUE")
#     core <- core[,-rm.column]
#   }
#   
#   return(core)
#   
# }
# 
# for(i in c("Mollusc", "Corallinales", "Lobophora", "Unknown")){
#   
#   raw_matrices <- lapply(raw_matrices, genus.remove, genus = i)
#   
# }

#### 2. taxonomic community classification ####

### time-sort the raw abundance matrices

# bin_width = 20
# 
# # check the limits of our time period across all Qld cores
# core_ages <- range(sapply(raw_matrices, function(x){as.numeric(rownames(x))}))
# 
# # create a sequence of our time bins that encompasses our time period
# binseq_qld <- seq(((core_ages[1] %/% bin_width) - 1) * bin_width,
#                   ((core_ages[2] %/% bin_width) + 1) * bin_width,
#                   bin_width)
# 
# # function for time-sorting the raw abundance matrices
# time.bin.sort <- function(core, binseq){
#   
#   # cut the core years into bins
#   age_char <- as.character(cut(as.numeric(rownames(core)), binseq))
#   
#   # calculate mean of bin range, string numbers
#   bin_centre <- rowMeans(cbind(lower = as.numeric(sub("\\((.+),.*", "\\1", age_char)),
#                                upper = as.numeric(sub("[^,]*,([^]]*)\\]", "\\1", age_char))))
#   
#   # make a matrix with coral genus relabundance for each bin year
#   # ifelse does time aggregation for bins with >1 sample
#   bin_mat <- sapply(split(core, f = bin_centre), function(x){
#     if (!is.null(ncol(x))) {return(colSums(x))} 
#     else {return(x)}
#   })
#   
#   # sort so we see observations of genera forward through time
#   bin_mat <- t(bin_mat)
#   bin_mat <- bin_mat[order(as.numeric(row.names(bin_mat)), decreasing = TRUE),]
#   
#   return(bin_mat)
#   
# }
# 
# # time sort the raw abundance matrices
# abund_matrices <- lapply(raw_matrices, time.bin.sort, binseq = binseq_qld)
# 
# # remove empty rows (with no coral abundance data)
# # in High 1.4, Big Woody 1.4, Four Mile 1.3, 1.4, Pt Vernon 1.1
# abund_matrices <- lapply(abund_matrices, function(core){
#   
#   missing <- which(rowSums(core) == 0)
#   if(length(missing) > 0){
#     core <- core[-c(missing),]
#   }
#   
#   return(core)
#   
# })
# 
# # remove cores with a single genus through time or a single time bin
# abund_matrices$Palm_Islands.Pandora.2.F <- NULL
# abund_matrices$Palm_Islands.Pandora.2.G <- NULL
# abund_matrices$Hervey_Bay.Pt_Vernon_West_Reef.1.1 <- NULL
# 
# # create relative abundance matrices
# sorted_matrices <- lapply(abund_matrices, function(raw_abund){
#   
#   rel_abund <- raw_abund/rowSums(raw_abund)
#   
#   return(rel_abund)
#   
# })
# 
# # select our final taxonomic compositional matrices
# # keeping in mind that the novelty detection framework requires ~5 obs.
# 
# # check how many time bins each core has
# sort(sapply(sorted_matrices, function(x){dim(x)[1]}))
# 
# # remove cores with <5 aggregated time points
# sorted_matrices <- sapply(sorted_matrices, function(core){
#   
#   if(dim(core)[1] > 4){return(core)}
#   else {core <- NA}
#   
# })
# sorted_matrices <- subset(sorted_matrices, !is.na(sorted_matrices))
# 
# # we now have 55 taxonomic matrices we can analyse!
# names(sorted_matrices)
# 
# # save for online publication and use:
# saveRDS(sorted_matrices, file = "output/sorted-matrices.RDS")

#### 3. wrangling the coral trait data ####

# import necessary datasets
species_list <- read.csv("data/species-list.csv")
trait_data_ctdb <- read.csv("data/trait-data-ctdb.csv")
trait_data_research <- read.csv("data/trait-data-research.csv")

# subset to keep species with known Queensland distributions
ctdb_species <- match(trait_data_ctdb$specie_name, #species with data
                      species_list$species) #species along Queensland

trait_data_ctdb <- trait_data_ctdb[-which(is.na(ctdb_species)),]

# create a column for genus
trait_data_ctdb$genus <- word(trait_data_ctdb$specie_name, 1)

# update some genus names (given "worms" package, and insider knowledge)
for(i in c("Danafungia", "Lithophyllon", "Lobactis", "Pleuractis")){
  
  trait_data_ctdb$genus <- str_replace(trait_data_ctdb$genus, i, "Fungia")
  
}

for(i in c("Astrea", "Paramontastraea")){
  
  trait_data_ctdb$genus <- str_replace(trait_data_ctdb$genus, i, "Montastraea")
  
}

for(i in c("Homophyllia", "Parascolymia")){
  
  trait_data_ctdb$genus <- str_replace(trait_data_ctdb$genus, i, "Scolymia")
  
}

# recreate the trait_data_ctdb but just with the useful columns
trait_dat <- data.frame(genus = trait_data_ctdb$genus,
                        species = trait_data_ctdb$specie_name,
                        trait = trait_data_ctdb$trait_name,
                        value = trait_data_ctdb$value,
                        units = trait_data_ctdb$standard_unit,
                        value_type = trait_data_ctdb$value_type,
                        method = trait_data_ctdb$methodology_name)

# add trait data from the species you had to research
trait_dat <- rbind(trait_dat, trait_data_research)

# keep only the desired traits
trait_vec <- c("Sexual system","Oocyte size at maturity",
              "Coloniality","Zooxanthellate","Growth form typical",
              "Life history strategy","Mode of larval development",
              "Colony maximum diameter","Growth rate","Calcification rate",
              "Tissue thickness","Skeletal density","Colony fecundity",
              "Corallite width maximum", "Corallite width minimum")

trait_dat <- trait_dat[!is.na(match(trait_dat$trait, trait_vec)),]

# how's she looking (she = representation of observations per trait)
sort(table(trait_dat$trait))

### sort the data into a workable list of trait data frames
### and clean up those with multiple units/annoying duplicates

# make a list of data per trait
trait_list <- split(trait_dat, f = as.factor(trait_dat$trait))

# are traits measured using multiple units?
lapply(trait_list, function(trait){
  
  return(unique(trait$units))
  
})

# select standard units for Calcification rate
table(trait_list[[1]]$units, trait_list[[1]]$genus)
rm1 <- which(trait_list[[1]]$units != "g cm ^-2 yr^-1")
trait_list[[1]] <- trait_list[[1]][-rm1,]

# remove duplicate rows for Coloniality
rm2 <- which(duplicated(trait_list[[2]]) == TRUE)
trait_list[[2]] <- trait_list[[2]][-rm2,]

# fix unit measurements for Colony maximum diameter
rm4 <- which(trait_list[[4]]$units == "m")
trait_list[[4]][rm4, "units"] <- "cm"
trait_list[[4]][rm4, "value"] <- 200

# treat measurements as numeric for Corallite Width Max & Min
trait_list[[5]]$value <- as.numeric(trait_list[[5]]$value)
trait_list[[6]]$value <- as.numeric(trait_list[[6]]$value)

# select standard units for Growth rate
table(trait_list[[8]]$genus, trait_list[[8]]$units)
rm8 <- which(trait_list[[8]]$units != "mm yr^-1")
trait_list[[8]] <- trait_list[[8]][-rm8,]

### calculate species-level trait means

# do any traits have more than one observation per species?
lapply(trait_list, function(x){
  
  obs <- length(x[,1])
  spp <- length(unique(x$species))
  if(obs==spp) {return(FALSE)} 
  else {return(TRUE)}
  
})

# there are some traits with multiple entries for same species
# but they're all continuous traits, so we can merge by mean-ing
spec.means <- function(ctdb_dat){
  
  obsv <- length(ctdb_dat[,1])
  spec <- length(unique(ctdb_dat$species))
  
  # return trait values as normal for species with one observation
  if (obsv==spec) {
    return(data.frame(genus=ctdb_dat$genus,
                      species=ctdb_dat$species,
                      value=ctdb_dat$value))
  }
  
  # or, return mean trait values for species with multiple observations
  else {
    new_vals <- sapply(split(ctdb_dat, f = as.factor(ctdb_dat$species)), function(x){
      # find the species who have more than one observation
      # calculate and return the mean for >1 entry
      if(!is.null(nrow(x))){
        meanval <- mean(as.numeric(x$value))
        return(meanval)
      } else {return(as.numeric(x$value))} # or give the single value
    })
    
    # put these new values into a nice data frame
    return(data.frame(genus = ctdb_dat$genus[match(names(new_vals), ctdb_dat$species)],
                      species = names(new_vals),
                      value = new_vals,
                      row.names = NULL))
    
  }
}

# find the species-level means for each functional trait
spec_dat <- lapply(trait_list, spec.means)

# calculate genus-level trait means

# write a function to calculate genus-level means
# yielding a single column for continuous traits (mean trait value),
# a single column for binary categorical traits (proportion max category),
# and multiple columns for multi-categorical traits (proportion each category)
genus.means <- function(species.dat, trait.name){
  
  # first let's deal with the continuous traits:
  if(is.numeric(species.dat$value)){
    
    # split by genus and calculate the mean if genera contain multiple species observations
    genus.vals <- sapply(split(species.dat, as.factor(species.dat$genus)), function(x){
      if(length(x$value > 1)){return(mean(as.numeric(x$value)))
      } else {return(as.numeric(x$value))} # single species per genus
    })
    
    # return a nice data frame
    output <- data.frame(genus=names(genus.vals),
                         value=genus.vals,
                         row.names=NULL)
    colnames(output)[2] <- trait.name
    return(output)
  }
  
  # then let's deal with the categorical traits:
  else {
    
    # tabulate the trait data per category
    gentab <- table(species.dat$genus, species.dat$value)
    
    # for binary categorical traits:
    if(length(colnames(gentab))==2){
      
      # find the most common category in the binary
      maincol <- which.max(colSums(gentab))
     
      # calculate genus mean - proportion of entries for the main category
      mainprop <- gentab[,maincol]/rowSums(gentab)
      
      # coerce these data to an output data frame
      output <- data.frame(genus=row.names(gentab),
                           value=mainprop,
                           row.names=NULL)
      colnames(output)[2] <- str_c(trait.name, "proportion", names(maincol), sep=" ")
    }
    
    # lastly consider multi-categorical traits:
    else {
      
      # calculate genus mean - proportion of entries per category
      proptab <- gentab/rowSums(gentab)
      
      # coerce these data to an output data frame
      propmat <- matrix(data=proptab, ncol=ncol(proptab))
      output <- as.data.frame(propmat, row.names(proptab))
      output <- cbind(row.names(proptab), output)
      
      colnames(output) <- c("genus", str_c(trait.name, colnames(proptab), sep=" "))
      row.names(output) <- NULL
      
    }
    
    return(output)
    
  }
}

# find the genus-level means for each functional trait
genus_dat <- lapply(names(spec_dat), function(x){
  output <- genus.means(spec_dat[[x]], x)
  return(output)
})

#### 4. trait-based community classification ####

### create the genus trait matrix

# create a vector of all our 43 genera
genus_vec <- unique(species_list$genus)

# equalise all trait data frames by including all genera (n = 43) as rows
genus_dat <- lapply(genus_dat, function(x){
  
  # for this trait, which genera are we missing data for?
  missing <- genus_vec[which(is.na(match(genus_vec, x$genus)))]
  
  # create the missing data
  add_dat <- matrix(NA, 
                    nrow = length(missing),
                    ncol = (ncol(x)-1))
  add_dat <- as.data.frame(cbind(missing, add_dat))
  colnames(add_dat) <- colnames(x)
  
  # bind the trait data with the missing data and clean up
  new_dat <- rbind(x, add_dat)
  new_dat <- new_dat[order(new_dat$genus),]
  row.names(new_dat) <- NULL
  
  # we now include all 43 genera in our trait datasets
  return(new_dat)
  
})

# trait data should have the same n(rows) for same n(genera)
lapply(genus_dat, function(x){length(x[,1])})

# isolate the trait data
value_dat <- lapply(genus_dat, function(x){
  
  values <- x[,-1]
  return(values)
  
})

# bind trait values
genus_trait_df <- as.data.frame(do.call("cbind", value_dat))

# genus observations (rows)
row.names(genus_trait_df) <- genus_dat[[1]]$genus

# trait variables (columns)
trait_cols <- unlist(lapply(genus_dat, names))
trait_cols <- trait_cols[-which(trait_cols == "genus")]
colnames(genus_trait_df) <- trait_cols

# make the trait values numeric
for(i in 1:dim(genus_trait_df)[2]){
  genus_trait_df[,i] <- as.numeric(genus_trait_df[,i])
}

### select a subset of traits

# how representative are our trait data?
# let's see the percentage of all genera with at least one data point
# across each of our 15 relevant traits (genera with data / total n genera)
names(value_dat) <- names(trait_list)

sort(unlist(lapply(value_dat, function(x){
  
  filled <- length(which(!is.na(x)))/
    (length(unlist(x))/43)
  
  return(filled/43*100)
  
})))

# remove under-represented traits
# thus leaving 10 traits containing data for the majority of genera (>60%)
genus_trait_mat <- genus_trait_df[,-c(1,5,22,30,31)]

### calculate trait-based community dissimilarity

# calculate Gower's dissimilarity among genera for all traits
gower_diss <- gowdis(genus_trait_mat,
                     w = c((1/3),(1/3),(1/3), 1, 1, 1,
                           (1/13),(1/13),(1/13),(1/13),(1/13),(1/13),
                           (1/13),(1/13),(1/13),(1/13),(1/13),(1/13),(1/13),
                           0.25,0.25,0.25,0.25, 1, 1, 1, 1))

# calculate Mean Pairwise Distance (MPD) between time points in time series
functional_dissims <- lapply(sorted_matrices, 
                             comdist,
                             as.matrix(gower_diss),
                             abundance.weighted = TRUE)

# now we have 55 community functional dissimilarity matrices!
names(functional_dissims)

#### 5. novelty detection framework ####

### use the framework to detect ecological novelty through each time series

# taxonomic novelty
framework_tax <- do.call("rbind", lapply(sorted_matrices, identify.novel.gam,
                                         alpha = 0.05, metric = "bray", site = "qld"))
output_tax <- as.data.frame(framework_tax)

# functional novelty
framework_fun <- do.call("rbind", lapply(functional_dissims, identify.novel.gam,
                                         alpha = 0.05, metric = "gower", site = "qld",
                                         reverse.mat = FALSE))
output_fun <- as.data.frame(framework_fun)

### create a dataset for comparing taxonomic and functional novelty

# establish our cores and reef sites across all observed coral communities
core_rows <- sapply(sorted_matrices, row.names)
core_length <- sapply(core_rows, length)
core_ID <- rep(names(sorted_matrices), core_length)
reef_site <- str_sub(core_ID, end = -3)

# write a function to convert the output's "category" column into binary classification
cat.convert <- function(cat.column, novelty.category){
  
  ifelse(cat.column == novelty.category, 1, 0)
  
}

# make the data frame
novelty <- data.frame(sample = c(paste(core_ID, output_tax$bins, sep = ".")),
                      core = core_ID,
                      bins = output_tax$bins,
                      reefsite = reef_site,
                      region = (c(rep("frankland", 3), rep("hervey", 2), rep("keppel", 2),
                                  rep("palm", 5)))[as.factor(reef_site)],
                      taxonomic = output_tax$cat,
                      functional = output_fun$cat,
                      tax_novel = cat.convert(output_tax$cat, "novel"),
                      tax_instant = cat.convert(output_tax$cat, "instant"),
                      tax_cumul = cat.convert(output_tax$cat, "cumul"),
                      fun_novel = cat.convert(output_fun$cat, "novel"),
                      fun_instant = cat.convert(output_fun$cat, "instant"),
                      fun_cumul = cat.convert(output_fun$cat, "cumul"))

# this "novelty" data frame will be heavily relied on for following analyses

# head to 2-statistical-analysis.R!

#