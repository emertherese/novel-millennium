##-----------------------------------------------##
##  NOVEL MILLENNIUM - TESTING                   ##
##  Ecological novelty in Queensland corals      ##
##  Emer Cunningham 2024                         ##
##-----------------------------------------------##

# if you didn't come from 1-data-preparation.R,
source("1-data-preparation.R")

# and bring in novel community emergence data for comparisons
aim1_dat <- readRDS("output/aim1-dat.RDS")

#

#### testing time bin size ####

# using 20-year time bins:

# what percentage of samples follow a large temporal gap (> 60 years)?
length(which(output_tax$bin.lag > 60)) / 1130 * 100
# 1.85%

# what percentage of samples are consecutive (20 years apart)?
length(which(output_tax$bin.lag == 20)) / 1130 * 100
# 80.18%

# what about other time bin sizes?

# a function to sort raw data into core matrices
bin.trial <- function(bin_size){
  
  # check the limits of our time period across all Qld cores
  core_ages <- range(sapply(raw_matrices, function(x){as.numeric(rownames(x))}))
  
  # create a sequence of our time bins that encompasses our time period
  binseq_qld <- seq(((core_ages[1] %/% bin_size) - 1) * bin_size,
                    ((core_ages[2] %/% bin_size) + 1) * bin_size,
                    bin_size)
  
  # function for time-sorting the raw abundance matrices
  time.bin.sort <- function(core, binseq){
    # cut the core years into bins
    age_char <- as.character(cut(as.numeric(rownames(core)), binseq))
    # calculate mean of bin range, string numbers
    bin_centre <- rowMeans(cbind(lower = as.numeric(sub("\\((.+),.*", "\\1", age_char)),
                                 upper = as.numeric(sub("[^,]*,([^]]*)\\]", "\\1", age_char))))
    # make a matrix with coral genus relative abundance for each bin year
    # ifelse does time aggregation for bins with >1 sample
    bin_mat <- sapply(split(core, f = bin_centre), function(x){
      if (!is.null(ncol(x))) {return(colSums(x))} 
      else {return(x)}
    })
    # sort so we see observations of genera forward through time
    bin_mat <- t(bin_mat)
    bin_mat <- bin_mat[order(as.numeric(row.names(bin_mat)), decreasing = TRUE),]
    return(bin_mat)
  }
  
  # time sort the raw abundance matrices
  abund_matrices <- lapply(raw_matrices, time.bin.sort, binseq = binseq_qld)
  
  # remove empty rows (with no coral abundance data)
  # in High 1.4, Big Woody 1.4, Four Mile 1.3, 1.4, Pt Vernon 1.1
  abund_matrices <- lapply(abund_matrices, function(core){
    missing <- which(rowSums(core) == 0)
    if(length(missing) > 0){
      core <- core[-c(missing),]
    }
    return(core)
  })
  
  # remove cores with a single genus through time or a single time bin
  abund_matrices$Palm_Islands.Pandora.2.F <- NULL
  abund_matrices$Palm_Islands.Pandora.2.G <- NULL
  abund_matrices$Hervey_Bay.Pt_Vernon_West_Reef.1.1 <- NULL
  
  # create relative abundance matrices
  sorted_matrices <- lapply(abund_matrices, function(raw_abund){
    rel_abund <- raw_abund/rowSums(raw_abund)
    return(rel_abund)
  })
  
  # select our final taxonomic compositional matrices
  # keeping in mind that the novelty detection framework requires ~5 obsv.
  
  # check how many time bins each core has
  sort(sapply(sorted_matrices, function(x){dim(x)[1]}))
  # remove cores with <5 aggregated time points
  sorted_matrices <- sapply(sorted_matrices, function(core){
    if(dim(core)[1] > 4){return(core)}
    else {core <- NA}
  })
  sorted_matrices <- subset(sorted_matrices, !is.na(sorted_matrices))
  
  # we now have 55 taxonomic matrices we can analyse!
  return(sorted_matrices)
  
}

# a function to return some summary values, like the percentages above
bin.summary <- function(core_matrices, bin_size){
  
  # run the framework
  output <- do.call("rbind", lapply(core_matrices, identify.novel.gam,
                                    alpha = 0.05, metric = "bray", site = "qld",
                                    plot = FALSE))
  
  # tests:
  
  # what percentage of samples follow a large temporal gap (> 3*time bin size)?
  result_gap <- length(which(output$bin.lag > (bin_size*3))) /
    length(output$bin.lag) * 100
  
  # what percentage of samples are consecutive (a single time bin apart)?
  result_seq <- length(which(output$bin.lag == bin_size)) / 
    length(output$bin.lag) * 100
  
  # novelty detection compared to the chosen 20-year time bins
  result_novel <- length(which(output$novel == TRUE))
  result_ratio <- length(which(output$novel == TRUE))/dim(output)[1]
    
  # summarise in a nice table
  dat <- data.frame(test = c("Time bin size",
                             "Samples following large temporal gap",
                             "Consecutive samples",
                             "True novel communities",
                             "Ratio of novelty"),
                    result = c(bin_size,
                               round(result_gap, 2),
                               round(result_seq, 2),
                               result_novel,
                               result_ratio))
  
  # and show this table
  return(dat)
  
}

# run using, for example, 50-year time bins

your_bin = 50

# 1. create genus matrices at the temporal grain you just set
core_data <- bin.trial(your_bin)
# 2. run summary statistics on the fullness of time series and rates of novelty
core_summary <- bin.summary(core_data, your_bin)
# 3. have a look!
core_summary

#

# how accurate are our time bin ratios across novelty detection studies?

# read in the data using a separate R script
source("data/supp-novelty-lengths.R")

# this script calculates the ratio of "time series length:bin size"
# i.e., gives X in X:1 of time scale:bin

# summary of time scale:bin for Pandolfi et al. 2020
summary(as.vector(grainToExtent_plankton)) # mean 316:1

# summary of time scale:bin for Staples et al. 2022
summary(as.vector(grainToExtent_pollen)) # mean 49:1

# let's calculate similar summary statistics for this study
grainToExtent_coral <- unlist(lapply(sorted_matrices, function(core){
  
  # extract the time bins from each core
  years <- as.numeric(rownames(core))
  
  # calculate the extent of this time series
  extent <- abs(max(years) - min(years))
  
  # divide time series length by the bin size to get the scale:bin ratio
  ratio <- extent/20
  
  # and save this ratio
  return(ratio)
  
}))

# summary of time scale:bin for this study
summary(grainToExtent_coral) # mean 24:1

#

#### testing trait methods ####

# do we see a difference in the number of functional novel communities
# when we choose mode instead of mean genus trait values?

# use mode instead of mean for categorical traits,
# re-using the original code from 1-data-preparation.R

#### 3. wrangling the coral trait data

# write a function to calculate genus-level modes
genus.modes <- function(species.dat, trait.name){
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
    # here's where we want the mode
    genlist <- split(species.dat, f = species.dat$genus)
    genmodes <- lapply(genlist, function(gen){
      mode <- which.max(table(gen$value))
      return(names(mode))
    })
    output <- data.frame(genus = names(genlist),
                         value = unlist(genmodes),
                         row.names = NULL)
  }
  return(output)
}

# find the genus-level means for each functional trait
genus_dat_mode <- lapply(names(spec_dat), function(x){
  output <- genus.modes(spec_dat[[x]], x)
  return(output)
})

#### 4. trait-based community classification

### create the genus trait matrix

# create a vector of all our 43 genera
genus_vec <- unique(species_list$genus)

# equalise all trait data frames by including all genera (n=43) as rows
genus_dat_mode <- lapply(genus_dat_mode, function(x){
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

# isolate the trait data
value_dat_mode <- lapply(genus_dat_mode, function(x){
  values <- x[,-1]
  return(values)
})

# bind trait values
genus_trait_mode <- as.data.frame(do.call("cbind", value_dat_mode))
# genus observations (rows)
row.names(genus_trait_mode) <- genus_dat_mode[[1]]$genus
# trait variables (columns)
trait_cols_m <- unlist(lapply(genus_dat_mode, names))
trait_cols_m <- trait_cols_m[-which(trait_cols_m == "genus")]
colnames(genus_trait_mode) <- trait_cols_m

### select a subset of traits

# how representative are our trait data?
# let's see the percentage of all genera with at least one data point
# across each of our 15 relevant traits (genera with data / total n genera)
names(value_dat_mode) <- names(trait_list)
sort(unlist(lapply(value_dat_mode, function(x){
  filled <- length(which(!is.na(x)))/(length(unlist(x))/43)
  return(filled/43*100)
})))

# remove under-represented traits
# thus leaving 10 traits containing data for >60% genera
genus_trait_mode_mat <- genus_trait_mode[,-c(1,5,22,30,31)]

### calculate trait-based community dissimilarity

# calculate Gower's dissimilarity among genera for all traits
gower_diss_mode <- gowdis(genus_trait_mode_mat)

# calculate Mean Pairwise Distance (MPD) between time points in time series
functional_diss_mode <- lapply(sorted_matrices, comdist,
                              as.matrix(gower_diss_mode),
                              abundance.weighted = TRUE)

# now we have 55 community functional dissimilarity matrices!
names(functional_diss_mode)

### framework

framework_mode <- do.call("rbind", lapply(functional_diss_mode, identify.novel.gam,
                                          alpha = 0.05, metric = "gower", site = "qld",
                                          reverse.mat = FALSE))
output_mode <- as.data.frame(framework_mode)

# make the data frame
novelty_mode <- data.frame(sample = c(paste(core_ID, output_tax$bins, sep = ".")),
                           core = core_ID,
                           bins = output_tax$bins,
                           reefsite = reef_site,
                           region = (c(rep("frankland", 3), rep("hervey", 2), rep("keppel", 2),
                                       rep("palm", 5)))[as.factor(reef_site)],
                           taxonomic = output_tax$cat,
                           functional = output_mode$cat,
                           tax_novel = cat.convert(output_tax$cat, "novel"),
                           tax_instant = cat.convert(output_tax$cat, "instant"),
                           tax_cumul = cat.convert(output_tax$cat, "cumul"),
                           fun_novel = cat.convert(output_mode$cat, "novel"),
                           fun_instant = cat.convert(output_mode$cat, "instant"),
                           fun_cumul = cat.convert(output_mode$cat, "cumul"))

#### model the probability of novel community emergence

# novel
mmodel1.t.n <- glmer(tax_novel ~ 1 + (1|reefsite), family = binomial, data = novelty_mode)
mmodel1.f.n <- glmer(fun_novel ~ 1 + (1|reefsite), family = binomial, data = novelty_mode)
# instant
mmodel1.t.i <- glmer(tax_instant ~ 1 + (1|reefsite), family = binomial, data = novelty_mode)
mmodel1.f.i <- glmer(fun_instant ~ 1 + (1|reefsite), family = binomial, data = novelty_mode)
# cumul
mmodel1.t.c <- glmer(tax_cumul ~ 1 + (1|reefsite), family = binomial, data = novelty_mode,
                     control = glmerControl(optimizer = "bobyqa"))
mmodel1.f.c <- glmer(fun_cumul ~ 1 + (1|reefsite), family = binomial, data = novelty_mode)

# listify the models
aim1_mmodels <- list(mmodel1.t.n, mmodel1.t.i, mmodel1.t.c,
                     mmodel1.f.n, mmodel1.f.i, mmodel1.f.c)
names(aim1_mmodels) <- c("tax_novel", "tax_instant", "tax_cumul",
                         "fun_novel", "fun_instant", "fun_cumul")

# extract their probabilities
aim1_mprobs <- lapply(aim1_mmodels, function(model){
  probs <- summary(model)$coefficients[1]
  return(plogis(probs))
})
aim1_mupper <- lapply(aim1_mmodels, function(model){
  upper <- summary(model)$coefficients[1] + (1.96 * summary(model)$coefficients[2])
  return(plogis(upper))
})
aim1_mlower <- lapply(aim1_mmodels, function(model){
  lower <- summary(model)$coefficients[1] - (1.96 * summary(model)$coefficients[2])
  return(plogis(lower))
})

# tabulate the output of the models
aim1_dat_mode <- data.frame(class = c(rep("tax", 3), rep("fun", 3)),
                            category = c("novel", "instant", "cumul",
                                         "novel", "instant", "cumul"),
                            probability = unlist(aim1_mprobs),
                            se_upper = unlist(aim1_mupper),
                            se_lower = unlist(aim1_mlower),
                            position = c(1, 2, 3,
                                         1, 2, 3),
                            row.names = NULL)

aim1_dat_mode

### plot side by side

# here, original results are shown in grey, while the coloured estimates
# are novelty emergence probabilities derived from genus mode trait values

dev.off()

par(mfrow = c(1, 2))

# A. taxonomic novel community emergence

plot(NULL, yaxs = "i",
     xlim = c(0.5, 3.5), ylim = c(0, 0.047),
     xaxt = "n", yaxt = "n",
     xlab = "Novelty classification", ylab = "Emergence probability",
     main = "taxonomic")

axis(1, at = c(1,2,3), labels = c("I", "C", "N"))
axis(2, at = c(0, 0.01, 0.02, 0.03, 0.04), las = 1)

with(aim1_dat[aim1_dat$class == "tax", ],
     arrows(x0 = position, y0 = se_upper,
            x1 = position, y1 = se_lower,
            length = 0, col = "grey"))

with(aim1_dat[aim1_dat$class == "tax", ],
     points(x = position, y = probability,
            pch = 19, col = "grey"))

# B. functional novel community emergence

plot(NULL, yaxs = "i",
     xlim = c(0.5, 3.5), ylim = c(0, 0.047),
     xaxt = "n", yaxt = "n",
     xlab = "Novelty classification", ylab = "Emergence probability",
     main = "functional")

axis(1, at = c(1,2,3), labels = c("I", "C", "N"))
axis(2, at = c(0, 0.01, 0.02, 0.03, 0.04), las = 1)

with(aim1_dat[aim1_dat$class == "fun", ],
     arrows(x0 = position - 0.1, y0 = se_upper,
            x1 = position - 0.1, y1 = se_lower,
            length = 0, col = "grey"))

with(aim1_dat[aim1_dat$class == "fun", ],
     points(x = position - 0.1, y = probability,
            pch = 19, col = "grey"))

# add the probabilities derived from genus mode trait values
with(aim1_dat_mode[aim1_dat_mode$class == "fun", ],
     arrows(x0 = c(3.1, 1.1, 2.1), y0 = se_upper,
            x1 = c(3.1, 1.1, 2.1), y1 = se_lower,
            length = 0, col = c("orange", "red", cumul.col)))

with(aim1_dat_mode[aim1_dat_mode$class == "fun", ],
     points(x = c(3.1, 1.1, 2.1), y = probability,
            pch = 21, col = "black",
            bg = c("orange", "red", cumul.col)))

#

#### testing the trait selection ####

# how does trait selection affects the number and identity of novel communities?

# create trait dissimilarity matrices for the following scenarios,
# considering a subset of coral genera that have data for growth rate

#  1. ten most data-rich traits
#  2. eleven most data-rich traits (add growth rate)
#  3. nine most data-rich traits (remove oocyte size)
#  4. eight most data-rich traits (remove oocyte size and life history strategy)

# first, specify the subset of genera with growth rate data
discard <- which(is.na(genus_trait_df$`Growth rate`))

#  1. ten most data-rich traits

diss10 <- gowdis(genus_trait_df[-discard, -c(1, 5, 22, 30, 31)],
                 w = c((1/3),(1/3),(1/3), 1, 1, 1,
                       (1/13),(1/13),(1/13),(1/13),(1/13),(1/13),
                       (1/13),(1/13),(1/13),(1/13),(1/13),(1/13),(1/13),
                       0.25,0.25,0.25,0.25, 1, 1, 1, 1))

diss10_all <- gowdis(genus_trait_df[, -c(1, 5, 22, 30, 31)],
                     w = c((1/3),(1/3),(1/3), 1, 1, 1,
                           (1/13),(1/13),(1/13),(1/13),(1/13),(1/13),
                           (1/13),(1/13),(1/13),(1/13),(1/13),(1/13),(1/13),
                           0.25,0.25,0.25,0.25, 1, 1, 1, 1))

mpd10 <- lapply(sorted_matrices, comdist,
                as.matrix(diss10_all),
                abundance.weighted = TRUE)

#  2. eleven most data rich traits (add growth rate)

diss11 <- gowdis(genus_trait_df[-discard, -c(1, 5, 30, 31)],
                 w = c((1/3),(1/3),(1/3), 1, 1, 1,
                       (1/13),(1/13),(1/13),(1/13),(1/13),(1/13),
                       (1/13),(1/13),(1/13),(1/13),(1/13),(1/13),(1/13),
                       1, 0.25,0.25,0.25,0.25, 1, 1, 1, 1))

diss11_all <- gowdis(genus_trait_df[, -c(1, 5, 30, 31)],
                     w = c((1/3),(1/3),(1/3), 1, 1, 1,
                           (1/13),(1/13),(1/13),(1/13),(1/13),(1/13),
                           (1/13),(1/13),(1/13),(1/13),(1/13),(1/13),(1/13),
                           1, 0.25,0.25,0.25,0.25, 1, 1, 1, 1))

mpd11 <- lapply(sorted_matrices, comdist,
                as.matrix(diss11_all),
                abundance.weighted = TRUE)

# compare genus dissimilarities

t.test(x = as.vector(diss10),
       y = as.vector(diss11),
       paired = TRUE)

#  3. nine most data-rich traits (remove oocyte size)

diss9 <- gowdis(genus_trait_df[-discard, -c(1, 5, 22, 28, 30, 31)],
                w = c((1/3),(1/3),(1/3), 1, 1, 1,
                      (1/13),(1/13),(1/13),(1/13),(1/13),(1/13),
                      (1/13),(1/13),(1/13),(1/13),(1/13),(1/13),(1/13),
                      0.25,0.25,0.25,0.25, 1, 1, 1))

diss9_all <- gowdis(genus_trait_df[, -c(1, 5, 22, 28, 30, 31)],
                    w = c((1/3),(1/3),(1/3), 1, 1, 1,
                          (1/13),(1/13),(1/13),(1/13),(1/13),(1/13),
                          (1/13),(1/13),(1/13),(1/13),(1/13),(1/13),(1/13),
                          0.25,0.25,0.25,0.25, 1, 1, 1))

mpd9 <- lapply(sorted_matrices, comdist,
               as.matrix(diss9_all),
               abundance.weighted = TRUE)

# compare genus dissimilarities

t.test(x = as.vector(diss10),
       y = as.vector(diss9),
       paired = TRUE)

#  4. eight most data-rich traits (remove oocyte size and life history strategy)

diss8 <- gowdis(genus_trait_df[-discard, -c(1, 5, 22, 23:26, 28, 30, 31)],
                w = c((1/3),(1/3),(1/3), 1, 1, 1,
                      (1/13),(1/13),(1/13),(1/13),(1/13),(1/13),
                      (1/13),(1/13),(1/13),(1/13),(1/13),(1/13),(1/13),
                      1, 1, 1))

diss8_all <- gowdis(genus_trait_df[, -c(1, 5, 22, 23:26, 28, 30, 31)],
                    w = c((1/3),(1/3),(1/3), 1, 1, 1,
                          (1/13),(1/13),(1/13),(1/13),(1/13),(1/13),
                          (1/13),(1/13),(1/13),(1/13),(1/13),(1/13),(1/13),
                          1, 1, 1))

mpd8 <- lapply(sorted_matrices, comdist,
               as.matrix(diss8_all),
               abundance.weighted = TRUE)

# compare genus dissimilarities

t.test(x = as.vector(diss10),
       y = as.vector(diss8),
       paired = TRUE)

# compare mean and standard errors for these trait dissimilarities

basic.stat <- function(dissim){
  
  av <- mean(as.vector(dissim))
  sd <- sd(as.vector(dissim))
  se <- sd(as.vector(dissim) / sqrt(length(as.vector(dissim))))
  
  return(data.frame(average = av,
                    std_dev = sd,
                    std_error = se))
  
}

basic.stat(diss10)
basic.stat(diss11)
basic.stat(diss9)
basic.stat(diss8)

# compare novelty framework output

novel.compare <- function(dissim_1, dissim_2){
  
  output1 <- as.data.frame(do.call("rbind", lapply(dissim_1, identify.novel.gam,
                                                   alpha = 0.05, metric = "gower", site = "qld",
                                                   reverse.mat = FALSE)))
  
  output2 <- as.data.frame(do.call("rbind", lapply(dissim_2, identify.novel.gam,
                                                   alpha = 0.05, metric = "gower", site = "qld",
                                                   reverse.mat = FALSE)))
  
  df <- data.frame(novel = c(length(which(output1$novel == TRUE)),
                             length(which(output2$novel == TRUE))),
                   instant = c(length(which(output1$instant == TRUE)),
                               length(which(output2$instant == TRUE))),
                   cumul = c(length(which(output1$cumul == TRUE)),
                             length(which(output2$cumul == TRUE))),
                   row.names = c("original", "altered"))
  
  novel1 <- which(output1$novel == TRUE)
  novel2 <- which(output2$novel == TRUE)
  
  final <- list(df, novel1, novel2)
  
  return(final)
  
}

# [[1]] shows the number of novel communities
# [[2]] shows the identity of original novel communities (original 10 traits)
# [[3]] shows the identity of altered trait selection (11, 9, or 8 traits)

novel.compare(mpd10, mpd11)

novel.compare(mpd10, mpd9)

novel.compare(mpd10, mpd8)

#