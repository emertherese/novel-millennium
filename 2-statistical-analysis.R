##-----------------------------------------------##
##  NOVEL MILLENNIUM - ANALYSIS                  ##
##  Ecological novelty in Queensland corals      ##
##  Emer Cunningham 2024                         ##
##-----------------------------------------------##

# load necessary packages
library(lme4)
library(vegan)

# if you didn't come from 1-data-preparation.R,
source("1-data-preparation.R")

#

#### 1. novel community emergence ####

### model the probability of emergence across all coral community samples (n = 1130)

# novel
model1.t.n <- glmer(tax_novel ~ 1 + (1|reefsite), family = binomial, data = novelty)
model1.f.n <- glmer(fun_novel ~ 1 + (1|reefsite), family = binomial, data = novelty)

# instant
model1.t.i <- glmer(tax_instant ~ 1 + (1|reefsite), family = binomial, data = novelty)
model1.f.i <- glmer(fun_instant ~ 1 + (1|reefsite), family = binomial, data = novelty)

# cumul
model1.t.c <- glmer(tax_cumul ~ 1 + (1|reefsite), family = binomial, data = novelty,
                    control = glmerControl(optimizer = "bobyqa"))
model1.f.c <- glmer(fun_cumul ~ 1 + (1|reefsite), family = binomial, data = novelty)

# list the models
aim1_models <- list(model1.t.n, model1.t.i, model1.t.c,
                    model1.f.n, model1.f.i, model1.f.c)

names(aim1_models) <- c("tax_novel", "tax_instant", "tax_cumul",
                        "fun_novel", "fun_instant", "fun_cumul")

# extract their probabilities
aim1_probs <- lapply(aim1_models, function(model){
  probs <- summary(model)$coefficients[1]
  return(plogis(probs))
})

aim1_upper <- lapply(aim1_models, function(model){
  upper <- summary(model)$coefficients[1] + (1.96 * summary(model)$coefficients[2])
  return(plogis(upper))
})

aim1_lower <- lapply(aim1_models, function(model){
  lower <- summary(model)$coefficients[1] - (1.96 * summary(model)$coefficients[2])
  return(plogis(lower))
})

# tabulate the output of the models
aim1_dat <- data.frame(class = c(rep("tax", 3), rep("fun", 3)),
                       category = c("novel", "instant", "cumul",
                                    "novel", "instant", "cumul"),
                       probability = unlist(aim1_probs),
                       se_upper = unlist(aim1_upper),
                       se_lower = unlist(aim1_lower),
                       position = c(3, 1, 2,
                                    3, 1, 2),
                       row.names = NULL)

aim1_dat

### model the probability of emergence per region

# iterate the above
region.emerge <- function(region.character){
  
  # data
  dat <- novelty[which(novelty$region == region.character),]
  
  # model
  tax.novel <- glmer(tax_novel ~ 1 + (1|reefsite), family = binomial, data = dat)
  fun.novel <- glmer(fun_novel ~ 1 + (1|reefsite), family = binomial, data = dat)
  tax.instant <- glmer(tax_instant ~ 1 + (1|reefsite), family = binomial, data = dat)
  fun.instant <- glmer(fun_instant ~ 1 + (1|reefsite), family = binomial, data = dat)
  tax.cumul <- glmer(tax_cumul ~ 1 + (1|reefsite), family = binomial, data = dat)
  fun.cumul <- glmer(fun_cumul ~ 1 + (1|reefsite), family = binomial, data = dat)
  
  # list
  region_models <- list(tax.novel, fun.novel, tax.instant, fun.instant, tax.cumul, fun.cumul)
  
  # extract
  region_probs <- lapply(region_models, function(model){
    probs <- summary(model)$coefficients[1]
    return(plogis(probs))
  })
  region_upper <- lapply(region_models, function(model){
    upper <- summary(model)$coefficients[1] + (1.96 * summary(model)$coefficients[2])
    return(plogis(upper))
  })
  region_lower <- lapply(region_models, function(model){
    lower <- summary(model)$coefficients[1] - (1.96 * summary(model)$coefficients[2])
    return(plogis(lower))
  })
  
  # summarise
  output <- data.frame(region = region.character,
                       class = rep(c("tax", "fun"),
                                   length(region_probs)/2),
                       category = c("novel", "novel", "instant",
                                    "instant", "cumul", "cumul"),
                       probability = unlist(region_probs),
                       se_upper = unlist(region_upper),
                       se_lower = unlist(region_lower))
  
  # return
  return(output)
  
}

# emergence data per region
aim1_region_dat <- rbind(region.emerge("frankland"),
                         region.emerge("palm"),
                         region.emerge("hervey"))

aim1_region_dat

### alter these data for plotting

# separate data sets per plot, to avoid later subsetting
region_dat <- split(aim1_region_dat, f = paste(aim1_region_dat$class, 
                                               aim1_region_dat$category,
                                               sep = ""))

# add a position column for regions: F > P > H
for(i in 1:length(region_dat)){
  
  region_dat[[i]]$position <- c(1, 2, 3)
  
}

#

#### 2. novel community co-occurrence ####

### model pairwise taxonomic-functional co-occurrence probabilities

# novel
model2.t.n <- glmer(tax_novel ~ fun_novel + fun_instant + fun_cumul + (1|reefsite),
                    family = binomial, data = novelty)
model2.f.n <- glmer(fun_novel ~ tax_novel + tax_instant + tax_cumul + (1|reefsite),
                    family = binomial, data = novelty)

# instant
model2.t.i <- glmer(tax_instant ~ fun_novel + fun_instant + fun_cumul + (1|reefsite),
                    family = binomial, data = novelty)
model2.f.i <- glmer(fun_instant ~ tax_novel + tax_instant + tax_cumul + (1|reefsite),
                    family = binomial, data = novelty)

# cumul
model2.t.c <- glmer(tax_cumul ~ fun_novel + fun_instant + fun_cumul + (1|reefsite),
                    control = glmerControl(optimizer = "bobyqa"),
                    family = binomial, data = novelty)
model2.f.c <- glmer(fun_cumul ~ tax_novel + tax_instant + tax_cumul + (1|reefsite),
                    family = binomial, data = novelty)

# list the models
aim2_models <- list(model2.t.n, model2.t.i, model2.t.c,
                    model2.f.n, model2.f.i, model2.f.c)

names(aim2_models) <- c("tax_novel", "tax_instant", "tax_cumul",
                        "fun_novel", "fun_instant", "fun_cumul")

# extract their probabilities
aim2_probs <- lapply(aim2_models, function(model){
  p1 <- summary(model)$coefficients[1] + summary(model)$coefficients[2]
  p2 <- summary(model)$coefficients[1] + summary(model)$coefficients[3]
  p3 <- summary(model)$coefficients[1] + summary(model)$coefficients[4]
  return(c(plogis(p1),
           plogis(p2),
           plogis(p3)))
})

aim2_upper <- lapply(aim2_models, function(model){
  e1 <- (summary(model)$coefficients[1] + summary(model)$coefficients[2]) + (1.96 * summary(model)$coefficients[6])
  e2 <- (summary(model)$coefficients[1] + summary(model)$coefficients[3]) + (1.96 * summary(model)$coefficients[7])
  e3 <- (summary(model)$coefficients[1] + summary(model)$coefficients[4]) + (1.96 * summary(model)$coefficients[8])
  return(c(plogis(e1),
           plogis(e2),
           plogis(e3)))
})

aim2_lower <- lapply(aim2_models, function(model){
  e1 <- (summary(model)$coefficients[1] + summary(model)$coefficients[2]) - (1.96 * summary(model)$coefficients[6])
  e2 <- (summary(model)$coefficients[1] + summary(model)$coefficients[3]) - (1.96 * summary(model)$coefficients[7])
  e3 <- (summary(model)$coefficients[1] + summary(model)$coefficients[4]) - (1.96 * summary(model)$coefficients[8])
  return(c(plogis(e1),
           plogis(e2),
           plogis(e3)))
})

# tabulate the output of the models
aim2_dat <- data.frame(resp_class = c(rep("tax", 9), rep("fun", 9)),
                       resp_category = rep(c(rep("novel", 3),
                                             rep("instant", 3),
                                             rep("cumul", 3)), times = 2),
                       pred_class = c(rep("fun", 9), rep("tax", 9)),
                       pred_category = rep(c("novel", 
                                             "instant", 
                                             "cumul"),
                                           times = 6),
                       probability = unlist(aim2_probs),
                       se_upper = unlist(aim2_upper),
                       se_lower = unlist(aim2_lower),
                       row.names = NULL)

aim2_dat

#

#### 3. novel community composition ####

### create multivariate taxonomic datasets
### across all coral community samples (n = 1130):
###    qld_dat = composition and metadata
###    qld_env = descriptive metadata
###    tax_dat = taxonomic composition data

# convert taxonomic compositional matrices to data frames
sorted_dfs <- lapply(sorted_matrices, as.data.frame)

# join them up
join.fn <- function(df1, df2){
  
  join(df1, df2,
       verbose = FALSE, na.replace = TRUE, split = FALSE)
  
}
qld_dat <- Reduce(join.fn, sorted_dfs)

# add metadata
qld_dat <- cbind(qld_dat, novelty)
row.names(qld_dat) <- qld_dat$sample

# create an associated descriptive dataset
qld_env <- qld_dat[, which(is.na(match(names(qld_dat), 
                                       genus_vec)))]
row.names(qld_env) <- qld_env$sample

# create just the taxonomic compositional matrix
tax_dat <- qld_dat[, which(!is.na(match(names(qld_dat),
                                        genus_vec)))]

# and one with an abundance cut-off
sub.comm.dat <- function(comm_dat, cutoff){
  
  # calculate the summed relative abundances of "species" in community
  relabunds <- colSums(comm_dat)
  
  # choose those species to include over a certain cutoff
  species <- relabunds[which(relabunds > cutoff)]
  
  # we hate Alveopora
  if(length(unique(match(names(species), "Alveopora"))) > 1){
    species <- species[-which(names(species) == "Alveopora")]
  }
  
  # create a new community data frame with species above our cutoff
  new_dat <- comm_dat[,names(species)]
  
  # calculate the new relative abundances
  new_dat <- new_dat/rowSums(new_dat)
  
  # remove communities that only featured rare species (below cutoff)
  new_dat <- new_dat[-which(is.nan(rowSums(new_dat)) == TRUE),]
  
  # return the new community data frame
  return(new_dat)
  
}

# exclude uncommon taxa
sub_tax_dat <- sub.comm.dat(tax_dat, 10)

# match this subset with their metadata
rm_comms <- which(is.na(match(row.names(qld_dat), row.names(sub_tax_dat))))
qld_dat[rm_comms,] #it's the Alveopora-dominated samples
sub_qld_dat <- qld_dat[-rm_comms,]

# we're now working with the same n(coral community samples)
dim(sub_qld_dat)
dim(sub_tax_dat)

### create multivariate trait-based datasets
### across all coral community samples (n = 1130):
###    fun_dat = trait-based composition dissimilarity matrix

fun_dat <- comdist(tax_dat,    #community matrix
                   gower_diss, #distance matrix
                   abundance.weighted = TRUE)

# or if you don't have approximately 11 minutes, run this line:
# fun_dat <- read.csv("output/functional-dissimilarities.csv", row.names = 1)

### test for variation in community composition over space and time

# taxonomic PERMANOVA
tax_perm <- adonis2(tax_dat ~ reefsite * core * as.numeric(bins),
                    method = "bray", data = qld_env)

# trait-based PERMANOVA
fun_perm <- adonis2(fun_dat ~ reefsite * core * as.numeric(bins),
                    method = "gower", data = qld_env)

# # or if you don't have approximately 50 minutes, run these lines:
# tax_perm <- read.csv("output/permanova-taxonomic.csv")
# fun_perm <- read.csv("output/permanova-functional.csv")

### variation in community composition by novelty category

# whole community ordination
tax_ord <- metaMDS(sub_tax_dat, distance = "bray")
fun_ord <- metaMDS(fun_dat, distance = "gower")

# check stress values
tax_ord$stress
fun_ord$stress

# calculate species vectors for taxonomic classifications
tax_vec <- envfit(ord = tax_ord, env = sub_tax_dat)

#

#### 4. genus compositional change ####

### data for model

# create long form dataset of genus relative abundance across Queensland
rela_dat <- sub_qld_dat %>%
  pivot_longer(cols = Acropora:Caryophyllia,
               names_to = "genus", values_to = "rel_abund") %>%
  select(!c(taxonomic, tax_instant, tax_cumul,
            functional, fun_instant, fun_cumul))

# write a function that calculates changes in relative abundance
delta.rel.abund <- function(tax.dat, env.dat, framework.output){
  
  # establish time t
  t <- 2:dim(tax.dat)[1]
  
  # calculate change in relative abundance
  dat <- tax.dat[t, ] - tax.dat[(t - 1), ]
  
  # remove each core's first bin
  first.bins <- env.dat$sample[which(is.na(framework.output$cat.bef))]
  rm.bins <- which(!is.na(match(row.names(dat),
                                first.bins)))
  dat <- dat[-rm.bins, ]
  
  # keep sample as a column
  meta_dat <- cbind(dat, env.dat[match(row.names(dat), env.dat$sample),
                                 c(2:8, 11)])
  
  # pop this into a long format, for modelling
  long_dat <- meta_dat %>%
    rownames_to_column(var = "sample") %>%
    pivot_longer(cols = Acropora:Millepora,
                 names_to = "genus", values_to = "rel_abund")
  
  # remove no-change observations
  long_dat <- long_dat[-which(long_dat$rel_abund == 0), ]
  
  # save modelling dataset
  return(long_dat)
  
}

# create long form dataset of relative abundance shifts across Queensland
aim3_dat <- delta.rel.abund(sub_tax_dat, 
                            qld_env, 
                            framework_tax)

### model genus relative abundance shifts

# model
model3.rel <- lmer(rel_abund ~ -1 + genus * tax_novel + (1|reefsite),
                   data = aim3_dat)

summary(model3.rel)

# significant genus shifting in novel or non-novel states
which(abs(rel(summary(model3.rel)$coefficients[,3])) > 1.96)

# predict changes in relative abundance
rel_dat <- expand.grid(genus = sort(unique(aim3_dat$genus)),
                       tax_novel = c(0, 1))

rel_dat$abund <- predict(model3.rel, 
                         rel_dat, 
                         re.form = NA)

rel_dat$colour <- c(rep("grey", 14), rep("orange", 14))

# add confidence intervals for model estimates
rel_dat <- add_ci(rel_dat, model3.rel, alpha = 0.05,
                  includeRanef = FALSE, type = "boot", nSims = 9999)

# or if you're short on time, run this line:
# rel_dat <- read.csv("output/rel-dat.csv")

#

#

#### 5. save all model outputs ####

### novel community emergence

# Queensland models

# extract model summary and arrange into table
table_s4 <- matrix(unlist(lapply(aim1_models, function(model){
  
  return(summary(model)$coefficients)
  
})), 

nrow = 6, ncol = 4, byrow = TRUE, 
dimnames = list(c(names(aim1_models)), 
                c("estimate", "std. error", "z", "p"))

)

# save as data frame
table_s4 <- as.data.frame(table_s4)

# add random effect variance
table_s4$reVar <- unlist(sapply(aim1_models, function(model){
  
  VarCorr(model)$reefsite[1]
  
}))

# save summary table
write.csv(table_s4, "output/table-S4.csv")

# regional models

# use function from above to save models, not model estimates
region.emerge.models <- function(region.character){
  
  # data
  dat <- novelty[which(novelty$region == region.character),]
  
  # model
  tax.novel <- glmer(tax_novel ~ 1 + (1|reefsite), family = binomial, data = dat)
  fun.novel <- glmer(fun_novel ~ 1 + (1|reefsite), family = binomial, data = dat)
  tax.instant <- glmer(tax_instant ~ 1 + (1|reefsite), family = binomial, data = dat)
  fun.instant <- glmer(fun_instant ~ 1 + (1|reefsite), family = binomial, data = dat)
  tax.cumul <- glmer(tax_cumul ~ 1 + (1|reefsite), family = binomial, data = dat)
  fun.cumul <- glmer(fun_cumul ~ 1 + (1|reefsite), family = binomial, data = dat)
  
  # list
  region_models <- list(tax.novel, fun.novel, tax.instant, fun.instant, tax.cumul, fun.cumul)
  
  # return
  return(region_models)
  
}

# save models as a list
frankland_models <- region.emerge.models("frankland")
palm_models <- region.emerge.models("palm")
hervey_models <- region.emerge.models("hervey")

# table s5

# extract model summary and arrange into table
table_s5 <- matrix(unlist(lapply(frankland_models, function(model){
  
  return(summary(model)$coefficients)
  
})), 

nrow = 6, ncol = 4, byrow = TRUE, 
dimnames = list(c(names(frankland_models)), 
                c("estimate", "std. error", "z", "p"))

)

# save as data frame
table_s5 <- as.data.frame(table_s5)

# add random effect variance
table_s5$reVar <- unlist(sapply(frankland_models, function(model){
  
  VarCorr(model)$reefsite[1]
  
}))

# save summary table
write.csv(table_s5, "output/table-S5.csv")

# table s6

# extract model summary and arrange into table
table_s6 <- matrix(unlist(lapply(palm_models, function(model){
  
  return(summary(model)$coefficients)
  
})), 

nrow = 6, ncol = 4, byrow = TRUE, 
dimnames = list(c(names(palm_models)), 
                c("estimate", "std. error", "z", "p"))

)

# save as data frame
table_s6 <- as.data.frame(table_s6)

# add random effect variance
table_s6$reVar <- unlist(sapply(palm_models, function(model){
  
  VarCorr(model)$reefsite[1]
  
}))

# save summary table
write.csv(table_s6, "output/table-S6.csv")

# table s7

# extract model summary and arrange into table
table_s7 <- matrix(unlist(lapply(hervey_models, function(model){
  
  return(summary(model)$coefficients)
  
})), 

nrow = 6, ncol = 4, byrow = TRUE, 
dimnames = list(c(names(hervey_models)), 
                c("estimate", "std. error", "z", "p"))

)

# save as data frame
table_s7 <- as.data.frame(table_s7)

# add random effect variance
table_s7$reVar <- unlist(sapply(hervey_models, function(model){
  
  VarCorr(model)$reefsite[1]
  
}))

# save summary table
write.csv(table_s7, "output/table-S7.csv")

### novel community co-occurrence

aim2_models

# extract model summary and arrange into table
table_s8 <- matrix(unlist(lapply(aim2_models, function(model){
  
  return(summary(model)$coefficients)
  
})), 

nrow = 24, ncol = 4, byrow = TRUE, 
dimnames = list(c(rep(names(aim2_models)[1], 4),
                  rep(names(aim2_models)[2], 4),
                  rep(names(aim2_models)[3], 4),
                  rep(names(aim2_models)[4], 4),
                  rep(names(aim2_models)[5], 4),
                  rep(names(aim2_models)[6], 4)), 
                c("estimate", "std. error", "z", "p"))

)

# save as data frame
table_s8 <- as.data.frame(table_s8)

# add random effect variance
table_s8$reVar <- as.vector(sapply(aim2_models, function(model){
  
  reVar <- VarCorr(model)$reefsite[1]
  
  return(rep(reVar, 4))
  
}))

# save summary table
write.csv(table_s8, "output/table-S8.csv")

### novel community composition

# save the functional dissimilarity matrix data
write.csv(as.matrix(fun_dat), file = "output/functional-dissimilarities.csv")

# save the PERMANOVA output
write.csv(as.data.frame(tax_perm), "output/permanova-taxonomic.csv")
write.csv(as.data.frame(fun_perm), "output/permanova-functional.csv")

### genus compositional change

# save model output plus confidence intervals
write.csv(rel_dat, "output/rel-dat.csv", row.names = FALSE)

# fixed effects
table_s9 <- as.data.frame(summary(model3.rel)$coefficients,
                          colnames = c("estimate", "std. error", "t"))

# random effects
table_s9_ranef <- as.data.frame(VarCorr(model3.rel))

# save both summary tables
write.csv(table_s9, "output/table-S9.csv")
write.csv(table_s9_ranef, "output/table-S9-ranef.csv")

#