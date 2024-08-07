##-----------------------------------------------##
##  NOVEL MILLENNIUM - FIGURES                   ##
##  Ecological novelty in Queensland corals      ##
##  Emer Cunningham 2024                         ##
##-----------------------------------------------##

# close and clear
close.screen(all.screens = TRUE)
dev.off()

# load necessary packages
library(cowplot)
library(ggspatial)
library(sf)
library(tidyverse)

# if you didn't come from 2-statistical-analysis.R,
source("2-statistical-analysis.R")

#

#### import data and load functions ####

# import co-ordinates data for core samples
map_dat_region <- read.csv("data/map-data-per-region.csv")
map_dat_reef <- read.csv("data/map-data-per-reef.csv")

# import Australia shapefile
# these are data from GADM, taken from their website on 27/07/2023
# https://gadm.org/download_country.html#google_vignette
aus_shp <- read_sf(dsn = "data/gadm-files/gadm41_AUS_0.shp")

# import GBR shapefile
# these are data from GBRMPA, taken from eAtlas on 18/08/2023
# https://eatlas.org.au/data/uuid/ac8e8e4f-fc0e-4a01-9c3d-f27e4a8fac3c
gbr_shp <- read_sf(dsn = "data/GBRMPA-files/Great_Barrier_Reef_Features.shp")

# function to plot individual cores' coral communities in nMDS space
core.nmds <- function(core.char, light.col, dark.col){
  
  # find the coral community samples from your core
  core_names <- qld_env$sample[qld_env$core == core.char]
  
  # find their ages and assign colours using time
  core_ages <- sort(as.numeric(qld_env$bins[qld_env$core == core.char]))
  
  cols <- colorRampPalette(c(light.col, dark.col))
  
  core_cols <- cols(length(core_names))[as.numeric(cut(core_ages,
                                                       breaks = length(core_ages)))]
  
  # age segments
  for(i in 1:(length(core_names) - 1)){
    segments(x0 = tax_ord$points[core_names[i], 1],
             y0 = tax_ord$points[core_names[i], 2],
             x1 = tax_ord$points[core_names[i+1], 1],
             y1 = tax_ord$points[core_names[i+1], 2],
             col = core_cols[i+1], lwd = 2)
  }
  
  # core communities
  points(tax_ord$points[core_names,],
         pch = 21, col = "black", bg = core_cols, cex = 1.2)
  
  # novel communities
  novels <- which(qld_env$tax_novel[qld_env$core == core.char] == 1)
  
  points(x = tax_ord$points[core_names,][novels, ][1],
         y = tax_ord$points[core_names,][novels, ][2],
         pch = 21, col = "black", bg = "orange", cex = 1.2)
  
}

# function to plot individual cores' genus relative abundance through time
core.barplot <- function(core.char){
  
  # extract the specified core
  dat <- sorted_matrices[[core.char]]
  
  # change time bin centres to calendar years
  years <- (1950 - as.numeric(row.names(dat)))
  row.names(dat) <- years
  
  # transpose and reorder data to stretch from past to present
  plot_dat <- t(dat)
  plot_dat <- plot_dat[order(row.names(plot_dat)), ]
  
  # create a grayscale colour scheme
  genus_all <- sort(unique(unlist(lapply(sorted_matrices[c("Frankland_Islands.Normanby_Island.3.2",
                                                          "Palm_Islands.Havannah.2.F")],
                                         colnames))))
  scheme <- colorRampPalette(c("grey40", "lightgrey"))
  genus_cols <- scheme(length(genus_all))
  names(genus_cols) <- genus_all
  
  plot_cols <- genus_cols[match(row.names(plot_dat),
                                names(genus_cols))]
  
  # add specific colours for featured taxa
  if(core.char == "Frankland_Islands.Normanby_Island.3.2"){
    
    plot_cols["Acropora"] <- "orchid1"
    plot_cols["Porites"] <- "orchid4"
    
  }
  
  # add specific colours for featured taxa
  if(core.char == "Palm_Islands.Havannah.2.F"){
    
    plot_cols["Pavona"] <- "lightpink"
    plot_cols["Goniopora"] <- "red"
    
  }
  
  # plot
  barplot(plot_dat,
          xlab = "Year", ylab = "Relative abundance", xpd = NA,
          col = plot_cols, border = "black")
  
}

# get plotting!

#

#### figure 1. map ####

# Queensland inset
qld_map <- ggplot(map_dat_region) +
  labs(x = NULL, y = NULL) + 
  xlim(140, 155) +
  ylim(-27, -10) +
  geom_sf(data = aus_shp, 
          col = "grey", fill = "grey") +
  # geom_point(aes(x = mid_lon,
  #                y = mid_lat),
  #            colour = "black") +
  geom_text(aes(x = mid_lon,
                y = mid_lat,
                label = c("A", "D", "C", "B"),
                fontface = "bold")) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))

#

# Frankland Islands
franklands_map <- ggplot(filter(map_dat_reef, 
                                region == "Frankland Islands")) +
  labs(x = NULL, y = NULL) +
  scale_x_continuous(limit = c(145.89, 146.14),
                     breaks = c(145.9, 146, 146.1)) +
  scale_y_continuous(limit = c(-17.3, -17.1)) +
  geom_sf(data = gbr_shp,
          col = "grey", fill = "grey") +
  geom_point(aes(x = mid_lon,
                 y = mid_lat),
             colour = "black") +
  geom_text(aes(x = mid_lon,
                y = mid_lat,
                label = c("High", "Normanby", "Russell"),
                vjust = c(-2, 0.5, 2),
                hjust = c(0.5, 1.1, 0.5))) +
  geom_text(aes(x = -Inf, hjust = -0.5,
                y = Inf, vjust = 1.5,
                label = "A", fontface = "bold")) +
  annotation_scale(location = "bl") +
  theme(panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))

# Palm Islands
palm_map <- ggplot(filter(map_dat_reef, 
                          region == "Palm Islands")) +
  labs(x = NULL, y = NULL) + 
  scale_x_continuous(limit = c(146.20, 146.70),
                     breaks = c(146.2, 146.4, 146.6)) +
  scale_y_continuous(limit = c(-19.00, -18.60)) +
  geom_sf(data = gbr_shp, 
          col = "grey", fill = "grey") +
  geom_point(aes(x = mid_lon,
                 y = mid_lat),
             colour = "black") +
  geom_text(aes(x = mid_lon,
                y = mid_lat,
                label = c("Havannah", "Pandora"),
                vjust = c(1.5, -1),
                hjust = c(0.5, 0.5))) +
  geom_text(aes(x = -Inf, hjust = -0.5,
                y = Inf, vjust = 1.5,
                label = "B", fontface = "bold")) +
  annotation_scale(location = "bl") +
  theme(panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))

# Keppel Islands
keppel_map <- ggplot(filter(map_dat_reef, 
                            region == "Keppel Islands")) +
  labs(x = NULL, y = NULL) + 
  scale_x_continuous(limit = c(150.675, 151.05),
                     breaks = c(150.7, 150.8, 150.9, 151)) +
  scale_y_continuous(limit = c(-23.30, -23.00),
                     breaks = c(-23.3, -23.2, -23.1, -23.0))+
  geom_sf(data = gbr_shp, 
          col = "grey", fill = "grey") +
  geom_point(aes(x = mid_lon,
                 y = mid_lat),
             colour = "black") +
  geom_text(aes(x = mid_lon,
                y = mid_lat,
                label = c("Middle", " North\nKeppel"),
                vjust = c(0.5, 0.5),
                hjust = c(1.3, -0.3))) +
  geom_text(aes(x = -Inf, hjust = -0.5,
                y = Inf, vjust = 1.5,
                label = "C", fontface = "bold")) +
  annotation_scale(location = "bl") +
  theme(panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))

# Hervey Bay
hervey_map <- ggplot(filter(map_dat_reef, 
                            region == "Hervey Bay")) +
  labs(x = NULL, y = NULL) + 
  scale_x_continuous(limit = c(152.15, 153.40),
                     breaks = c(152.3, 152.6, 152.9, 153.2)) +
  scale_y_continuous(limit = c(-25.70, -24.70)) +
  geom_sf(data = gbr_shp, 
          col = "grey", fill = "grey") +
  geom_point(aes(x = mid_lon,
                 y = mid_lat),
             colour = "black") +
  geom_text(aes(x = mid_lon,
                y = mid_lat,
                label = c("Big\n         Woody", "Four Mile"),
                vjust = c(-0.3, 0.5),
                hjust = c(0.9, -0.2))) +
  geom_text(aes(x = -Inf, hjust = -0.5,
                y = Inf, vjust = 1.5,
                label = "D", fontface = "bold")) +
  annotation_scale(location = "bl") +
  theme(panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))

# plot
pdf(paste("figures/", Sys.Date(), "-figure-1.pdf", sep = ""),
    width = 8, height = 5, useDingbats = FALSE)

# arrange
plot_grid(qld_map, franklands_map, palm_map,
          NULL, keppel_map, hervey_map,
          nrow = 2, ncol = 3, rel_widths = c(1, 2, 2))

# close
dev.off()

#

#

#### figure 2. emergence ####

pdf(paste("figures/", Sys.Date(), "-figure-2.pdf", sep = ""),
    width = 6, height = 4, useDingbats = FALSE)

par(mar = c(0,0,0,0), ps = 10, tcl = -0.25, las = 1, xpd = NA)

split.screen(rbind(c(0.15, 0.525, 0.15, 0.95),
                   c(0.575, 0.95, 0.15, 0.95)))

# A. taxonomic novel community emergence

screen(1)

plot(NULL, yaxs = "i",
     xlim = c(0.5, 3.5), ylim = c(0, 0.047),
     xaxt = "n", yaxt = "n",
     xlab = "", ylab = "Emergence probability")

axis(1, at = c(1,2,3), labels = c("I", "C", "N"), mgp = c(2,0.2,0))
axis(2, at = c(0, 0.01, 0.02, 0.03, 0.04), 
     labels = c("0", "0.01", "0.02", "0.03", "0.04"), mgp = c(2,0.5,0))

with(aim1_region_dat[aim1_region_dat$class == "tax",],
     points(x = c(2, 1, 3)[as.factor(category)], y = probability,
            col = alpha(c(cumul.col, "red", "orange")[as.factor(category)], 0.5),
            pch = c(15, 17, 19)[as.factor(region)]))

with(aim1_dat[aim1_dat$class == "tax", ],
     arrows(x0 = position, y0 = se_upper,
            x1 = position, y1 = se_lower,
            length = 0, col = "black"))

with(aim1_dat[aim1_dat$class == "tax", ],
     points(x = position, y = probability,
            cex = 1.8, pch = 21, col = "black",
            bg = c("orange", "red", cumul.col)))

# text(x = 0.6, y = 0.0447, cex = 1.5,
#      labels = substitute(paste(bold("A"))))

text(x = 1, y = 0.0447, cex = 1,
     labels = substitute(paste(bold("taxonomic"))))

text(x = 3.75, y = -0.006, "Novelty classification", xpd = NA)

close.screen(1)

# B. functional novel community emergence

screen(2)

plot(NULL, yaxs = "i",
     xlim = c(0.5, 3.5), ylim = c(0, 0.047),
     xlab = "", ylab = "",
     xaxt = "n", yaxt = "n")

axis(1, at = c(1,2,3), labels = c("I", "C", "N"), mgp = c(2,0.2,0))
axis(2, at = c(0, 0.01, 0.02, 0.03, 0.04), 
     labels = FALSE)

with(aim1_region_dat[aim1_region_dat$class == "fun",],
     points(x = c(2, 1, 3)[as.factor(category)], y = probability,
            col = alpha(c(cumul.col, "red", "orange")[as.factor(category)], 0.5),
            pch = c(15, 17, 19)[as.factor(region)]))

with(aim1_dat[aim1_dat$class == "fun", ],
     arrows(x0 = position, y0 = se_upper,
            x1 = position, y1 = se_lower,
            length = 0, col = "black"))

with(aim1_dat[aim1_dat$class == "fun", ],
     points(x = position, y = probability,
            cex = 1.8, pch = 21, col = "black",
            bg = c("orange", "red", cumul.col)))

# legend(x = 1.95, y = 0.047, cex = 0.8, box.lty = 0,
#        legend = c("I = instantaneous", "C = cumulative", "T = true"),
#        text.col = c("red", cumul.col, "orange"))

# text(x = 0.6, y = 0.0447, cex = 1.5,
#      labels = substitute(paste(bold("B"))))

text(x = 1, y = 0.0447, cex = 1,
     labels = substitute(paste(bold("functional"))))

legend("topright", cex = 1, box.lty = 0,
       legend = c(" I = instantaneous", "C = cumulative", "N = true"),
       text.col = c("red", cumul.col, "orange"))

close.screen(2)

# close and clear
close.screen(all.screens = TRUE)
dev.off()

#

#### figure 3. co-occurrence ####

pdf(paste("figures/", Sys.Date(), "-figure-3.pdf", sep = ""),
    width = 5, height = 5.5, useDingbats = FALSE)

par(mar = c(5,4,4,2), mgp = c(2,0.4,0), ps = 10, tcl = -0.25, las = 1)

plot(NULL, #xaxs = "i", yaxs = "i",
     xlim = c(0, 0.6), ylim = c(0, 0.6),
     xlab = "P(taxonomic|functional)", ylab = "P(functional|taxonomic)",
     xaxt = "n", yaxt = "n")

axis(1)
axis(2)

abline(a = 0, b = 1, lty = 2)

with(aim2_dat,
     segments(x0 = se_lower[1:9], y0 = probability[10:18],
              x1 = se_upper[1:9], y1 = probability[10:18],
              col = alpha("grey", 0.3), lwd = 3))
with(aim2_dat,
     segments(x0 = probability[1:9], y0 = se_lower[10:18],
              x1 = probability[1:9], y1 = se_upper[10:18],
              col = alpha("grey", 0.3), lwd = 3))

points(x = aim2_dat$probability[1:9],
       y = aim2_dat$probability[10:18],
       pch = 19, cex = 1.2)

# labels, for adding outside of R given visual clutter
with(aim2_dat, text(x = probability[1:9],
                    y = probability[10:18],
                    labels = paste(toupper(str_sub(resp_category, end = 1)),
                                   toupper(str_sub(pred_category, end = 1)),
                                   sep = "~"), pos = 4))

dev.off()

#

#### figure 4. nMDS plots ####

pdf(paste("figures/", Sys.Date(), "-figure-4.pdf", sep = ""),
    width = 8, height = 4, useDingbats = FALSE)

par(mar = c(0,0,0,0), mgp = c(1.5,0.2,0), ps = 10, tcl = -0.25, las = 1, xpd = NA)

split.screen(rbind(c(0.10, 0.475, 0.15, 0.90),
                   c(0.575, 0.95, 0.15, 0.90)))

# firstly find the identity of the novel communities
tax_nov <- qld_env$sample[which(qld_env$tax_novel == 1)]
fun_nov <- qld_env$sample[which(qld_env$fun_novel == 1)]

# and parse apart their regions
tax_nov_region <- as.factor(str_sub(tax_nov, end = 1))
fun_nov_region <- as.factor(str_sub(fun_nov, end = 1))

# A. taxonomic novelty

screen(1)

plot(tax_ord, 
     type = "n", 
     xaxt = "n", yaxt = "n",
     xlab = "nMDS1", ylab = "nMDS2")

axis(1, mgp = c(1,0.2,0))
axis(2, mgp = c(1,0.4,0))

points(tax_ord,
       pch = 19, col = alpha("grey", alpha = 0.5))

points(tax_ord$points[!is.na(match(row.names(tax_ord$points), 
                                   tax_nov)), ],
       pch = 21, col = "black", bg = "orange")

sapply(which(tax_vec$vectors$pvals == 0.001), function(x){
  arrows(x0 = 0, y0 = 0,
         x1 = tax_ord$species[x, ][1], y1 = tax_ord$species[x, ][2],
         length = 0.08)
})

# sapply(which(tax_vec$vectors$pvals == 0.001), function(x){
#   text(x = tax_ord$species[x, ][1], y = tax_ord$species[x, ][2],
#        labels = row.names(tax_ord$species)[x], cex = 0.75)
# })

text(x = -1.8, y = 2.8, cex = 1.5,
     labels = substitute(paste(bold("A"))))

close.screen(1)

# B. functional novelty

screen(2)

plot(fun_ord, 
     type = "n", 
     xaxt = "n", yaxt = "n",
     xlab = "nMDS1", ylab = "nMDS2")

axis(1, mgp = c(1,0.2,0))
axis(2, mgp = c(1,0.4,0))

points(fun_ord,
       pch = 19, col = alpha("grey", alpha = 0.5))

points(fun_ord$points[!is.na(match(row.names(fun_ord$points), 
                                   fun_nov)), ],
       pch = 21, col = "black", bg = "orange")

# no genus centroids --- but trait value groupings

text(x = -0.43, y = 0.26, cex = 1.5,
     labels = substitute(paste(bold("B"))))

close.screen(2)

# close and clear
close.screen(all.screens = TRUE)
dev.off()

#

#### figure 5. relative abundance shifts ####

pdf(paste("figures/", Sys.Date(), "-figure-5.pdf", sep = ""),
    width = 3.5, height = 4.5, useDingbats = FALSE)

par(mar = c(4,5,2,1), ps = 10, tcl = -0.25, las = 1, xpd = FALSE)

plot(NULL,
     xlim = c(-0.8, 1.2), ylim = c(1, 14),
     xlab = "Change in relative abundance", ylab = "",
     xaxs = "i", yaxt = "n", mgp = c(2,0.2,0))

# axis(2, at = c(14:1), las = 1, mgp = c(2,0.5,0),
#      labels = rel_dat$genus[1:14])

axis(2, at = c(14:1), las = 1, mgp = c(2,0.5,0),
     labels = c(expression(italic("Acropora")), 
                expression(italic("Anacropora")),
                expression(italic("Cyphastrea")),
                expression(italic("Echinopora")), 
                expression(italic("Galaxea")),
                expression(italic("Goniopora")),
                expression(italic("Heliopora")), 
                expression(italic("Millepora")),
                expression(italic("Montipora")),
                expression(italic("Pavona")),
                expression(italic("Pocillopora")), 
                expression(italic("Porites")), 
                expression(italic("Seriatopora")),
                expression(italic("Turbinaria"))))

polygon(x = c(-0.03, 0.03, 0.03, -0.03),
        y = c(0, 0, 15, 15),
        col = alpha("grey", alpha = 0.5), border = NA)

segments(x0 = rel_dat$LCB0.025, y0 = length(rel_dat$genus):1,
         x1 = rel_dat$UCB0.975, y1 = length(rel_dat$genus):1,
         col = rel_dat$colour)

points(rel_dat$abund, rev(as.factor(rel_dat$genus)),
       pch = 21, bg = rel_dat$colour)

dev.off()

#

#### figure 6. nMDS core examples ####

pdf(paste("figures/", Sys.Date(), "-figure-6.pdf", sep = ""),
    width = 8, height = 4, useDingbats = FALSE)

par(mar = c(0,0,0,0), mgp = c(1.5,0.5,0), ps = 10, tcl = -0.25, las = 1, xpd = NA)

split.screen(rbind(c(0.06, 0.46, 0.15, 0.9),
                   c(0.54, 0.86, 0.6, 0.9),
                   c(0.54, 0.86, 0.15, 0.45),
                   c(0.85, 0.99, 0.15, 0.9),
                   c(0, 1, 0, 1)))

# 1. nMDS plot with highlighted examples of coral communities through time

screen(1)

plot(tax_ord, type = "n",
     xlab = "nMDS1", ylab = "nMDS2")

points(tax_ord,
       pch = 19, col = alpha("grey", alpha = 0.5))

core.nmds("Frankland_Islands.Normanby_Island.3.2",
          "orchid1", "orchid4")

core.nmds("Palm_Islands.Havannah.2.F",
          "lightpink", "red")

tax_ord$points[str_detect(row.names(tax_ord$points),
                          "Frankland_Islands.Normanby_Island.3.2"), ]

text(x = -0.6185995, y = -1.715284922,
     labels = substitute(paste(bold("B"))), pos = 1,
     col = "orchid4")

tax_ord$points[str_detect(row.names(tax_ord$points), 
                          "Palm_Islands.Havannah.2.F"), ]

text(x = 2.1263273, y = 0.33383536,
     labels = substitute(paste(bold("C"))), pos = 4,
     col = "red")

close.screen(1)

# 2. relative abundance example - top

screen(2)

core.barplot("Frankland_Islands.Normanby_Island.3.2")

close.screen(2)

# 3. relative abundance example - bottom

screen(3)

core.barplot("Palm_Islands.Havannah.2.F")

close.screen(3)

# 4. legend

screen(4)

legend("center",
       legend = c(expression(italic("Acropora")),
                  expression(italic("Porites")),
                  expression(italic("Pavona")),
                  expression(italic("Goniopora"))),
       col = c("orchid1", 
               "orchid4",
               "lightpink",
               "red"),
       pch = 19, bty = "n")

close.screen(4)

# 5. labels

screen(5)

text(x = 0.038, y = 0.97, cex = 1.5,
     labels = substitute(paste(bold("A"))))

text(x = 0.569, y = 0.97, cex = 1.5,
     labels = substitute(paste(bold("B"))))

text(x = 0.569, y = 0.48, cex = 1.5,
     labels = substitute(paste(bold("C"))))

close.screen(5)

# close and clear
close.screen(all.screens = TRUE)
dev.off()

#

#### figure S2. regional differences ####

pdf(paste("figures/", Sys.Date(), "-figure-S2.pdf", sep = ""),
    width = 8, height = 4, useDingbats = FALSE)

par(mar = c(0,0,0,0), mgp = c(1.5,0.2,0), ps = 10, tcl = -0.25, las = 1, xpd = NA)

split.screen(rbind(c(0.10, 0.475, 0.15, 0.90),
                   c(0.575, 0.95, 0.15, 0.90)))

# firstly find the identity of the novel communities
tax_nov <- qld_env$sample[which(qld_env$tax_novel == 1)]
fun_nov <- qld_env$sample[which(qld_env$fun_novel == 1)]

# and parse apart their regions
tax_nov_region <- as.factor(str_sub(tax_nov, end = 1))
fun_nov_region <- as.factor(str_sub(fun_nov, end = 1))

# A. taxonomic novelty

screen(1)

plot(tax_ord, 
     type = "n", 
     xaxt = "n", yaxt = "n",
     xlab = "nMDS1", ylab = "nMDS2")

axis(1, mgp = c(1,0.2,0))
axis(2, mgp = c(1,0.4,0))

# I'm using colour-blind safe colours from IBM Design
# accessed via: https://www.color-hex.com/color-palette/1044488

points(tax_ord,
       pch = 19, col = alpha(c("#ffb000", "#785ef0",
                               "#648fff", "#fe6100")[as.factor(sub_qld_dat$region)],
                             alpha = 0.2))

points(tax_ord$points[!is.na(match(row.names(tax_ord$points), 
                                   tax_nov)), ],
       pch = 21, col = "black", bg = c("#ffb000", "#785ef0",
                                       "#fe6100")[as.factor(str_sub(tax_nov, end = 1))])

sapply(which(tax_vec$vectors$pvals == 0.001), function(x){
  arrows(x0 = 0, y0 = 0,
         x1 = tax_ord$species[x, ][1], y1 = tax_ord$species[x, ][2],
         length = 0.08)
})

# sapply(which(tax_vec$vectors$pvals == 0.001), function(x){
#   text(x = tax_ord$species[x, ][1], y = tax_ord$species[x, ][2],
#        labels = row.names(tax_ord$species)[x], cex = 0.75)
# })

legend("topleft", bty = "n", 
       legend = substitute(paste(bold("taxonomic"))))

close.screen(1)

# B. functional novelty

screen(2)

plot(fun_ord, 
     type = "n", 
     xaxt = "n", yaxt = "n",
     xlab = "nMDS1", ylab = "nMDS2")

axis(1, mgp = c(1,0.2,0))
axis(2, mgp = c(1,0.4,0))

points(fun_ord,
       pch = 19, col = alpha(c("#ffb000", "#785ef0",
                               "#648fff", "#fe6100")[as.factor(sub_qld_dat$region)],
                             alpha = 0.2))

points(fun_ord$points[!is.na(match(row.names(fun_ord$points), 
                                   fun_nov)), ],
       pch = 21, col = "black", bg = c("#ffb000", "#785ef0",
                                       "#648fff", "#fe6100")[as.factor(str_sub(fun_nov, end = 1))])

legend("topleft", bty = "n", 
       legend = substitute(paste(bold("functional"))))

close.screen(2)

# close and clear
close.screen(all.screens = TRUE)
dev.off()

#

#