pand2020 = readRDS(url("https://github.com/TimothyStaples/novelty-cenozoic-microplankton/raw/master/outputs/processedNovelData.rds"))

pandLong = do.call("rbind", lapply(1:length(pand2020), function(x){
  temp = do.call('rbind', pand2020[[x]]$novel)
  temp$taxa = c("nano", "foram", "radio", "diatom")[x]
  return(temp)
}))

# time series ratio (length to bin count)
grainToExtent_plankton = tapply(pandLong$bins, list(pandLong$site, pandLong$taxa),
                                function(x){diff(range(as.numeric(as.character(x)))) / 0.1})
hist(grainToExtent_plankton)
summary(as.vector(grainToExtent_plankton))

# Staples 2022
Stap2022 <- readRDS(url("https://github.com/TimothyStaples/novel_comm_quaternary/raw/master/outputs/all%20neotoma%20novelty%20(sub-sampled).rds"))

StapLong <- do.call("rbind", Stap2022$novel)

grainToExtent_pollen = tapply(StapLong$bins, StapLong$site,
                              function(x){diff(range(as.numeric(as.character(x)))) / 200})
hist(grainToExtent_pollen)
summary(as.vector(grainToExtent_pollen))
