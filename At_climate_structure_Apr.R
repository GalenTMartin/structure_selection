library(data.table)
library(tidyverse)
#install.packages(c("sf", "rnaturalearth", "rnaturalearthdata"))
library("rnaturalearth")
library("rnaturalearthdata")
library("sf")
library(ggpubr)
library(ggpmisc)

vcf <- fread("unpaired_5UTR.vcf.gz") # Change diff SNP types

coords <- fread("../Downloads/A. thaliana Master Accession List - Sheet1.tsv")
tots <- fread("accession_totals.txt")
colnames(vcf) <- c("#CHROM",  
                   "POS",     
                   "ID",
                   "REF",
                   "ALT",
                   "QUAL",
                   "FILTER",
                   "INFO",
                   "FORMAT",
                   tots$V1)
colnames(tots) <- c("ID", "num_paired_snps")

coords <- merge(coords, tots, by.x = "id", by.y = "ID")

maxlat = max(coords$latitude, na.rm = T)
minlat = min(coords$latitude, na.rm = T)
maxlong = max(coords$longitude, na.rm = T)
minlong = min(coords$longitude, na.rm = T)

latdis <- maxlat - minlat
longdis <- maxlong - minlong

res <- longdis / latdis

nlat <- 15
nlong <- floor(nlat*5)

#latbreaks <- seq(minlat, maxlat, length.out = nlat)
#longbreaks <- seq(minlong, maxlong, length.out = nlong)
#latbreaks <- rep(latbreaks, each=5)

#breaks <- data.frame(longbreaks, latbreaks)
#colnames(breaks) <- c("longitude", "latitude")
library(raster)
library(sp)
library(geodata)
#r <- geodata::getdata_path("worldclim",var="bio",res=10)
r <- worldclim_global(path = "worldclim", var="bio",res=10)
coords <- filter(coords, !is.na(coords$latitude))
points <- SpatialPoints(coords[,5:6])

#clim<-raster(choose.files())

#Get WorldClim
values <- list()
i <- 1
for (i in i:nrow(coords)) {
  skip_to_next <- FALSE
  
  # Note that print(b) fails since b doesn't exist
  
  values[[i]] <- tryCatch(extract(r,coords[,6:5][i], 
                   method = 'simple'), error = function(e) { skip_to_next <<- TRUE})
  
  if(skip_to_next) { next }     
  
}


test <- lapply(values, as.data.frame)
test <- bind_rows(test)
coords$BIO1  <- test$wc2.1_10m_bio_1 
coords$BIO2  <- test$wc2.1_10m_bio_2 
coords$BIO3  <- test$wc2.1_10m_bio_3 
coords$BIO4  <- test$wc2.1_10m_bio_4 
coords$BIO5  <- test$wc2.1_10m_bio_5 
coords$BIO6  <- test$wc2.1_10m_bio_6 
coords$BIO7  <- test$wc2.1_10m_bio_7 
coords$BIO8  <- test$wc2.1_10m_bio_8 
coords$BIO9  <- test$wc2.1_10m_bio_9 
coords$BIO10 <- test$wc2.1_10m_bio_10
coords$BIO11 <- test$wc2.1_10m_bio_11
coords$BIO12 <- test$wc2.1_10m_bio_12
coords$BIO13 <- test$wc2.1_10m_bio_13
coords$BIO14 <- test$wc2.1_10m_bio_14
coords$BIO15 <- test$wc2.1_10m_bio_15
coords$BIO16 <- test$wc2.1_10m_bio_16
coords$BIO17 <- test$wc2.1_10m_bio_17
coords$BIO18 <- test$wc2.1_10m_bio_18
coords$BIO19 <- test$wc2.1_10m_bio_19



griddf <- expand.grid(latcoords = seq(minlong, maxlong, length.out = nlong),
                      lngcoords = seq(minlong, maxlong, length.out = nlong))

changelat <- griddf$latcoords[[2]] - griddf$latcoords[[1]]

changelng <- 1.0423

griddf$latnext <- griddf$latcoords+ changelat
griddf$lngnext <- griddf$lngcoords+ changelat

griddf <- filter(griddf, between(griddf$latcoords, left = minlat, right = maxlat) )


snpeff <- fread("../Downloads/snpeff_simpler.txt")

colnames(snpeff) <- c("CHR", "POS", "REF", "ALT", "EFF")

snpeff <- filter(snpeff, snpeff$EFF == "synonymous_variant"
                 )
vcf <- filter(vcf, paste(vcf$`#CHROM`, vcf$POS, sep = "_") %in% paste(snpeff$CHR, snpeff$POS, sep = "_"))

vcf <- sapply(vcf[,!1:9], function(x) as.numeric(substr(x, 1,1)) )
gc()

#Remove extras
coords <- filter(coords, coords$id %in% colnames(vcf))

listy <- list()
griddf$freq <- NA
griddf$num_in <- NA
griddf$BIO1 <- NA
griddf$BIO2 <- NA
griddf$BIO3 <- NA
griddf$BIO4 <- NA
griddf$BIO5 <- NA
griddf$BIO6 <- NA
i <- 1
for (i in i:nrow(griddf)) {
  subcoords <- filter(coords, 
                     coords$latitude  >= griddf$latcoords[[i]] &
                     coords$latitude  <= griddf$latnext[[i]] &
                     coords$longitude >= griddf$lngcoords[[i]] &
                     coords$longitude <= griddf$lngnext[[i]]
                   )
  subpop <- as.character(subcoords$id)
  num_in <- length(subpop)
  if (num_in > 1) {
  vcf_sub <- subset(vcf_alt,select = c(subpop)) %>% as.data.frame()
  vcf_sub$sum <- NA
  for (c in num_in-1) {
    vcf_sub$sum <- (vcf_sub[,c] + vcf_sub[,(c+1)]) / num_in
  }
  griddf$freq[[i]] <- sum(vcf_sub$sum, na.rm = TRUE)  / nrow(vcf_sub)
  griddf$num_in[[i]] <- num_in
  griddf$BIO1[[i]] <- mean(subcoords$BIO1)
  griddf$BIO2[[i]] <- mean(subcoords$BIO2)
  griddf$BIO3[[i]] <- mean(subcoords$BIO3)
  griddf$BIO4[[i]] <- mean(subcoords$BIO4)
  griddf$BIO5[[i]] <- mean(subcoords$BIO5)
  griddf$BIO6[[i]] <- mean(subcoords$BIO6)
  } else {
    griddf$freq[[i]] <- NA
    griddf$num_in[[i]] <- num_in
    griddf$BIO1[[i]] <- NA
    griddf$BIO2[[i]] <- NA
    griddf$BIO3[[i]] <- NA
    griddf$BIO4[[i]] <- NA
    griddf$BIO5[[i]] <- NA
    griddf$BIO6[[i]] <- NA
  }
}

#write_delim(griddf, "At_pairedSNPs_n15_mapdata_Aug14_summed.txt", delim = "\t")
#griddf <- fread("At_pairedSNPs_n15_mapdata_june12_summed.txt")
ggplot(data = na.omit(griddf), aes(
  #x = latcoords, 
   #                       y = lngcoords,
                          fill = freq
                          )) +
  geom_rect(aes(xmin=lngcoords, 
                ymin=latcoords,
                xmax=lngnext,
                ymax=latnext)) +
  theme_classic() +
  #geom_raster() +
  scale_fill_viridis_c()
  #geom_map(map = )

world <- ne_countries(scale = "medium", returnclass = "sf")
ggplot(data = world) + 
  geom_sf(fill = "white",
          col = "lightgrey",
          ) + 
  theme_classic() +
  #geom_point(data = df, aes(long, lat)) +
  geom_rect(data = na.omit(griddf[griddf$num_in > 2,]), 
            aes(xmin=lngcoords, 
                ymin=latcoords,
                xmax=lngnext,
                ymax=latnext, fill = freq),
            alpha = 0.9) +
  scale_fill_viridis_c() +
  xlim(-25, 100) +
  ylim(20, 90) +
  geom_sf(fill = alpha("#2C77BF", 0), col = "grey",
  )
ggsave("map_paired_SNP_n15_frequency_summed_Aug14_Aug14.png", 
       device = "png", 
       scale = 3)

ggplot(data = world) + 
  geom_sf(fill = "white",
          col = "lightgrey",
  ) + 
  theme_classic() +
  #geom_point(data = df, aes(long, lat)) +
  geom_rect(data = na.omit(griddf[griddf$num_in > 2,]), 
            aes(xmin=lngcoords, 
                ymin=latcoords,
                xmax=lngnext,
                ymax=latnext, fill =freq),
            alpha = 0.9) +
  scale_fill_viridis_c() +
  xlim(-120, 130) +
  ylim(20, 90) +
  geom_sf(fill = alpha("#2C77BF", 0), col = "black",
  )

has <- na.omit(griddf)
has <- filter(has,has$num_in >= 5)

formula2 <- y ~ log(x)
formula1 <- y ~ poly(x, 2, raw = TRUE)
formula_line <- y ~ x


# By admixutre

groups <- fread("1001genomes_admixturegroups.txt")
groups <- groups[,c("V1", "V11")]
colnames(groups) <- c("id", "group")

coords <- merge(coords, groups, by = "id")

ad <- unique(coords$group)
ad <- data.frame(ad)



listy <- list()
ad$freq <- NA
ad$num_in <- NA




ad$BIO1 <- NA
ad$BIO2 <- NA
ad$BIO3 <- NA
ad$BIO4 <- NA
ad$BIO5 <- NA
ad$BIO6 <- NA
ad$BIO7 <- NA
ad$BIO8 <- NA
ad$BIO9 <- NA
ad$BIO10 <- NA
ad$BIO11 <- NA
ad$BIO12 <- NA
ad$BIO13 <- NA
ad$BIO14 <- NA
ad$BIO15 <- NA
ad$BIO16 <- NA
ad$BIO17 <- NA
ad$BIO18 <- NA
ad$BIO19 <- NA

i <- 1
for (i in i:nrow(ad)) {
  gc()
  subcoords <- filter(coords, coords$group == ad$ad[[i]])
  subpop <- as.character(subcoords$id)
  num_in <- length(subpop)
  if (num_in > 1) {
    vcf_sub <- subset(vcf,select = c(subpop)) %>% as.data.frame()
    vcf_sub$sum <- NA
    for (c in num_in-1) {
      vcf_sub$sum <- (vcf_sub[,c] + vcf_sub[,(c+1)]) / num_in
    }
    ad$freq[[i]] <- mean(vcf_sub$sum, na.rm = TRUE) 
    ad$num_in[[i]] <- num_in
    ad$BIO1[[i]] <- mean(na.omit(subcoords$BIO1 ))
    ad$BIO2[[i]] <- mean(na.omit(subcoords$BIO2 ))
    ad$BIO3[[i]] <- mean(na.omit(subcoords$BIO3 ))
    ad$BIO4[[i]] <- mean(na.omit(subcoords$BIO4 ))
    ad$BIO5[[i]] <- mean(na.omit(subcoords$BIO5 ))
    ad$BIO6[[i]] <- mean(na.omit(subcoords$BIO6 ))
    ad$BIO7[[i]] <- mean(na.omit(subcoords$BIO7 ))
    ad$BIO8[[i]] <- mean(na.omit(subcoords$BIO8 ))
    ad$BIO9[[i]] <- mean(na.omit(subcoords$BIO9 ))
    ad$BIO10[[i]] <- mean(na.omit(subcoords$BIO10))
    ad$BIO11[[i]] <- mean(na.omit(subcoords$BIO11))
    ad$BIO12[[i]] <- mean(na.omit(subcoords$BIO12))
    ad$BIO13[[i]] <- mean(na.omit(subcoords$BIO13))
    ad$BIO14[[i]] <- mean(na.omit(subcoords$BIO14))
    ad$BIO15[[i]] <- mean(na.omit(subcoords$BIO15))
    ad$BIO16[[i]] <- mean(na.omit(subcoords$BIO16))
    ad$BIO17[[i]] <- mean(na.omit(subcoords$BIO17))
    ad$BIO18[[i]] <- mean(na.omit(subcoords$BIO18))
    ad$BIO19[[i]] <- mean(na.omit(subcoords$BIO19))
  } else {
    ad$freq[[i]] <- NA
    ad$num_in[[i]] <- num_in
    
  }
}

library(ggpmisc)
ggplot(data = ad, aes(x = BIO1, y = freq)) +
  geom_smooth(method = "lm", col = "black") + stat_poly_eq(formula = y~ x,
                                            aes(label = paste(..eq.label..,..rr.label..,..p.value.., sep = "*`,`~")),
                                            parse = TRUE,
                                            label.x.npc = "middle",
                                            size =2,
                                            vstep =0.05)+
  scale_color_manual(values=c("#E69F00FF", "#F0E442FF", "#009E73FF", "#56B4E9FF", "#D55E00FF", "#CC79A7FF", "#0072B2FF", "#7570B3FF", "#000000FF")) +
  geom_point(col = "white", size = 2.5) +
  geom_point(aes(col = ad), size = 2) +
  theme_bw()
ggsave("freqbyBIO1_paired_SNP_admixgroup.png", device = "png", width = 5, height = 3) 

ad_long <- pivot_longer(ad, cols = colnames(ad[,4:22]), names_to = "BIO")
ad_long$BIO <- factor(ad_long$BIO, levels = colnames(ad[,4:22]))
ggplot(data = ad_long, aes(x = value, y = freq)) +
  geom_smooth(method = "lm", col = "black") + stat_poly_eq(formula = y~ x,
                                                           aes(label = paste(..eq.label..,..rr.label..,..p.value.., sep = "*`,`~")),
                                                           parse = TRUE,
                                                           label.x.npc = "middle",
                                                           size =2,
                                                           vstep =0.05)+
  scale_color_manual(values=c("#E69F00FF", "#F0E442FF", "#009E73FF", "#56B4E9FF", "#D55E00FF", "#CC79A7FF", "#0072B2FF", "#7570B3FF", "#000000FF")) +
  geom_point(col = "white", size = 2.5) +
  geom_point(aes(col = ad), size = 2) +
  theme_bw() +
  facet_wrap(BIO ~ ., ncol = 4, scales = "free")
ggsave("allBIO_lm_upM_UTR_Apr.png", device = "png", width = 6, height = 5, scale = 2)


ggplot(data = world) + 
  geom_sf(fill = "white",
          col = "lightgrey",
  ) + 
  theme_classic() +
  #geom_point(data = df, aes(long, lat)) +
  geom_point(data = new, 
            aes(x=longitude, 
                y=latitude,
                col = freq
                #alpha = freq
                ),
            size = 0.5) +
  scale_color_viridis_c() +
  xlim(-25, 100) +
  ylim(20, 90)
ggsave("world_admixture_sN_freq.png", device = "png", width = 7, height = 4)

