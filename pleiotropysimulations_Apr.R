## DFE pleiotropy sims

library(data.table)
library(tidyverse)
library(IRanges)
library(GenomicRanges)

## Assign DFE gamma parameters

#two_epoch
## paired
two_epoch_shape1 = 1.08
two_epoch_scale1 = 0.23

## missense
two_epoch_shape2 = 1.94329796e-01
two_epoch_scale2 = 1.21489631e+03

#growth
## paired
growth_shape1 = 0.74
growth_scale1 = 0.31

## missense
growth_shape2 = 2.35134414e-01
growth_scale2 = 5.27494215e+01

#bottlegrowth
## paired
bottlegrowth_shape1 = 0.64
bottlegrowth_scale1 = 0.38

## missense
bottlegrowth_shape2 = 2.19512301e-01
bottlegrowth_scale2 = 1.41138573e+02

#snm
## paired
snm_shape1 = 0.18
snm_scale1 = 714.6

## missense
snm_shape2 = 3.23913796e-01
snm_scale2 = 3.86390225e+03




dat <- fread("../Documents/At_cds.txt")


dat$randomfactor_ns <- runif(nrow(dat), min=0.75, max=.85) # Manipulate by Mu params
dat$adj_cdslen <- dat$`CDS length` * dat$randomfactor_ns
dat$adj_cdslen <- ceiling(dat$adj_cdslen)
dat$`CDS length` <- dat$adj_cdslen

# Find # coding nucleotides
cds_len <- sum(dat$`CDS length`)


# old, non-gene specific way
#gene <- 1:(cds_len*(2/3)) %>% as.data.frame()

num_permutes = 100

# Assign each CDS nucleotide a gene identifier
listy <- list()
i=1
for (i in i:nrow(dat)) {
  id <- dat$gene[i]
  len <- dat$`CDS length`[i]
  #vec <- rep(id, times = len)
  listy[[i]] <- rep(id, times = len)
}
gene <- unlist(listy) %>% as.data.frame()

# assign structure changing or not
#num_paired <- ()
frac_sN <- 0.0043 # From paper
num_sN <-( nrow(gene) * frac_sN ) %>% round()
num_sS <- nrow(gene) - num_sN




# Permutation loop

listy2 <- list()
sample_vec1 <- 1:num_permutes
sample_vec2 <- 1:num_permutes
sample_vec3 <- 1:num_permutes
sample_vec4 <- 1:num_permutes
sample_vec5 <- 1:num_permutes
sample_vec6 <- 1:num_permutes

j <- 1

for (j in j:num_permutes) {
#gene$paired <- c(rep("paired", times = num_sN), rep("unpaired", times = num_sS)) %>% sample()
# Structure effects
## Randomly distribute paired states
gene$paired <- c(rep("paired", times = num_sN), rep("unpaired", times = num_sS)) %>% sample()
gene_save <- gene
gene <- filter(gene,gene$paired == "paired")
v <- rep(-1.0, num_sN)
k <- 1
while (TRUE) {
  q <- rgamma(1, shape=growth_shape1, scale=growth_scale1)
  if (q > 0.0 && q < 1000) {
    v[k] <- q
    k<-k+1
    if (k>num_sN)
      break
  }
}

samples <- v

# Protein effects
v <- rep(-1.0, nrow(gene))
k <- 1
while (TRUE) {
  q <- rgamma(1, shape=growth_shape2, scale=growth_scale2)
  if (q > 0.0 && q < 1000) {
    v[k] <- q
    k<-k+1
    if (k>nrow(gene))
      break
  }
}

samples_ms <- v

gene$struc_eff <- ifelse(gene$paired == "paired",
                         samples,
                         NA)
gene$aa_eff <- samples_ms

#gene$aa_s <- gene$aa_eff / (1500*2)
#gene$struc_s <- gene$struc_eff / (1500*2)

gene$eff_diff <- gene$aa_eff - gene$struc_eff

listy2[[j]] <- gene

aa_over <- filter(gene,gene$eff_diff >0)
st_over <- filter(gene,gene$eff_diff <0)

#aa_over <- filter(aa_over,aa_over$struc_eff > 1)
#st_over <- filter(st_over,st_over$aa_eff > 1)

sample_vec1[j] <- nrow(aa_over)
sample_vec2[j] <- nrow(st_over)
sample_vec3[j] <- filter(gene,gene$eff_diff > 1)[,1] %>% sort() %>% unique() %>% length()
sample_vec4[j] <- filter(gene,gene$eff_diff < 0)[,1] %>% sort() %>% unique() %>% length()
sample_vec5[j] <- filter(gene, gene$struc_eff > 1) %>% nrow()
sample_vec6[j] <- filter(gene, gene$struc_eff > 1)[,1] %>% sort() %>% unique() %>% length()


# st_over$. %>% sort() %>% unique() %>% length()
# aa_over$. %>% sort() %>% unique() %>% length()

} ## Big permutation loop over here ##

growth_df <- cbind(rep("growth", times =num_permutes),sample_vec1,sample_vec2,sample_vec3,sample_vec4,sample_vec5,sample_vec6) %>% as.data.frame()
colnames(growth_df) <-  c("model", "aa_over", "st_over","genes_higheraa", "genes_higherstruc",  "num_nonneutral", "genes_nonneutral")

#growth_df$mean_gamma <- growth_scale1*growth_scale2

#################
### two_epoch ###
#################

# Permutation loop

listy2_two_epoch <- list()
two_epoch_vec1 <- 1:num_permutes
two_epoch_vec2 <- 1:num_permutes
two_epoch_vec3 <- 1:num_permutes
two_epoch_vec4 <- 1:num_permutes
two_epoch_vec5 <- 1:num_permutes
two_epoch_vec6 <- 1:num_permutes

j <- 1

for (j in j:num_permutes) {
  ## Randomly distribute paired states
  gene$paired <- c(rep("paired", times = num_sN), rep("unpaired", times = num_sS)) %>% sample()
  gene_save <- gene
  gene <- filter(gene,gene$paired == "paired")
  v <- rep(-1.0, num_sN)
  k <- 1
  while (TRUE) {
    q <- rgamma(1, shape=two_epoch_shape1, scale=two_epoch_scale1)
    if (q > 0.0 && q < 1000) {
      v[k] <- q
      k<-k+1
      if (k>num_sN)
        break
    }
  }
  
  samples <- v
  
  # Protein effects
  v <- rep(-1.0, nrow(gene))
  k <- 1
  while (TRUE) {
    q <- rgamma(1, shape=two_epoch_shape2, scale=two_epoch_scale2)
    if (q > 0.0 && q < 1000) {
      v[k] <- q
      k<-k+1
      if (k>nrow(gene))
        break
    }
  }
  
  samples_ms <- v
  
  gene$struc_eff <- ifelse(gene$paired == "paired",
                           samples,
                           NA)
  gene$aa_eff <- samples_ms
  
  #gene$aa_s <- gene$aa_eff / (1500*2)
  #gene$struc_s <- gene$struc_eff / (1500*2)
  
  gene$eff_diff <- gene$aa_eff - gene$struc_eff
  
  listy2_two_epoch[[j]] <- gene
  
  aa_over <- filter(gene,gene$eff_diff >0)
  st_over <- filter(gene,gene$eff_diff <0)
  
  #aa_over <- filter(aa_over,aa_over$struc_eff > 1)
  #st_over <- filter(st_over,st_over$aa_eff > 1)
  
  two_epoch_vec1[j] <- nrow(aa_over)
  two_epoch_vec2[j] <- nrow(st_over)
  two_epoch_vec3[j] <- filter(gene,gene$eff_diff > 1 )[,1] %>% sort() %>% unique() %>% length()
  two_epoch_vec4[j] <- filter(gene,gene$eff_diff < 0 )[,1] %>% sort() %>% unique() %>% length()
  two_epoch_vec5[j] <- filter(gene, gene$struc_eff > 1) %>% nrow()
  two_epoch_vec6[j] <- filter(gene, gene$struc_eff > 1)[,1] %>% sort() %>% unique() %>% length()
  
  
  # st_over$. %>% sort() %>% unique() %>% length()
  # aa_over$. %>% sort() %>% unique() %>% length()
  
} ## Big permutation loop over here ##

two_epoch_df <- cbind(rep("two_epoch", times =num_permutes),two_epoch_vec1,two_epoch_vec2,two_epoch_vec3,two_epoch_vec4, two_epoch_vec5, two_epoch_vec6) %>% as.data.frame()
colnames(two_epoch_df) <-  c("model", "aa_over", "st_over","genes_higheraa", "genes_higherstruc",  "num_nonneutral", "genes_nonneutral")


#################
### bottlegrowth ###
#################

# Permutation loop

listy2_bottlegrowth <- list()
bottlegrowth_vec1 <- 1:num_permutes
bottlegrowth_vec2 <- 1:num_permutes
bottlegrowth_vec3 <- 1:num_permutes
bottlegrowth_vec4 <- 1:num_permutes
bottlegrowth_vec5 <- 1:num_permutes
bottlegrowth_vec6 <- 1:num_permutes

j <- 1

for (j in j:num_permutes) {
  ## Randomly distribute paired states
  gene$paired <- c(rep("paired", times = num_sN), rep("unpaired", times = num_sS)) %>% sample()
  gene_save <- gene
  gene <- filter(gene,gene$paired == "paired")
  v <- rep(-1.0, num_sN)
  k <- 1
  while (TRUE) {
    q <- rgamma(1, shape=bottlegrowth_shape1, scale=bottlegrowth_scale1)
    if (q > 0.0 && q < 1000) {
      v[k] <- q
      k<-k+1
      if (k>num_sN)
        break
    }
  }
  
  samples <- v
  
  # Protein effects
  v <- rep(-1.0, nrow(gene))
  k <- 1
  while (TRUE) {
    q <- rgamma(1, shape=bottlegrowth_shape2, scale=bottlegrowth_scale2)
    if (q > 0.0 && q < 1000) {
      v[k] <- q
      k<-k+1
      if (k>nrow(gene))
        break
    }
  }
  
  samples_ms <- v
  
  gene$struc_eff <- ifelse(gene$paired == "paired",
                           samples,
                           NA)
  gene$aa_eff <- samples_ms
  
  #gene$aa_s <- gene$aa_eff / (1500*2)
  #gene$struc_s <- gene$struc_eff / (1500*2)
  
  gene$eff_diff <- gene$aa_eff - gene$struc_eff
  
  listy2_bottlegrowth[[j]] <- gene
  
  aa_over <- filter(gene,gene$eff_diff >0)
  st_over <- filter(gene,gene$eff_diff <0)
  
  #aa_over <- filter(aa_over,aa_over$struc_eff > 1)
  #st_over <- filter(st_over,st_over$aa_eff > 1)
  
  bottlegrowth_vec1[j] <- nrow(aa_over)
  bottlegrowth_vec2[j] <- nrow(st_over)
  bottlegrowth_vec3[j] <- filter(gene,gene$eff_diff > 1)[,1] %>% sort() %>% unique() %>% length()
  bottlegrowth_vec4[j] <- filter(gene,gene$eff_diff < 0)[,1] %>% sort() %>% unique() %>% length()
  bottlegrowth_vec5[j] <- filter(gene, gene$struc_eff > 1) %>% nrow()
  bottlegrowth_vec6[j] <- filter(gene, gene$struc_eff > 1)[,1] %>% sort() %>% unique() %>% length()
  
  
  # st_over$. %>% sort() %>% unique() %>% length()
  # aa_over$. %>% sort() %>% unique() %>% length()
  
} ## Big permutation loop over here ##

bottlegrowth_df <- cbind(rep("bottlegrowth", times =num_permutes),bottlegrowth_vec1,bottlegrowth_vec2,bottlegrowth_vec3,bottlegrowth_vec4, bottlegrowth_vec5, bottlegrowth_vec6) %>% as.data.frame()
colnames(bottlegrowth_df) <-  c("model", "aa_over", "st_over","genes_higheraa", "genes_higherstruc",  "num_nonneutral", "genes_nonneutral")


#################
### snm ###
#################

# Permutation loop

listy2_snm <- list()
snm_vec1 <- 1:num_permutes
snm_vec2 <- 1:num_permutes
snm_vec3 <- 1:num_permutes
snm_vec4 <- 1:num_permutes
snm_vec5 <- 1:num_permutes
snm_vec6 <- 1:num_permutes

j <- 1

for (j in j:num_permutes) {
  ## Randomly distribute paired states
  gene$paired <- c(rep("paired", times = num_sN), rep("unpaired", times = num_sS)) %>% sample()
  gene_save <- gene
  gene <- filter(gene,gene$paired == "paired")
  v <- rep(-1.0, num_sN)
  k <- 1
  while (TRUE) {
    q <- rgamma(1, shape=snm_shape1, scale=snm_scale1)
    if (q > 0.0 && q < 1000) {
      v[k] <- q
      k<-k+1
      if (k>num_sN)
        break
    }
  }
  
  samples <- v
  
  # Protein effects
  v <- rep(-1.0, nrow(gene))
  k <- 1
  while (TRUE) {
    q <- rgamma(1, shape=snm_shape2, scale=snm_scale2)
    if (q > 0.0 && q < 1000) {
      v[k] <- q
      k<-k+1
      if (k>nrow(gene))
        break
    }
  }
  
  samples_ms <- v
  
  gene$struc_eff <- ifelse(gene$paired == "paired",
                           samples,
                           NA)
  gene$aa_eff <- samples_ms
  
  #gene$aa_s <- gene$aa_eff / (1500*2)
  #gene$struc_s <- gene$struc_eff / (1500*2)
  
  gene$eff_diff <- gene$aa_eff - gene$struc_eff
  
  listy2_snm[[j]] <- gene
  
  aa_over <- filter(gene,gene$eff_diff >0)
  st_over <- filter(gene,gene$eff_diff <0)
  
  #aa_over <- filter(aa_over,aa_over$struc_eff > 1)
  #st_over <- filter(st_over,st_over$aa_eff > 1)
  
  snm_vec1[j] <- nrow(aa_over)
  snm_vec2[j] <- nrow(st_over)
  snm_vec3[j] <- filter(gene,gene$eff_diff > 1)[,1] %>% sort() %>% unique() %>% length()
  snm_vec4[j] <- filter(gene,gene$eff_diff < 0)[,1] %>% sort() %>% unique() %>% length()
  snm_vec5[j] <- filter(gene, gene$struc_eff > 1) %>% nrow()
  snm_vec6[j] <- filter(gene, gene$struc_eff > 1)[,1] %>% sort() %>% unique() %>% length()
  
  
  # st_over$. %>% sort() %>% unique() %>% length()
  # aa_over$. %>% sort() %>% unique() %>% length()
  
} ## Big permutation loop over here ##

snm_df <- cbind(rep("snm", times =num_permutes),snm_vec1,snm_vec2,snm_vec3,snm_vec4, snm_vec5, snm_vec6) %>% as.data.frame()
colnames(snm_df) <-  c("model", "aa_over", "st_over","genes_higheraa", "genes_higherstruc",  "num_nonneutral", "genes_nonneutral")




all_df <- rbind(growth_df, two_epoch_df, bottlegrowth_df, snm_df)
all_df$genes_higherstruc <- all_df$genes_higherstruc %>% as.numeric()
all_df$genes_higheraa <- all_df$genes_higheraa %>% as.numeric()
all_df$num_nonneutral <- as.numeric(all_df$num_nonneutral)
all_df$st_over <- as.numeric(all_df$st_over)

all_df <- filter(all_df, all_df$model != "snm")

pl_al <- ggplot(data = all_df, aes(x = model, y = (st_over), col = model)) +
  geom_boxplot(outlier.shape = NA) +
  theme_bw() +
  xlab("Demographic model") +
  ylab("Number of pleitropic alleles") + theme(legend.title=element_blank()) +
  theme(legend.position="none")
  #theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggsave("numnplei_boxplot.png", device = "png", width = 5, height =4)

ggplot(data = all_df, aes(x = model, y = (num_nonneutral), col = model)) +
  geom_boxplot() +
  theme_bw() +
  xlab("Demographic model") +
  ylab("Number of non-neutral alleles")
ggsave("numnoneutral_boxplot_ten.png", device = "png", width = 5, height =4)

all_df_longer <- pivot_longer(all_df, cols = c("genes_higherstruc", "genes_higheraa"), names_to = "genetype")
all_df_longer$'Pleiotropy type' <- ifelse(all_df_longer$genetype == "genes_higheraa", "Stronger protein effect", "Stronger RNA effect")
pl_gene <- ggplot(data = all_df_longer, aes(x = model, y = value, col = all_df_longer$'Pleiotropy type')) +
  geom_boxplot(outlier.shape = NA) +
  theme_bw() +
  xlab("Demographic model") +
  ylab("Number of genes") +
  scale_color_brewer(type = "qual", palette = 3) + 
  theme(legend.title=element_blank())  +
  theme(legend.position="none") + theme(axis.text = element_text(size = 14), 
                                      axis.title = element_text(size = 14))

pl_al <- ggplot(data = all_df, aes(x = model, y = (st_over), col = model)) +
  geom_boxplot(outlier.shape = NA) +
  theme_bw() +
  xlab("Demographic model") +
  ylab("Number of pleitropic alleles") + theme(legend.title=element_blank()) +
  theme(legend.position="none") +theme(axis.text = element_text(size = 14), 
                                       axis.title = element_text(size = 14))

ggsave("pl_al.png", pl_al, device = "png", dpi = 300, height = 3, width = 3)
ggsave("pl_al.png", pl_al, device = "png", dpi = 300, height = 5, width = 5)


#ggsave("numgeneswithpleiotropy1_boxplot_ten.png", device = "png", width = 5)

library(gridExtra)
gridded <- grid.arrange(pl_al, pl_gene, ncol = 1)
ggsave(plot = gridded,filename = "nsimsum_pl_gp.png", 
       device = "png", 
       height = 5,width = 3, dpi = 300, scale = 1.5)

library(grid)
grid.newpage()
grid.draw(rbind(ggplotGrob(pl_al), ggplotGrob(pl_gene), size = "last"))
ggsave("pl_al_gene.png", device = "png", scale = 4)

png("pl_al_gene.png", height = 1440, width = 1440)

grid.draw(rbind(ggplotGrob(pl_al), ggplotGrob(pl_gene), size = "last")) # print it

dev.off()

