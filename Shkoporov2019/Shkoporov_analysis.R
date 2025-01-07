library(tidyverse)
library(phyloseq)
library(vegan)
library(reshape2)

library(microshades)
library(cowplot)
library(gridExtra)
library(grid)
library(viridis)
library(ggpubr)


ps <- readRDS("phyloseq_vir_raw.RDS") # import phyloseq object

ps <- ps %>% subset_samples(Subject != "GEND") # GEND samples don't have metadata, remove them
ps.relabund <- transform_sample_counts(ps, function(x) x / sum(x) )
ps.melted <- ps.relabund %>% psmelt() %>% 
  select(OTU, Sample, Abundance, Subject, Time_point) # melt ps object, keeping some metadata only, takes a few mins

iphop_out <- read_csv("Host_prediction_to_genome_m90-merged.csv") # import iphop genome-level results

iphop_out.dedup <- iphop_out %>% group_by(Virus) %>% arrange(desc(`Confidence score`)) %>% slice(1) # pick most confident hit
iphop_out.dedup <- iphop_out.dedup[, c("Virus", "Host taxonomy", "Host genome", "Confidence score")] # keep relevant columns

# split and cleanup taxonomic ranks
iphop_out.split <- iphop_out.dedup %>% 
  separate(col = `Host taxonomy`, into = c("domain", "phylum", "class", "order", "family", "genus", "species"), sep = ";")

iphop_out.split$domain <- gsub("d__", "", iphop_out.split$domain)
iphop_out.split$phylum <- gsub("p__", "", iphop_out.split$phylum)
iphop_out.split$class <- gsub("c__", "", iphop_out.split$class)
iphop_out.split$order <- gsub("o__", "", iphop_out.split$order)
iphop_out.split$family <- gsub("f__", "", iphop_out.split$family)
iphop_out.split$genus <- gsub("g__", "", iphop_out.split$genus)
iphop_out.split$species <- gsub("s__", "", iphop_out.split$species)

iphop_out.split$phylum <- sapply(strsplit(iphop_out.split$phylum, "_"), `[`, 1)
iphop_out.split$class <- sapply(strsplit(iphop_out.split$class, "_"), `[`, 1)
iphop_out.split$order <- sapply(strsplit(iphop_out.split$order, "_"), `[`, 1)
iphop_out.split$genus <- sapply(strsplit(iphop_out.split$genus, "_"), `[`, 1)
iphop_out.split$species <- sapply(strsplit(iphop_out.split$species, "_"), `[`, 1)

ps.melted.host <- ps.melted %>% 
  left_join(iphop_out.split, by = c("OTU" = "Virus")) # merge iphop data with melted ps data

#FIGURE 2A - plot relative abundance of phages at PHF level
color_objs <- create_color_dfs(as.data.frame(ps.melted.host),
                                   selected_groups = c("Desulfobacterota", "Proteobacteria", "Actinobacteriota",
                                                       "Bacteroidota", "Firmicutes"),
                                   group_level = "phylum",
                                   subgroup_level = "family",
                                   top_orientation = TRUE)
mdf <- color_objs$mdf
cdf <- color_objs$cdf
legend <- custom_legend(mdf, cdf, group_level = "phylum", subgroup_level = "family", legend_key_size = 0.6, legend_text_size = 14)
plot_diff <- plot_microshades(mdf, cdf, x = "Sample") + 
  scale_y_continuous(labels = scales::percent, expand = expansion(0)) +
  facet_grid(. ~ Subject, scales = "free", space = "free") +
  theme(axis.text.x = element_text(size= 10)) +
  theme(legend.position = "none") +
  theme(axis.text.x=element_blank())
fig2a <- plot_grid(plot_diff, arrangeGrob(legend, nullGrob(), ncol = 1, heights = c(0.80, 0.20)),  rel_widths = c(1, .25))


# calculate percent unknown host per sample
ps.melted.host.summary <- ps.melted.host %>% 
  group_by(Sample, domain) %>%
  summarize(relabund = 100*sum(Abundance))

ps.melted.host.summary[is.na(ps.melted.host.summary)] <- "Unknown"

ps.melted.host.summary <- ps.melted.host.summary %>%
  filter(domain == "Unknown") %>%
  select(Sample, relabund)

colnames(ps.melted.host.summary) <- c("Sample", "unknown")

# make OTU table for phyloseq at predicted host genome level
ps.melt.host.filt <- ps.melted.host %>% 
  ungroup() %>% 
  group_by(Sample, `Host genome`) %>% 
  summarize(relabund = sum(Abundance)) %>% 
  filter(!is.na(`Host genome`)) %>% 
  filter(relabund > 0)

votu.table.genome <- as.data.frame(pivot_wider(ps.melt.host.filt[, c("Sample", "relabund", "Host genome")], id_cols = "Host genome", names_from = "Sample", values_from = "relabund"))

rownames(votu.table.genome) <- votu.table.genome$`Host genome`
votu.table.genome <- votu.table.genome[, -1]
votu.table.genome[is.na(votu.table.genome)] <- 0

PHAGE_OTU.genome <- otu_table(votu.table.genome, taxa_are_rows = T)

# update metadata with % unknown host info for each sample
metadata.orig <- sample_data(ps.relabund)
class(metadata.orig) <- "data.frame"
metadata.orig$Sample <- rownames(metadata.orig)
metadata.joined <- left_join(metadata.orig, ps.melted.host.summary, by = c("Sample" = "Sample"))
rownames(metadata.joined) <- metadata.joined$Sample
sample_data(ps.relabund) <- sample_data(metadata.joined)

PHAGE_METADATA <- sample_data(metadata.joined)

# prepare phyloseq taxtable from iphop output
taxtable <- iphop_out.split %>% ungroup() %>% select(domain, phylum, class, order, family, genus, species, `Host genome`) %>% distinct() %>% as.data.frame()
rownames(taxtable) <- taxtable$`Host genome`
taxtable <- taxtable %>% select(domain, phylum, class, order, family, genus, species)

TAXONOMY <- tax_table(as.matrix(taxtable))

# import gtdb tree
TREE.gtdb.r202 <- read_tree("bac120_r202.tree")

# make new host taxonomy aware phyloseq object
ps.relabund.gtdb <- phyloseq(PHAGE_OTU.genome, PHAGE_METADATA, TREE.gtdb.r202, TAXONOMY)


# filter for samples with % unknown <30%
ps.relabund.gtdb.filt <- ps.relabund.gtdb %>% subset_samples(unknown < 30)

# calculate distances at votu level
dist.bray.filt <- as.matrix(vegdist(t(otu_table(ps.relabund.gtdb.filt)), method = "bray")) %>% melt() %>% filter(Var1 != Var2)
colnames(dist.bray.filt) <- c("sample1", "sample2", "dist.bray.filt")

#tax glom at family level
ps.relabund.gtdb.filt.family <- ps.relabund.gtdb.filt %>% tax_glom(taxrank = "family")

# calculate distanced at phf level
dist.bray.filt.phf <- as.matrix(distance(ps.relabund.gtdb.filt.family, method = "bray")) %>% melt() %>% filter(Var1 != Var2)
colnames(dist.bray.filt.phf) <- c("sample1", "sample2", "dist.bray.filt.phf")

dist.wuni.filt.phf <- as.matrix(distance(ps.relabund.gtdb.filt.family, method = "wunifrac")) %>% melt() %>% filter(Var1 != Var2)
colnames(dist.wuni.filt.phf) <- c("sample1", "sample2", "dist.wuni.filt.phf")

dist.merged.all.filt <- left_join(dist.bray.filt %>% select(sample1, sample2, dist.bray.filt), 
                                  dist.wuni.filt.phf %>% select(sample1, sample2, dist.wuni.filt.phf)) %>% 
  left_join(dist.bray.filt.phf %>% select(sample1, sample2, dist.bray.filt.phf)) %>% 
  left_join(metadata.joined %>% select(Subject, Sample), by = c("sample1" = "Sample")) %>% 
  left_join(metadata.joined %>% select(Subject, Sample), by = c("sample2" = "Sample")) %>% 
  melt()

colnames(dist.merged.all.filt) <- c("sample1", "sample2", "subject1", "subject2", "metric", "value")

dist.merged.all.filt$metric <- factor(dist.merged.all.filt$metric, levels = c("dist.bray.filt", "dist.bray.filt.phf", "dist.wuni.filt.phf"))

dist.merged.all.filt <- dist.merged.all.filt %>% mutate(same.subject = subject1 == subject2)

# FIGURE 2B - distances grouped by comparison type for all 3 metrics
labels <- c("Bray-Curtis (contig)", "Bray-Curtis (PHF)", "Weighted UniFrac (PHF)")
fig2b <- dist.merged.all.filt %>% filter(!metric %in% c("dist.aitch.filt", "dist.aitch.filt.pbf")) %>% ggplot(aes(x = same.subject, y = value, color = metric)) +
  geom_violin(scale = "width", aes(fill=metric), alpha = 0.3, position = position_dodge(width = 0.9)) +
  geom_boxplot(outlier.shape = NA, width = 0.3, aes(fill=metric), color="black", position = position_dodge(width = 0.9)) +
  scale_fill_viridis(discrete = T, begin = 0.1, end = 0.7, labels = labels) +
  scale_colour_viridis(discrete = T, begin = 0.1, end = 0.7, labels = labels) +
  scale_x_discrete(labels = c("FALSE" = "inter-individual (across)", "TRUE"  = "intra-individual (within)")) +
  scale_y_continuous(breaks = c(0.25, 0.5, 0.75, 1.0), limits = c(0, 1.1)) +
  theme_bw() +
  theme(legend.position = "none") +
  labs(x = "Comparison", y = "Distance", color = "Distance metric", fill = "Distance metric") +
  annotate("segment", x = 0.66, xend = 0.96, y = 1.03, yend = 1.03, 
           linetype = "solid", color = "black", size = 0.75) +
  annotate("text", x = 0.82, y = 1.04, label = "****", size = 5, color = "black") +
  annotate("segment", x = 0.99, xend = 1.29, y = 1.03, yend = 1.03, 
           linetype = "solid", color = "black", size = 0.75) +
  annotate("text", x = 1.15, y = 1.04, label = "****", size = 5, color = "black") +
  annotate("segment", x = 0.66, xend = 1.29, y = 1.08, yend = 1.08, 
           linetype = "solid", color = "black", size = 0.75) +
  annotate("text", x = 0.98, y = 1.09, label = "****", size = 5, color = "black") +

  annotate("segment", x = 1.66, xend = 1.96, y = 1.03, yend = 1.03, 
           linetype = "solid", color = "black", size = 0.75) +
  annotate("text", x = 1.82, y = 1.04, label = "***", size = 5, color = "black") +
  annotate("segment", x = 1.99, xend = 2.29, y = 1.03, yend = 1.03, 
           linetype = "solid", color = "black", size = 0.75) +
  annotate("text", x = 2.15, y = 1.04, label = "****", size = 5, color = "black") +
  annotate("segment", x = 1.66, xend = 2.29, y = 1.08, yend = 1.08, 
           linetype = "solid", color = "black", size = 0.75) +
  annotate("text", x = 1.98, y = 1.09, label = "****", size = 5, color = "black")

# STATS FOR FIGURE 2B
# diff subject (inter-individual)
dist.merged.all.stats.diffsubject <- dist.merged.all.filt %>% 
  mutate(comparison = paste(sample1, "-", sample2, sep = "")) %>% filter(!same.subject)
dist.merged.all.stats.diffsubject$metric <- factor(dist.merged.all.stats.diffsubject$metric)
dist.merged.all.test.samesubject <- friedman.test(y = dist.merged.all.stats.diffsubject$value, groups = dist.merged.all.stats.diffsubject$metric, blocks = dist.merged.all.stats.diffsubject$comparison)
pairwise.wilcox.test(dist.merged.all.stats.diffsubject$value, dist.merged.all.stats.diffsubject$metric, p.adj = "bonf")

# same subject (intra-individual)
dist.merged.all.stats.samesubject <- dist.merged.all.filt %>% 
  mutate(comparison = paste(sample1, "-", sample2, sep = "")) %>% filter(same.subject)
dist.merged.all.stats.samesubject$metric <- factor(dist.merged.all.stats.samesubject$metric)
dist.merged.all.test.samesubject <- friedman.test(y = dist.merged.all.stats.samesubject$value, groups = dist.merged.all.stats.samesubject$metric, blocks = dist.merged.all.stats.samesubject$comparison)
pairwise.wilcox.test(dist.merged.all.stats.samesubject$value, dist.merged.all.stats.samesubject$metric, p.adj = "bonf")






# stability over time

# overview of number of samples remaining per subject
metadata.joined %>% 
  select(Subject, Time_point) %>% 
  group_by(Subject) %>% 
  summarize(nrows = n())

dist.merged.all.filt.timepoints <- dist.merged.all.filt %>% 
  filter(same.subject) %>% 
  filter(!metric %in% c("dist.aitch.filt", "dist.aitch.filt.phf")) %>% 
  select(sample1, sample2, metric, value) %>% 
  left_join(metadata.joined %>% select(Date, Sample, Subject), by = c("sample1" = "Sample")) %>% 
  left_join(metadata.joined %>% select(Date, Sample, Subject), by = c("sample2" = "Sample"), suffix = c(".sample1", ".sample2")) %>% 
  mutate()


dist.merged.all.filt.timepoints$Date.sample1 <- as.character(dist.merged.all.filt.timepoints$Date.sample1)
dist.merged.all.filt.timepoints$Date.sample2 <- as.character(dist.merged.all.filt.timepoints$Date.sample2)

dist.merged.all.filt.timepoints$Date.sample1 <- as.Date(dist.merged.all.filt.timepoints$Date.sample1, "%d/%m/%y")
dist.merged.all.filt.timepoints$Date.sample2 <- as.Date(dist.merged.all.filt.timepoints$Date.sample2, "%d/%m/%y")

dist.merged.all.filt.timepoints <- aggregate(value ~ Subject.sample1 + Date.sample1 + Date.sample2 + metric, 
                                             data = dist.merged.all.filt.timepoints, 
                                             FUN = mean)


dist.merged.all.filt.timepoints.final <- dist.merged.all.filt.timepoints %>% 
  mutate(timepoint_after = Date.sample2 > Date.sample1) %>% 
  filter(timepoint_after) %>% 
  mutate(timepoint_diff = Date.sample2 - Date.sample1) %>% 
  group_by(Subject.sample1, metric, Date.sample1) %>% 
  arrange(timepoint_diff) %>% 
  slice(1)

# FIGUERE 2C
fig2c <- dist.merged.all.filt.timepoints.final %>% 
  filter(!metric %in% c("dist.wuni.filt.phf")) %>% 
  ggplot(aes(x = Subject.sample1, y = 1-value, group = interaction(Subject.sample1, metric), color = metric)) +
    geom_boxplot(outlier.shape = NA, width = 0.75, aes(fill = metric), alpha = 0.4, color = "black", position = position_dodge(width = 1)) +
    geom_point(position = position_jitterdodge()) +
    scale_fill_viridis(discrete = T, begin = 0.1, end = 0.3, labels = labels) +  
    scale_colour_viridis(discrete = T, begin = 0.1, end = 0.3, labels = labels) +
    theme_bw() +
    theme(legend.position = "none") +
    labs(x = "Individual", y = "Stability (1 - distance)", color = "Distance metric", fill = "Distance metric") +
    annotate("segment", x = 0.75, xend = 1.25, y = 0.86, yend = 0.86, 
             linetype = "solid", color = "black", size = 0.75) +
    annotate("text", x = 1, y = 0.87, label = "*", size = 6, color = "black") +
    annotate("segment", x = 2.75, xend = 3.25, y = 0.86, yend = 0.86, 
             linetype = "solid", color = "black", size = 0.75) +
    annotate("text", x = 3, y = 0.87, label = "**", size = 6, color = "black")+
    annotate("segment", x = 8.75, xend = 9.25, y = 0.86, yend = 0.86, 
             linetype = "solid", color = "black", size = 0.75) +
    annotate("text", x = 9, y = 0.87, label = "***", size = 6, color = "black")

# STATS FOR FIGURE 2C
dist.merged.all.filt.timepoints.final %>% 
  filter(!metric %in% c("dist.wuni.filt.phf")) %>% 
  ungroup() %>% 
  group_by(Subject.sample1) %>% 
  do(w = wilcox.test((1-value) ~ metric, data = ., paired = T)) %>% summarize(Subject.sample1, Wilcox = w$p.value)

### COMBINE ALL FIGURES
fig2b.2row <- fig2b + guides(fill = guide_legend(nrow = 2, ncol = 2))
fig2c.2row <- fig2c + guides(fill = guide_legend(nrow = 2, ncol = 2))
fig2 <- ggarrange(fig2a, 
          ggarrange(fig2b.2row, fig2c.2row, ncol = 1, common.legend = T, legend = "bottom", labels = c("B", "C")),
          ncol = 2, labels = "A", widths = c(2, 0.8))

ggsave(file="Figure2.pdf", plot=fig2, width=16, height=8)




# SUPPLEMENTARY FIGURE 1 - CHECKV FILTERING
checkv.out <- read_tsv("quality_summary.tsv")
checkv.filt <- checkv.out %>% filter(completeness > 50)
dim(checkv.filt)


ps.melted.host.filt <- ps.melted.host %>% filter(OTU %in% checkv.filt$contig_id)

#scale abundances considering we removed contigs <=50% complete
ps.melted.host.filt <- ps.melted.host.filt %>% group_by(Sample) %>% mutate(Abundance = Abundance/sum(Abundance))

ps.melted.host.filt$`Host genome`[is.na(ps.melted.host.filt$`Host genome`)] <- "Unknown"


color_objs.hic <- create_color_dfs(as.data.frame(ps.melted.host.filt),
                                   selected_groups = c("Desulfobacterota", "Proteobacteria", "Actinobacteriota",
                                                       "Bacteroidota", "Firmicutes"),
                                   group_level = "phylum",
                                   subgroup_level = "family",
                                   top_orientation = TRUE)
mdf.hic <- color_objs.hic$mdf
cdf.hic <- color_objs.hic$cdf
legend.hic <- custom_legend(mdf.hic, cdf.hic, group_level = "phylum", subgroup_level = "family", legend_key_size = 0.6, legend_text_size = 14)


plot_diff.hic <- plot_microshades(mdf.hic, cdf.hic, x = "Sample") + 
  scale_y_continuous(labels = scales::percent, expand = expansion(0)) +
  facet_grid(. ~ Subject, scales = "free", space = "free") +
  theme(axis.text.x = element_text(size= 10)) +
  theme(legend.position = "none") +
  theme(axis.text.x=element_blank())

suppl.fig1 <-  plot_grid(plot_diff.hic, arrangeGrob(legend, nullGrob(), ncol = 1, heights = c(0.80, 0.20)),  rel_widths = c(1, .25))

ggsave(file="SupplFigure1.pdf", plot=suppl.fig1, width=12.8, height=8)
