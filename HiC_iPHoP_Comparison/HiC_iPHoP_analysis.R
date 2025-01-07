library(tidyverse)
library(data.table)

bin.tax.summary <- read_tsv("gtdbtk.bac120.summary.tsv") %>% select(user_genome, classification) %>% 
  separate(col = classification, into = c("domain", "phylum", "class", "order", "family", "genus", "species"), sep = ";")

bin.tax.summary$domain <- gsub("d__", "", bin.tax.summary$domain)
bin.tax.summary$phylum <- gsub("p__", "", bin.tax.summary$phylum)
bin.tax.summary$class <- gsub("c__", "", bin.tax.summary$class)
bin.tax.summary$order <- gsub("o__", "", bin.tax.summary$order)
bin.tax.summary$family <- gsub("f__", "", bin.tax.summary$family)
bin.tax.summary$genus <- gsub("g__", "", bin.tax.summary$genus)
bin.tax.summary$species <- gsub("s__", "", bin.tax.summary$species)

bin.tax.summary$phylum <- sapply(strsplit(bin.tax.summary$phylum, "_"), `[`, 1)

viralhost.summary <- read_tsv("summary-viral_host_associations_filtered.tsv")

viralhost.summary.merged <- left_join(viralhost.summary, bin.tax.summary, by = c("cluster_name" = "user_genome"))

viralhost.summary.merged$sample <- sapply(str_split(viralhost.summary.merged$mobile_contig_name, "-k141"), `[`, 1)


length(viralhost.summary.merged$mobile_contig_name) # 1,577 phage-host pairings
length(unique(viralhost.summary.merged$mobile_contig_name)) # 1,547 unique phage contigs with hosts
length(unique(viralhost.summary.merged$genus)) # 77 unique genus-level hosts

### import iphop data
iphop_out.hic <- read_csv("iphop-genus.csv")

length(iphop_out.hic$Virus) # 1,990 phage-host pairings
length(unique(iphop_out.hic$Virus)) # 1,543 unique phage contigs with hosts

iphop_out.hic.dedup <- iphop_out.hic %>% group_by(Virus)

iphop_out.hic.dedup <- iphop_out.hic.dedup[, c("Virus", "Host genus", "Confidence score")]

iphop_out.hic.split <- iphop_out.hic.dedup %>% 
  separate(col = `Host genus`, into = c("domain", "phylum", "class", "order", "family", "genus"), sep = ";")

iphop_out.hic.split$domain <- gsub("d__", "", iphop_out.hic.split$domain)
iphop_out.hic.split$phylum <- gsub("p__", "", iphop_out.hic.split$phylum)
iphop_out.hic.split$class <- gsub("c__", "", iphop_out.hic.split$class)
iphop_out.hic.split$order <- gsub("o__", "", iphop_out.hic.split$order)
iphop_out.hic.split$family <- gsub("f__", "", iphop_out.hic.split$family)
iphop_out.hic.split$genus <- gsub("g__", "", iphop_out.hic.split$genus)

iphop_out.hic.split$phylum <- sapply(strsplit(iphop_out.hic.split$phylum, "_"), `[`, 1)
iphop_out.hic.split$class <- sapply(strsplit(iphop_out.hic.split$class, "_"), `[`, 1)
iphop_out.hic.split$order <- sapply(strsplit(iphop_out.hic.split$order, "_"), `[`, 1)
iphop_out.hic.split$genus <- sapply(strsplit(iphop_out.hic.split$genus, "_"), `[`, 1)


iphop_out.hic.split %>% filter(Virus %in% viralhost.summary.merged$mobile_contig_name) %>% 
  dim() # 1,587 phage-host pairings for contigs shared with Hi-C

iphop_out.hic.split %>% filter(Virus %in% viralhost.summary.merged$mobile_contig_name) %>% 
  pull(Virus) %>% unique() %>% length() # 1,243 unique phage contigs making up these pairings

iphop_out.hic.split %>% filter(Virus %in% viralhost.summary.merged$mobile_contig_name) %>% 
  pull(genus) %>% unique() %>% length() # 108 unique genus-level hosts

### Reviewer 2 comment

iphop_out.hic.split %>% 
  filter(Virus %in% viralhost.summary.merged$mobile_contig_name) %>% 
  group_by(Virus) %>% 
  arrange(desc(`Confidence score`)) %>% 
  filter(n() > 1)
#252 phage contigs have >1 genus-level predicted host
#(596 phage-host pairings total, so average of 2.37 hosts per phage)

iphop_out.hic.split %>% 
  filter(Virus %in% viralhost.summary.merged$mobile_contig_name) %>% 
  group_by(Virus) %>% 
  arrange(desc(`Confidence score`)) %>% 
  filter(n() > 1) %>% 
  select(Virus, family, `Confidence score`) %>% 
  summarise(n_family = as.factor(length(unique(family)))) %>% 
  summary()
#of the 252 contigs, 231 (91.7%) of them had the same family-level prediction,
#20 (7.9%) of them had different 2 family-level predictions,
#and 1 (0.4%) of them had 3 different family-level predictions


multi.family.phages <- iphop_out.hic.split %>% 
  filter(Virus %in% viralhost.summary.merged$mobile_contig_name) %>% 
  group_by(Virus) %>% 
  arrange(desc(`Confidence score`)) %>% 
  filter(n() > 1) %>% 
  select(Virus, family, `Confidence score`) %>% 
  summarise(n_family = length(unique(family))) %>% 
  filter(n_family > 1) %>% pull(Virus)

iphop_out.hic.split %>% 
  filter(Virus %in% viralhost.summary.merged$mobile_contig_name) %>% 
  group_by(Virus) %>% 
  arrange(desc(`Confidence score`)) %>% 
  filter(n() > 1) %>% 
  select(Virus, order, `Confidence score`) %>% 
  filter(Virus %in% multi.family.phages) %>% 
  summarise(n_order = as.factor(length(unique(order)))) %>% 
  summary()
#of the 21 phages with more than 1 family hit, 15 (71%) of them were within the same predicted order,
#the remaining 6 (29%) of them had a different order predicted

###


# Comparison between HiC and iPHoP - TOP HIT (HI-C) VS TOP HIT (IPHOP)
comparison <- inner_join(iphop_out.hic.split %>% 
                           filter(!is.na(genus)) %>% 
                           group_by(Virus) %>% 
                           arrange(desc(`Confidence score`)) %>% 
                           slice(1),
                         viralhost.summary.merged %>% 
                           filter(domain != "Unclassified") %>% 
                           filter(domain != "Unclassified Bacteria") %>% 
                           filter(!is.na(genus)) %>% 
                           group_by(mobile_contig_name) %>% 
                           arrange(desc(`adjusted_inter_connective_linkage_density (reads/kbp^2)`)) %>% 
                           slice(1) %>% 
                           select("mobile_contig_name", "domain", "phylum", "class", "order", "family", "genus"),
                         by = c("Virus" = "mobile_contig_name"),
                         suffix = c(".iphop", ".hic")
)

comparison <- comparison %>% 
  mutate(
    same.phylum = (phylum.iphop == phylum.hic),
    same.class = (class.iphop == class.hic),
    same.order = (order.iphop == order.hic),
    same.family = (family.iphop == family.hic),
    same.genus = (genus.iphop == genus.hic)
  ) %>% ungroup()

comparison.tophit.results.summary <- as.data.frame(t(apply(comparison %>% select(same.phylum, same.class, same.order, same.family, same.genus), 2, function(x) sum(x == TRUE) / length(x))))
comparison.tophit.results.summary$comparison <- c("top.hit")
colnames(comparison.tophit.results.summary) <- c("Phylum", "Class", "Order", "Family", "Genus", "comparison")


# Comparison between HiC and iPHoP - ALL HITS (HI-C) VS ALL HITS (IPHOP)
hic.summary <- setDT(viralhost.summary.merged %>% filter(!is.na(genus)) 
                     %>% filter(domain != "Unclassified") %>% 
                       filter(domain != "Unclassified Bacteria"))[, list(hic.phyla = list(phylum), 
                        hic.class = list(class),
                        hic.order = list(order),
                        hic.family = list(family),
                        hic.genus = list(genus)),
                        by = c('mobile_contig_name')] %>%
  select(mobile_contig_name, hic.phyla, hic.class, hic.order, hic.family, hic.genus)


iphop.summary <- setDT(iphop_out.hic.split %>% filter(!is.na(genus)))[, list(iphop.phyla = list(phylum), 
                                                                             iphop.class = list(class),
                                                                             iphop.order = list(order),
                                                                             iphop.family = list(family),
                                                                             iphop.genus = list(genus)),
                                                                      by = c('Virus')] %>%
  select(Virus, iphop.phyla, iphop.class, iphop.order, iphop.family, iphop.genus)


comparison.allhits <- inner_join(iphop.summary, hic.summary, by = c("Virus" = "mobile_contig_name"))

comparison.allhits.results <- comparison.allhits %>% group_by(Virus) %>% 
  mutate(
    same.phylum = ifelse(length(intersect(unlist(hic.phyla), unlist(iphop.phyla))) > 0, 1, 0),
    same.class = ifelse(length(intersect(unlist(hic.class), unlist(iphop.class))) > 0, 1, 0),
    same.order = ifelse(length(intersect(unlist(hic.order), unlist(iphop.order))) > 0, 1, 0),
    same.family = ifelse(length(intersect(unlist(hic.family), unlist(iphop.family))) > 0, 1, 0),
    same.genus = ifelse(length(intersect(unlist(hic.genus), unlist(iphop.genus))) > 0, 1, 0)
  ) %>% ungroup()



comparison.allhits.results.summary <- as.data.frame(t(apply(comparison.allhits.results %>% select(same.phylum, same.class, same.order, same.family, same.genus), 2, function(x) sum(x == 1) / length(x))))
comparison.allhits.results.summary$comparison <- c("all.hits")
colnames(comparison.allhits.results.summary) <- c("Phylum", "Class", "Order", "Family", "Genus", "comparison")


# Comparison between HiC and iPHoP - ALL HITS (HI-C) VS TOP HIT (IPHOP)
iphop.tophit <- iphop_out.hic.split %>% 
  filter(!is.na(genus)) %>% 
  group_by(Virus) %>% 
  arrange(desc(`Confidence score`)) %>% 
  slice(1)

iphop.summary.tophit <- setDT(iphop.tophit %>% filter(!is.na(genus)))[, list(iphop.phyla = list(phylum), 
                                                                             iphop.class = list(class),
                                                                             iphop.order = list(order),
                                                                             iphop.family = list(family),
                                                                             iphop.genus = list(genus)),
                                                                      by = c('Virus')] %>%
  select(Virus, iphop.phyla, iphop.class, iphop.order, iphop.family, iphop.genus)

comparison.topiphop.allhic <- inner_join(iphop.summary.tophit, hic.summary, by = c("Virus" = "mobile_contig_name"))


comparison.topiphop.allhic.results <- comparison.topiphop.allhic %>% group_by(Virus) %>% 
  mutate(
    same.phylum = ifelse(length(intersect(unlist(hic.phyla), unlist(iphop.phyla))) > 0, 1, 0),
    same.class = ifelse(length(intersect(unlist(hic.class), unlist(iphop.class))) > 0, 1, 0),
    same.order = ifelse(length(intersect(unlist(hic.order), unlist(iphop.order))) > 0, 1, 0),
    same.family = ifelse(length(intersect(unlist(hic.family), unlist(iphop.family))) > 0, 1, 0),
    same.genus = ifelse(length(intersect(unlist(hic.genus), unlist(iphop.genus))) > 0, 1, 0)
  ) %>% ungroup()

comparison.topiphop.allhic.results.summary <- as.data.frame(t(apply(comparison.topiphop.allhic.results %>% select(same.phylum, same.class, same.order, same.family, same.genus), 2, function(x) sum(x == 1) / length(x))))
comparison.topiphop.allhic.results.summary$comparison <- c("top.iphop")
colnames(comparison.topiphop.allhic.results.summary) <- c("Phylum", "Class", "Order", "Family", "Genus", "comparison")


comparison.merged <- rbind(comparison.tophit.results.summary, comparison.allhits.results.summary, comparison.topiphop.allhic.results.summary) %>% melt(value.name = "concordance", id = "comparison")

comparison.merged$comparison <- factor(comparison.merged$comparison, levels = c("top.hit", "top.iphop", "all.hits"))

comparison.plot <- comparison.merged %>% 
  ggplot(aes(x = variable, y = 100*concordance, fill = comparison)) +
  geom_bar(stat = "identity", position = "dodge", colour = "black") +
  ylab("% concordance between iPHoP and Hi-C") +
  xlab("Taxonomic rank") +
  scale_fill_manual(values = c("#CC6666", "#9999CC", "#66CC99"), name = "Comparison", labels = c("Top hit (iPHoP) vs top hit (Hi-C)", "Top hit (iPHoP) vs all hits (Hi-C)", "All hits (iPHoP) vs all hits (Hi-C)")) +
  theme_bw() +
  ggbreak::scale_y_break(c(75,90)) +
  ggbreak::scale_y_break(c(0,60)) +
  scale_y_continuous(limits = c(0, 100), breaks = c(60, 65, 70, 75, 90, 92, 94, 96, 98), position = "left") +
  theme(legend.position = "bottom", legend.title = element_blank(), axis.title.x = element_blank()) +
  guides(fill=guide_legend(nrow=2, byrow=TRUE))

comparison.plot

ggsave(file="Figure1.svg", plot=comparison.plot, width=7, height=6)

