

## Example R code

Example R code to look at PHF membership from iPHoP output, generating a plot of # of vOTUs per PHF:

```
library(tidyverse)

# import iPHoP output file (genome level)
iphop_out <- read_csv("Host_prediction_to_genome_m90-merged.csv")

# group iPHoP data by Virus (ie. vOTU/contig), arrange in descending order of prediction confidence, and keep only most confident hit
iphop_out.dedup <- iphop_out %>% 
	select(c("Virus", "Host taxonomy", "Host genome", "Confidence score")) %>%
	group_by(Virus) %>% 
	arrange(desc(`Confidence score`)) %>% 
	slice(1)

# split iPHoP host prediction into proper taxonomic ranks, cleanup rank names
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

# calculate number of vOTUs/contigs per PHF
iphop_out.split.summary <- iphop_out.split %>%
	group_by(family) %>%
	summarize(membership_count = n())

# make plot of # of vOTUs per PHF
ggplot(iphop_out.split.summary, aes(x = reorder(family, -membership_count) , y = membership_count)) +
    geom_point(size = 2) +  
    labs(x = "PHF", y = "# of vOTUs", title = "vOTUs per PHF")  +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust=1))
```

