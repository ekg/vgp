#!/usr/bin/env Rscript

library(dplyr)

# fai file is the first argument
if (length(commandArgs(trailingOnly = TRUE)) != 1) {
  stop("Usage: ", commandArgs()[1], " <fai file>")
}
fai_file <- commandArgs(trailingOnly = TRUE)[1]

x <- read.delim(fai_file, header = FALSE, sep = "\t", col.names = c("seq", "length", "offset", "line.chars", "line.chars1"))

# create a new variable, "genome" which contains the genome prefix
# where a hash '#' is used to separate the genome prefix from the rest of the sequence name
# group by genome and summarize the total length
x$genome <- sapply(strsplit(x$seq, "#"), `[`, 1)

result <- x %>%
  group_by(genome) %>%                               # Group by genome
  summarize(genome_length = sum(length))              # Sum lengths per group

print(result)

require(ggplot2)

# plot a density plot
ggplot(result, aes(x = genome_length)) +
  geom_density(fill = "blue", alpha = 0.5) +
  ggtitle("Density plot of total genome length") +
  xlab("Total genome length") +
  ylab("Density")
# this should concatenate the filenames
# and save to the file
ggsave(paste(fai_file, ".genome_length_density.pdf", sep = ""), width=8, height=4)

# plot a histogram
ggplot(result, aes(x = genome_length)) +
  geom_histogram(binwidth = 1e8, fill = "blue", alpha = 0.5) +
  ggtitle("Histogram of total genome length") +
  xlab("Total genome length") +
  ylab("Frequency")
# save the file
ggsave(paste(fai_file, ".genome_length_histogram.pdf", sep = ""), width=8, height=4)


# Sort the data by genome_length in decreasing order
result <- result[order(result$genome_length, decreasing = TRUE), ]

# Calculate the maximum width of the labels
#max_label_width <- max(nchar(result$genome)) * 0.02
#print(max_label_width)

#  geom_text(aes( label = genome), size = 2, hjust = -0.2, vjust = 0.5) +

# Create the plot
plot <- ggplot(result, aes(x = genome_length, y = reorder(genome, genome_length))) +
  geom_point(size = 1.5) +
  geom_segment(aes(x = 0, xend = genome_length, yend = genome), linetype = "solid", color = "black", size = 0.2) +
  geom_text(aes(label = paste0(round(genome_length / 1e6, 2), " Mb   ", genome)), size = 2, hjust = -0.1, vjust = 0.5) +
  labs(x = "Genome Length", y = "Assembly", title = "Genome Length Distribution") +
  theme_minimal() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  ) +
  scale_x_continuous(breaks = seq(0, max(result$genome_length), by = 500e6), labels = function(x) paste0(x / 1e9, " Gb"), expand = c(0, 0.05))

# Save the plot
ggsave(paste(fai_file, ".genome_length_species.pdf", sep = ""), plot = plot, width = 20, height = 40)


# write a new version of the plot that uses a cool ggplot compatible R library to display all the species names on it
#require(ggplot2)
#require(ggrepel)
# avoid ggrepel: 475 unlabeled data points (too many overlaps). Consider increasing max.overlaps
# set max.overlaps
# make the text really small
#ggplot(result, aes(x = genome_length, y = 0, label = genome)) +
#  geom_text_repel(max.overlaps = Inf) +
#  ggtitle("Genome length") +
#  xlab("Total genome length") +
#  ylab("Genome") +
#  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
#  theme(axis.text.x = element_text(size = 4))
#ggsave(paste(fai_file, ".genome_length_text.png", sep = ""), width=20, height=40)

# write the genome list out sorted by size from longest to shortest
write.table(result[order(-result$genome_length),], paste(fai_file, "genome_list.txt", sep = ""), quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")


# yay, we're done