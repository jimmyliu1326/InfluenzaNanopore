#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))

# parse arguments
args = commandArgs(trailingOnly = T)
segment <- args[3]

# read data
centrifuge_res <- fread(args[1], header = T, sep = "\t")
segment_ids <- fread(args[2], header = F, sep = ",") %>% filter(V2 == segment) %>% pull(V1)
segment_mapping <- fread(args[2], header = F, sep = ",") %>% rename(ID = V1, segment = V2)

# segment expected length
# Segment 1: 2341
# Segment 2: 2341
# Segment 3: 2233
# Segment 4: 1778
# Segment 5: 1565
# Segment 6: 1413
# Segment 7: 1027
# Segment 8: 890
expected_length <- c(2341, 2341, 2233, 1778, 1565, 1413, 1027, 890)[as.double(gsub("segment_", "", segment))]

# filter by taxonomy and alignment score
centrifuge_res <- centrifuge_res %>%
  filter(seqID %in% segment_mapping$ID,
         score >= 100)

# map segment to seqID
centrifuge_res <- centrifuge_res %>%
  left_join(segment_mapping, by = c("seqID" = "ID"))

# identify chimeric reads
chimeric_reads <- centrifuge_res %>%
  filter(numMatches >= 2) %>%
  group_by(readID, segment) %>%
  tally() %>%
  ungroup() %>%
  select(readID) %>%
  group_by(readID) %>%
  tally() %>%
  filter(n >= 2) %>%
  pull(readID)
print(paste0(segment, " chimeric: ", length(chimeric_reads)))

# identify target based on sequence length
target <- centrifuge_res %>%
  filter(!(readID %in% chimeric_reads),
         seqID %in% segment_ids,
         queryLength <= (expected_length+50)) %>%
  group_by(queryLength) %>%
  arrange(desc(score)) %>%
  slice(1) %>%
  ungroup() %>%
  arrange(desc(queryLength)) %>%
  slice(1) %>%
  pull(readID)

# identify fastq sequence IDs for polishing
fastq_ids <- centrifuge_res %>%
  filter(!(readID %in% chimeric_reads),
         seqID %in% segment_ids) %>%
         #queryLength <= expected_length+200) %>%
  filter(readID != target) %>%
  pull(readID) %>%
  unique()

# write output
write.table(target, file = args[4], quote = F, row.names = F, col.names = F)
write.table(fastq_ids, file = args[5], quote = F, row.names = F, col.names = F)