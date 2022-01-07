#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(pastecs))

# parse arguments
args = commandArgs(trailingOnly = T)
segment <- args[3]

# read data
centrifuge_res <- fread(args[1], header = T, sep = "\t")
segment_ids <- fread(args[2], header = F, sep = ",") %>% filter(V2 == segment) %>% pull(V1)
segment_mapping <- fread(args[2], header = F, sep = ",") %>% rename(ID = V1, segment = V2)
mode <- args[7]

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

# select a suitable read as the initial template for error correction
# run read length distribution peak detection if dynamic mode
if (mode == "dynamic") {
## calculate read lengths density
tryCatch( {
  length_density <- centrifuge_res %>%
  filter(!(readID %in% chimeric_reads),
         seqID %in% segment_ids) %>%
  pull(queryLength) %>%
  density() 
  },
  error=function(e) {
    message(paste0("influenza_consensus: Cannot calculate the read length distribution for Segment ",
    gsub("segment_", "", segment),
    ". Falling back to using --mode static"))
    best_peak <- expected_length
  }
)
## find peaks in read lengths density
tryCatch( {
  tps <- turnpoints(length_density$y)
  peaks <- data.frame(x = length_density$x[which(extract(tps) == 1)],
                    y = length_density$y[which(extract(tps) == 1)],
                    prob = tps$prob[which(tps$tppos %in% which(extract(tps) == 1))]) %>%
            filter(prob <= 0.00005)
}, error=function(e) {
  message(paste0("influenza_consensus: Failed to detect any peak signals for Segment ",
    gsub("segment_", "", segment),
    ". Falling back to using --mode static"))
    best_peak <- expected_length
})

if ( exists("peaks") == F ) {
    best_peak <- expected_length
  } else if (nrow(peaks) == 0) {
    best_peak <- expected_length
  } else {
    ## select the peak closest to the expected length
    best_peak <- peaks %>% 
      mutate(difference = abs(x - expected_length)) %>%
      arrange(-desc(difference)) %>%
      slice(1) %>%
      pull(x)
  }

} else {
  best_peak <- expected_length
}
print(paste0(gsub("_", " ", segment), " peak found at:", best_peak, " bps"))

# identify fastq sequence IDs for consensus building
fastq_ids <- centrifuge_res %>%
  filter(!(readID %in% chimeric_reads),
         seqID %in% segment_ids,
         queryLength <= best_peak+50,
         queryLength >= best_peak-50) %>%
  #select(readID, queryLength) %>%
  #distinct() %>%
  #mutate(difference = abs(queryLength - best_peak)) %>%
  #arrange(-desc(difference)) %>%
  #slice(1:subsample_n) %>%
  rename(lengths = queryLength)

# identify fastq sequence IDs for all reads mapped to target segment
binned_fastq_ids <- centrifuge_res %>%
  filter(!(readID %in% chimeric_reads),
         seqID %in% segment_ids) %>%
  rename(lengths = queryLength)

# write output
write.table(unique(binned_fastq_ids$readID), file = args[4], quote = F, row.names = F, col.names = F)
write.table(unique(fastq_ids$readID), file = args[5], quote = F, row.names = F, col.names = F)
write.table(select(binned_fastq_ids, lengths), file = args[6], quote = F, row.names = F, col.names = T)