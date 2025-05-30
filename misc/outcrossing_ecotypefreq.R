#R script to infer the lower bound to rate of outcrossing within GrENE-net experiment
setwd("/global/scratch/users/ruthkepstein/grenet_outcrossing")
setwd("/global/scratch/users/ruthkepstein/grenet_outcrossing")
setwd("~/Shareddrives//MOI-LAB/PROJECTS/grenenet/GRENE-net_PHASE1/grenephase1-analyses/") # GOOGLE DRVIE

library(ggplot2)
library(dplyr)
library(data.table)

#file about each sample, 2145 samples across 3 generations & 32 sites
meta = read.table("samples_data_fix57.csv", sep = ",", header = T)
meta = meta[!is.na(meta$sampleid),]

#example with a single sample
sample.id = "MLFH090320200606"
sample_ex = read.table(paste("/global/home/users/xingwu/moi_lab/PROJECTS/grenenet-phase1/frequency/hapFIRE_frequencies/samples/ecotype_frequency/",sample.id,"_ecotype_frequency.txt",sep=""),row.names = NULL, header = F)
sample_ex = read.table(paste("data/ecotype_frequency/",sample.id,"_ecotype_frequency.txt",sep=""),row.names = NULL, header = F)
colnames(sample_ex) = c("Ecotype", "Frequency")

#consider only ecotypes with a frequency above 10%
sample_ex_filt = sample_ex %>% filter(Frequency >= 0.1)

#finding the row in the meta file with the sample name
p = meta[meta$sampleid %like% sample.id,]

#creating a function to do all samples at once
outcrossing_df = data.frame(matrix(NA, nrow = 2415, ncol = 7))
colnames(outcrossing_df) = c("sample.id", "num_flowers", "num_ecotypes", "outcross", "year", "site", "half")
outcross = function(meta,num_flow,freq){
  for(i in 1:nrow(meta)){
    sample.id = meta$sampleid[i]
    # sample_curr = read.table(paste("/global/home/users/xingwu/moi_lab/PROJECTS/grenenet-phase1/frequency/hapFIRE_frequencies/samples/ecotype_frequency/",sample.id,"_ecotype_frequency.txt",sep=""),row.names = NULL, header = F)
    sample_curr = read.table(paste("data/ecotype_frequency/",sample.id,"_ecotype_frequency.txt",sep=""),row.names = NULL, header = F)
    colnames(sample_curr) = c("Ecotype", "Frequency")
    sample_curr_filt = sample_curr %>% filter(Frequency >= freq)
    
    p = meta[meta$sampleid %like% sample.id,]
    
    #checking if ecotypes are 50/50 suggesting recent outcrossing of flowers = 1
    top_2 = sample_curr_filt %>% top_n(2)
    if(isTRUE(top_2$Frequency[1] > .45 & top_2$Frequency[1] < .55 & top_2$Frequency[2] > .45 & top_2$Frequency[2] < .55)){
      outcrossing_df$half[i] = "recent"
    } else{
      outcrossing_df$half[i] = "no"
    }
    
    outcrossing_df$sample.id[i] = sample.id
    outcrossing_df$num_flowers[i] = p$flowerscollected
    outcrossing_df$num_ecotypes[i] = nrow(sample_curr_filt)
    outcrossing_df$year[i] = p$year
    outcrossing_df$site[i] = p$site
    outcrossing_df$month[i] = p$month
  }
  outcrossing_filt = outcrossing_df %>% filter(num_flowers <= num_flow)
  for(k in 1:nrow(outcrossing_filt)){
    if(isTRUE(outcrossing_filt$num_ecotypes[k] > outcrossing_filt$num_flowers[k])){
      outcrossing_filt$outcross[k] = "yes"
    } else{
      outcrossing_filt$outcross[k] = "no"
      }
  }
  return(outcrossing_filt)
}

outcrossprop = outcross(meta,2,0.1)

#plotting outcrossing rate per year
library(scales)
n_years <- length(unique(all_num$year))
palette <- gradient_n_pal(c("lightgreen", "darkgreen"))(seq(0, 1, length.out = n_years))
rate = ggplot(all_freq, aes(x = num_flowers, y = outcrossing_freq, color = year, group = year)) +
  geom_line(size = 1) + 
  geom_point(size = 3) +    
  scale_color_manual(values = palette, name = "Year")+
  scale_x_continuous(breaks = seq(1, 10, 1)) +
  scale_y_continuous(breaks = seq(1, 18, 1)) +
  labs(
    x = "# of flowers per sample",
    y = "outcrossing rate",
    color = "Year"            
  ) + xlim(0,8) + ylim(0,18) +
  theme_minimal()

#plotting the number of outcrossing events as well
count = ggplot(all_num, aes(x = num_flowers, y = outcrossing_freq, color = year, group = year)) +
  geom_line(size = 1) +        
  geom_point(size = 3) +       
  scale_color_manual(values = palette, name = "Year")+
  scale_x_continuous(breaks = seq(1, 10, 1)) +
  scale_y_continuous(breaks = seq(1, 25, 1)) +
  labs(
    x = "# of flowers per sample",
    y = "outcrossing count",
    color = "Year"             # Legend title
  ) + xlim(0,8) + ylim(0,25)+
  theme_minimal()

#combining the 2 plots
library(cowplot)
combined_plot <- plot_grid(count, rate, labels = c("A", "B"), ncol = 1, nrow = 2)

print(combined_plot)

all_freq = rbind(first,second,third)
all_freq$year = as.factor(all_freq$year)

all_num = rbind(first_num,second_num,third_num)
all_num$year = as.factor(all_num$year)

write.table(first_num, file = "2018_outcrossingnum_numoflowers.txt")
write.table(second_num, file = "2019_outcrossingnum_numoflowers.txt")
write.table(third_num, file = "2020_outcrossingnum_numoflowers.txt")

#plotting 50/50 F1 outcrossed individual
top5_freq_df <- sample_ex %>%
  slice_max(Frequency, n = 5)
top5_freq_df$Frequency = top5_freq_df$Frequency*100

ggplot(top5_freq_df, aes(x = reorder(Ecotype, -Frequency), y = Frequency)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  geom_text(
    aes(label = ifelse(Frequency > 10, paste0(round(Frequency, 1), "%"), "")),
    vjust = -0.3,
    size = 3
  ) + ylim(0, 50) + 
  labs(
    title = "Frequency of Each Ecotype for sample MLFH090320200606",
    x = "Ecotypes detected",
    y = "Frequency of each ecotype (%)"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6),
    plot.title = element_text(hjust = 0.5)
  )

##factorial design for 2 flower samples to predict most likely crosses found within each sample
library(combinat)

is_ecotype_present <- function(frequency) {
  return(frequency > 0.1)
}

count_present_ecotypes <- function(sample) {
  return(sum(sapply(sample, is_ecotype_present)))
}

generate_ecotype_combinations <- function() {
  #all combinations of crosses with 2 flowers
  combinations <- list(
    two_inbreds = list(list(e1 = 1, e2 = 1), list(e1 = 2)),
    inbred_F1 = list(list(e1 = 1, e2 = 0.5, e3 = 0.5), list(e1 = 1, e2 = 1)),
    two_F1s = list(
      list(e1 = 1),
      list(e1 = 0.75, e2 = 0.25),
      list(e1 = 0.5, e2 = 0.5),
      list(e1 = 0.5, e2 = 0.25, e3 = 0.25),
      list(e1 = 0.25, e2 = 0.25, e3 = 0.25, e4 = 0.25)
    ),
    inbred_F2 = list(
      list(e1 = 0.75, e2 = 0.25),
      list(e1 = 0.5, e2 = 0.375, e3 = 0.125),
      list(e1 = 0.5, e2 = 0.25, e3 = 0.25)
    ),
    inbred_F1_backcross = list(
      list(e1 = 0.875, e2 = 0.125),
      list(e1 = 0.625, e2 = 0.375),
      list(e1 = 0.5, e2 = 0.375, e3 = 0.125)
    ),
    inbred_F3 = list(
      list(e1 = 0.75, e2 = 0.25),
      list(e1 = 0.5, e2 = 0.375, e3 = 0.125),
      list(e1 = 0.5, e2 = 0.25, e3 = 0.25)
    ),
    two_F2s = list(
      list(e1 = 1),
      list(e1 = 0.75, e2 = 0.25),
      list(e1 = 0.5, e2 = 0.5),
      list(e1 = 0.5, e2 = 0.25, e3 = 0.25),
      list(e1 = 0.25, e2 = 0.25, e3 = 0.25, e4 = 0.25)
    ),
    F1_F2 = list(
      list(e1 = 0.625, e2 = 0.375),
      list(e1 = 0.5, e2 = 0.3125, e3 = 0.1875),
      list(e1 = 0.375, e2 = 0.3125, e3 = 0.3125),
      list(e1 = 0.375, e2 = 0.3125, e3 = 0.1875, e4 = 0.125)
    ),
    F2_F3 = list(
      list(e1 = 0.625, e2 = 0.375),
      list(e1 = 0.5, e2 = 0.3125, e3 = 0.1875),
      list(e1 = 0.375, e2 = 0.3125, e3 = 0.3125),
      list(e1 = 0.375, e2 = 0.3125, e3 = 0.1875, e4 = 0.125)
    ),
    F1_F3 = list(
      list(e1 = 0.625, e2 = 0.375),
      list(e1 = 0.5, e2 = 0.3125, e3 = 0.1875),
      list(e1 = 0.375, e2 = 0.3125, e3 = 0.3125),
      list(e1 = 0.375, e2 = 0.3125, e3 = 0.1875, e4 = 0.125)
    )
  )
  return(combinations)
}

calculate_category_probabilities <- function(sample, num_flowers, ecotype_combinations, outcrossing_rate) {
  present_ecotypes <- sum(sample > 0.1)
  cat("Present ecotypes:", present_ecotypes, "\n")
  category_probabilities <- list()
  
  outcrossing_events <- list(
    two_inbreds = 0, inbred_F1 = 1, two_F1s = 2, inbred_F2 = 2, 
    inbred_F1_backcross = 2, inbred_F3 = 3, two_F2s = 4, 
    F1_F2 = 3, F2_F3 = 5, F1_F3 = 4
  )
  
  for (category in names(ecotype_combinations)) {
    category_prob <- 0
    cat("\nCategory:", category, "\n")
    for (combo in ecotype_combinations[[category]]) {
      combo_ecotypes <- length(combo)
      if (combo_ecotypes == present_ecotypes) {
        similarity <- exp(-sum(abs(sort(sample[sample > 0.1], decreasing = TRUE) - sort(unlist(combo), decreasing = TRUE))))
        category_prob <- category_prob + similarity
        cat("Similarity:", similarity, "\n")
      }
    }
    
    events <- outcrossing_events[[category]]
    category_prob <- category_prob * (outcrossing_rate ^ events) * ((1 - outcrossing_rate) ^ (num_flowers - events))
    cat("Category prob after weighting:", category_prob, "\n")
    
    category_probabilities[[category]] <- category_prob
  }
  
  total_prob <- sum(unlist(category_probabilities))
  cat("Total probability before normalization:", total_prob, "\n")
  if (total_prob > 0) {
    category_probabilities <- lapply(category_probabilities, function(x) x / total_prob)
  }
  
  return(category_probabilities)
}

#process a single sample
process_sample <- function(sample_frequencies, num_flowers, outcrossing_rate) {
  ecotype_combinations <- generate_ecotype_combinations()
  probabilities <- calculate_category_probabilities(sample_frequencies, num_flowers, ecotype_combinations, outcrossing_rate)
  return(probabilities)
}

library(dplyr)
library(data.table)

#process all samples at once
process_all_samples <- function(meta, num_flowers, outcrossing_rate, freq_threshold) {
  results <- data.frame(sampleid = character(), highest_prob_category = character(), probability = numeric(), stringsAsFactors = FALSE)
  
  for(i in 1:nrow(meta)) {
    sample.id <- meta$sampleid[i]
    cat("Processing sample:", sample.id, "\n")
    
    sample_file <- paste0("/global/home/users/xingwu/moi_lab/PROJECTS/grenenet-phase1/frequency/hapFIRE_frequencies/samples/ecotype_frequency/", sample.id, "_ecotype_frequency.txt")
    sample_curr <- read.table(sample_file, row.names = NULL, header = FALSE)
    colnames(sample_curr) <- c("Ecotype", "Frequency")
    sample_curr_filt <- sample_curr %>% filter(Frequency >= freq_threshold)
    
    sample_frequencies <- sample_curr_filt$Frequency
    names(sample_frequencies) <- sample_curr_filt$Ecotype
    
    result <- process_sample(sample_frequencies, num_flowers, outcrossing_rate)
    
    #highest probability cross found in sample
    highest_prob <- which.max(unlist(result))
    highest_prob_category <- names(result)[highest_prob]
    highest_prob_value <- unlist(result)[highest_prob]
    
    # Add to results dataframe
    results <- rbind(results, data.frame(sampleid = sample.id, 
                                         highest_prob_category = highest_prob_category, 
                                         probability = highest_prob_value))
  }
  
  return(results)
}

#meta is metadata on all samples
meta = read.table("samples_data_fix57.csv", sep = ",", header = T)
meta = meta[!is.na(meta$sampleid),]
#need samples with only 2 flowers
meta_2 = meta %>% filter(flowerscollected == 2)
num_flowers <- 2
outcrossing_rate <- 0.1
freq_threshold <- 0.1

results <- process_all_samples(meta_2, num_flowers, outcrossing_rate, freq_threshold)

#summary of results
print(table(results$highest_prob_category))

plot_results <- function(results) {
  category_counts <- results %>%
    group_by(highest_prob_category) %>%
    summarise(count = n()) %>%
    mutate(percentage = count / sum(count) * 100)
  
  p1 <- ggplot(category_counts, aes(x = reorder(highest_prob_category, -count), y = count)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    geom_text(aes(label = sprintf("%.1f%%", percentage)), vjust = -0.5) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Distribution of Highest Probability Categories",
         x = "Category", y = "Count")
  
  print(p1)
}

plot_results(results)