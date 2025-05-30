
.libPaths("/home/tbellagio/miniforge3/envs/r-environment/lib/R/library")
library(gradientForest)
library(data.table)

# Get the command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Extract split and batch numbers
split_num <- as.integer(args[1])
#batch_num <- as.integer(args[2])

#cat(sprintf("Starting Gradient Forest training for Split %d, Batch %d\n", split_num, batch_num))

# Read environmental data
env <- fread("env_scaled_bio1.csv")

# Construct the split folder path
split_dir <- paste0("splits_samples/split_", split_num, "/")
train_file <- paste0(split_dir, "train_samples.csv")
test_file <- paste0(split_dir, "test_samples.csv")

# Read train samples
train_samples <- fread(train_file)$sample

# Filter environmental data for train samples
env_train <- env[samples %in% train_samples]

# Construct the batch folder path
#batch_dir <- paste0("batch_", batch_num, "/")
#snp_filename <- paste0(batch_dir, "snp_batch_", batch_num, ".csv")

# Load the SNP batch data
delta_p <- fread('deltap_ldprunned_gen1.csv')

# Select only the columns that match the train samples
delta_p_train <- delta_p[, ..train_samples]

delta_p_train_t = as.data.frame(t(delta_p_train))

# Reorder delta_p_train_t to match the order of env_samples
delta_p_train_t <- delta_p_train_t[match(env_train$samples, rownames(delta_p_train_t)), ]

# Check alignment again after reordering
if (identical(env_train$samples, rownames(delta_p_train_t))) {
  cat("Sample names are now aligned after reordering.\n")
} else {
  cat("Alignment issue persists!\n")
}

colnames(delta_p_train_t) <- paste("snp", 1:ncol(delta_p_train_t), sep = "_")

# Remove the 'samples' column from the environmental data
env_train <- env_train[, -1, with = FALSE]
#delta_p_train_t <- delta_p_train_t[, -1, with = FALSE]

# Bind the aligned data
gf_data <- cbind(env_train, delta_p_train_t)

## collect test samples
test_samples <- fread(test_file)$sample

# Filter environmental data for train samples
env_test <- env[samples %in% test_samples]
# Remove the 'samples' column from the test data
env_test <- env_test[, -1, with = FALSE]

# Remove duplicate rows from env_test
env_test <- unique(env_test)


# Train the Gradient Forest model
gf_model <- gradientForest(
  data = gf_data, 
  predictor.vars = colnames(env_train),  # Exclude 'site'
  response.vars = colnames(delta_p_train_t), 
  ntree = 500, 
  maxLevel = floor(log2(0.368 * nrow(env_train) / 2)), 
  trace = TRUE, 
  corr.threshold = 0.50
)

# Save the model inside the split and batch folder
model_filename <- paste0(split_dir, "bio1_gf_model_split_", split_num, ".rds")
saveRDS(gf_model, model_filename)

# make trian data unique 
env_train_unique <- unique(env_train)

# 1. Transform environments to "genetic space"
natural_trans <- predict(gf_model, env_train_unique)  # Should be 30 x n
garden_trans <- predict(gf_model, env_test)          # Should be 1 x n

# 2. Convert to matrices if they aren't already
natural_trans <- as.matrix(natural_trans)
garden_trans <- as.matrix(garden_trans)

# 3. Verify dimensions
stopifnot(ncol(natural_trans) == ncol(garden_trans))  # Must have same # of GF axes

# 4. Calculate genomic offset (Euclidean distance)
offsets <- sqrt(rowSums((natural_trans - garden_trans[rep(1, nrow(natural_trans)), ])^2))

# Save transformed environments and offsets as CSV files
natural_filename_csv <- paste0(split_dir, "bio1_natural_trans_split_", split_num, ".csv")
garden_filename_csv <- paste0(split_dir, "bio1_garden_trans_split_", split_num, ".csv")
offsets_filename_csv <- paste0(split_dir, "bio1_offsets_split_", split_num, ".csv")

# Convert to data frames (if they are matrices)
natural_trans_df <- as.data.frame(natural_trans)
garden_trans_df <- as.data.frame(garden_trans)
offsets_df <- data.frame(offsets = offsets)

# Save as CSV
write.csv(natural_trans_df, natural_filename_csv, row.names = FALSE)
write.csv(garden_trans_df, garden_filename_csv, row.names = FALSE)
write.csv(offsets_df, offsets_filename_csv, row.names = FALSE)

cat(sprintf("Natural transformed data saved at %s\n", natural_filename_csv))
cat(sprintf("Garden transformed data saved at %s\n", garden_filename_csv))
cat(sprintf("Genomic offsets saved at %s\n", offsets_filename_csv))


