#!/usr/bin/env python
# coding: utf-8

"""
Weighted Z-score Analysis (WZA) Script with Polynomial Correction

This script implements the Weighted Z-score Analysis (WZA) method to combine evidence
across closely linked SNPs within defined genomic windows (e.g., blocks or genes).
It calculates a composite Z-score for each window based on SNP-level association
statistics (e.g., p-values or Z-scores) from GEA or GWAS results. It also performs
a binomial top-candidate test and can apply a spline-based correction to the WZA
scores based on the number of SNPs per window.

Command-line Arguments:
--correlations / -c (required): Path to the input file containing SNP-level association results.
--summary_stat / -s (required): Name of the column in the input file containing the summary statistic (e.g., 'p_wald', 'z_score').
--window / -w (required): Name of the column in the input file containing the window/block identifier (e.g., 'blocks', 'gene').
--output (required): Path for the output file where the WZA results will be saved.
--sample_snps (optional, default=0): Number of SNPs to downsample to within each window if the window has more SNPs than this value. If -1, uses the 75th percentile of SNP counts. If 0, no downsampling.
--resamples (optional, default=100): Number of times to resample SNPs when downsampling is applied (sample_snps > 0).
--min_snps (optional, default=2): Minimum number of SNPs required for a window to be included in the analysis.
--large_i_small_p (optional, flag): Set this flag if large values of the summary statistic indicate stronger association (e.g., Z-scores or test statistics), otherwise small values are assumed to indicate stronger association (e.g., p-values).
--top_candidate_threshold (optional, default=99): Percentile threshold for the top-candidate binomial test (e.g., 99 means top 1%).
--verbose / -v (optional, flag): Enable verbose output.
--MAF (optional): Name of the column containing Minor Allele Frequency (MAF) if it's not named 'MAF'.
--sep (optional, default="\t"): Separator character used in the input file.
--retain (optional, nargs='+'): List of additional column names from the input file to include in the output. For string columns, the first value in the window is taken; for numeric columns, the mean is taken.
--no_SNP_number_correction (optional, flag): Skip the spline-based correction for SNP number.
--empiricalP (optional, flag): Treat the summary statistic as empirical p-values (already ranked/scaled between 0 and 1).

Inputs:
- Input file specified by --correlations: Contains SNP-level data including window identifiers, summary statistics, and MAF.

Outputs:
- Output file specified by --output: A CSV file containing WZA results per window, including the window ID, number of SNPs, raw WZA score (Z), top candidate p-value, and potentially adjusted WZA p-value (Z_pVal) and retained columns.
"""

# --- Import Libraries ---
import pandas as pd
import scipy.stats
import numpy as np
import sys, argparse
from scipy.stats import norm
from pandas.api.types import is_string_dtype, is_numeric_dtype

# --- WZA Calculation Function ---
def WZA(gea_df, statistic, MAF_filter=0.01):
    """
    Calculates the Weighted Z-score (WZA) for a given DataFrame and statistic.

    Args:
        gea_df (pd.DataFrame): DataFrame subset for a specific window/block.
        statistic (str): Name of the column containing the p-values or transformed statistic.
        MAF_filter (float): Minimum Minor Allele Frequency to include a SNP.

    Returns:
        float: The calculated Weighted Z-score, or NaN if the window is filtered out
               or results in a zero denominator.
    """
    # Ensure p-values are within a valid range for z-score conversion
    gea_df[statistic] = gea_df[statistic].clip(lower=1e-15).replace(1, 1-1e-15) # Use 1e-15 for upper bound too to avoid inf

    # Convert p-values to 1-sided Z scores (assuming small p indicates significance)
    # Using sf (survival function) is more robust for very small p-values
    gea_df["z_score"] = scipy.stats.norm.ppf(1 - gea_df[statistic]) # Using 1 - p for significance in upper tail

    # Calculate pbar_qbar (2*p*q) which is proportional to variance of allele frequency
    # Assuming 'MAF' column exists after input processing
    gea_df["pbar_qbar"] = gea_df["MAF"] * (1 - gea_df["MAF"])

    # Apply the MAF filter
    gea_filt = gea_df[gea_df["MAF"] >= MAF_filter].copy()
    
    if gea_filt.shape[0] != gea_df.shape[0]:
        # print(f"Window filtered out due to MAF. Initial SNPs: {gea_df.shape[0]}, After MAF filter: {gea_filt.shape[0]}") # Verbose optional
        pass

    # If no SNPs remain after filtering, return NaN
    if gea_filt.shape[0] == 0:
        # print("Warning: No SNPs remaining after MAF filter for this window.") # Verbose optional
        return np.nan

    ## Calculate the numerator and the denominator for the WZA
    # Numerator: Sum of (pbar_qbar * z_score)
    gea_filt["weiZ_num"] = gea_filt["pbar_qbar"] * gea_filt["z_score"]
    numerator = gea_filt["weiZ_num"].sum()

    # Denominator: Square root of the sum of (pbar_qbar^2)
    gea_filt["weiZ_den"] = gea_filt["pbar_qbar"] ** 2
    denominator = np.sqrt(gea_filt["weiZ_den"].sum())

    # Avoid division by zero
    if denominator == 0:
        # print("Warning: Denominator is zero for window. Returning NaN.") # Verbose optional
        return np.nan

    ## Calculate the Weighted Z-score (ratio of numerator to denominator)
    weiZ = numerator / denominator

    return weiZ

# --- Top Candidate Test Function ---
def top_candidate(gea_df, thresh):
    """
    Performs a binomial test based on the count of SNPs exceeding a significance threshold.

    Args:
        gea_df (pd.DataFrame): DataFrame subset for a specific window/block.
        thresh (float): The p-value threshold for defining a 'hit'.

    Returns:
        tuple: A tuple containing the binomial test p-value and the number of hits.
    """
    # Count the number of SNPs in the window with p-values below the threshold
    hits = (gea_df["pVal"] < thresh).sum()
    snps = gea_df.shape[0] # Total number of SNPs in the window
    
    # Perform a binomial test: Is the number of hits significantly higher than expected by chance?
    # Null hypothesis: probability of a hit is 'thresh'
    # Alternative hypothesis: probability of a hit is greater than 'thresh'
    # binomtest requires integer k (hits) and n (snps), and float p (thresh)
    if snps == 0:
         return 1.0, 0 # Return p=1 and 0 hits if no SNPs

    try:
        top_candidate_p = scipy.stats.binomtest(k=hits, n=snps, p=thresh, alternative="greater").pvalue
    except ValueError as e:
        print(f"Error in binomial test for hits={hits}, snps={snps}, thresh={thresh}: {e}")
        return np.nan, hits # Return NaN p-value on error

    return top_candidate_p, hits

# --- Spline-based WZA Correction Function ---
def adjust_WZA_with_spline(wza_df, roller=50, minEntries=10, chosen_degree=7):
    """
    Adjusts raw WZA Z-scores using a spline based on SNP count to calculate adjusted p-values.

    It calculates a rolling mean and variance of Z-scores with respect to SNP count,
    fits a polynomial (spline) to these, and uses the fitted models to predict the
    expected mean and standard deviation of Z-scores for each window's SNP count.
    Adjusted p-values are then calculated using the normal distribution with these
    predicted parameters.

    Args:
        wza_df (pd.DataFrame): DataFrame containing raw WZA results per window, including 'gene', 'SNPs', and 'Z'.
        roller (int): Window size for rolling mean and variance calculation.
        minEntries (int): Minimum number of windows required in a rolling window.
        chosen_degree (int): Degree of the polynomial for spline fitting.

    Returns:
        pd.DataFrame: Original DataFrame with an added 'Z_pVal' column for adjusted p-values.
    """
    # Remove rows with NaN Z-scores as they cannot be used for fitting or adjustment
    wza_t = wza_df[~wza_df.Z.isnull()].reset_index(drop=True) # reset index after filtering

    # Sort by SNP count for rolling window calculation
    wza_s = wza_t.sort_values('SNPs')

    # Identify windows with fewer SNPs than the minimum required for rolling calculation
    problematic_windows = wza_df[wza_df["SNPs"] < minEntries]
    if not problematic_windows.empty:
        print(f"Warning: {problematic_windows.shape[0]} windows have fewer than {minEntries} SNPs and may affect rolling calculations or receive NaN adjusted p-values.")
        # problematic_windows.to_csv('problematic_windows_for_spline.csv') # Optional: save these

    # Calculate rolling mean and variance of Z-scores based on sorted SNP count
    # min_periods ensures that rolling stats are only calculated when enough data points are available
    rolled_Z_vars = wza_s.Z.rolling(window=roller, min_periods=minEntries).var()
    rolled_Z_means = wza_s.Z.rolling(window=roller, min_periods=minEntries).mean()
    rolled_mean_SNP_number = wza_s.SNPs.rolling(window=roller, min_periods=minEntries).mean()

    # Filter out NaN values that result from the rolling calculation (especially at the beginning)
    masking_array = ~rolled_Z_vars.isnull()
    rolled_Z_sd = np.sqrt(rolled_Z_vars[masking_array])
    rolled_Z_means_masked = rolled_Z_means[masking_array]
    rolled_mean_SNP_number_masked = rolled_mean_SNP_number[masking_array]

    if rolled_Z_sd.empty or rolled_Z_means_masked.empty or rolled_mean_SNP_number_masked.empty:
         print("Error: Rolling window calculation resulted in no valid data points for spline fitting. Skipping adjustment.")
         wza_df["Z_pVal"] = np.nan # Add Z_pVal column with NaNs
         return wza_df


    # Fit polynomial models to predict mean and standard deviation of Z-scores based on SNP number
    try:
        sd_weights = np.polyfit(rolled_mean_SNP_number_masked, rolled_Z_sd, deg=chosen_degree)
        sd_polynomial_model = np.poly1d(sd_weights)

        mean_weights = np.polyfit(rolled_mean_SNP_number_masked, rolled_Z_means_masked, deg=chosen_degree)
        mean_polynomial_model = np.poly1d(mean_weights)
    except np.linalg.LinAlgError:
        print(f"Error: Polynomial fitting failed (LinAlgError) with degree {chosen_degree}. Skipping adjustment.")
        wza_df["Z_pVal"] = np.nan
        return wza_df

    # Predict expected standard deviation and mean for each window's SNP count
    # Use the original wza_df to ensure predictions are made for all windows
    sd_predictions = sd_polynomial_model(wza_df["SNPs"])
    mean_predictions = mean_polynomial_model(wza_df["SNPs"])

    # Ensure predicted standard deviations are non-negative
    sd_predictions = np.clip(sd_predictions, a_min=1e-9, a_max=None) # Clip at a small positive value to avoid division by zero or negative scale

    # Calculate adjusted p-values using the normal distribution with predicted mean and scale
    # Use the survival function (sf) to get P(X > obs) for upper tail p-values
    # Handle potential NaNs in Z scores or predictions
    adjusted_p_values = []
    for i in range(wza_df.shape[0]):
        z_score = wza_df["Z"].iloc[i]
        mean_pred = mean_predictions.iloc[i] # Use .iloc[i] for Series
        sd_pred = sd_predictions.iloc[i]     # Use .iloc[i] for Series

        if pd.isna(z_score) or pd.isna(mean_pred) or pd.isna(sd_pred) or sd_pred <= 0:
            adjusted_p_values.append(np.nan)
        else:
            try:
                # Calculate p-value using the cumulative distribution function (CDF)
                # norm.sf(z, loc=mean, scale=sd) calculates P(X > z)
                p_val = norm.sf(z_score, loc=mean_pred, scale=sd_pred)
                adjusted_p_values.append(p_val)
            except Exception as e:
                print(f"Error calculating adjusted p-value for index {i}: {e}. Z={z_score}, Mean={mean_pred}, SD={sd_pred}")
                adjusted_p_values.append(np.nan)

    # Add the adjusted p-values as a new column
    wza_df["Z_pVal"] = adjusted_p_values

    return wza_df

# --- Main Execution Block ---
def main():
    # Define command line arguments
    parser = argparse.ArgumentParser(description="A script that implements the WZA, a method for combining evidence across closely linked SNPs in GEA studies.")

    parser.add_argument("--correlations", "-c", required=True, dest="correlations", type=str, help="The file containing the correlations")
    parser.add_argument("--summary_stat", "-s", required=True, dest="summary_stat", type=str, help="The name of the column you are analysing")
    parser.add_argument("--window", "-w", required=True, dest="window", type=str, help="The name of column containing the windows you want to analyse")
    parser.add_argument("--output", required=True, dest="output", type=str, help="The name of the output file")
    parser.add_argument("--sample_snps", required=False, dest="sample_snps", type=int, default=0, help="[OPTIONAL] Give the number of SNPs you want to downsample to.")
    parser.add_argument("--resamples", required=False, dest="resamples", type=int, default=100, help="[OPTIONAL] Number of times to resample WZA scores")
    parser.add_argument("--min_snps", required=False, dest="min_snps", type=int, default=2, help="[OPTIONAL] Minimum number of SNPs per window")
    parser.add_argument("--large_i_small_p", required=False, action="store_true", help="[OPTIONAL] Extreme values of the summary stat you're using are large values.")
    parser.add_argument("--top_candidate_threshold", required=False, dest="top_candidate_threshold", type=float, default=99, help="[OPTIONAL] Percentile threshold for top-candidate test")
    parser.add_argument("--verbose", "-v", required=False, action="store_true", help="[OPTIONAL] Verbose mode")
    parser.add_argument("--MAF", required=False, dest="maf_column", type=str, default='MAF', help="[OPTIONAL] MAF column name (default: MAF).")
    parser.add_argument("--sep", required=False, dest="sep", type=str, default="\t", help="What separator do you use in your file?")
    parser.add_argument("--retain", required=False, dest="retain", nargs="+", type=str, help="Columns to add to output file")
    parser.add_argument("--no_SNP_number_correction", required=False, action="store_true", help="Provide this flag if you just want the raw WZA scores")
    parser.add_argument("--empiricalP", required=False, action="store_true", help="Flag for empirical p-values")
    parser.add_argument("--spline_degree", required=False, dest="spline_degree", type=int, default=7, help="[OPTIONAL] Degree of the polynomial for spline correction (default: 7).")
    parser.add_argument("--rolling_window_size", required=False, dest="rolling_window_size", type=int, default=50, help="[OPTIONAL] Window size for rolling mean/variance in spline correction (default: 50).")
    parser.add_argument("--min_rolling_entries", required=False, dest="min_rolling_entries", type=int, default=10, help="[OPTIONAL] Minimum entries for rolling calculation in spline correction (default: 10).")

    args = parser.parse_args()

    # Print arguments if verbose mode is enabled
    if args.verbose:
        print("\nArguments passed:")
        for arg in vars(args):
            print(f"{arg}: {getattr(args, arg)}")

    # --- Load and Prepare Data ---
    print("Loading and preparing data...")
    try:
        csv = pd.read_csv(args.correlations, sep=args.sep, engine="python")
    except FileNotFoundError:
        print(f"Error: Input file not found at {args.correlations}")
        sys.exit(1)
    except Exception as e:
        print(f"Error reading input file {args.correlations}: {e}")
        sys.exit(1)

    # Validate required columns
    required_cols = [args.summary_stat, args.window]
    if args.maf_column != 'MAF':
        required_cols.append(args.maf_column)
    if args.retain:
        required_cols.extend(args.retain)

    for col in required_cols:
        if col not in csv.columns:
            print(f"Error: Required column '{col}' not found in the input file.")
            print(f"Available columns: {list(csv.columns)}")
            sys.exit(1)

    # Ensure MAF column is correctly named 'MAF' for the WZA function
    if args.maf_column != 'MAF':
        csv['MAF'] = csv[args.maf_column].copy()

    # Apply MAF filter
    maf_filter_value = 0.05 # Hardcoded MAF filter as seen in original script
    initial_rows = csv.shape[0]
    csv = csv[csv['MAF'] > maf_filter_value]
    if args.verbose and csv.shape[0] < initial_rows:
        print(f"Applied MAF filter > {maf_filter_value}. Removed {initial_rows - csv.shape[0]} SNPs.")

    # Prepare 'pVal' column based on summary statistic and flags
    if args.empiricalP:
        # Assuming empirical p-values are already scaled between 0 and 1
        csv["pVal"] = csv[args.summary_stat].copy()
    else:
        # Rank-based p-values if not empirical
        if args.large_i_small_p:
            # Rank such that larger values get smaller p-values (rank from largest to smallest)
            if args.summary_stat == "RDA": # Special handling for RDA as seen in original script
                 # Assuming RDA statistic, rank (value^2) from largest to smallest
                 csv["pVal"] = 1 - (csv[args.summary_stat] ** 2).rank(method='average', ascending=True) / csv.shape[0]
            else:
                 # Rank other statistics (larger value = smaller p-value)
                 csv["pVal"] = 1 - csv[args.summary_stat].rank(method='average', ascending=True) / csv.shape[0]
        else:
            # Rank such that smaller values get smaller p-values (rank from smallest to largest)
            csv["pVal"] = csv[args.summary_stat].rank(method='average', ascending=True) / csv.shape[0]

    # Filter out rows where the window identifier is 'None'
    initial_rows = csv.shape[0]
    csv_genes = csv[csv[args.window].astype(str) != "None"].copy() # Ensure window column is string for comparison
    if args.verbose and csv_genes.shape[0] < initial_rows:
         print(f"Removed {initial_rows - csv_genes.shape[0]} rows with window ID 'None'.")

    # Group data by window/block identifier
    csv_gb_gene = csv_genes.groupby(args.window)

    # Determine the maximum number of SNPs to use per window (for optional downsampling)
    max_SNP_count = args.sample_snps
    if max_SNP_count == -1:
        # Calculate 75th percentile of SNP counts among windows with at least min_snps
        num_SNP_list = np.array([g[1].shape[0] for g in csv_gb_gene if g[1].shape[0] >= args.min_snps])
        if num_SNP_list.size > 0:
             max_SNP_count = int(np.percentile(num_SNP_list, 75))
             if args.verbose:
                print(f"Using the 75th percentile number of SNPs ({max_SNP_count}) as the maximum in each gene.")
        else:
             max_SNP_count = 1e6 # Default to a large number if no windows meet the criteria
             if args.verbose:
                 print("No windows meet the minimum SNP count for percentile calculation. No downsampling will be applied.")
    elif max_SNP_count == 0:
        max_SNP_count = 1e9 # Use a very large number to effectively disable downsampling
        if args.verbose:
            print("Downsampling disabled.")
    else:
         if args.verbose:
            print(f"Downsampling to a maximum of {max_SNP_count} SNPs per window.")

    # --- Perform WZA and Top Candidate Test per Window ---
    print("Calculating WZA and Top Candidate p-values per window...")
    all_genes_results = []
    
    total_windows = len(csv_gb_gene)
    for i, (gene, gene_df) in enumerate(csv_gb_gene):
        # Filter windows with fewer than the minimum required SNPs
        if gene_df.shape[0] < args.min_snps:
            if args.verbose:
                 print(f"Skipping window '{gene}' with {gene_df.shape[0]} SNPs (less than min_snps={args.min_snps}).")
            continue

        original_snp_count = gene_df.shape[0]  # Track the original number of SNPs

        # Perform WZA, potentially with downsampling and resampling
        if original_snp_count <= max_SNP_count:
            # No downsampling needed for this window
            wza_score = WZA(gene_df, "pVal", MAF_filter=maf_filter_value)
            snp_count_used = original_snp_count
        else:
            # Downsampling and resampling required
            if args.verbose:
                print(f"Downsampling window '{gene}' ({original_snp_count} SNPs) to {max_SNP_count} SNPs with {args.resamples} resamples.")
            # Perform WZA multiple times on resampled SNPs and take the mean score
            wza_scores_resampled = []
            for _ in range(args.resamples):
                # Sample max_SNP_count SNPs without replacement
                gene_df_sampled = gene_df.sample(max_SNP_count, replace=False)
                wza_score_resample = WZA(gene_df_sampled, "pVal", MAF_filter=maf_filter_value)
                if not pd.isna(wza_score_resample):
                    wza_scores_resampled.append(wza_score_resample)
            
            if wza_scores_resampled:
                 wza_score = np.mean(wza_scores_resampled)
            else:
                 wza_score = np.nan # If all resamples resulted in NaN

            snp_count_used = max_SNP_count

        if pd.isna(wza_score) and args.verbose:
            print(f"Warning: WZA calculation resulted in NaN for window: {gene}")

        # Perform Top Candidate test
        # Calculate the p-value threshold for the top candidate test based on the percentile
        # We need the actual p-value corresponding to the percentile across ALL SNPs, not just in this window.
        # This requires having access to all p-values before grouping.
        # A simpler approach, as implied by the original code structure, might be to use a fixed small p-value threshold for the top candidate test.
        # However, the original code uses a percentile threshold and applies it to the p-values *within* the window.
        # Let's reinterpret the original logic: thresh is the percentile of p-values *within* the window.
        # This doesn't seem right for a standard top-candidate test where the threshold is typically fixed (e.g., 0.001).
        # Let's assume the intended logic for `top_candidate_threshold` (e.g., 99) means the top 1% of *all* SNP p-values.
        # We need the p-value corresponding to the (100 - threshold) percentile of the full dataset's p-values.

        # Recalculate the actual p-value threshold for the top candidate test based on the percentile across the full dataset
        if 0 <= args.top_candidate_threshold <= 100:
            # Calculate the p-value value at the specified percentile (e.g., 1st percentile for threshold=99)
            top_candidate_pvalue_threshold = np.percentile(csv['pVal'], 100 - args.top_candidate_threshold)
            top_candidate_p_val, hits = top_candidate(gene_df, top_candidate_pvalue_threshold)
        else:
            print(f"Warning: Invalid top_candidate_threshold ({args.top_candidate_threshold}). Skipping top candidate test for window {gene}.")
            top_candidate_p_val = np.nan
            hits = np.nan

        # Store results for the current window
        output = {
            args.window: gene, # Use the original window column name
            "SNPs": snp_count_used,
            "original_SNPs": original_snp_count, # Keep original count
            "hits": hits,
            "Z": wza_score,
            "top_candidate_p": top_candidate_p_val
        }

        # Add retained columns
        if args.retain:
            for r in args.retain:
                if r in gene_df.columns:
                    if is_string_dtype(gene_df[r]):
                        # Take the first value for string columns
                        output[r] = gene_df[r].iloc[0] if not gene_df.empty else None
                    elif is_numeric_dtype(gene_df[r]):
                        # Take the mean for numeric columns
                        output[r] = gene_df[r].mean() if not gene_df.empty else np.nan
                else:
                     output[r] = np.nan # Column not found in window df (shouldn't happen if validation passed)

        all_genes_results.append(output)

        if args.verbose and (i + 1) % 100 == 0:
             print(f"Processed {i + 1}/{total_windows} windows...")

    print("Finished calculating WZA and Top Candidate p-values.")

    # --- Create DataFrame and Apply Correction ---
    wza_results_df = pd.DataFrame(all_genes_results)

    # Apply SNP number correction unless disabled or no variation in SNP count
    if not args.no_SNP_number_correction:
        if wza_results_df.empty or wza_results_df["SNPs"].isnull().all() or wza_results_df["SNPs"].var() == 0:
            if args.verbose:
                print("Skipping SNP number correction: No variation in SNP count among windows or empty results.")
            wza_results_df["Z_pVal"] = np.nan # Add Z_pVal column with NaNs if correction skipped
        else:
            print("Applying spline-based SNP number correction...")
            try:
                 # Use the specified spline parameters
                 wza_results_df = adjust_WZA_with_spline(wza_results_df, roller=args.rolling_window_size, minEntries=args.min_rolling_entries, chosen_degree=args.spline_degree)
            except Exception as e:
                 print(f"Error during spline correction: {e}. Adjusted p-values may be missing.")
                 wza_results_df["Z_pVal"] = np.nan # Add Z_pVal column with NaNs on error
    else:
         if args.verbose:
            print("SNP number correction disabled by flag.")
         wza_results_df["Z_pVal"] = np.nan # Add Z_pVal column with NaNs if correction skipped

    # --- Save Results ---
    print(f"Saving results to {args.output}...")
    try:
        wza_results_df.to_csv(args.output, index=False)
        print("Script finished successfully.")
    except Exception as e:
        print(f"Error saving output file {args.output}: {e}")
        sys.exit(1)

# --- Script Entry Point ---
if __name__ == "__main__":
    main()