
# Step 0: Example dataset of site #24
d <- data.frame(
  freq = c(0.005342900, 0.063700882, 0.004166750, 0.001089000, 0.007331319, 
           0.058988000, 0.058988000, 0.058988000, 0.002227000, 0.058988000, 
           0.058988000, 0.011098368, 0.009820083, 0.058988000, 0.007058000, 
           0.007612279, 0.058988000, 0.014579571, 0.001838000, 0.058988000, 
           0.058988000, 0.058988000),
  startfreq = rep(0.058988, 22),  # Same initial frequency for all rows
  flowers = c(20, 51, 4, 7, 75, NA, NA, NA, 11, NA, NA, 38, 48, NA, 13, 43, NA, 56, 17, NA, NA, NA),
  year = c(1, 1, 2, 2, 1, 0, 0, 0, 1, 0, 0, 1, 1, 0, 1, 1, 0, 1, 1, 0, 0, 0)
)

# Load necessary library
library(dplyr)

# Fit the binomial regression model
mymod <- 
  d %>%
  # Step 1: Pre-process the dataset
  # - Calculate "success" as the number of successful events (alleles) by multiplying `freq` and `flowers`
  # - Define "total" as the total number of trials (flowers)
  # - Convert `year` to numeric for regression analysis
  dplyr::mutate(
    success = round(freq * flowers),  # Calculate the number of successful events
    total = flowers,                  # Define the total number of trials
    year = as.numeric(year)           # Convert year to numeric
  ) %>%
  # Step 2: Fit the binomial regression model
  # - Use the glm function to model the success rate as a function of year
  # - Use a quasibinomial family to account for overdispersion
  glm(
    data = .,                         # Use the modified dataset
    cbind(success, total) ~ year,     # Formula: success and failures modeled against year
    family = quasibinomial(link = logit)  # Logistic regression with quasibinomial link
  )

# Step 3: Summarize the fitted model
# - This provides the coefficients, standard errors, and statistical significance of the model
summary(mymod)

# Define a function to calculate the change in allele frequency
calculate_frequency_change <- function(p0, beta) {
  # Step 1: Convert initial frequency to odds
  odds_0 <- p0 / (1 - p0)
  
  # Step 2: Calculate new odds after applying the beta (regression coefficient)
  odds_1 <- odds_0 * exp(beta)
  
  # Step 3: Convert new odds back to frequency
  p1 <- odds_1 / (1 + odds_1)
  
  # Step 4: Calculate frequency change
  delta_p <- p1 - p0
  
  # Return the frequency change
  return(delta_p) #
  # return(odds_1) # uncomment to instead return odds_1
  # return(exp(beta)) # uncomment to instead return odds_1
}

# Step 4: Extract the regression coefficient for `year`
# - This coefficient represents the change in log-odds of allele frequency per year
beta <- mymod$coefficients[2]

# Step 5: Define the initial frequency
# - Use the first value of `startfreq` as the initial frequency (p0)
p0 <- d$startfreq[1]

# Step 6: Calculate the frequency change using the custom function
frequency_change <- calculate_frequency_change(p0, beta)

# Step 7: Print the frequency change
# - This represents the expected change in allele frequency per year
print(frequency_change)




