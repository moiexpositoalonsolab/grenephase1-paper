# Define the custom theme function
theme_grenenet <- function() {
  theme_minimal() +
    theme(
      axis.title = element_text(size = 12, color = "black"),
      axis.text = element_text(size = 12, color = "black"),
      legend.text = element_text(size = 12, color = "black"),
      legend.title = element_text(size = 12, color = "black")
    )
}
