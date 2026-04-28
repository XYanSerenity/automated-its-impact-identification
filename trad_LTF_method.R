library(lmtest)
library(ggplot2)

# Parameter settings
true_b <- 0
n <- 100
x_p <- c(rep(0, 50), 1, rep(0, 49))  # Pulse intervention at t = 51

# Read the simulated data
# Each column represents one simulated time series, with 100 simulations in total
df <- read.csv("It=xtb=0_h=1_r=1_w0=0.5_r1=0.6_data.csv")
num_sim <- ncol(df)

# Function for fitting the LTF model
LTFxpFun <- function(y, x_p, rlist) {
  lagged_xp <- embed(x_p, 11)  # Construct lagged intervention variables from lag 0 to lag 10
  colnames(lagged_xp) <- paste0("xpv", 0:10)
  
  y_clean <- tail(y, nrow(lagged_xp))
  
  TFpt <- try(
    arima(y_clean, order = rlist, xreg = lagged_xp, method = "ML"),
    silent = TRUE
  )
  
  if (inherits(TFpt, "try-error")) {
    return(NULL)
  } else {
    return(TFpt)
  }
}

# Function for identifying the delay time b
# The delay is determined based on the first coefficient exceeding the peak-based threshold
get_b_from_model <- function(modelxp) {
  coefpt <- coeftest(modelxp)
  ptvcoef <- tail(coefpt, 11)  # Extract the last 11 coefficients corresponding to intervention variables
  
  coefs <- ptvcoef[, 1]
  p_vals <- ptvcoef[, 4]
  
  peak <- max(abs(coefs), na.rm = TRUE)
  threshold <- peak / 3
  
  is_insig_or_small <- abs(coefs) < threshold
  
  b_est <- 0
  
  for (i in seq_along(is_insig_or_small)) {
    if (!is_insig_or_small[i]) break
    b_est <- b_est + 1
  }
  
  return(b_est)
}

# Run the delay identification procedure
b_results <- rep(NA, num_sim)
rlistp0 <- c(1, 0, 0)

for (j in 1:num_sim) {
  yj <- df[[j]]  # The j-th column represents one simulated time series
  
  TFptMod <- LTFxpFun(yj, x_p, rlistp0)
  
  if (!is.null(TFptMod)) {
    b_results[j] <- get_b_from_model(TFptMod)
  }
}

# Summary statistics
b_results_clean <- na.omit(b_results)
table_b <- table(b_results_clean)
accuracy <- mean(b_results_clean == true_b)

print(table_b)
print(accuracy)

# Visualization of the identified delay time
df_b <- data.frame(
  Estimated_b = factor(b_results_clean, levels = 0:10)
)

ggplot(df_b, aes(x = Estimated_b)) +
  geom_bar(fill = "steelblue") +
  geom_vline(
    xintercept = true_b + 1,
    linetype = "dashed",
    color = "red",
    size = 1
  ) +
  labs(
    title = "Distribution of Estimated Delay Time b Based on the LTF Model",
    subtitle = "Red dashed line indicates the true delay time b = 0",
    x = "Estimated Delay Time b",
    y = "Frequency"
  ) +
  theme_minimal()