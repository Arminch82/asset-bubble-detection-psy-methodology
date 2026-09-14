#===============================================================================
# Apple (AAPL)
# BUbble Detection and Regime Behavior Analysis
#===============================================================================

# 1. Packages ------------------------------------------------------------------
library(tidyquant)
library(exuber)
library(dplyr)
library(ggplot2)
library(tidyr)
library(xts)
library(psymonitor)
library(stargazer)
library(tidyverse)
library(quantmod)
library(tseries)
library(readr)
library(cointReg)
library(frenchdata)
library(lubridate)
library(rlang)

# 2. Severity function ---------------------------------------------------------
calculate_bubble_severity <- function(data, date_col, value_col, crisis_dates = NULL, conf_level = 0.90) {
  
  # TOGGLE: If no crisis dates are provided (e.g., for full period), use the entire dataset's date range
  if (is.null(crisis_dates)) {
    crisis_dates <- data.frame(
      start = min(data[[date_col]], na.rm = TRUE),
      end = max(data[[date_col]], na.rm = TRUE)
    )
  }
  
  alpha_level <- 1 - conf_level
  
  # Initialize storage vectors
  theta_n_list <- numeric(nrow(crisis_dates))
  gamma_n_list <- numeric(nrow(crisis_dates))
  ci_lower_list <- numeric(nrow(crisis_dates))
  ci_upper_list <- numeric(nrow(crisis_dates))
  n_obs_list <- numeric(nrow(crisis_dates))
  
  # Loop through either the multiple crisis dates or the single full period
  for (i in seq_len(nrow(crisis_dates))) {
    start_date <- crisis_dates$start[i]
    end_date <- crisis_dates$end[i]
    
    # Filter data for the specific window
    bubble_data <- data %>%
      filter(!!sym(date_col) >= start_date & !!sym(date_col) <= end_date)
    
    n <- nrow(bubble_data)
    n_obs_list[i] <- n
    
    # Need at least 3 observations to calculate lags and have a valid T_obs > 1
    if (n > 2) {
      y <- bubble_data[[value_col]]
      y <- y[!is.na(y)]
      
      # T_obs is the sample size minus 1 (number of transitions)
      T_obs <- length(y) - 1
      
      y_t <- y[2:(T_obs + 1)]
      y_t_minus_1 <- y[1:T_obs]
      
      # 1. Point Estimates (Phillips 2023)
      theta_hat <- sum(y_t * y_t_minus_1) / sum(y_t_minus_1^2)
      gamma_n <- -log(abs(theta_hat - 1)) / log(T_obs)
      
      # 2. Confidence Interval Calculation
      if (theta_hat > 1) {
        # Mildly Explosive Root (MER) -> Cauchy distribution
        cv <- qcauchy(1 - alpha_level / 2)
        margin_error <- cv * (2 / ((1 + 1 / (T_obs^gamma_n))^T_obs * log(T_obs)))
      } else {
        # Mildly Integrated Root (MIR) -> Normal distribution
        cv <- qnorm(1 - alpha_level / 2)
        margin_error <- cv * (sqrt(2) / (T_obs^((1 - gamma_n) / 2) * log(T_obs)))
      }
      
      theta_n_list[i] <- theta_hat
      gamma_n_list[i] <- gamma_n
      
      # Rate boundary defaults to 0 on the lower end
      ci_lower_list[i] <- max(0, gamma_n - margin_error) 
      ci_upper_list[i] <- gamma_n + margin_error
      
    } else {
      theta_n_list[i] <- NA
      gamma_n_list[i] <- NA
      ci_lower_list[i] <- NA
      ci_upper_list[i] <- NA
    }
  }
  
  # Compile final results
  results <- crisis_dates %>%
    mutate(
      Duration_Months = n_obs_list,
      Theta_n = round(theta_n_list, 4),
      Severity_Gamma_n = round(gamma_n_list, 4),
      CI_Lower = round(ci_lower_list, 4),
      CI_Upper = round(ci_upper_list, 4),
      CI_Length = round(ci_upper_list - ci_lower_list, 4)
    )
  
  return(results)
}

# 3. Data collection & Cleaning ------------------------------------------------
prices <- tq_get(c("AAPL", "^GSPC", "^IXIC"), get = "stock.prices",
                 from = "2010-01-01", to = Sys.Date(),
                 periodicity = "monthly") |>
  select(symbol, date, adjusted)

prices_wide <- prices|>
  pivot_wider(names_from = symbol, values_from = adjusted) |>
  rename(apple = AAPL, sp500 = `^GSPC`, nasdaq = `^IXIC`) |>
  drop_na()

prices_clean <- prices_wide

# Revenue
revenue <- read.csv("apple_revenue.csv") |>
  mutate(date = as.Date(date))

prices_wide <- prices_wide |>
  left_join(revenue, by = "date") |>
  arrange(date) |>
  fill(sales, .direction = "down") |>
  drop_na()

# Sample split
period_1 <- prices_wide |>
  filter(date <= as.Date("2022-11-30")) |>
  mutate(
    apple_norm = (apple / first(apple)) * 100,
    sp500_norm = (sp500 / first(sp500)) * 100,
    price_ratio = apple_norm / sp500_norm,
    sales = sales / 1000,
    price_sales = apple / sales) 

period_2 <- prices_wide |>
  filter(date >= as.Date("2022-12-01")) |>
  mutate(
    apple_norm = (apple / first(apple)) * 100,
    sp500_norm = (sp500 / first(sp500)) * 100,
    price_ratio = apple_norm / sp500_norm,
    sales = sales / 1000,
    price_sales = apple / sales)

# 4. PSY Bubble Detection & Regime Severity on Price/Sales ---------------------
# PSY test on period 1 for P/S
y <- period_1$price_sales
obs <- length(y)
swindow0 <- floor(obs * (0.01 + 1.8 / sqrt(obs)))
yr <- 2
Tb <- 12*yr + swindow0 - 1

bsadf <- PSY(y, swindow0 = swindow0, IC = 2, adflag = 6)

quantilesBsadf <- cvPSYwmboot(y, swindow0 = swindow0, IC = 2, adflag = 6, Tb = Tb, nboot = 2000)

dim <- obs - swindow0 + 1
monitorDates <- period_1$date[swindow0:obs]
quantile95 <- quantilesBsadf %*% matrix(1, nrow = 1, ncol = dim)
ind95 <- (bsadf > t(quantile95[1, ])) * 1
periods <- locate(ind95, monitorDates)
# only if errors
bubble_dates <- monitorDates[which(ind95 == 1)]
bubble_summary <- data.frame(
  start = min(bubble_dates),
  end = max(bubble_dates),
  Duration_Months = length(bubble_dates)
)
print(bubble_summary)

apple_m1_p1 <- ggplot() +
  geom_rect(data = bubble_summary, aes(xmin = start, xmax = end,
                                       ymin = -Inf, ymax = Inf), alpha = 0.5) +
  geom_line(data = period_1, aes(date, price_sales), lwd = 1) +
  geom_vline(data = filter(bubble_summary, start == end), 
             aes(xintercept = start), 
             alpha = 0.5, color = "grey50", linewidth = 1) +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  labs(x = NULL, y = NULL) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(face = "bold", size = 14),
    axis.text.y = element_text(face = "bold", size = 14),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black")
  )
apple_m1_p1
###
crisisDates <- disp(periods, obs)
print(crisisDates)

# Calculate Severity for Price/Sales
ps_severity_summary <- calculate_bubble_severity(
  data = period_1, 
  date_col = "date", 
  value_col = "price_sales", 
  crisis_dates = bubble_dates
)
print("Price/Sales Bubble Severity:")
print(ps_severity_summary)

# Plot 1
apple_m1_p1<- ggplot() + 
  geom_rect(data = crisisDates, aes(xmin = start, xmax = end,
                                    ymin = -Inf, ymax = Inf), alpha = 0.5) +
  geom_vline(data = filter(crisisDates, start == end), 
             aes(xintercept = start), 
             alpha = 0.5, color = "grey50", linewidth = 1) +
  geom_line(data = period_1, aes(date, price_sales), lwd = 1) +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  labs(x = NULL, y = NULL) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(face = "bold", size = 14),
    axis.text.y = element_text(face = "bold", size = 14),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black")
  )
apple_m1_p1
ggsave("apple_price_sales_p1.pdf", plot = apple_m1_p1, width = 6, height = 4)

# Period 2
y_p2_ps <- as.numeric(period_2$price_sales)
radf_result <- radf(y_p2_ps, minw = 18)
cv_gsadf    <- radf_wb_cv(y_p2_ps, nboot = 2000, minw = 18)
summary(radf_result, cv_gsadf)

# Calculate Severity for Price/Sales Period 2
ps_severity_summary_p2 <- calculate_bubble_severity(
  data = period_2, 
  date_col = "date", 
  value_col = "price_sales"
)

print("Period 2: Price/Sales Bubble Severity (Whole Sample):")
print(ps_severity_summary_p2)

# 5. PSY Bubble Detection & Regime Severity on Fama-French 3 Factor IC ---------
prices_clean <- prices_clean |>
  mutate(
    delta_p = log(apple) - lag(log(apple)),
    delta_n = log(nasdaq) - lag(log(nasdaq))
  ) |>
  drop_na(delta_p)

french_data <- get_french_data_list()
ff_3_factors <- download_french_data("Fama/French 3 Factors")
monthly_ff_3_factors <- ff_3_factors$subsets$data[[1]]
monthly_ff_3_factors <- monthly_ff_3_factors |>
  filter(date > 201001) |>
  mutate(
    date = ymd(paste0(date, "01"))
  )
monthly_ff_3_factors <- monthly_ff_3_factors |>
  left_join(prices_clean, by = "date")
monthly_ff_3_factors <- monthly_ff_3_factors |>
  mutate(
    RF_dec = RF / 100,
    SMB_dec = SMB / 100,
    HML_dec = HML / 100,
    RF_log = log(1 + RF_dec),
    SMB_log = log(1 + SMB_dec),
    HML_log = log(1 + HML_dec),
    excess_return = delta_n - RF_log
  )
training_data <- monthly_ff_3_factors |>
  filter(date >= as.Date("2010-01-01") & date <= as.Date("2014-12-01"))

ff_3_reg <- lm(delta_p ~ 0 + excess_return + SMB_log + HML_log, data = training_data)
monthly_ff_3_factors$fitted_delta_p <- predict(ff_3_reg, newdata = monthly_ff_3_factors)

monthly_ff_3_factors <- monthly_ff_3_factors |>
  arrange(date) |>
  mutate(
    fitted_delta_p_adj = if_else(row_number() == 1, 0, fitted_delta_p),
    M_hat = first(log(apple)) + cumsum(fitted_delta_p_adj),
    IC = log(apple) - M_hat
  )
monthly_ff_3_factors_p1 <- monthly_ff_3_factors |>
  filter(date <= as.Date("2022-11-01"))
monthly_ff_3_factors_p2 <- monthly_ff_3_factors |> filter(date >= as.Date("2022-12-01"))

# PSY test on the IC components
y <- as.numeric(monthly_ff_3_factors_p1$IC)
obs <- length(y)
swindow0 <- floor(0.01 * obs + 1.8 * sqrt(obs))
yr <- 2
Tb <- 12*yr + swindow0 - 1

bsadf <- PSY(y, swindow0 = swindow0, IC = 2, adflag = 6)

quantilesBsadf <- cvPSYwmboot(y, swindow0 = swindow0, IC = 2, adflag = 6, Tb = Tb, nboot = 2000)

dim <- obs - swindow0 + 1
monitorDates <- monthly_ff_3_factors_p1$date[swindow0:obs]
quantile95 <- quantilesBsadf %*% matrix(1, nrow = 1, ncol = dim)
ind95 <- (bsadf > t(quantile95[1, ])) * 1
periods <- locate(ind95, monitorDates)
crisisDates <- disp(periods, obs)
print(crisisDates)

# Calculate Severity for FF3 ICs
ff3_severity_summary <- calculate_bubble_severity(
  data = monthly_ff_3_factors_p1, 
  date_col = "date", 
  value_col = "IC", 
  crisis_dates = crisisDates
)
print("Fama-French 3-Factor IC Bubble Severity:")
print(ff3_severity_summary)

# Plot the bubbles
apple_ic_p1<- ggplot() + 
  geom_rect(data = crisisDates, aes(xmin = start, xmax = end,
                                    ymin = -Inf, ymax = Inf), alpha = 0.5) +
  geom_vline(data = filter(crisisDates, start == end), 
             aes(xintercept = start), 
             alpha = 0.5, color = "grey50", linewidth = 1) +
  geom_line(data = monthly_ff_3_factors_p1, aes(date, IC), lwd = 1) +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  labs(x = NULL, y = NULL) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(face = "bold", size = 14),
    axis.text.y = element_text(face = "bold", size = 14),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black")
  )
apple_ic_p1
ggsave("apple_ic_p1.pdf", plot = apple_ic_p1, width = 6, height = 4)

# Period 2
monthly_ff_3_factors_p2 <- monthly_ff_3_factors |> filter(date >= as.Date("2022-12-01"))
y_p2_ff3 <- as.numeric(monthly_ff_3_factors_p2$IC)
radf_result <- radf(y_p2_ff3, minw = 18)
cv_gsadf    <- radf_wb_cv(y_p2_ff3, nboot = 2000, minw = 18)
summary(radf_result, cv_gsadf)

# Calculate Severity for Fama-French 3 factor IC Period 2
ff3_severity_p2 <- calculate_bubble_severity(
  data = monthly_ff_3_factors_p2, 
  date_col = "date", 
  value_col = "IC"
)

print("Period 2: Fama-French 3-Factor IC Bubble Severity (Whole Sample):")
print(ff3_severity_p2)

# 6. PSY Bubble Detection & Regime Severity on Fama-French 5 Factor IC ---------
prices_clean <- prices_clean |>
  mutate(
    delta_p = log(apple) - lag(log(apple)),
    delta_n = log(nasdaq) - lag(log(nasdaq))
  ) |>
  drop_na(delta_p)

french_data <- get_french_data_list()
ff_5_factors <- download_french_data("Fama/French 5 Factors (2x3)")
monthly_ff_5_factors <- ff_5_factors$subsets$data[[1]]
monthly_ff_5_factors <- monthly_ff_5_factors |>
  filter(date > 201001) |>
  mutate(
    date = ymd(paste0(date, "01"))
  )
monthly_ff_5_factors <- monthly_ff_5_factors |>
  left_join(prices_clean, by = "date")
monthly_ff_5_factors <- monthly_ff_5_factors |>
  mutate(
    RF_dec = RF / 100,
    SMB_dec = SMB / 100,
    HML_dec = HML / 100,
    RMW_dec = RMW / 100,
    CMA_dec = CMA / 100,
    RF_log = log(1 + RF_dec),
    excess_return = delta_n - RF_log
  )
monthly_ff_5_factors <- monthly_ff_5_factors |>
  arrange(date) |>
  drop_na(apple, nasdaq, delta_p, delta_n, RF_log, SMB_dec, HML_dec, excess_return, 
          RMW_dec, CMA_dec)
training_data <- monthly_ff_5_factors |>
  filter(date >= as.Date("2010-01-01") & date <= as.Date("2014-12-01"))

ff_5_reg <- lm(delta_p ~ 0 + excess_return + SMB_dec + HML_dec + RMW_dec + CMA_dec, data = training_data)
monthly_ff_5_factors$fitted_delta_p <- predict(ff_5_reg, newdata = monthly_ff_5_factors)

monthly_ff_5_factors <- monthly_ff_5_factors |>
  arrange(date) |>
  mutate(
    fitted_delta_p_adj = if_else(row_number() == 1, 0, fitted_delta_p),
    M_hat = first(log(apple)) + cumsum(fitted_delta_p_adj),
    IC = log(apple) - M_hat
  )

monthly_ff_5_factors_p1 <- monthly_ff_5_factors |>
  filter(date <= as.Date("2022-11-01"))
monthly_ff_5_factors_p2 <- monthly_ff_5_factors |> filter(date >= as.Date("2022-12-01"))

# PSY test on the IC components
y <- as.numeric(monthly_ff_5_factors_p1$IC)
obs <- length(y)
swindow0 <- floor(0.01 * obs + 1.8 * sqrt(obs))
yr <- 2
Tb <- 12*yr + swindow0 - 1

bsadf <- PSY(y, swindow0 = swindow0, IC = 2, adflag = 6)

quantilesBsadf <- cvPSYwmboot(y, swindow0 = swindow0, IC = 2, adflag = 6, Tb = Tb, nboot = 2000)

dim <- obs - swindow0 + 1
monitorDates <- monthly_ff_5_factors_p1$date[swindow0:obs]
quantile95 <- quantilesBsadf %*% matrix(1, nrow = 1, ncol = dim)
ind95 <- (bsadf > t(quantile95[1, ])) * 1
periods <- locate(ind95, monitorDates)
crisisDates <- disp(periods, obs)
print(crisisDates)

# Calculate Severity for FF5 ICs
ff5_severity_summary <- calculate_bubble_severity(
  data = monthly_ff_5_factors_p1, 
  date_col = "date", 
  value_col = "IC", 
  crisis_dates = crisisDates
)
print("Fama-French 5-Factor IC Bubble Severity:")
print(ff5_severity_summary)

# Plot the bubbles
apple_ic2_p1<- ggplot() + 
  geom_line(data = monthly_ff_5_factors_p1, aes(date, IC), lwd = 1) +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  labs(x = NULL, y = NULL) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(face = "bold", size = 14),
    axis.text.y = element_text(face = "bold", size = 14),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black")
  )
apple_ic2_p1
ggsave("apple_ic2_p1.pdf", plot = apple_ic2_p1, width = 6, height = 4)

# Period 2
monthly_ff_5_factors_p2 <- monthly_ff_5_factors |> filter(date >= as.Date("2022-12-01"))
y_p2_ff5 <- as.numeric(monthly_ff_5_factors_p2$IC)
radf_result <- radf(y_p2_ff5, minw = 18)
cv_gsadf    <- radf_wb_cv(y_p2_ff5, nboot = 2000, minw = 18)
summary(radf_result, cv_gsadf)

# severity period 2
ff5_severity_p2 <- calculate_bubble_severity(
  data = monthly_ff_5_factors_p2, 
  date_col = "date", 
  value_col = "IC"
)

print("Period 2: Fama-French 5-Factor IC Bubble Severity (Whole Sample):")
print(ff5_severity_p2)

# 7. PSY Bubble Detection & Regime Severity on P/I -----------------------------
y <- period_1$price_ratio
obs <- length(y)
swindow0 <- floor(obs * (0.01 + 1.8 / sqrt(obs)))
yr <- 2
Tb <- 12*yr + swindow0 - 1

bsadf <- PSY(y, swindow0 = swindow0, IC = 2, adflag = 6)

quantilesBsadf <- cvPSYwmboot(y, swindow0 = swindow0, IC = 2, adflag = 6, Tb = Tb, nboot = 2000)

dim <- obs - swindow0 + 1
monitorDates <- period_1$date[swindow0:obs]
quantile95 <- quantilesBsadf %*% matrix(1, nrow = 1, ncol = dim)
ind95 <- (bsadf > t(quantile95[1, ])) * 1
periods <- locate(ind95, monitorDates)
# only if errors
bubble_dates <- monitorDates[which(ind95 == 1)]
bubble_summary <- data.frame(
  start = min(bubble_dates),
  end = max(bubble_dates),
  Duration_Months = length(bubble_dates)
)
print(bubble_summary)

apple_m2_p1 <- ggplot() +
  geom_rect(data = bubble_summary, aes(xmin = start, xmax = end,
                                       ymin = -Inf, ymax = Inf), alpha = 0.5) +
  geom_line(data = period_1, aes(date, price_ratio), lwd = 1) +
  geom_vline(data = filter(bubble_summary, start == end), 
             aes(xintercept = start), 
             alpha = 0.5, color = "grey50", linewidth = 1) +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  labs(x = NULL, y = NULL) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(face = "bold", size = 14),
    axis.text.y = element_text(face = "bold", size = 14),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black")
  )
apple_m2_p1
###
crisisDates <- disp(periods, obs)
print(crisisDates)
# Plot 1
apple_m2_p1<- ggplot() + 
  geom_rect(data = crisisDates, aes(xmin = start, xmax = end,
                                    ymin = -Inf, ymax = Inf), alpha = 0.5) +
  geom_vline(data = filter(crisisDates, start == end), 
             aes(xintercept = start), 
             alpha = 0.5, color = "grey50", linewidth = 1) +
  geom_line(data = period_1, aes(date, price_ratio), lwd = 1) +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  labs(x = NULL, y = NULL) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(face = "bold", size = 14),
    axis.text.y = element_text(face = "bold", size = 14),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black")
  )
apple_m2_p1
ggsave("apple_price_ratio_p1.pdf", plot = apple_m2_p1, width = 6, height = 4)

# Calculate Severity for Price/Sales
ps_severity_summary <- calculate_bubble_severity(
  data = period_1, 
  date_col = "date", 
  value_col = "price_ratio", 
  crisis_dates = crisisDates
)
print("Price/Index Bubble Severity:")
print(ps_severity_summary)

# Period 2
y_p2_pi <- as.numeric(period_2$price_ratio)
radf_result <- radf(y_p2_pi, minw = 18)
cv_gsadf    <- radf_wb_cv(y_p2_pi, nboot = 2000, minw = 18)
summary(radf_result, cv_gsadf)

# period 2 severity
pi_severity_p2 <- calculate_bubble_severity(
  data = period_2, 
  date_col = "date", 
  value_col = "price_ratio"
)

print("Period 2: Price/Index Bubble Severity (Whole Sample):")
print(pi_severity_p2)