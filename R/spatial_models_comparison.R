# ==============================================================================
# Spatial Model Comparison: OLS, SAR (Spatial Lag), and SER (Spatial Error)
# ==============================================================================

# --- 1. Load Required Libraries -------------------------------------------
library(spdep)
library(spatialreg)
library(xtable)

# --- 2. Compute Residuals -------------------------------------------------
ols_resid <- residuals(ols1)
sar_resid <- residuals(sarr)  # Spatial Lag Model
ser_resid <- residuals(serr)  # Spatial Error Model

# --- 3. Moran's I Test on Residuals ---------------------------------------
# Test for remaining spatial autocorrelation in residuals
moran_ols <- moran.test(ols_resid, yieldlist)
moran_sar <- moran.test(sar_resid, yieldlist)
moran_ser <- moran.test(ser_resid, yieldlist)

# --- 4. Compute Performance Metrics ---------------------------------------
# Define helper functions
rmse <- function(y, yhat) sqrt(mean((y - yhat)^2))
r2_score <- function(y, yhat) 1 - sum((y - yhat)^2) / sum((y - mean(y))^2)

# Get actual yield values (original scale)
y_original <- cocoa_pos$yld_pr_h

# Get log-transformed yield (if that's what was modeled)
y_log <- log(cocoa_pos$yld_pr_h)

# --- 5. Generate Predictions with Proper Back-transformation -------------

# OLS predictions (if modeling log yield)
ols_pred_log <- predict(ols1)
# Back-transform with smearing estimator
smearing_factor_ols <- mean(exp(ols_resid))
ols_pred <- exp(ols_pred_log) * smearing_factor_ols

# SAR predictions (if modeling log yield)
sar_pred_log <- y_log - sar_resid
# Back-transform with smearing estimator
smearing_factor_sar <- mean(exp(sar_resid))
sar_pred <- exp(sar_pred_log) * smearing_factor_sar

# SER predictions (if modeling log yield)
ser_pred_log <- y_log - ser_resid
# Back-transform with smearing estimator for SER
smearing_factor_ser <- mean(exp(ser_resid))
ser_pred <- exp(ser_pred_log) * smearing_factor_ser

# --- Alternative: If NOT using log transformation ------------------------
# Uncomment these lines if your models are NOT on log scale:
# ols_pred <- predict(ols1)
# sar_pred <- y_original - sar_resid
# ser_pred <- y_original - ser_resid

# --- 6. Create Summary Table ----------------------------------------------
model_comparison <- data.frame(
  Model = c("OLS", "SAR (Spatial Lag)", "SER (Spatial Error)"),
  AIC = c(AIC(ols1), AIC(sarr), AIC(serr)),
  BIC = c(BIC(ols1), BIC(sarr), BIC(serr)),
  LogLik = c(logLik(ols1)[1], logLik(sarr)[1], logLik(serr)[1]),
  RMSE = c(
    rmse(y_original, ols_pred),
    rmse(y_original, sar_pred),
    rmse(y_original, ser_pred)
  ),
  R2 = c(
    r2_score(y_original, ols_pred),
    r2_score(y_original, sar_pred),
    r2_score(y_original, ser_pred)
  ),
  Moran_I = c(
    moran_ols$estimate["Moran I statistic"],
    moran_sar$estimate["Moran I statistic"],
    moran_ser$estimate["Moran I statistic"]
  ),
  Moran_pvalue = c(
    moran_ols$p.value,
    moran_sar$p.value,
    moran_ser$p.value
  ),
  N_obs = c(
    length(ols_resid),
    length(sar_resid),
    length(ser_resid)
  ),
  Smearing_Factor = c(
    smearing_factor_ols,
    smearing_factor_sar,
    smearing_factor_ser
  )
)

# --- 7. Round Numeric Columns ---------------------------------------------
model_comparison_rounded <- model_comparison
numeric_cols <- sapply(model_comparison_rounded, is.numeric)
model_comparison_rounded[numeric_cols] <- lapply(
  model_comparison_rounded[numeric_cols], 
  round, 
  4
)

# --- 8. Print Model Comparison Summary ------------------------------------
cat("\n=== 📊 Spatial Model Comparison Summary ===\n\n")
print(model_comparison_rounded, row.names = FALSE)

# --- 9. Rank Models by AIC -----------------------------------------------
cat("\n=== Models Ranked by AIC (lower is better) ===\n\n")
model_ranked <- model_comparison_rounded[order(model_comparison_rounded$AIC), ]
print(model_ranked, row.names = FALSE)

# --- 10. Calculate RMSE Improvement --------------------------------------
rmse_ols <- model_comparison$RMSE[1]
rmse_sar <- model_comparison$RMSE[2]
rmse_ser <- model_comparison$RMSE[3]

cat("\n=== RMSE Improvement Over OLS ===\n")
cat(sprintf("SAR: %.2f%% reduction\n", 
            (rmse_ols - rmse_sar) / rmse_ols * 100))
cat(sprintf("SER: %.2f%% reduction\n", 
            (rmse_ols - rmse_ser) / rmse_ols * 100))

# --- 11. Display Smearing Factors -----------------------------------------
cat("\n=== Smearing Factors (Back-transformation Adjustment) ===\n")
cat(sprintf("OLS: %.4f\n", smearing_factor_ols))
cat(sprintf("SAR: %.4f\n", smearing_factor_sar))
cat(sprintf("SER: %.4f\n", smearing_factor_ser))
cat("\nNote: Smearing factors correct for retransformation bias when\n")
cat("      converting predictions from log scale to original scale.\n")

# --- 12. Interpret Moran's I Results --------------------------------------
cat("\n=== Interpretation ===\n")
cat("Moran's I on Residuals:\n")
cat("  - Tests for remaining spatial autocorrelation in model residuals\n")
cat("  - p-value > 0.05: No significant spatial dependence (good)\n")
cat("  - p-value < 0.05: Spatial dependence remains (model may be misspecified)\n\n")

cat("Model Selection Criteria:\n")
cat("  - Lower AIC/BIC: Better model fit with parsimony penalty\n")
cat("  - Lower RMSE: Better predictive accuracy (on original scale)\n")
cat("  - Higher R²: Better explained variance\n")
cat("  - Moran's I near 0: Spatial dependence properly captured\n")

# --- 13. Statistical Significance of Spatial Parameters -------------------
cat("\n=== Spatial Parameters ===\n")

# Extract spatial parameters safely
if("rho" %in% names(sarr)) {
  cat(sprintf("SAR Rho (ρ): %.4f\n", sarr$rho))
}

if("lambda" %in% names(serr)) {
  cat(sprintf("SER Lambda (λ): %.4f\n", serr$lambda))
}

# --- 14. Export to LaTeX Table (without Smearing Factor) ------------------
# Prepare table for export (exclude smearing factor from publication table)
model_comparison_tex <- model_comparison_rounded[order(model_comparison_rounded$AIC), 
                                                 c("Model", "AIC", "BIC", "LogLik", 
                                                   "RMSE", "R2", "Moran_I", "Moran_pvalue", "N_obs")]

# Convert to LaTeX
latex_table <- xtable(
  model_comparison_tex,
  caption = "Model Performance Comparison: OLS, Spatial Autoregressive (SAR), and Spatial Error (SER) Models. RMSE calculated on original scale with smearing adjustment.",
  label = "tab:model_comparison",
  align = c("l", "l", rep("r", ncol(model_comparison_tex) - 1))
)

# Print to console (preview)
cat("\n=== LaTeX Table Preview ===\n")
print(latex_table,
      include.rownames = FALSE,
      caption.placement = "top",
      booktabs = TRUE,
      sanitize.text.function = identity)

# Export to file
print(latex_table,
      include.rownames = FALSE,
      caption.placement = "top",
      booktabs = TRUE,
      sanitize.text.function = identity,
      file = "model_comparison_table.tex")

cat("\n✓ LaTeX table exported to: model_comparison_table.tex\n")

# --- 15. Create Publication-Ready Summary Table ---------------------------
# Select key metrics for main text
summary_table <- model_comparison_tex[, c("Model", "AIC", "RMSE", "R2", "Moran_I", "Moran_pvalue")]

# Add significance stars for Moran's I
summary_table$Spatial_Autocorr <- ifelse(
  summary_table$Moran_pvalue < 0.001, "***",
  ifelse(summary_table$Moran_pvalue < 0.01, "**",
         ifelse(summary_table$Moran_pvalue < 0.05, "*", "ns"))
)

# Clean column names
colnames(summary_table) <- c("Model", "AIC", "RMSE", "R²", "Moran's I", "p-value", "Sig.")

# Export simplified version
latex_summary <- xtable(
  summary_table,
  caption = "Spatial Model Comparison Summary",
  label = "tab:model_summary",
  align = c("l", "l", rep("r", 5), "c")
)

print(latex_summary,
      include.rownames = FALSE,
      caption.placement = "top",
      booktabs = TRUE,
      file = "model_summary_table.tex")

cat("✓ Summary table exported to: model_summary_table.tex\n")

# --- 16. Likelihood Ratio Test -------------------------------------------
# Compare nested models
cat("\n=== Likelihood Ratio Tests ===\n")

# LR test: OLS vs SAR
lr_sar <- -2 * (logLik(ols1)[1] - logLik(sarr)[1])
df_sar <- 1  # One additional parameter (rho)
p_sar <- 1 - pchisq(lr_sar, df_sar)
cat(sprintf("OLS vs SAR: LR = %.4f, df = %d, p-value = %.4f\n", 
            lr_sar, df_sar, p_sar))

# LR test: OLS vs SER
lr_ser <- -2 * (logLik(ols1)[1] - logLik(serr)[1])
df_ser <- 1  # One additional parameter (lambda)
p_ser <- 1 - pchisq(lr_ser, df_ser)
cat(sprintf("OLS vs SER: LR = %.4f, df = %d, p-value = %.4f\n", 
            lr_ser, df_ser, p_ser))

cat("\n=== Conclusion ===\n")
best_model <- model_ranked$Model[1]
cat(sprintf("Best performing model: %s\n", best_model))
cat(sprintf("  - Lowest AIC: %.4f\n", model_ranked$AIC[1]))
cat(sprintf("  - Lowest RMSE (original scale): %.4f\n", model_ranked$RMSE[1]))
cat(sprintf("  - Moran's I p-value: %.4f\n", model_ranked$Moran_pvalue[1]))