# ===================================================================
# INFLUENTIAL POINTS ANALYSIS FOR SPATIAL ERROR MODEL (SER)
# ===================================================================

library(spatialreg)
library(car)
library(dplyr)
library(sf)

# Note: Spatial models don't have Cook's distance equivalent
# We use standardized residuals as the primary diagnostic

# 1. EXTRACT RESIDUALS AND FITTED VALUES (LOG SCALE)
# ===================================================================
serr_residuals <- residuals(serr)
serr_fitted <- fitted(serr)
n <- length(serr_residuals)

# 2. CALCULATE DIAGNOSTIC MEASURES
# ===================================================================
# For spatial error models, standardized residuals are the 
# primary diagnostic tool (Anselin 1988, LeSage & Pace 2009)

# Standardized residuals
std_residuals <- serr_residuals / sd(serr_residuals)

# Extract spatial parameter
lambda <- serr$lambda

# 3. IDENTIFY INFLUENTIAL OBSERVATIONS
# ===================================================================
# Use multiple criteria based on standardized residuals:

# Criterion 1: Large standardized residuals (>3 SD)
cutoff_resid <- 3
large_residuals <- which(abs(std_residuals) > cutoff_resid)

# Criterion 2: Top 5% by absolute standardized residual
cutoff_pct <- quantile(abs(std_residuals), 0.95)
top_residuals <- which(abs(std_residuals) > cutoff_pct)

# Combine all criteria
influential_points <- unique(c(large_residuals, top_residuals))

cat("\n===== INFLUENTIAL OBSERVATIONS IDENTIFIED =====\n")
cat("By residual > 3 SD:", length(large_residuals), "\n")
cat("By top 5% residuals:", length(top_residuals), "\n")
cat("Total unique influential points:", length(influential_points), 
    "(", round(100*length(influential_points)/n, 1), "%)\n\n")

# 4. CREATE DIAGNOSTIC TABLE
# ===================================================================
diagnostic_table <- data.frame(
  Obs_ID = influential_points,
  Std_Resid = round(std_residuals[influential_points], 3),
  Log_Yield = round(log(cocoa_pos$yld_pr_h[influential_points]), 3),
  Log_Fitted = round(serr_fitted[influential_points], 3),
  Residual = round(serr_residuals[influential_points], 3),
  Farm_Size = round(cocoa_pos$farm_sz[influential_points], 2),
  Tree_Density = round(cocoa_pos$tr_dnst[influential_points], 2),
  Yield_Actual = round(cocoa_pos$yld_pr_h[influential_points], 2)
) %>% arrange(desc(abs(Std_Resid)))

cat("===== TOP 15 INFLUENTIAL OBSERVATIONS (SER MODEL) =====\n")
print(head(diagnostic_table, 15))



# 7. CREATE REDUCED DATASET AND REFIT MODEL
# ===================================================================

# Remove most influential observation
cocoa_reduced <- cocoa_pos[-most_influential, ]

# Extract coordinates and recreate spatial weights for reduced dataset
coords_full <- st_coordinates(cocoa_pos)
coords_reduced <- coords_full[-most_influential, ]
coords_reduced <- as.matrix(coords_reduced)

# Recreate k-nearest neighbor list for reduced data (using same k=7)
yieldlist_reduced <- knn2nb(knearneigh(coords_reduced, k = 7))
yieldlist_reduced <- nb2listw(yieldlist_reduced, style = "W", zero.policy = TRUE)

# Refit SER model on reduced dataset
cat("Refitting model without observation", most_influential, "...\n")
serr_reduced <- errorsarlm(
  log(yld_pr_h) ~ farm_sz + srfc_mp + nb_f_2122 + rcvd_rm + c_21_22 + 
    d_frmr_ + fdp_prv + map_frm + gender + hshld_s + nb_frms + 
    log(yld_pr_t) + pds_pr_ + log(tr_dnst + 1),
  data = cocoa_reduced,
  listw = yieldlist_reduced,
  zero.policy = TRUE,
  na.action = na.omit
)

summary(serr_reduced)

# 8. COMPARE COEFFICIENTS
# ===================================================================

cat("\n===== COEFFICIENT COMPARISON =====\n")
coef_comparison <- data.frame(
  Variable = names(coef(serr)),
  Full_Model = round(coef(serr), 4),
  Without_Top = round(coef(serr_reduced), 4),
  Abs_Change = round(coef(serr_reduced) - coef(serr), 4),
  Pct_Change = round(((coef(serr_reduced) - coef(serr)) / abs(coef(serr))) * 100, 1)
)

# Sort by absolute percentage change
coef_comparison <- coef_comparison %>% 
  arrange(desc(abs(Pct_Change)))

print(coef_comparison)

# 9. COMPARE SPATIAL PARAMETERS
# ===================================================================

cat("\n===== SPATIAL PARAMETER COMPARISON =====\n")
cat("Lambda (full model):   ", round(serr$lambda, 4), "\n")
cat("Lambda (reduced):      ", round(serr_reduced$lambda, 4), "\n")
cat("Absolute change:       ", round(serr_reduced$lambda - serr$lambda, 4), "\n")
cat("Percentage change:     ", 
    round(100*(serr_reduced$lambda - serr$lambda)/serr$lambda, 2), "%\n\n")

# 10. COMPARE MODEL FIT
# ===================================================================

cat("===== MODEL FIT COMPARISON =====\n")
fit_comparison <- data.frame(
  Metric = c("Log-likelihood", "AIC", "Lambda", "N"),
  Full_Model = c(
    round(serr$LL, 2),
    round(AIC(serr), 2),
    round(serr$lambda, 4),
    nrow(cocoa_pos)
  ),
  Without_Top = c(
    round(serr_reduced$LL, 2),
    round(AIC(serr_reduced), 2),
    round(serr_reduced$lambda, 4),
    nrow(cocoa_reduced)
  ),
  Change = c(
    round(serr_reduced$LL - serr$LL, 2),
    round(AIC(serr_reduced) - AIC(serr), 2),
    round(serr_reduced$lambda - serr$lambda, 4),
    -1
  )
)
print(fit_comparison)

# 11. COMPARE R² ON ORIGINAL SCALE
# ===================================================================

# Back-transform predictions with smearing estimator for reduced model
serr_reduced_residuals <- residuals(serr_reduced)
serr_reduced_fitted <- fitted(serr_reduced)
smearing_factor_reduced <- mean(exp(serr_reduced_residuals))
predicted_reduced_original <- exp(serr_reduced_fitted) * smearing_factor_reduced

# Calculate R² on original scale
R2_reduced <- 1 - sum((cocoa_reduced$yld_pr_h - predicted_reduced_original)^2) / 
  sum((cocoa_reduced$yld_pr_h - mean(cocoa_reduced$yld_pr_h))^2)

RMSE_reduced <- sqrt(mean((cocoa_reduced$yld_pr_h - predicted_reduced_original)^2))
MAE_reduced <- mean(abs(cocoa_reduced$yld_pr_h - predicted_reduced_original))

cat("\n===== PREDICTION PERFORMANCE COMPARISON (ORIGINAL SCALE) =====\n")
performance_comparison <- data.frame(
  Metric = c("R²", "RMSE (kg/ha)", "MAE (kg/ha)", "Smearing Factor"),
  Full_Model = c(
    round(R2_original, 4),
    round(RMSE_original, 2),
    round(MAE_original, 2),
    round(smearing_factor, 4)
  ),
  Without_Top = c(
    round(R2_reduced, 4),
    round(RMSE_reduced, 2),
    round(MAE_reduced, 2),
    round(smearing_factor_reduced, 4)
  ),
  Change = c(
    round(R2_reduced - R2_original, 4),
    round(RMSE_reduced - RMSE_original, 2),
    round(MAE_reduced - MAE_original, 2),
    round(smearing_factor_reduced - smearing_factor, 4)
  )
)
print(performance_comparison)

# 12. ROBUSTNESS ASSESSMENT
# ===================================================================

# Identify maximum coefficient change
max_change_idx <- which.max(abs(coef_comparison$Pct_Change))
max_change_var <- coef_comparison$Variable[max_change_idx]
max_change_pct <- abs(coef_comparison$Pct_Change[max_change_idx])

# Identify variables with >10% change
large_changes <- coef_comparison %>% 
  filter(abs(Pct_Change) > 10)

cat("\n===== ROBUSTNESS ASSESSMENT =====\n")
cat("Most influential observation:", most_influential, "\n")
cat("Variable with maximum change:", max_change_var, "\n")
cat("Maximum coefficient change:", round(max_change_pct, 1), "%\n")
cat("Lambda change:", 
    round(100*abs(serr_reduced$lambda - serr$lambda)/serr$lambda, 1), "%\n")
cat("R² change:", round(abs(R2_reduced - R2_original), 4), "\n\n")

if(nrow(large_changes) > 0) {
  cat("Variables with >10% coefficient change:\n")
  print(large_changes[, c("Variable", "Full_Model", "Without_Top", "Pct_Change")])
  cat("\n")
}

# Overall robustness conclusion
if(max_change_pct < 5 & abs(serr_reduced$lambda - serr$lambda)/serr$lambda < 0.05) {
  cat("✓ CONCLUSION: Model is HIGHLY ROBUST to influential observations\n")
  cat("  - All coefficients change <5%\n")
  cat("  - Lambda change <5%\n")
  cat("  - Results are stable and reliable\n")
} else if(max_change_pct < 10 & abs(serr_reduced$lambda - serr$lambda)/serr$lambda < 0.10) {
  cat("✓ CONCLUSION: Model is ROBUST (changes <10%)\n")
  cat("  - Key findings remain stable\n")
  cat("  - Minor sensitivity to influential observations\n")
} else if(max_change_pct < 20) {
  cat("⚠ CONCLUSION: Model shows MODERATE sensitivity\n")
  cat("  - Recommend: Report sensitivity analysis in paper\n")
  cat("  - Consider investigating observation", most_influential, "\n")
} else {
  cat("⚠ CONCLUSION: Model shows HIGH sensitivity (change >20%)\n")
  cat("  - ACTION REQUIRED: Investigate observation", most_influential, "\n")
  cat("  - Check for data quality issues or outliers\n")
  cat("  - Consider exclusion or discuss in limitations section\n")
}

# 13. EXPORT RESULTS
# ===================================================================

write.csv(diagnostic_table, "ser_influential_observations.csv", row.names = FALSE)
write.csv(coef_comparison, "ser_sensitivity_coefficients.csv", row.names = FALSE)
write.csv(fit_comparison, "ser_model_fit_comparison.csv", row.names = FALSE)
write.csv(performance_comparison, "ser_performance_comparison.csv", row.names = FALSE)

cat("\n✓ Files exported:\n")
cat("  - ser_influential_observations.csv\n")
cat("  - ser_sensitivity_coefficients.csv\n")
cat("  - ser_model_fit_comparison.csv\n")
cat("  - ser_performance_comparison.csv\n")

cat("\n===== ANALYSIS COMPLETE =====\n")



# ===================================================================
# REFIT SPATIAL ERROR MODEL EXCLUDING OBSERVATION 95
# ===================================================================

library(spatialreg)
library(spdep)
library(sf)
library(dplyr)

# 1. CREATE FINAL DATASET WITHOUT OBSERVATION 95
# ===================================================================
cat("===== FINAL MODEL: EXCLUDING OBSERVATION 95 =====\n\n")

# Remove observation 95
cocoa_final <- cocoa_pos[-95, ]

# Extract coordinates and recreate spatial weights
coords_final <- st_coordinates(cocoa_final)
yieldlist_final <- knn2nb(knearneigh(coords_final, k = 7))
yieldlist_final <- nb2listw(yieldlist_final, style = "W", zero.policy = TRUE)

cat("Sample size (final model):", nrow(cocoa_final), "\n")
cat("Observations excluded: 1 (Obs 95)\n\n")

# 2. REFIT SPATIAL ERROR MODEL
# ===================================================================
serr_final <- errorsarlm(
  log(yld_pr_h) ~ farm_sz + srfc_mp + nb_f_2122 + rcvd_rm + c_21_22 + 
    d_frmr_ + fdp_prv + map_frm + gender + hshld_s + nb_frms + 
    log(yld_pr_t) + pds_pr_ + log(tr_dnst + 1),
  data = cocoa_final,
  listw = yieldlist_final,
  zero.policy = TRUE,
  na.action = na.omit
)

# 3. MODEL SUMMARY
# ===================================================================
cat("===== FINAL MODEL SUMMARY =====\n")
summary(serr_final)

# 4. BACK-TRANSFORMATION WITH SMEARING
# ===================================================================
serr_final_residuals <- residuals(serr_final)
serr_final_fitted <- fitted(serr_final)
smearing_factor_final <- mean(exp(serr_final_residuals))

predicted_final_original <- exp(serr_final_fitted) * smearing_factor_final
actual_final_original <- cocoa_final$yld_pr_h

# Calculate performance metrics on original scale
R2_final <- 1 - sum((actual_final_original - predicted_final_original)^2) / 
  sum((actual_final_original - mean(actual_final_original))^2)
RMSE_final <- sqrt(mean((actual_final_original - predicted_final_original)^2))
MAE_final <- mean(abs(actual_final_original - predicted_final_original))

cat("\n===== FINAL MODEL PERFORMANCE (ORIGINAL SCALE) =====\n")
cat("R²:", round(R2_final, 4), "\n")
cat("RMSE (kg/ha):", round(RMSE_final, 2), "\n")
cat("MAE (kg/ha):", round(MAE_final, 2), "\n")
cat("Smearing factor:", round(smearing_factor_final, 4), "\n")
cat("Lambda (spatial error):", round(serr_final$lambda, 4), "\n")
cat("Log-likelihood:", round(serr_final$LL, 2), "\n")
cat("AIC:", round(AIC(serr_final), 2), "\n\n")

# 5. CHECK RESIDUAL SPATIAL AUTOCORRELATION
# ===================================================================
moran_final <- moran.test(serr_final_residuals, yieldlist_final, zero.policy = TRUE)
cat("===== RESIDUAL SPATIAL AUTOCORRELATION TEST =====\n")
cat("Moran's I:", round(moran_final$estimate[1], 4), "\n")
cat("p-value:", round(moran_final$p.value, 4), "\n")
if(moran_final$p.value > 0.05) {
  cat("✓ No significant spatial autocorrelation in residuals\n\n")
} else {
  cat("⚠ Residuals still show spatial autocorrelation\n\n")
}

# 6. COEFFICIENT TABLE FOR PAPER
# ===================================================================
coef_final <- summary(serr_final)$Coef
coef_table <- data.frame(
  Variable = rownames(coef_final),
  Coefficient = round(coef_final[, 1], 4),
  Std_Error = round(coef_final[, 2], 4),
  z_value = round(coef_final[, 3], 4),
  p_value = round(coef_final[, 4], 4)
)

cat("===== COEFFICIENT TABLE (FINAL MODEL) =====\n")
print(coef_table)

# 7. KEY FINDINGS SUMMARY
# ===================================================================
tree_density_coef <- coef(serr_final)["log(tr_dnst + 1)"]
farm_size_coef <- coef(serr_final)["farm_sz"]
yield_per_tree_coef <- coef(serr_final)["log(yld_pr_t)"]

cat("\n===== KEY PRODUCTIVITY DETERMINANTS =====\n")
cat("Tree density elasticity:", round(tree_density_coef, 4), "\n")
cat("  → 10% increase in tree density → ", 
    round(tree_density_coef * 10, 2), "% yield increase\n")
cat("Farm size elasticity:", round(farm_size_coef, 4), "\n")
cat("  → 1 ha increase → ", round(farm_size_coef * 100, 2), "% yield change\n")
cat("Yield per tree elasticity:", round(yield_per_tree_coef, 4), "\n\n")

# 8. DIAGNOSTIC PLOTS
# ===================================================================

# Plot 1: Residuals vs Fitted

plot(serr_final_fitted, serr_final_residuals,
     xlab = "Fitted Values (log scale)",
     ylab = "Residuals",
     main = "Final SER Model: Residuals vs Fitted Values",
     pch = 19, col = rgb(0, 0, 0, 0.5))
abline(h = 0, col = "red", lty = 2, lwd = 2)
abline(h = c(-2, 2) * sd(serr_final_residuals), col = "blue", lty = 2)
dev.off()

# Plot 2: Actual vs Predicted (Original Scale)
plot(actual_final_original, predicted_final_original,
     xlab = "Actual Yield (kg/ha)",
     ylab = "Predicted Yield (kg/ha)",
     main = "Final SER Model: Actual vs Predicted",
     pch = 19, col = rgb(0, 0, 0, 0.5),
     xlim = c(0, max(actual_final_original)),
     ylim = c(0, max(actual_final_original)))
abline(0, 1, col = "red", lty = 2, lwd = 2)
text(max(actual_final_original) * 0.1, max(actual_final_original) * 0.9,
     labels = paste0("R² = ", round(R2_final, 3)), cex = 1.2)
dev.off()

# 9. EXPORT FINAL MODEL RESULTS
# ===================================================================
write.csv(coef_table, "ser_final_coefficients_for_removed_obs95.csv", row.names = FALSE)

# Save final model object
saveRDS(serr_final, "serr_final_model.rds")

# Save predictions
predictions_final <- data.frame(
  Actual = actual_final_original,
  Predicted = predicted_final_original,
  Residual = actual_final_original - predicted_final_original,
  Fitted_Log = serr_final_fitted,
  Residual_Log = serr_final_residuals
)
write.csv(predictions_final, "ser_final_predictions.csv", row.names = FALSE)

cat("✓ Final model results exported:\n")
cat("  - ser_final_coefficients.csv\n")
cat("  - ser_final_predictions.csv\n")
cat("  - serr_final_model.rds\n")
cat("  - Diagnostic plots in figures/\n")
