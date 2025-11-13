# ===================================================================
#SPATIAL ANALYSIS AND SPATIAL PROBIT MODEL
# ===================================================================
library(sf)            # spatial data handling
library(spdep)         # spatial weights, tests
library(spatialreg)    # errorsarlm, lagsarlm
# spatialprobit provides spatial probit/EM; optional but recommended
# If not installed, the script falls back to non-spatial glm + spaMM suggestion
suppressWarnings(suppressMessages(library(spatialprobit)))
library(ggplot2)
library(tmap)
library(dplyr)
library(texreg)
library(broom)




##################################################################
# SHAPEFILE CREATION AND VISUALIZATION
mydata <- read.csv("reg.csv")
str(mydata)

mydata_sf <- st_as_sf(mydata, coords = c("longitude", "latitude"))
str(mydata_sf)


#ASSIGN WGS84 PROJECTION
st_crs(mydata_sf) <- 4326
#CHECK
st_crs(mydata_sf)
mydata_sf
head(mydata_sf)
sf::st_geometry_type(mydata_sf)

# 
# 
#BASEMAP FOR MAP CREATION
library(ggmap)
citation("ggmap")
#REGISTER  STADIA MAP FOR API
register_stadiamaps(key = "aa4a6990-826c-4e2c-8525-862b0318eb09")

mydata_sf %>% st_bbox() %>% as.vector() %>%
 ggmap::get_stadiamap(zoom = 11, messaging = FALSE) -> baseMap;

p <- ggmap(baseMap) +
  geom_point(aes(x=X, y=Y),
             data = mydata_sf %>% st_coordinates() %>% tibble :: as_tibble(),
             color = "brown",
             size = 1,
             alpha = 0.5
             )
ggsave(
  filename = "high_res_map.png",
  plot = p,
  width = 10,        # in inches
  height = 8,        # in inches
  dpi = 400          # typical for print-quality
)


#To convert sf to sp just use, sf::st_as_sf
#Save the sf
st_write(mydata_sf, "cocoa.shp", delete_layer = TRUE)


#############################################################################

#READ NEW DATA
cocoa_sf <- st_read("cocoa.shp", quiet=TRUE)
head(cocoa_sf)
par(mar=c(1,1,1,1))
plot(cocoa_sf)

dev.off()
#Plot of yield per hectare
plot(cocoa_sf$yld_pr_h)
hist(cocoa_sf$yld_pr_h, axes = TRUE, ann = TRUE)

#BETTER VISUALIZE/DISTINGUISH YIELDS WITH LOW VALUES
ggplot(cocoa_sf) +
  geom_sf(aes(color = yld_pr_h)) +
  scale_color_viridis_c(option = "plasma") +
  theme_minimal() +
  labs(
    title = "Cocoa Yield Distribution",
    subtitle = "Yield values across sampled farm locations",
    color = "Yield (kg/ha)"
  )


#RECHECK THE GEOMETRY AGAIN TO CONFIRM
sf::st_geometry_type(cocoa_sf)


#SPATIAL AUTOCORRELATION (GLOBAL AND LISA)
head(cocoa_sf)
str(cocoa_sf)
names(cocoa_sf)


# 
# ggplot(cocoa_sf) +
#   geom_sf(aes(color = distrct)) +
#   scale_color_viridis_d(option = "plasma") +  # 👈 discrete scale
#   theme_minimal() +
#   labs(
#     title = "Cocoa Yield Distribution by District",
#     subtitle = "Each color represents a district",
#     color = "District"
#   )


#Weight based on neighboring relationships
library(magrittr)

# 1. Extract coordinates from your sf object
coords <- st_coordinates(cocoa_sf)

# 2. Find k nearest neighbors (e.g., k = 4)
knn_obj <- knearneigh(coords, k = 7)

# 3. Convert to neighbors list
nb_obj <- knn2nb(knn_obj)

# 4. Convert to weights list
cocoaNbList <- nb2listw(nb_obj, zero.policy = TRUE)

#Check how neighborhood works
summary(nb_obj)
n.comp.nb(nb_obj)   # number of connected components

# Global Moran's I test
options(scipen = 5)
cocoaNbList %>% 
  spdep::moran.test(cocoa_sf$yld_pr_h, ., zero.policy = TRUE)

# Global Moran's I plot
spdep::moran.plot(cocoa_sf$yld_pr_h,
                  cocoaNbList,
                  zero.policy = TRUE,
                  xlab = 'Cocoa yield per hectare',
                  ylab = 'Lagged cocoa yield per hectare (of Neighbors)',
                  pch = 20)

#Local Moran's I plot
lisarslt <- spdep::localmoran(cocoa_sf$yld_pr_h,
                              cocoaNbList,
                              zero.policy = TRUE,
                              na.action = na.omit
)
dim(lisarslt); dim(cocoa_sf);
head(lisarslt)

#Derive the cluster/outlier types (CoType in ArcGIS) for each spatial feature in the data
significanceLevel <- 0.05  # 95% confidence
meanVal <- mean(cocoa_sf$yld_pr_h)

#For Debugging
#lisarslt <- lisarslt %>% dplyr::select(-coType)
#names(lisarslt)

lisarslt <- lisarslt %>%
  tibble::as_tibble() %>%
  magrittr::set_colnames(c("Ii", "E.Ii", "Var.Ii", "Z.Ii", "Pr")) %>%   # ✅ renamed here
  dplyr::mutate(
    coType = dplyr::case_when(
      Pr > significanceLevel ~ "Insignificant",
      Pr <= significanceLevel & Ii >= 0 & cocoa_sf$yld_pr_h >= meanVal ~ "HH",
      Pr <= significanceLevel & Ii >= 0 & cocoa_sf$yld_pr_h <  meanVal ~ "LL",
      Pr <= significanceLevel & Ii <  0 & cocoa_sf$yld_pr_h >= meanVal ~ "HL",
      Pr <= significanceLevel & Ii <  0 & cocoa_sf$yld_pr_h <  meanVal ~ "LH",
      TRUE ~ "Undefined"
    )
  )

# Attach classification
cocoa_sf$lis_cluster <- lisarslt$coType

# Plot map
ggplot(cocoa_sf) +
  geom_sf(aes(color = lis_cluster)) +
  scale_color_manual(
    values = c(
      "HH" = "red",
      "LL" = "blue",
      "HL" = "orange",
      "LH" = "green",
      "Insignificant" = "grey80"
    )
  ) +
  theme_minimal() +
  labs(
    title = "Local Moran's I Cluster Map",
    subtitle = "Cocoa Yield Spatial Clusters",
    color = "Cluster Type",
    dpi = 400
  )




#LOGISTIC REGRESSION TO TEST PRODUCTIVITY
cocoa_sf <- cocoa_sf %>%
  mutate(yield_bin = ifelse(yld_pr_h > 0, 1L, 0L))

#Number of zeros(no yields), number of ones(yields)
table(cocoa_sf$yield_bin)

pa = yield_bin ~ farm_sz + srfc_mp + nb_f_2122 + rcvd_rm +  c_21_22 + d_frmr_ + 
  fdp_prv + map_frm + gender + hshld_s + nb_frms + yld_pr_t + pds_pr_ + tr_dnst



#Convert the weight list to matrix as part of spatial probit requirement
cocoaNbList_mat <- nb2mat(nb_obj, style = "W", zero.policy = TRUE)
cocoaNbList_mat <- Matrix(cocoaNbList_mat, sparse = TRUE)
diag(cocoaNbList_mat) <- 0


#Model1: SAR Spatial Probit model using MCMC
sar_pm <- spatialprobit::sarprobit(
  formula = pa,
  data    = cocoa_sf,
  W       = cocoaNbList_mat,
  ndraw   = 5000,        # start with 1000 draws; increase to 5000 for final model
  showProgress = TRUE
)

summary(sar_pm)











#PART B (NOW MODELLING YIELD WITH FARMS WITH POSITIVE YIELDS)
#Extract only yields > 0
# Define predictor names (vector of variable names)
predictors <- c("farm_sz", "srfc_mp", "nb_f_2122", "rcvd_rm", 
                "c_21_22", "d_frmr_", "fdp_prv", "map_frm", 
                "gender", "hshld_s", "nb_frms", "yld_pr_t", 
                "pds_pr_", "tr_dnst")

# Identify rows with positive yield and complete predictor data
head(cocoa_sf)
# Drop geometry temporarily
cocoa_df <- st_drop_geometry(cocoa_sf)
head(cocoa_df)

# Now safely subset
pos_yield <- which(cocoa_df$yld_pr_h > 0 & 
                     complete.cases(cocoa_df[, predictors]))

# Subset the original spatial data
cocoa_pos <- cocoa_sf[pos_yield, ]
head(cocoa_pos)

# Check
nrow(cocoa_pos)

#SPATIAL REGRESSION
names(cocoa_pos)

#Model1
library(car)

ols1 = lm(log(yld_pr_h) ~ farm_sz + srfc_mp + nb_f_2122 + rcvd_rm +  c_21_22 + d_frmr_ + 
            fdp_prv + map_frm + gender + hshld_s + nb_frms + log(yld_pr_t) + pds_pr_ + log(tr_dnst + 1),
          data = cocoa_pos)

summary(ols1)

#CHECK RESIDUALS
plot(ols1)
crPlots(ols1)
res <- residuals(ols1)
hist(res)

library(lmtest)
bptest(ols1)
#BP TEST IS SIGNFICANT, FURTHER NECESSITATING THE NEED FOR SPATIAL ERROR MODEL OR MORAN'S I RESIDUAL CHECK ON THE TRANSFORMED OLS


dev.off()

# Cook’s Distance plot
plot(ols1, which = 4, cook.levels = c(4/length(ols1$fitted.values)),
     main = "Cook's Distance (Influence of Observations)")


influencePlot(ols1, 
              id.method = "identify", 
              main = "Influence Plot: Studentized Residuals vs Leverage",
              sub = "Circle size ∝ Cook’s Distance")

plot(hatvalues(ols1), cooks.distance(ols1),
     xlab = "Leverage", ylab = "Cook's Distance",
     main = "Leverage vs. Cook's Distance")
abline(h = 4/length(ols1$fitted.values), col = "red", lty = 2)
abline(v = 2*mean(hatvalues(ols1)), col = "blue", lty = 2)


# Define Cook's distance threshold
n <- nrow(cocoa_pos)
cutoff <- 4 / n

# Find influential points
influential_points <- which(cooks.distance(ols1) > cutoff)

# Convert to regular data frame
infl_df <- as.data.frame(influence.measures(ols1)$infmat)
colnames(infl_df)


# View influence stats for those points
infl_df[influential_points, c("cook.d", "hat", "dffit", "cov.r")]
cocoa_pos[influential_points, ]




#Moran's residual test for spatial autocorrelation
#Redfine the neighborhood list
coords <- st_coordinates(cocoa_pos)

# 2. Find k nearest neighbors (e.g., k = 4)
knn_obj <- knearneigh(coords, k = 7)

# 3. Convert to neighbors list
nb_obj <- knn2nb(knn_obj)

# 4. Convert to weights list
yieldlist <- nb2listw(nb_obj, zero.policy = TRUE)


# --- Moran's I Test on OLS Residuals ---
ols_residuals <- residuals(ols1)
cat("\n=== Moran's I Test on OLS Residuals ===\n")
moran_test_ols <- moran.test(ols_residuals, yieldlist)
print(moran_test_ols)

# --- Rao's Score Tests for Spatial Dependence (recommended replacement for lm.LMtests) ---
cat("\n=== Rao's Score Tests (lm.RStests) ===\n")
rs_tests <- lm.RStests(ols1, yieldlist)
print(rs_tests)


#running SER model since moran's residual test is significant
serr <- spatialreg::errorsarlm(log(yld_pr_h) ~ farm_sz + srfc_mp + nb_f_2122 + rcvd_rm +  c_21_22 + d_frmr_ + 
                                 fdp_prv + map_frm + gender + hshld_s + nb_frms + log(yld_pr_t) + pds_pr_ + log(tr_dnst + 1),
                               data = cocoa_pos,
                               listw = yieldlist,
                               zero.policy = TRUE,
                               na.action = na.omit);


summary(serr)


#DIAGNOSTICS
serres <- residuals(serr)
hist(res)
shapiro_test <- shapiro.test(serres)
print("--- Shapiro-Wilk Normality Test ---")
print(shapiro_test)
# Simpler version: test against fitted values only
resids <- residuals(serr)
fitted_vals <- fitted(serr)

# Breusch-Pagan test
bp_simple <- lm(resids^2 ~ fitted_vals)
bp_result <- lmtest::bptest(bp_simple)
print(bp_result)

# White test (more general heteroskedasticity test)
white_model <- lm(resids^2 ~ fitted_vals + I(fitted_vals^2))
n <- length(resids)
R2_white <- summary(white_model)$r.squared
white_stat <- n * R2_white
white_p <- 1 - pchisq(white_stat, df = 2)

cat("White Test:\n")
cat("Chi-squared =", white_stat, "\n")
cat("P-value =", white_p, "\n")
if(white_p > 0.05) cat("✓ No heteroskedasticity (White test)\n")

# 5. BACK-TRANSFORMATION WITH SMEARING ESTIMATOR
# ===================================================================
# Duan's smearing factor
smearing_factor <- mean(exp(serr_residuals))

# Predicted yields on original scale
predicted_original <- exp(serr_fitted) * smearing_factor
actual_original <- cocoa_pos$yld_pr_h

# Model performance on original scale
R2_original <- 1 - sum((actual_original - predicted_original)^2) / 
  sum((actual_original - mean(actual_original))^2)
RMSE_original <- sqrt(mean((actual_original - predicted_original)^2))
MAE_original <- mean(abs(actual_original - predicted_original))

cat("\n===== MODEL FIT (ORIGINAL SCALE WITH SMEARING) =====\n")
cat("R² on original scale:", round(R2_original, 4), "\n")
cat("RMSE (kg/ha):", round(RMSE_original, 2), "\n")
cat("MAE (kg/ha):", round(MAE_original, 2), "\n")
cat("Mean observed yield:", round(mean(actual_original), 2), "kg/ha\n")
cat("Smearing factor:", round(smearing_factor, 4), "\n")
cat("Prediction error (% of mean):", 
    round(100 * RMSE_original / mean(actual_original), 2), "%\n\n")

# ===================================================================
# SENSITIVITY ANALYSIS FOR SER MODEL
# ===================================================================

# 6. IDENTIFY MOST INFLUENTIAL OBSERVATION
# ===================================================================

# Find most influential observation (largest absolute standardized residual)
most_influential <- influential_points[which.max(abs(std_residuals[influential_points]))]

cat("===== SENSITIVITY ANALYSIS =====\n")
cat("Most influential observation:", most_influential, "\n")
cat("Standardized residual:", round(std_residuals[most_influential], 3), "\n")
cat("Actual yield:", round(cocoa_pos$yld_pr_h[most_influential], 2), "kg/ha\n")
cat("Predicted yield (original scale):", 
    round(predicted_original[most_influential], 2), "kg/ha\n")
cat("Residual (original scale):", 
    round(actual_original[most_influential] - predicted_original[most_influential], 2), 
    "kg/ha\n\n")

#Test if the residuals have been taken care of
cat("\n=== Moran’s I on SEM Residuals ===\n")
moran_sem <- moran.test(residuals(serr), yieldlist)
print(moran_sem)

#SAR Model
sarr <- spatialreg::lagsarlm(log(yld_pr_h) ~ farm_sz + srfc_mp + nb_f_2122 + rcvd_rm +  c_21_22 + d_frmr_ + 
                               fdp_prv + map_frm + gender + hshld_s + nb_frms + log(yld_pr_t) + pds_pr_ + log(tr_dnst + 1),
                             data = cocoa_pos,
                             listw = yieldlist,
                             zero.policy = TRUE,
                             na.action = na.omit);
summary(sarr)











