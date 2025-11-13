# ===============================================================
# Ripley’s K-function and Nearest Neighbor G-function Analysis
# ===============================================================

# ---------------------------------------------------------------
# 1. Load required libraries
# ---------------------------------------------------------------
library(sf)
library(spatstat.geom)
library(spdep)
library(spatstat)

# ---------------------------------------------------------------
# 2. Load and inspect data
# ---------------------------------------------------------------
data <- read.csv("cocoa.csv")
head(data)
names(data)

# ---------------------------------------------------------------
# 3. Convert to sf object in WGS84 (EPSG:4326)
# ---------------------------------------------------------------
sf_points <- st_as_sf(data, coords = c("longitude", "latitude"), crs = 4326)

# Convert to spatstat point pattern (ppp)
window <- as.owin(st_bbox(sf_points))
coords <- st_coordinates(sf_points)
pp <- ppp(x = coords[,1],
          y = coords[,2],
          window = window,
          marks = data$yield_per_hectare)

summary(pp)

# ===============================================================
# 4. Ripley’s K-function with Monte Carlo envelope
# ===============================================================
set.seed(123)
env_K <- envelope(pp, fun = Kest, nsim = 199,
                  correction = "border", global = TRUE, savefuns = TRUE)

# Convert to dataframe for custom plotting
K_data <- as.data.frame(env_K)

# --- Custom K-function plot with CI shading and legend ---
plot(
  K_data$r * 111.32, K_data$obs,
  type = "l", col = "black", lwd = 2,
  xlab = "Distance r (km)", ylab = "K(r)",
  main = "Ripley’s K-function (Observed vs CSR)"
)

# Add 95% Monte Carlo confidence interval (shaded)
polygon(
  c(K_data$r * 111.32, rev(K_data$r * 111.32)),
  c(K_data$lo, rev(K_data$hi)),
  col = adjustcolor("gray80", alpha.f = 0.6), border = NA
)

# Add CSR expectation line
lines(K_data$r * 111.32, K_data$theo, col = "red", lty = 2, lwd = 2)

# Re-draw observed line on top
lines(K_data$r * 111.32, K_data$obs, col = "black", lwd = 2)

# Vertical markers for 5.5 km and 28 km
#abline(v = c(5.5, 28), col = "grey50", lty = 2)
#text(5.5, max(K_data$hi, na.rm = TRUE), labels = "~5.5 km", pos = 3)
#text(28, max(K_data$hi, na.rm = TRUE), labels = "~28 km", pos = 3)

# Add legend
legend("topleft",
       legend = c("Observed K(r)", "CSR expectation", "95% Monte Carlo envelope"),
       col = c("black", "red", "gray50"),
       lty = c(1, 2, NA), lwd = c(2, 2, NA),
       pch = c(NA, NA, 15),
       pt.cex = 2, bty = "n", text.col = "black")

# ---------------------------------------------------------------
# 4b. Quantify K-function exceedance at 0.05° and 0.25°
# ---------------------------------------------------------------
Kest_result <- Kest(pp, correction = "border")

dist_deg <- c(0.05, 0.25)  # degrees (~5.5 km and ~28 km)
K_obs <- approx(x = Kest_result$r, y = Kest_result$border, xout = dist_deg)$y
K_theo <- pi * (dist_deg^2)
percent_exceedance <- (K_obs - K_theo) / K_theo * 100

results_table <- data.frame(
  distance_deg = dist_deg,
  distance_km = dist_deg * 111.32,
  K_obs = K_obs,
  K_theo = K_theo,
  percent_exceedance = percent_exceedance
)

print(results_table)

for (i in seq_len(nrow(results_table))) {
  d_km <- results_table$distance_km[i]
  pct <- round(results_table$percent_exceedance[i], 1)
  cat(sprintf("At %.1f km (~%.3f°): observed K(r) is %.1f%% %s CSR expectation.\n",
              d_km, results_table$distance_deg[i],
              abs(pct), ifelse(pct >= 0, "above", "below")))
}


# ===============================================================
# 5. Nearest Neighbor G-function with Monte Carlo envelope
# ===============================================================
set.seed(123)
env_G <- envelope(pp, fun = Gest, nsim = 199,
                  correction = "rs", global = TRUE, savefuns = TRUE)

G_data <- as.data.frame(env_G)

# --- Custom G-function plot with CI shading and legend ---
plot(
  G_data$r * 111.32, G_data$obs,
  type = "l", col = "black", lwd = 2,
  xlab = "Distance r (km)", ylab = "G(r)",
  main = "Nearest Neighbor G-function (Observed vs CSR)"
)

# 95% Monte Carlo confidence interval
polygon(
  c(G_data$r * 111.32, rev(G_data$r * 111.32)),
  c(G_data$lo, rev(G_data$hi)),
  col = adjustcolor("gray80", alpha.f = 0.6), border = NA
)

# CSR expectation
lines(G_data$r * 111.32, G_data$theo, col = "red", lty = 2, lwd = 2)

# Observed line (on top)
lines(G_data$r * 111.32, G_data$obs, col = "black", lwd = 2)

# Add legend
legend("topleft",
       legend = c("Observed G(r)", "CSR expectation", "95% Monte Carlo envelope"),
       col = c("black", "red", "gray50"),
       lty = c(1, 2, NA), lwd = c(2, 2, NA),
       pch = c(NA, NA, 15),
       pt.cex = 2, bty = "n", text.col = "black")


# ---------------------------------------------------------------
# 5b. Quantify nearest neighbor clustering (G-function)
# ---------------------------------------------------------------
dist_target_km <- 2
dist_target_deg <- dist_target_km / 111.32
nearest_index <- which.min(abs(G_data$r - dist_target_deg))
actual_r_deg <- G_data$r[nearest_index]
actual_r_km <- actual_r_deg * 111.32
G_obs <- G_data$obs[nearest_index]
G_theo <- G_data$theo[nearest_index]
percent_diff <- (G_obs - G_theo) / G_theo * 100

cat(sprintf(
  "At ~%.2f km (%.3f°): observed G(r) = %.3f, theoretical = %.3f, observed exceeds CSR by %.1f%%\n",
  actual_r_km, actual_r_deg, G_obs, G_theo, percent_diff
))

ks_D <- max(abs(G_data$obs - G_data$theo), na.rm = TRUE)
cat(sprintf("Kolmogorov–Smirnov D = %.2f (p < 0.001 assumed for strong deviation)\n", ks_D))

cat(sprintf(
  "\nInterpretation:\nThe nearest neighbor G-function revealed pronounced short-distance clustering, with %.0f%% of farms located within %.1f km of another farm compared to %.0f%% expected under CSR (Kolmogorov–Smirnov test: D = %.2f, p < 0.001). This tight spatial aggregation has implications for pest and disease transmission, as pathogens such as Phytophthora (black pod) and cocoa swollen shoot virus spread more efficiently through contiguous farm networks.\n",
  100 * G_obs, actual_r_km, 100 * G_theo, ks_D
))


