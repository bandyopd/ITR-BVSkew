#####################################################
### R Code from Fan Y and Bandyopadhyay D. (2026+). Individualized treatment rules for bivariate outcomes with 
### informative follow-up: Application to periodontal health
### Authors: Yiwei Fan and Dipankar Bandyopadhyay 
### Last edited: 02/23/2026
### "analysis.R" calls the list of functions in "functions_submit.R", and the application dataset (example_data.csv)
#######################################################


rm(list = ls())

setwd("G:\\My Drive\\OneDrive\\Work\\YiweiFan_Renmin\\ITRProject_Code\\Latest\\Code\\")

# Vector of required packages
required_packages <- c(
  "haven", "dplyr", "ggplot2", "gridExtra", "lmerTest", "tidyr", "RColorBrewer", "scales"
)



# Identify packages not yet installed
new_packages <- required_packages[!(
  required_packages %in% rownames(installed.packages())
)]

# Install missing packages
if (length(new_packages) > 0) {
  install.packages(new_packages, dependencies = TRUE)
}

# Load all required packages
invisible(lapply(required_packages, library, character.only = TRUE))



###############################################
### The "example_data.csv" dataset contains a random sample of 2464 subjects from the full HealthPartners dataset
### Since the full dataset can't be shared due to ethical constraints, a random sample was generated for illustration
##############################################

subdf=read.csv("example_data.csv")
dim(subdf)

source("myfunctions.R")
## Usage
## Run
################################# mixed model
res1 <- fit_stage1_empirical(subdf, y = "COMP_next", id = "STUDY_ID", scale_continuous = FALSE)
summary(res1$fit_lmm)
head(res1$W_hat_df)


library(lmerTest)

fit_p <- lmerTest::as_lmerModLmerTest(res1$fit_lmm)

# Fixed effects table (df and p-values)
fe <- coef(summary(fit_p))  # matrix
terms <- rownames(fe)

# 95% CI (fixed effects only; Wald)
ci <- suppressMessages(confint(fit_p, parm = terms, method = "Wald"))

fixed_effect_table <- data.frame(
  term     = terms,
  estimate = fe[, "Estimate"],
  se       = fe[, "Std. Error"],
  df       = fe[, "df"],
  t        = fe[, "t value"],
  p        = fe[, "Pr(>|t|)"],
  ci_low   = ci[, 1],
  ci_high  = ci[, 2],
  row.names = NULL
)

fixed_effect_table

#write.csv(fixed_effect_table, "mixedmodel_fixed_effects_COMP.csv", row.names = FALSE)

#####################################
######### Building the intensity model
#####################################
# 1) Build intensity data
subdf$STUDY_ID <- as.character(subdf$STUDY_ID)
res1$W_hat_df$STUDY_ID <- as.character(res1$W_hat_df$STUDY_ID)

# Merge back to visits
subdf$W_hat <- res1$W_hat_df$W_hat[ match(subdf$STUDY_ID, res1$W_hat_df$STUDY_ID) ]

# Quick check
head(res1$W_hat_df)
head(subdf[, c("STUDY_ID","W_hat")])

int_dat <- build_intensity_data_from_visits_emp(subdf,
                                                id_col = "STUDY_ID",
                                                time_col = "cumdelta",  # remove if absent; will use cumulative delta
                                                delta_col = "delta",
                                                a_col = "TX_IND1",
                                                w_col = "W_hat")

# 2) Fit Cox intensity model
res_int <- fit_intensity_cox_emp(int_dat, a_col = "TX_IND1", w_col = "W_hat")
summary(res_int$fit_cox)
s <- summary(res_int$fit_cox)

cox_table <- data.frame(
  term    = rownames(s$coefficients),
  coef    = s$coefficients[, "coef"],
  se      = s$coefficients[, "se(coef)"],
  z       = s$coefficients[, "z"],
  p       = s$coefficients[, "Pr(>|z|)"],
  ci_low  = s$coefficients[, "coef"] - qnorm(0.975) * s$coefficients[, "se(coef)"],
  ci_high = s$coefficients[, "coef"] + qnorm(0.975) * s$coefficients[, "se(coef)"],
  row.names = NULL
)


cox_table
#write.csv(cox_table, "visit_intensity_cox_coef_p.csv", row.names = FALSE)

# 3) Add eta_rate / hr back to subdf (one row per visit)
subdf2 <- add_eta_rate_to_visits_emp(subdf, res_int$intensity_data, res_int$fit_cox)
head(subdf2)
subdf2[is.na(subdf2$hr),]
head(int_dat)
dim(res_int$intensity_data)
length(is.na(subdf2$hr))
range(subdf2$hr,na.rm=T)


################################### treatment model
res_trt <- fit_treatment(subdf2)          # 用 subdf2 拟合（或用 subdf 拟合也行，但贴回要贴到 subdf2）
summary(res_trt$fit_trt)

# attach pi_hat to subdf2
subdf2$pi_hat <- predict(res_trt$fit_trt, newdata = subdf2, type = "response")

fit <- res_trt$fit_trt
s <- summary(fit)

# Coefficient table (Estimate/SE/z/p)
tab <- as.data.frame(s$coefficients)
tab$term <- rownames(tab)
rownames(tab) <- NULL
names(tab) <- c("estimate","se","z","p","term")

# beta 95% CI (Wald)
tab$CI_low  <- tab$estimate - qnorm(0.975) * tab$se
tab$CI_high <- tab$estimate + qnorm(0.975) * tab$se

tab <- tab[, c("term","estimate","se","z","p","CI_low","CI_high")]
tab
#write.csv(tab, "treatment_model_coef_p.csv", row.names = FALSE)


## Usage:
## subdf2 must already contain pi_hat (treatment model) and hr (intensity model)
subdf2 <- build_weights_simple(subdf2, trim_q = c(0.05, 0.95), stabilized = FALSE)
range(subdf2$w_all_trim,na.rm=T)




################## joint spline model (robust runnable block)

## -------- 0) setup --------
w_col <- "w_all_trim"

X_vars <- c("GENDER1","RACE1","RACE2","AGE1","DM_IND1",
            "DENT_COMMERCIAL1","TOBACCO_IND1","TOBACCO_IND2","COMP",
            "BASE_CAL","BASE_PPD")

need_cols <- c("STUDY_ID","cumdelta","TX_IND1","CAL","PD",
               "hr","pi_hat", w_col, X_vars)

missing_cols <- setdiff(need_cols, names(subdf2))
if (length(missing_cols) > 0) {
  stop(paste0("subdf2 is missing columns: ", paste(missing_cols, collapse = ", ")))
}

## -------- 1) define variables required by fit_blip_spline() --------
subdf2$A      <- subdf2$TX_IND1
subdf2$Y1_obs <- subdf2$CAL
subdf2$Y2_obs <- subdf2$PD

## -------- 2) sort within subject --------
subdf2 <- subdf2[order(subdf2$STUDY_ID, subdf2$cumdelta), ]

## -------- 3) create visit index (0,1,2,...) within subject --------
subdf2$visit <- ave(seq_len(nrow(subdf2)), subdf2$STUDY_ID, FUN = seq_along) - 1

## -------- 4) keep complete/finite rows (NO visit filtering) --------
idx_keep <- rep(TRUE, nrow(subdf2))

# finite weight
idx_keep <- idx_keep & is.finite(subdf2[[w_col]])

# finite treatment/outcomes/covariates
idx_keep <- idx_keep & complete.cases(subdf2[, c("A","Y1_obs","Y2_obs", X_vars), drop = FALSE])

# also require finite hr/pi_hat (since weight components come from them)
idx_keep <- idx_keep & is.finite(subdf2$hr) & is.finite(subdf2$pi_hat)

subdf2_clean <- subdf2[idx_keep, ]

# stop early if nothing left
if (nrow(subdf2_clean) < 10) stop("Too few rows after cleaning. Check hr/pi_hat/weights and covariate missingness.")

## -------- 5) rescale weight (numerically safer; constant scaling) --------
mw <- mean(subdf2_clean[[w_col]])
if (!is.finite(mw) || mw <= 0) stop("Mean of weights is not finite/positive after cleaning.")
subdf2_clean[[w_col]] <- subdf2_clean[[w_col]] / mw

## quick checks
print(c(
  n_before = nrow(subdf2),
  n_after  = nrow(subdf2_clean),
  drop_n   = nrow(subdf2) - nrow(subdf2_clean)
))
print(summary(subdf2_clean[[w_col]]))
print(dim(subdf2_clean))

## -------- 6) define grids and dat_raw_emp (MUST exist before fitting) --------
tau1_vec <- as.numeric(quantile(subdf2_clean$CAL, probs = c(0.15,0.3,0.5,0.7,0.85), na.rm = TRUE))
tau2_vec <- as.numeric(quantile(subdf2_clean$PD,  probs = c(0.15,0.3,0.5,0.7,0.85), na.rm = TRUE))

# ensure strictly increasing unique grid (avoid accidental duplicates)
tau1_vec <- sort(unique(tau1_vec))
tau2_vec <- sort(unique(tau2_vec))
if (length(tau1_vec) < 3 || length(tau2_vec) < 3) stop("tau grids have too few unique values; choose different probs.")

dat_raw_emp <- list(tau1_vec = tau1_vec, tau2_vec = tau2_vec)

## choose (tau1, tau2) INSIDE the grid to avoid mismatch
tau1 <- tau1_vec[3]   # median grid point
tau2 <- tau2_vec[3]


##############################
myw=matrix(subdf2_clean$w_all_trim, 
           nrow = length(subdf2_clean$w_all_trim), 
           ncol = 25)
dim(myw)
fit_fast <- fit_blip_spline_fast(subdf2_clean, dat_raw_emp,
                                 myw, 
                                 tau1 = tau1, tau2 = tau2,
                                 X_vars = X_vars)

kappa = fit_fast$kappa_grid
dim(kappa)
kappa = kappa/max(kappa)
kappa = pmin(pmax(kappa, quantile(kappa,0.01)), quantile(kappa,0.99))
dim(kappa)
fit_fast <- fit_blip_spline_fast(subdf2_clean, dat_raw_emp,
                                 myw*(1/kappa), 
                                 tau1 = tau1, tau2 = tau2,
                                 X_vars = X_vars)



fit_blip <- fit_fast


library(dplyr)
library(tidyr)
library(ggplot2)
library(RColorBrewer)
library(scales)

## 1) 25×(2+q_x) long table
blip_df <- get_blip_grid(fit_blip)

blip_long <- blip_df %>%
  pivot_longer(cols = -c(tau1, tau2), names_to = "term", values_to = "coef") %>%
  mutate(
    tau1_f = factor(formatC(tau1, format = "f", digits = 1),
                    levels = formatC(sort(unique(tau1)), format = "f", digits = 1)),
    tau2_f = factor(formatC(tau2, format = "f", digits = 1),
                    levels = formatC(sort(unique(tau2)), format = "f", digits = 1)),
    lab = sprintf("%.3f", coef)
  )

pal <- rev(brewer.pal(11, "RdYlBu"))  # same “nice” diverging palette
## 0. 颜色：用稍微浅一点的 YlGnBu（去掉最深的那格）
ylgnbu_pal   <- RColorBrewer::brewer.pal(9, "YlGnBu")
ylgnbu_light <- ylgnbu_pal[2:7]  # 可以改成 [2:7] 更浅一点
pal=ylgnbu_light


plot_blip_one <- function(tt, mytit, digits = 3) {
  dd <- blip_long %>%
    filter(term == tt) %>%
    mutate(lab = sprintf(paste0("%.", digits, "f"), coef))
  
  lim_t <- max(abs(dd$coef), na.rm = TRUE)
  
  ggplot(dd, aes(x = tau1_f, y = tau2_f, fill = coef)) +
    geom_tile(color = "white", linewidth = 0.7) +
    geom_text(aes(label = lab), color = "black", size = 5) +
    coord_equal() +
    scale_fill_gradientn(colors = pal, limits = c(-lim_t, lim_t), oob = squish) +
    labs(title = mytit, x = expression(tau[1]), y = expression(tau[2])) +
    theme_minimal(base_size = 16) +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 16, color = "black"),
      axis.text  = element_text(size = 14, color = "black"),
      panel.grid = element_blank(),
      plot.margin = margin(3, 3, 3, 3, unit = "mm")
    )
}

terms <- unique(blip_long$term)
titlename<-c("Intercept","Sex","Race1", "Race2",
             "Age","Diabetes","Insurance","Smoking1",
             "Smoking2","Compliance","Base CAL",
             "Base PPD")
for (i in 1:12) {
  p <- plot_blip_one(terms[i], titlename[i], digits = 3)
  p
}





