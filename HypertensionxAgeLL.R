# =============================================================================
# Ischemic stroke incidence by age and hypertension
#
# Cohort of people followed for a first ischemic stroke. Hypertension is the
# exposure of interest. Age is modelled both in bands and as a spline, and
# the hypertension effect is allowed to change with age.
#
# Two follow-up definitions are analysed, in this order:
#
#   1. Full follow-up. Oral anticoagulation (OAC) is a time-varying state.
#      The OAC state starts at purchase_date and stops at
#      last_anticoagulation_date_plus_120_days.
#   2. Follow-up ends at exit_date_purchase_stroke_death_or_end: purchase,
#      or ischemic stroke, or death, or end of follow-up. OAC is not a
#      state in this analysis.
#
# Outcome on each row of follow-up: ischemic stroke (Lexis exit state "IS").
# Death without a stroke, and the administrative end of follow-up, are
# censoring (exit state "EOF"). They are not modelled as competing events.
#
# Input is one row per person in an SPSS .sav file. Column names used
# below are the names in that file.
#
# Rates are events per person-year. Tables and absolute differences are
# scaled per 100 person-years. Model coefficients are log incidence rates
# (Epi family "poisreg"); exponentiated coefficients are incidence rate
# ratios.
#
# Knots, calendar breaks, and the age range of the plots are the values used
# in the original analysis. They are set in the configuration section.
# =============================================================================

# ---- Packages ---------------------------------------------------------------
# Epi:        Lexis objects, cutLexis, splitLexis, ci.exp, Relevel, poisreg
# splines:    natural splines ns() for age
# cohorttools: mkratetable(); not on CRAN, install it before running
# haven:      read SPSS .sav
# survival:   Surv() in the rate-table formula
# knitr:      HTML tables
# ggplot2, RColorBrewer, grid, gridExtra: figures

library(cohorttools)
library(haven)
library(ggplot2)
library(RColorBrewer)
library(Epi)
library(splines)
library(survival)
library(knitr)
library(grid)
library(gridExtra)

# ---- Configuration ----------------------------------------------------------

# Person-level incident cohort. One row per person.
cohort_file <- "incident_cohort.sav"

# Figures and HTML tables are written here (the working directory).
output_dir <- "."

# Calendar-time split points (decimal years). Follow-up is cut at each
# break so calendar period can be used as a covariate.
calendar_breaks <- c(2009, 2011, 2013, 2015, 2017)

# Natural-spline knots for age (years).
#
# Adjusted models, full follow-up. Other knot sets were tried and not used:
#   c(63, 71, 78, 84)
#   c(65, 75, 83)
age_knots_adjusted_with_oac <- c(65, 80, 83)

# Crude age curves, full follow-up.
age_knots_crude_with_oac <- c(65, 74, 82)

# Adjusted models, follow-up censored at OAC.
# The original script never assigned these knots directly. In a
# top-to-bottom run they were whatever was last stored in `kn.age`,
# which was the crude full-follow-up set c(65, 74, 82). They are written
# here so this section does not depend on run order.
age_knots_adjusted_without_oac <- c(65, 74, 82)

# Crude age curves, follow-up censored at OAC.
age_knots_crude_without_oac <- c(65, 74)

# Ages at which curves are predicted. The axes of the figures are narrower
# than this grid (see coord_cartesian in each plot).
prediction_ages <- seq(20, 107, length.out = 87)
prediction_ages_fine <- seq(20, 107, length.out = 200)

# Full-follow-up rate-ratio curve and the single crude age curve use this
# pastel pair (line, ribbon). Rate-difference curves, and the curves for
# follow-up censored at OAC, use ColorBrewer Set1 (red, blue, green).
pastel_colors <- c("#AEC6CF", "#FFB347", "#B39EB5", "#77DD77", "#FF6961")
set1_colors <- brewer.pal(n = 3, name = "Set1")

# ---- Small helpers ----------------------------------------------------------

# Binary 0/1 fields arrive from SPSS as labelled values. `unclass` keeps the
# numeric code. Reference level is "no".
as_no_yes <- function(x) {
  factor(unclass(x), levels = 0:1, labels = c("no", "yes"))
}

# Print a one-way frequency table, including missing values.
print_table <- function(x, label) {
  cat("\n", label, "\n", sep = "")
  print(table(x, useNA = "always"))
  invisible(x)
}

# Baseline covariates. These do not change over follow-up, so they are
# coded once on the person-level data and then carried into each Lexis object.
#
# Clinical fields were assessed before oral anticoagulation.
# `cancer` is cancer recorded before or at cohort entry.
# `alcohol` is a baseline alcohol-related diagnosis.
# Income and diabetes medication are coded for description. They are not
# covariates in the regression models below. Income tertile 1 = low,
# 2 = mid, 3 = high. Diabetes medication codes: 1 = insulin only,
# 2 = insulin plus oral therapy, 3 = oral therapy only,
# 4 = no antidiabetic medication. Code 0 is not labelled.
recode_baseline_covariates <- function(data) {
  data$sex <- factor(
    unclass(data$sex_code),
    levels = 0:1,
    labels = c("female", "male")
  )
  data$hyperlipidemia <- as_no_yes(data$hyperlipidemia_code)
  data$hypertension <- as_no_yes(data$hypertension_code)
  data$diabetes <- as_no_yes(data$diabetes_code)
  data$heart_failure <- as_no_yes(data$heart_failure_code)
  data$abnormal_liver_function <- as_no_yes(data$abnormal_liver_code)
  data$abnormal_renal_function <- as_no_yes(data$abnormal_renal_code)
  data$prior_bleeding <- as_no_yes(data$prior_bleeding_code)
  data$alcohol <- as_no_yes(data$alcohol_code)
  data$vascular_disease <- as_no_yes(data$vascular_disease_code)
  data$cancer <- as_no_yes(data$cancer_code)
  data$dementia <- as_no_yes(data$dementia_code)
  data$psychiatric_disease <- as_no_yes(data$psychiatric_disease_code)
  data$prior_ischemic_stroke <- as_no_yes(data$prior_ischemic_stroke_code)
  data$income_tertile <- factor(
    unclass(data$income_tertile_code),
    levels = 1:3,
    labels = c("low", "mid", "high")
  )
  data$diabetes_medication <- factor(
    unclass(data$diabetes_medication_code),
    levels = 0:4
  )

  print_table(data$sex, "Sex (0 = female, 1 = male)")
  print_table(data$hyperlipidemia, "Hyperlipidemia")
  print_table(data$hypertension, "Hypertension")
  print_table(data$diabetes, "Diabetes")
  print_table(data$heart_failure, "Heart failure")
  print_table(data$abnormal_liver_function, "Abnormal liver function")
  print_table(data$abnormal_renal_function, "Abnormal renal function")
  print_table(data$prior_bleeding, "Prior bleeding")
  print_table(data$alcohol, "Alcohol-related diagnosis")
  print_table(data$vascular_disease, "Vascular disease")
  print_table(data$cancer, "Cancer before or at cohort entry")
  print_table(data$dementia, "Dementia")
  print_table(data$psychiatric_disease, "Psychiatric disease")
  print_table(data$prior_ischemic_stroke, "Prior ischemic stroke")
  print_table(data$income_tertile, "Income tertile (not used in the models)")
  print_table(data$diabetes_medication, "Diabetes medication code (not used in the models)")
  data
}

# Lexis object on three time scales (package Epi):
#   age   age in years at the start of the interval
#   fu    years since cohort entry (0 at entry)
#   per   calendar time in decimal years (Epi::cal.yr)
#
# `exit_date` and `event_flag` are column names in `data`.
# Entry state is SOF (start of follow-up) for every person.
# Exit state is IS when the event flag is 1, otherwise EOF.
# Lexis drops rows with follow-up length below its tolerance and prints
# how many were dropped.
build_lexis <- function(data, exit_date, event_flag) {
  entry_status <- factor(
    rep("SOF", nrow(data)),
    levels = c("SOF", "EOF", "IS")
  )
  exit_status <- factor(
    ifelse(data[[event_flag]] == 1, "IS", "EOF"),
    levels = c("SOF", "EOF", "IS")
  )
  print_table(entry_status, "Entry state")
  print_table(exit_status, "Exit state")

  entry_calendar <- cal.yr(data$CohortEntryDate, format = "%Y-%m-%d")
  exit_calendar <- cal.yr(data[[exit_date]], format = "%Y-%m-%d")

  Lexis(
    entry = list(age = data$Age, fu = 0, per = entry_calendar),
    duration = exit_calendar - entry_calendar,
    entry.status = entry_status,
    exit.status = exit_status,
    id = data$SID,
    data = data
  )
}

# Split follow-up by calendar year and attach the time variables used
# in the models.
#
# After splitLexis, `age` is age at the left edge of each slice and
# `lex.dur` is the slice length in years.
# `age` is then moved to the middle of the slice. Continuous-age models
# use that midpoint.
# `age_group` is cut from midpoint + lex.dur/2, which is age at the
# right edge of the slice. That is what the original script did; it is
# not the same age as the continuous-age models.
# `calendar_period` is the calendar band of the slice (a factor).
# `calendar_time_mid` is decimal calendar time at the middle of the slice.
# It is created for description and is not used in the models.
prepare_split_followup <- function(lexis, calendar_breaks) {
  lexis <- splitLexis(lex = lexis, time.scale = "per", breaks = calendar_breaks)
  lexis$calendar_period <- timeBand(lex = lexis, time.scale = "per", type = "factor")
  lexis$calendar_time_mid <- with(lexis, per + lex.dur / 2)
  lexis$age <- lexis$age + (lexis$lex.dur / 2)
  lexis$age_group <- cut(
    lexis$age + lexis$lex.dur / 2,
    breaks = c(0, 50, 60, 70, 80, 90, Inf)
  )
  lexis
}

# First level of a factor after sorting. For the no/yes factors this is
# "no". For sex it is "female". For calendar period it is the earliest
# band. Predictions use this profile. It is not the most common value
# and not the sample mean.
reference_level <- function(x) sort(unique(x))[1]

# Write one or more knitr tables to a single HTML file.
write_html <- function(path, expr) {
  sink(path)
  on.exit(sink(), add = TRUE)
  force(expr)
  invisible(path)
}

# Shared figure style.
theme_incidence <- function() {
  theme_minimal(base_family = "sans") +
    theme(
      text = element_text(color = "gray10"),
      plot.title = element_text(size = 15, face = "bold"),
      axis.title = element_text(size = 15),
      axis.text = element_text(size = 13),
      panel.grid.major = element_line(color = "gray90"),
      panel.grid.minor = element_blank()
    )
}

# PDF and JPEG of the same figure.
save_plot <- function(plot, filename, width, height) {
  ggsave(
    file.path(output_dir, paste0(filename, ".pdf")),
    plot, width = width, height = height, dpi = 300
  )
  ggsave(
    file.path(output_dir, paste0(filename, ".jpg")),
    plot, width = width, height = height, dpi = 300
  )
}

# Absolute rate difference, hypertension minus no hypertension,
# in events per 100 person-years.
#
# `pred_yes` and `pred_no` are lists with $fit and $se.fit on the
# log-rate scale (predict(..., type = "link")).
# The standard error uses the delta method on each rate and then adds
# the variances. Covariance between the two linear predictors is not
# included. That is the interval from the original script.
# The 1.96 multiplier is a normal 95% interval.
absolute_rate_difference <- function(pred_yes, pred_no, age) {
  rate_yes <- exp(pred_yes$fit)
  rate_no <- exp(pred_no$fit)
  rate_diff <- 100 * (rate_yes - rate_no)
  rate_diff_se <- 100 * sqrt(
    (rate_yes * pred_yes$se.fit)^2 +
      (rate_no * pred_no$se.fit)^2
  )
  data.frame(
    age = age,
    rate_diff = rate_diff,
    lo = rate_diff - 1.96 * rate_diff_se,
    hi = rate_diff + 1.96 * rate_diff_se
  )
}

# ---- Read the cohort --------------------------------------------------------

if (!file.exists(cohort_file)) {
  stop(
    "Cohort file not found: ", cohort_file,
    "\nSet `cohort_file` at the top of the script.",
    call. = FALSE
  )
}

cohort <- data.frame(read_sav(cohort_file))
cohort <- recode_baseline_covariates(cohort)

cat("\nPersons in the cohort:", nrow(cohort), "\n")

# Calendar span of the exit dates, as a check that the date variables
# parsed. Death date is not the Lexis exit time; follow-up length uses
# exit_date_full_followup or exit_date_purchase_stroke_death_or_end.
cat("\nFull follow-up exit date, decimal calendar years:\n")
print(range(cal.yr(cohort$exit_date_full_followup, format = "%Y-%m-%d")))
cat("\nDeath date, decimal calendar years:\n")
print(range(cal.yr(cohort$death_date, format = "%Y-%m-%d"), na.rm = TRUE))
cat("\nExit date (purchase or stroke or death or end), decimal calendar years:\n")
print(range(cal.yr(cohort$exit_date_purchase_stroke_death_or_end, format = "%Y-%m-%d")))

# =============================================================================
# 1. Full follow-up, with OAC as a time-varying state
# =============================================================================
#
# States after the two cuts:
#   SOF           follow-up before the OAC state (start of follow-up)
#   OAC           from purchase_date
#   OAC_stopped   from last_anticoagulation_date_plus_120_days
#   IS            ischemic stroke, absorbing
#   EOF           end of follow-up without a stroke, absorbing
#
# A missing purchase_date is not used as a cut.

cat("\n\n==== Full follow-up, OAC as a time-varying state ====\n")

lexis_entry <- build_lexis(
  data = cohort,
  exit_date = "exit_date_full_followup",
  event_flag = "ischemic_stroke_in_followup"
)
summary(lexis_entry)
timeScales(lexis_entry)

lexis_oac <- cutLexis(
  data = lexis_entry,
  cut = cal.yr(lexis_entry$purchase_date, format = "%Y-%m-%d"),
  timescale = "per",
  new.state = "OAC",
  new.scale = "OAC_time"
)
summary(lexis_oac)

lexis_oac <- cutLexis(
  data = lexis_oac,
  cut = cal.yr(lexis_oac$last_anticoagulation_date_plus_120_days, format = "%Y-%m-%d"),
  timescale = "per",
  new.state = "OAC_stopped",
  new.scale = "OAC_stopped_time"
)
summary(lexis_oac)

# Multistate diagram. show.persons = FALSE prints transition counts only.
boxesLx(lexis_oac, show.persons = FALSE)

cat("\nCalendar time at the start of each slice, before the year split:\n")
print(range(lexis_oac$per))

lexis_with_oac <- prepare_split_followup(lexis_oac, calendar_breaks)

cat("\nAge group by exit state of the slice:\n")
print(table(lexis_with_oac$age_group, lexis_with_oac$lex.Xst, useNA = "always"))

# ---- Rate table, full follow-up ---------------------------------------------
#
# mkratetable: events, person-time, rate, and rate ratios.
# scale = 100 reports rates per 100 person-years.
# lex.Xst == "IS" is 1 only on the slice that ends in a stroke.
# Surv(time, event) here is the person-time and the event indicator;
# it is not a Cox model.
#
# Relevel(lex.Cst, list(1, 2, 3, 4)) keeps the first four current-state
# levels and sets any later level to missing (Epi::Relevel). The levels
# are printed first so it is clear which four enter the table.
# Under treatment contrasts the first of those four is the reference.

cat("\nCurrent Lexis state, all levels, before Relevel:\n")
print(table(lexis_with_oac$lex.Cst, useNA = "always"))
cat("Level order:\n")
print(levels(lexis_with_oac$lex.Cst))

rate_table_with_oac <- mkratetable(
  Surv(lex.dur, lex.Xst == "IS") ~
    calendar_period + age_group +
    Relevel(lex.Cst, list(1, 2, 3, 4)) +
    sex + hyperlipidemia + hypertension + diabetes +
    psychiatric_disease + dementia + cancer + vascular_disease +
    diabetes_medication + alcohol + prior_bleeding +
    prior_ischemic_stroke + abnormal_liver_function +
    abnormal_renal_function + heart_failure,
  data = lexis_with_oac,
  scale = 100,
  add.RR = TRUE
)

write_html(file.path(output_dir, "rate_table_full_followup.html"), {
  print(kable(
    rate_table_with_oac,
    caption = "Rate table (1/100 person years)",
    format = "html",
    digits = 3
  ))
})

# ---- Categorical Poisson models, full follow-up -----------------------------
#
# family = "poisreg" (Epi): the response cbind(event, person-time) is a
# rate. Coefficients are log rate ratios.
# The same four Lexis states as in the rate table are used here.
# The second model adds a diabetes x calendar-period interaction and is
# a check only; the hypertension figures do not use it.

model_categorical_with_oac <- glm(
  cbind(lex.Xst == "IS", lex.dur) ~
    Relevel(lex.Cst, list(1, 2, 3, 4)) + sex + age_group + hypertension +
    hyperlipidemia + psychiatric_disease + dementia + cancer +
    vascular_disease + alcohol + prior_bleeding + prior_ischemic_stroke +
    abnormal_liver_function + abnormal_renal_function + heart_failure +
    diabetes + calendar_period,
  data = lexis_with_oac,
  family = "poisreg"
)

model_categorical_diabetes_period_with_oac <- update(
  model_categorical_with_oac,
  ~ . + diabetes:calendar_period
)

cat("\nIncidence rate ratios, age in bands, full follow-up:\n")
print(round(ci.exp(model_categorical_with_oac), 3))
cat("\nWith a diabetes x calendar-period interaction:\n")
print(round(ci.exp(model_categorical_diabetes_period_with_oac), 3))

write_html(file.path(output_dir, "irr_full_followup.html"), {
  print(kable(
    round(ci.exp(model_categorical_with_oac), 3),
    caption = "IRR from Poisson regression model",
    format = "html"
  ))
  print(kable(
    round(ci.exp(model_categorical_diabetes_period_with_oac), 3),
    caption = "IRR from Poisson regression model",
    format = "html"
  ))
})

# ---- Continuous age, full follow-up -----------------------------------------
# Added 21 October 2024.
# Age in the models below is the midpoint of each slice.

# Diagram scaled to transitions per 1000 persons.
boxesLx(
  lexis_with_oac,
  prop.penwidth = TRUE,
  scale.Y = 1000,
  show.loop = FALSE,
  show.persons = TRUE
)

# Age distribution used when the knots were chosen. Not used by the models.
cat("\nAge at the middle of each slice, full follow-up:\n")
print(summary(lexis_with_oac$age))
cat("\nQuantiles of slice-midpoint age (full follow-up):\n")
print(trunc(quantile(lexis_with_oac$age, c(0:5) / 5)))
print(trunc(quantile(lexis_with_oac$age, c(0:3) / 3)))
print(trunc(quantile(lexis_with_oac$age, c(0:4) / 4)))
print(trunc(quantile(lexis_with_oac$age, c(0:10) / 10)))
cat("\nQuantiles of age at cohort entry, before OAC cuts and splitting:\n")
print(trunc(quantile(lexis_entry$age, c(0:3) / 3)))
print(trunc(quantile(lexis_entry$age, c(0:4) / 4)))

# Drop current-state levels that have no rows, so the spline models
# do not keep empty levels. The categorical models above were fitted
# before this drop, as in the original script.
lexis_with_oac$lex.Cst <- droplevels(lexis_with_oac$lex.Cst)

# Diabetes x age is a likelihood-ratio check. The figures below are the
# hypertension models, not this diabetes model.
model_age_spline_with_oac <- glm(
  cbind(lex.Xst == "IS", lex.dur) ~
    lex.Cst + sex + ns(age, knots = age_knots_adjusted_with_oac) +
    hypertension + hyperlipidemia + psychiatric_disease + dementia + cancer +
    vascular_disease + alcohol + prior_bleeding + prior_ischemic_stroke +
    abnormal_liver_function + abnormal_renal_function + heart_failure +
    diabetes + calendar_period,
  data = lexis_with_oac,
  family = "poisreg"
)

model_age_spline_diabetes_with_oac <- update(
  model_age_spline_with_oac,
  ~ . + diabetes:ns(age, knots = age_knots_adjusted_with_oac)
)

cat("\nLikelihood-ratio test, diabetes x age spline, full follow-up:\n")
print(anova(
  model_age_spline_with_oac,
  model_age_spline_diabetes_with_oac,
  test = "Chisq"
))

# Hypertension, no age interaction, then hypertension x age spline.
# Unlike the categorical model, lex.Cst is not passed through Relevel,
# so every remaining current state is in the model.
model_hypertension_with_oac <- glm(
  cbind(lex.Xst == "IS", lex.dur) ~
    lex.Cst + sex + ns(age, knots = age_knots_adjusted_with_oac) +
    hypertension + diabetes + hyperlipidemia + psychiatric_disease +
    dementia + cancer + alcohol + prior_bleeding + prior_ischemic_stroke +
    abnormal_liver_function + abnormal_renal_function + heart_failure +
    vascular_disease + calendar_period,
  data = lexis_with_oac,
  family = "poisreg"
)

model_hypertension_age_with_oac <- update(
  model_hypertension_with_oac,
  ~ . + hypertension:ns(age, knots = age_knots_adjusted_with_oac)
)

cat("\nLikelihood-ratio test, hypertension x age spline, full follow-up:\n")
print(anova(
  model_hypertension_with_oac,
  model_hypertension_age_with_oac,
  test = "Chisq"
))

# Prediction grid for the hypertension rate ratio.
# lex.Cst varies over every current state. There is no state x
# hypertension interaction, so the rate ratio is the same in every
# state and the curves are drawn on top of each other.
# Every other covariate is fixed at reference_level().
prediction_grid_irr_with_oac <- with(
  lexis_with_oac,
  expand.grid(
    lex.Cst = sort(unique(lex.Cst)),
    sex = reference_level(sex),
    calendar_period = reference_level(calendar_period),
    hypertension = c("no", "yes"),
    diabetes = reference_level(diabetes),
    hyperlipidemia = reference_level(hyperlipidemia),
    psychiatric_disease = reference_level(psychiatric_disease),
    dementia = reference_level(dementia),
    cancer = reference_level(cancer),
    alcohol = reference_level(alcohol),
    prior_bleeding = reference_level(prior_bleeding),
    prior_ischemic_stroke = reference_level(prior_ischemic_stroke),
    abnormal_liver_function = reference_level(abnormal_liver_function),
    abnormal_renal_function = reference_level(abnormal_renal_function),
    vascular_disease = reference_level(vascular_disease),
    age = prediction_ages,
    heart_failure = reference_level(heart_failure)
  )
)

# Design matrix in the same term order as the fitted interaction model.
# The interaction term is written ns(age):hypertension here and
# hypertension:ns(age) in the model formula, matching the original script.
design_irr_with_oac <- model.matrix(
  ~ lex.Cst + sex + ns(age, knots = age_knots_adjusted_with_oac) +
    hypertension + diabetes + hyperlipidemia + psychiatric_disease +
    dementia + cancer + alcohol + prior_bleeding + prior_ischemic_stroke +
    abnormal_liver_function + abnormal_renal_function + heart_failure +
    vascular_disease + calendar_period +
    ns(age, knots = age_knots_adjusted_with_oac):hypertension,
  data = prediction_grid_irr_with_oac
)

# Contrast: log rate with hypertension minus log rate without it.
# Rows line up because both subsets come from the same expand.grid.
contrast_irr_with_oac <-
  design_irr_with_oac[prediction_grid_irr_with_oac$hypertension == "yes", ] -
  design_irr_with_oac[prediction_grid_irr_with_oac$hypertension == "no", ]

irr_by_age_with_oac <- data.frame(
  ci.exp(model_hypertension_age_with_oac, ctr.mat = contrast_irr_with_oac)
)
colnames(irr_by_age_with_oac)[1:3] <- c(
  "rate_ratio", "rate_ratio_low", "rate_ratio_high"
)
irr_by_age_with_oac$age <- prediction_grid_irr_with_oac$age[
  prediction_grid_irr_with_oac$hypertension == "yes"
]

# Adjusted incidence rate ratio, hypertension versus no hypertension.
plot_irr_with_oac <- ggplot(irr_by_age_with_oac, aes(x = age, y = rate_ratio)) +
  geom_ribbon(
    aes(ymin = rate_ratio_low, ymax = rate_ratio_high),
    fill = pastel_colors[2],
    alpha = 0.3
  ) +
  geom_line(color = pastel_colors[1], linewidth = 1.2) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray40", linewidth = 0.6) +
  scale_y_continuous(breaks = seq(0, 5, 1)) +
  coord_cartesian(ylim = c(0, 5), xlim = c(40, 100)) +
  labs(
    x = "Age (years)",
    y = "Adjusted incidence rate ratio",
    title = "Hypertension vs. no hypertension"
  ) +
  theme_incidence()

# ---- Absolute rate difference, full follow-up -------------------------------
#
# Same covariates and the same age x hypertension interaction as
# model_hypertension_age_with_oac. Refitted here, as in the original
# script, and then used for predicted rates.
# The grid fixes lex.Cst at its first level, so this is one absolute
# difference, not one curve per state.

model_rate_difference_with_oac <- glm(
  cbind(lex.Xst == "IS", lex.dur) ~
    ns(age, knots = age_knots_adjusted_with_oac) * hypertension +
    diabetes + hyperlipidemia + psychiatric_disease + dementia + cancer +
    alcohol + prior_bleeding + prior_ischemic_stroke +
    abnormal_liver_function + abnormal_renal_function + heart_failure +
    vascular_disease + calendar_period + lex.Cst + sex,
  data = lexis_with_oac,
  family = "poisreg"
)

prediction_grid_ard_with_oac <- with(
  lexis_with_oac,
  expand.grid(
    lex.Cst = reference_level(lex.Cst),
    sex = reference_level(sex),
    calendar_period = reference_level(calendar_period),
    hypertension = c("no", "yes"),
    diabetes = reference_level(diabetes),
    hyperlipidemia = reference_level(hyperlipidemia),
    psychiatric_disease = reference_level(psychiatric_disease),
    dementia = reference_level(dementia),
    cancer = reference_level(cancer),
    alcohol = reference_level(alcohol),
    prior_bleeding = reference_level(prior_bleeding),
    prior_ischemic_stroke = reference_level(prior_ischemic_stroke),
    abnormal_liver_function = reference_level(abnormal_liver_function),
    abnormal_renal_function = reference_level(abnormal_renal_function),
    heart_failure = reference_level(heart_failure),
    vascular_disease = reference_level(vascular_disease),
    age = prediction_ages
  )
)

predicted_rates_with_oac <- predict(
  model_rate_difference_with_oac,
  newdata = prediction_grid_ard_with_oac,
  type = "link",
  se.fit = TRUE
)

yes_with_oac <- prediction_grid_ard_with_oac$hypertension == "yes"
no_with_oac <- prediction_grid_ard_with_oac$hypertension == "no"

rate_difference_with_oac <- absolute_rate_difference(
  pred_yes = list(
    fit = predicted_rates_with_oac$fit[yes_with_oac],
    se.fit = predicted_rates_with_oac$se.fit[yes_with_oac]
  ),
  pred_no = list(
    fit = predicted_rates_with_oac$fit[no_with_oac],
    se.fit = predicted_rates_with_oac$se.fit[no_with_oac]
  ),
  age = prediction_grid_ard_with_oac$age[yes_with_oac]
)

plot_rate_difference_with_oac <- ggplot(
  rate_difference_with_oac,
  aes(x = age, y = rate_diff)
) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = set1_colors[2], alpha = 0.3) +
  geom_line(color = set1_colors[1], linewidth = 1.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.6) +
  coord_cartesian(xlim = c(40, 100), ylim = c(-2, 4)) +
  scale_y_continuous(breaks = seq(-2, 4, 1)) +
  labs(
    x = "Age (years)",
    y = "Rate difference (events per 100 patient-years)",
    title = "Hypertension vs. no hypertension"
  ) +
  theme_incidence()

print(plot_irr_with_oac)
print(plot_rate_difference_with_oac)
save_plot(plot_irr_with_oac, "hypertension_irr_full_followup", width = 8, height = 6)
save_plot(
  plot_rate_difference_with_oac,
  "hypertension_rate_difference_full_followup",
  width = 8,
  height = 6
)

# ---- Crude rates by age, full follow-up -------------------------------------
# No covariate adjustment. Separate knots from the adjusted models.
# The first curve is the whole cohort. The second fits a separate
# age spline in each hypertension group (age x hypertension, no other
# covariates). These are rates, not rate differences.

model_crude_age_with_oac <- glm(
  cbind(lex.Xst == "IS", lex.dur) ~ ns(age, knots = age_knots_crude_with_oac),
  data = lexis_with_oac,
  family = "poisreg"
)

crude_age_grid_with_oac <- data.frame(age = prediction_ages_fine)
crude_age_pred_with_oac <- predict(
  model_crude_age_with_oac,
  newdata = crude_age_grid_with_oac,
  type = "link",
  se.fit = TRUE
)
crude_age_grid_with_oac$rate <- 100 * exp(crude_age_pred_with_oac$fit)
crude_age_grid_with_oac$lo <- 100 * exp(
  crude_age_pred_with_oac$fit - 1.96 * crude_age_pred_with_oac$se.fit
)
crude_age_grid_with_oac$hi <- 100 * exp(
  crude_age_pred_with_oac$fit + 1.96 * crude_age_pred_with_oac$se.fit
)

plot_crude_rate_with_oac <- ggplot(crude_age_grid_with_oac, aes(x = age, y = rate)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = pastel_colors[2], alpha = 0.3) +
  geom_line(color = pastel_colors[1], linewidth = 1.2) +
  coord_cartesian(xlim = c(40, 100), ylim = c(0, 6)) +
  labs(
    title = "Entire follow-up",
    x = "Age (years)",
    y = "Crude rate per 100 patient-years"
  ) +
  theme_incidence()

model_crude_hypertension_with_oac <- glm(
  cbind(lex.Xst == "IS", lex.dur) ~
    ns(age, knots = age_knots_crude_with_oac) * hypertension,
  data = lexis_with_oac,
  family = "poisreg"
)

crude_hypertension_grid_with_oac <- expand.grid(
  age = prediction_ages_fine,
  hypertension = c("yes", "no")
)
crude_hypertension_pred_with_oac <- predict(
  model_crude_hypertension_with_oac,
  newdata = crude_hypertension_grid_with_oac,
  type = "link",
  se.fit = TRUE
)
crude_hypertension_grid_with_oac$rate <- 100 * exp(crude_hypertension_pred_with_oac$fit)
crude_hypertension_grid_with_oac$lo <- 100 * exp(
  crude_hypertension_pred_with_oac$fit - 1.96 * crude_hypertension_pred_with_oac$se.fit
)
crude_hypertension_grid_with_oac$hi <- 100 * exp(
  crude_hypertension_pred_with_oac$fit + 1.96 * crude_hypertension_pred_with_oac$se.fit
)

plot_crude_by_hypertension_with_oac <- ggplot(
  crude_hypertension_grid_with_oac,
  aes(x = age, y = rate, color = hypertension, fill = hypertension)
) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.20, color = NA) +
  geom_line(linewidth = 1.2) +
  coord_cartesian(xlim = c(40, 100), ylim = c(0, 8)) +
  labs(
    title = "Entire follow-up",
    x = "Age (years)",
    y = "Rate per 100 patient-years",
    color = "Hypertension",
    fill = "Hypertension"
  ) +
  theme_incidence()

print(plot_crude_rate_with_oac)
print(plot_crude_by_hypertension_with_oac)

# =============================================================================
# 2. Follow-up ending at purchase, stroke, death, or end
# =============================================================================
#
# Exit date is exit_date_purchase_stroke_death_or_end:
# purchase, or ischemic stroke, or death, or end of follow-up.
# The exit state is IS only when ischemic_stroke_before_oac == 1.
# There is no OAC state and no cutLexis.
# The cohort is not re-read. Baseline covariates were attached above.

cat("\n\n==== Follow-up ending at purchase, stroke, death, or end ====\n")

lexis_before_oac <- build_lexis(
  data = cohort,
  exit_date = "exit_date_purchase_stroke_death_or_end",
  event_flag = "ischemic_stroke_before_oac"
)
summary(lexis_before_oac)
timeScales(lexis_before_oac)
boxesLx(lexis_before_oac, show.persons = FALSE)

cat("\nCalendar time at the start of each slice, before the year split:\n")
print(range(lexis_before_oac$per))

lexis_without_oac <- prepare_split_followup(lexis_before_oac, calendar_breaks)

cat("\nAge group by exit state of the slice:\n")
print(table(lexis_without_oac$age_group, lexis_without_oac$lex.Xst, useNA = "always"))

# ---- Rate table, censored at OAC --------------------------------------------
# Current Lexis state is omitted. Before the stroke or the end of
# follow-up, person-time is in SOF, so lex.Cst does not vary across
# the rows that contribute time.

rate_table_without_oac <- mkratetable(
  Surv(lex.dur, lex.Xst == "IS") ~
    calendar_period + age_group +
    sex + hyperlipidemia + hypertension + diabetes +
    psychiatric_disease + dementia + cancer + vascular_disease +
    diabetes_medication + alcohol + prior_bleeding +
    prior_ischemic_stroke + abnormal_liver_function +
    abnormal_renal_function + heart_failure,
  data = lexis_without_oac,
  scale = 100,
  add.RR = TRUE
)

write_html(file.path(output_dir, "rate_table_censored_at_oac.html"), {
  print(kable(
    rate_table_without_oac,
    caption = "Rate table (1/100 person years)",
    format = "html",
    digits = 3
  ))
})

# ---- Categorical Poisson models, censored at OAC ----------------------------
# The original script wrote these ratios into the same HTML file as the
# full-follow-up model, and so replaced that file. They are written to
# irr_censored_at_oac.html instead.

model_categorical_without_oac <- glm(
  cbind(lex.Xst == "IS", lex.dur) ~
    sex + age_group + hypertension +
    hyperlipidemia + psychiatric_disease + dementia + cancer +
    vascular_disease + alcohol + prior_bleeding + prior_ischemic_stroke +
    abnormal_liver_function + abnormal_renal_function + heart_failure +
    diabetes + calendar_period,
  data = lexis_without_oac,
  family = "poisreg"
)

model_categorical_diabetes_period_without_oac <- update(
  model_categorical_without_oac,
  ~ . + diabetes:calendar_period
)

cat("\nIncidence rate ratios, age in bands, follow-up censored at OAC:\n")
print(round(ci.exp(model_categorical_without_oac), 3))
cat("\nWith a diabetes x calendar-period interaction:\n")
print(round(ci.exp(model_categorical_diabetes_period_without_oac), 3))

write_html(file.path(output_dir, "irr_censored_at_oac.html"), {
  print(kable(
    round(ci.exp(model_categorical_without_oac), 3),
    caption = "IRR from Poisson regression model",
    format = "html"
  ))
  print(kable(
    round(ci.exp(model_categorical_diabetes_period_without_oac), 3),
    caption = "IRR from Poisson regression model",
    format = "html"
  ))
})

# ---- Continuous age, censored at OAC ----------------------------------------

boxesLx(
  lexis_without_oac,
  prop.penwidth = TRUE,
  scale.Y = 1000,
  show.loop = FALSE,
  show.persons = TRUE
)

cat("\nAge at the middle of each slice, follow-up censored at OAC:\n")
print(summary(lexis_without_oac$age))
cat("\nQuantiles of slice-midpoint age (censored at OAC):\n")
print(trunc(quantile(lexis_without_oac$age, c(0:4) / 4)))
print(trunc(quantile(lexis_without_oac$age, c(0:3) / 3)))
print(trunc(quantile(lexis_without_oac$age, c(0:10) / 10)))

model_age_spline_without_oac <- glm(
  cbind(lex.Xst == "IS", lex.dur) ~
    sex + ns(age, knots = age_knots_adjusted_without_oac) +
    hypertension + hyperlipidemia + psychiatric_disease + dementia + cancer +
    vascular_disease + alcohol + prior_bleeding + prior_ischemic_stroke +
    abnormal_liver_function + abnormal_renal_function + heart_failure +
    diabetes + calendar_period,
  data = lexis_without_oac,
  family = "poisreg"
)

model_age_spline_diabetes_without_oac <- update(
  model_age_spline_without_oac,
  ~ . + diabetes:ns(age, knots = age_knots_adjusted_without_oac)
)

cat("\nLikelihood-ratio test, diabetes x age spline, censored at OAC:\n")
print(anova(
  model_age_spline_without_oac,
  model_age_spline_diabetes_without_oac,
  test = "Chisq"
))

model_hypertension_without_oac <- glm(
  cbind(lex.Xst == "IS", lex.dur) ~
    sex + ns(age, knots = age_knots_adjusted_without_oac) +
    hypertension + diabetes + hyperlipidemia + psychiatric_disease +
    dementia + cancer + alcohol + prior_bleeding + prior_ischemic_stroke +
    abnormal_liver_function + abnormal_renal_function + heart_failure +
    vascular_disease + calendar_period,
  data = lexis_without_oac,
  family = "poisreg"
)

model_hypertension_age_without_oac <- update(
  model_hypertension_without_oac,
  ~ . + hypertension:ns(age, knots = age_knots_adjusted_without_oac)
)

cat("\nLikelihood-ratio test, hypertension x age spline, censored at OAC:\n")
print(anova(
  model_hypertension_without_oac,
  model_hypertension_age_without_oac,
  test = "Chisq"
))

# Same contrast as in the full follow-up, without a Lexis-state term.
# sort(unique(hypertension)) keeps the factor levels of the fitted model.
prediction_grid_irr_without_oac <- with(
  lexis_without_oac,
  expand.grid(
    sex = reference_level(sex),
    calendar_period = reference_level(calendar_period),
    hypertension = sort(unique(hypertension)),
    diabetes = reference_level(diabetes),
    hyperlipidemia = reference_level(hyperlipidemia),
    psychiatric_disease = reference_level(psychiatric_disease),
    dementia = reference_level(dementia),
    cancer = reference_level(cancer),
    alcohol = reference_level(alcohol),
    prior_bleeding = reference_level(prior_bleeding),
    prior_ischemic_stroke = reference_level(prior_ischemic_stroke),
    abnormal_liver_function = reference_level(abnormal_liver_function),
    abnormal_renal_function = reference_level(abnormal_renal_function),
    heart_failure = reference_level(heart_failure),
    vascular_disease = reference_level(vascular_disease),
    age = prediction_ages
  )
)

design_irr_without_oac <- model.matrix(
  ~ sex + ns(age, knots = age_knots_adjusted_without_oac) +
    hypertension + diabetes + hyperlipidemia + psychiatric_disease +
    dementia + cancer + alcohol + prior_bleeding + prior_ischemic_stroke +
    abnormal_liver_function + abnormal_renal_function + heart_failure +
    vascular_disease + calendar_period +
    ns(age, knots = age_knots_adjusted_without_oac):hypertension,
  data = prediction_grid_irr_without_oac
)

contrast_irr_without_oac <-
  design_irr_without_oac[prediction_grid_irr_without_oac$hypertension == "yes", ] -
  design_irr_without_oac[prediction_grid_irr_without_oac$hypertension == "no", ]

irr_by_age_without_oac <- data.frame(
  ci.exp(model_hypertension_age_without_oac, ctr.mat = contrast_irr_without_oac)
)
colnames(irr_by_age_without_oac)[1:3] <- c(
  "rate_ratio", "rate_ratio_low", "rate_ratio_high"
)
irr_by_age_without_oac$age <- prediction_grid_irr_without_oac$age[
  prediction_grid_irr_without_oac$hypertension == "yes"
]

plot_irr_without_oac <- ggplot(irr_by_age_without_oac, aes(x = age, y = rate_ratio)) +
  geom_ribbon(
    aes(ymin = rate_ratio_low, ymax = rate_ratio_high),
    fill = set1_colors[2],
    alpha = 0.3
  ) +
  geom_line(color = set1_colors[1], linewidth = 1.2) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray40", linewidth = 0.6) +
  scale_y_continuous(breaks = seq(0, 5, by = 1)) +
  coord_cartesian(ylim = c(0, 5), xlim = c(40, 100)) +
  labs(
    x = "Age (years)",
    y = "Adjusted incidence rate ratio",
    title = "Hypertension vs. no hypertension (without OAC)"
  ) +
  theme_incidence()

print(plot_irr_without_oac)

# ---- Absolute rate difference, censored at OAC ------------------------------
# Uses the interaction model already fitted (model_hypertension_age_without_oac).
# The full-follow-up section refitted an equivalent model instead.

prediction_grid_ard_without_oac <- with(
  lexis_without_oac,
  expand.grid(
    sex = reference_level(sex),
    calendar_period = reference_level(calendar_period),
    hypertension = sort(unique(hypertension)),
    diabetes = reference_level(diabetes),
    hyperlipidemia = reference_level(hyperlipidemia),
    psychiatric_disease = reference_level(psychiatric_disease),
    dementia = reference_level(dementia),
    cancer = reference_level(cancer),
    alcohol = reference_level(alcohol),
    prior_bleeding = reference_level(prior_bleeding),
    prior_ischemic_stroke = reference_level(prior_ischemic_stroke),
    abnormal_liver_function = reference_level(abnormal_liver_function),
    abnormal_renal_function = reference_level(abnormal_renal_function),
    heart_failure = reference_level(heart_failure),
    vascular_disease = reference_level(vascular_disease),
    age = prediction_ages
  )
)

predicted_rates_without_oac <- predict(
  model_hypertension_age_without_oac,
  newdata = prediction_grid_ard_without_oac,
  type = "link",
  se.fit = TRUE
)

yes_without_oac <- prediction_grid_ard_without_oac$hypertension == "yes"
no_without_oac <- prediction_grid_ard_without_oac$hypertension == "no"

rate_difference_without_oac <- absolute_rate_difference(
  pred_yes = list(
    fit = predicted_rates_without_oac$fit[yes_without_oac],
    se.fit = predicted_rates_without_oac$se.fit[yes_without_oac]
  ),
  pred_no = list(
    fit = predicted_rates_without_oac$fit[no_without_oac],
    se.fit = predicted_rates_without_oac$se.fit[no_without_oac]
  ),
  age = prediction_grid_ard_without_oac$age[yes_without_oac]
)

plot_rate_difference_without_oac <- ggplot(
  rate_difference_without_oac,
  aes(x = age, y = rate_diff)
) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = set1_colors[2], alpha = 0.3) +
  geom_line(color = set1_colors[1], linewidth = 1.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.6) +
  coord_cartesian(xlim = c(40, 100), ylim = c(-2, 4)) +
  scale_y_continuous(breaks = seq(-2, 4, by = 1)) +
  labs(
    x = "Age (years)",
    y = "Rate difference (events per 100 patient-years)",
    title = "Hypertension vs. no hypertension (without OAC)"
  ) +
  theme_incidence()

print(plot_rate_difference_without_oac)
save_plot(plot_irr_without_oac, "hypertension_irr_censored_at_oac", width = 8, height = 6)
save_plot(
  plot_rate_difference_without_oac,
  "hypertension_rate_difference_censored_at_oac",
  width = 8,
  height = 6
)

# ---- Crude rates by age, censored at OAC ------------------------------------

model_crude_age_without_oac <- glm(
  cbind(lex.Xst == "IS", lex.dur) ~ ns(age, knots = age_knots_crude_without_oac),
  data = lexis_without_oac,
  family = "poisreg"
)

# This overall curve uses 87 ages. The by-hypertension curve below uses 200.
crude_age_grid_without_oac <- data.frame(age = prediction_ages)
crude_age_pred_without_oac <- predict(
  model_crude_age_without_oac,
  newdata = crude_age_grid_without_oac,
  type = "link",
  se.fit = TRUE
)
crude_age_grid_without_oac$rate <- 100 * exp(crude_age_pred_without_oac$fit)
crude_age_grid_without_oac$lo <- 100 * exp(
  crude_age_pred_without_oac$fit - 1.96 * crude_age_pred_without_oac$se.fit
)
crude_age_grid_without_oac$hi <- 100 * exp(
  crude_age_pred_without_oac$fit + 1.96 * crude_age_pred_without_oac$se.fit
)

plot_crude_rate_without_oac <- ggplot(
  crude_age_grid_without_oac,
  aes(x = age, y = rate)
) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = set1_colors[2], alpha = 0.3) +
  geom_line(color = set1_colors[1], linewidth = 1.2) +
  coord_cartesian(xlim = c(40, 100), ylim = c(0, 6)) +
  labs(
    title = "Follow-up without anticoagulation",
    x = "Age (years)",
    y = "Crude rate per 100 patient-years"
  ) +
  theme_incidence()

print(plot_crude_rate_without_oac)

model_crude_hypertension_without_oac <- glm(
  cbind(lex.Xst == "IS", lex.dur) ~
    ns(age, knots = age_knots_crude_without_oac) * hypertension,
  data = lexis_without_oac,
  family = "poisreg"
)

crude_hypertension_grid_without_oac <- expand.grid(
  age = prediction_ages_fine,
  hypertension = c("yes", "no")
)
crude_hypertension_pred_without_oac <- predict(
  model_crude_hypertension_without_oac,
  newdata = crude_hypertension_grid_without_oac,
  type = "link",
  se.fit = TRUE
)
crude_hypertension_grid_without_oac$rate <- 100 * exp(
  crude_hypertension_pred_without_oac$fit
)
crude_hypertension_grid_without_oac$lo <- 100 * exp(
  crude_hypertension_pred_without_oac$fit -
    1.96 * crude_hypertension_pred_without_oac$se.fit
)
crude_hypertension_grid_without_oac$hi <- 100 * exp(
  crude_hypertension_pred_without_oac$fit +
    1.96 * crude_hypertension_pred_without_oac$se.fit
)

plot_crude_by_hypertension_without_oac <- ggplot(
  crude_hypertension_grid_without_oac,
  aes(x = age, y = rate, color = hypertension, fill = hypertension)
) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.20, color = NA) +
  geom_line(linewidth = 1.2) +
  coord_cartesian(xlim = c(40, 100), ylim = c(0, 8)) +
  labs(
    title = "Follow-up without anticoagulation",
    x = "Age (years)",
    y = "Rate per 100 patient-years",
    color = "Hypertension",
    fill = "Hypertension"
  ) +
  theme_incidence()

print(plot_crude_by_hypertension_without_oac)

# Side-by-side crude rates by hypertension: full follow-up (left) and
# follow-up censored at OAC (right).
crude_rates_panel <- grid.arrange(
  plot_crude_by_hypertension_with_oac,
  plot_crude_by_hypertension_without_oac,
  ncol = 2,
  top = textGrob(
    "Crude ischemic stroke rate by age",
    gp = gpar(fontsize = 18, fontface = "bold")
  )
)
save_plot(crude_rates_panel, "hypertension_crude_rates", width = 12, height = 6)

# Package versions, for the deposit.
writeLines(
  capture.output(sessionInfo()),
  file.path(output_dir, "session_info.txt")
)
