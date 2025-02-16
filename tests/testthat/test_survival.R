context("baggr() calls with survival regression model")
library(baggr)
library(testthat)


# prepare inputs ----------------------------------------------------------
set.seed(1990)

N = 400
K = 6
treatment = rbinom(N, 1, 0.5)
x = rnorm(N)
intervals = matrix(rweibull(N*2, 1, scale = exp(1 + 0.2*x + 5*treatment)), ncol = 2)
df_surv = data.frame(
  treatment = treatment, 
  interval_left = apply(intervals, 1, min),
  interval_right  = apply(intervals, 1, max), 
  censoring = rbinom(N, 1, 0.5), 
  site = sample(1:K, size = N, replace = TRUE), 
  x = x
)
df_surv[df_surv[, "censoring"] == 1, "interval_right"] = Inf

# tests ----------------------------------------------------------
# library(baggr)
bg5_n <- expect_warning(baggr(df_surv,
                              "survival", 
                              outcome = c("interval_left", "interval_right"),
                              group = "site", 
                              pooling = "none",
                              iter = 150, chains = 2, refresh = 0,
                              show_messages = F))
bg5_p <- expect_warning(baggr(df_surv,
                              "survival", 
                              outcome = c("interval_left", "interval_right"),
                              group = "site", 
                              pooling = "partial",
                              iter = 150, chains = 2, refresh = 0,
                              show_messages = F))

bg5_f <- expect_warning(baggr(df_surv,
                              "survival", 
                              outcome = c("interval_left", "interval_right"),
                              group = "site", 
                              pooling = "full",
                              iter = 150, chains = 2, refresh = 0,
                              show_messages = F))

bg5_ppd <- expect_warning(baggr(df_surv,
                              "survival", 
                              outcome = c("interval_left", "interval_right"),
                              group = "site", 
                              ppd = TRUE,
                              iter = 150, chains = 2, refresh = 0,
                              show_messages = F))


test_that("Estimated TEs approx correct", {

  te5_p = colMeans(group_effects(bg5_p))

  lapply(te5_p, expect_lte, 6)
  lapply(te5_p, expect_gte, 4)
  te5_f = unique(colMeans(group_effects(bg5_f))) 
  expect_lte(te5_f, 6)
  expect_gte(te5_f, 4)
})




test_that("Different pooling methods work for Rubin model", {
  expect_is(bg5_n, "baggr")
  expect_is(bg5_p, "baggr")
  expect_is(bg5_f, "baggr")
})


test_that("Various attr of baggr object are correct", {
  expect_equal(bg5_n$pooling, "none")
  expect_equal(bg5_p$pooling, "partial")
  expect_equal(bg5_f$pooling, "full")
  expect_equal(bg5_p$n_parameters, 1)
  expect_equal(bg5_p$n_groups, K)
  expect_equal(bg5_p$effects, "mean")
  expect_equal(bg5_p$model, "survival")
  expect_is(bg5_p$fit, "stanfit")
})

test_that("Data are available in baggr object", {
  expect_is(bg5_n$data, "data.frame")
  expect_is(bg5_p$data, "data.frame")
  expect_is(bg5_f$data, "data.frame")
  expect_is(bg5_n$summary_data, "data.frame")
  expect_is(bg5_p$summary_data, "data.frame")
  expect_is(bg5_f$summary_data, "data.frame")
})

test_that("Pooling metrics", {
  # all pooling stats are 0 if no pooling
  expect_equal(unique(as.numeric(bg5_n$pooling_metric)), 0)
  # full pooling means 1's everywhere
  expect_equal(unique(as.numeric(bg5_f$pooling_metric)), 1)

  # pp <- pooling(bg5_p)
  # expect_is(pp, "array")
  # expect_gt(min(pp), 0)
  # expect_lt(max(pp), 1)
  # expect_identical(bg5_p$pooling_metric, pooling(bg5_p))

  # since all SEs are the same, pooling should be the same for all sites
  # capture_output(print(pp))
  # expect_equal(pp[2,,1], .75, tolerance = .1) #YUGE tolerance as we only do 150 iter
  # expect_equal(length(unique(pp[2,,1])), 1)
  # expect_equal(as.numeric(pp[2,1,1]), .75, tolerance = .1)
})

test_that("extra pooling stats work", {
  # Extra pooling checks
  # Calculation of I^2 and H^2
  i2 <- pooling(bg5_p, metric = "isq")
  expect_is(i2, "array")
  expect_gte(min(i2), 0)
  expect_lte(max(i2), 1)
  h2 <- pooling(bg5_p, metric = "hsq")
  expect_is(h2, "array")
  expect_gte(min(h2), 1)
  # Calculation of weights makes sense
  wt <- weights(bg5_p)
  expect_is(wt, "array")
  expect_equal(dim(wt), c(3,K,1))
  expect_equal(sum(wt[2,,1]), 1)
  expect_lte(sum(wt[1,,1]), sum(wt[2,,1]))
  expect_gte(sum(wt[3,,1]), sum(wt[2,,1]))
  expect_gte(sum(wt[1,,1]), 0)
  wt2 <- pooling(bg5_p, metric = "weights")
  expect_identical(wt, wt2)
})

test_that("Calculation of effects works", {
  expect_is(group_effects(bg5_p), "array")
  expect_is(treatment_effect(bg5_p), "list")
  expect_length(treatment_effect(bg5_p, summary = TRUE)$tau, 5)
  expect_length(treatment_effect(bg5_p, summary = TRUE)$sigma_tau, 5)

  expect_identical(dim(group_effects(bg5_n)), as.integer(c(150, K , 1)))
  expect_identical(dim(group_effects(bg5_p)), as.integer(c(150, K , 1)))
  expect_identical(dim(group_effects(bg5_f)), as.integer(c(150, K , 1)))
  expect_identical(names(treatment_effect(bg5_p)), c("tau", "sigma_tau"))
})


test_that("Plotting and printing works", {
  expect_is(plot(bg5_n), "gg")
  expect_is(plot(bg5_p, transform = exp), "gg")
  expect_is(plot(bg5_p, transform = exp, hyper = TRUE), "gg")
  expect_is(plot(bg5_p, hyper = TRUE), "gg")
  expect_is(plot(bg5_p, order = TRUE), "gg")
  expect_is(plot(bg5_f, order = FALSE), "gg")
  skip("Not working for ed")
  # This has been changed in forestplot 2.0:
  expect_is(forest_plot(bg5_n), "gforge_forestplot")
  expect_is(forest_plot(bg5_p), "gforge_forestplot")
  expect_is(forest_plot(bg5_f), "gforge_forestplot")
  expect_is(forest_plot(bg5_f, graph.pos = 1), "gforge_forestplot")
  # but we can crash it easily if
  expect_error(plot(bg5_n, style = "rubbish"), "argument must be one of")

  capture_output(print(bg5_n))
  capture_output(print(bg5_p))
  capture_output(print(bg5_f))
  capture_output(print(bg5_p, exponent = TRUE))
})

# test_that("Test data can be used in the Rubin model", {
#   bg_lpd <- expect_warning(baggr(df_binary[1:6,], test_data = df_binary[7:8,],
#                                  iter = 500, refresh = 0))
#   expect_is(bg_lpd, "baggr")
#   # make sure that we have 6 sites, not 8:
#   expect_equal(dim(group_effects(bg_lpd)), c(1000, 6, 1))
#   # make sure it's not 0
#   expect_equal(mean(rstan::extract(bg_lpd$fit, "logpd[1]")[[1]]), -3.6, tolerance = 1)
#
#   # wrong test_data
#   df_na <- df_binary[7:8,]; df_na$tau <- NULL
#   expect_error(baggr(df_binary[1:6,], test_data = df_na), "must be of the same format as input")
# })


# test helpers -----

test_that("Extracting treatment/study effects works", {
  expect_error(treatment_effect(df_surv))
  expect_is(treatment_effect(bg5_p), "list")
  expect_identical(names(treatment_effect(bg5_p)), c("tau", "sigma_tau"))
  expect_is(treatment_effect(bg5_p)$tau, "numeric") #this might change to accommodate more dim's
  expect_message(treatment_effect(bg5_n), "no treatment effect estimated when")

  # Drawing values of tau:
  expect_error(effect_draw(cars))
  expect_is(effect_draw(bg5_p), "numeric")
  expect_is(effect_draw(bg5_p, transform = exp), "numeric")
  expect_length(effect_draw(bg5_p), 150)
  expect_length(effect_draw(bg5_p,7), 7)

  # Plotting tau:
  expect_is(effect_plot(bg5_p), "gg")
  expect_is(effect_plot(bg5_p, bg5_f), "gg")
  expect_is(effect_plot("Model A" = bg5_p, "Model B" = bg5_f), "gg")
  # Crashes when passing nonsense
  expect_error(effect_plot(cars), "baggr class")
  expect_error(effect_plot(cars, cars, bg5_f), "baggr class")
})



# covariates ------
sa <- df_surv
sa$a <- rnorm(nrow(df_surv))
sa$b <- rnorm(nrow(df_surv))
sb <- sa
sb$b <- NULL

library(baggr)
test_that("Model with covariates works fine", {
  bg_cov <- expect_warning(
    baggr(sa, 
    outcome = c("interval_left", "interval_right"),
    group = "site",
    covariates = c("a", "b"), iter = 150, chains = 1, refresh = 0))

  expect_is(bg_cov, "baggr")
  expect_error(baggr(sa,
    outcome = c("interval_left", "interval_right"),
    group = "site",
   covariates = c("made_up_covariates")), "are not columns")
  expect_error(baggr(sa, 
    outcome = c("interval_left", "interval_right"),
    group = "site",
    covariates = c("a", "b", "made_up_covariates")))


  expect_length(bg5_p$covariates, 0)
  expect_length(bg_cov$covariates, 2)
  expect_null(bg_cov$mean_lpd)

  # Fixed effects extraction
  expect_is(fixed_effects(bg_cov), "matrix")
  expect_is(fixed_effects(bg_cov, transform = exp), "matrix")
  expect_equal(dim(fixed_effects(bg_cov, summary = TRUE)), c(2,5,1))
  expect_equal(dim(fixed_effects(bg_cov, summary = FALSE))[2], 2)
})



### TODO: Generate Quantities for loo etc.

# tests for helper functions -----


# test_that("baggr_compare basic cases work with logit models", {
#   # try to make nonexistant comparison:
#   expect_error(baggr_compare(bg5_p, bg5_n, bg5_f, compare = "sreffects"))
#   # Compare existing models:
#   expect_is(plot(baggr_compare(bg5_p, bg5_n, bg5_f)), "gg")
# })

# test_that("loocv", {
#   # Rubbish model
#   # expect_error(loocv(schools, model = "mutau"))
#   # Can't do pooling none
#   expect_error(loocv(df_surv, pooling = "none"))

#   skip_on_cran()
#   skip_on_travis()

#   loo_model <- expect_warning(loocv(df_binary, model = "logit",
#                                     return_models = TRUE, iter = 150, chains = 1, refresh = 0))
#   expect_is(loo_model, "baggr_cv")
#   capture_output(print(loo_model))

#   loo_full <- expect_warning(loocv(df_binary, model = "logit", pooling = "full",
#                                     return_models = TRUE, iter = 150, chains = 1, refresh = 0))
#   expect_is(loo_full, "baggr_cv")
#   capture_output(print(loo_full))
# })

# comp_pl <- expect_warning(baggr_compare(
#   df_binary, model = "logit", iter = 150, what = "pooling"
# ))

# comp_pr <- expect_warning(baggr_compare(
#   df_binary, model = "logit", iter = 150, what = "prior"
# ))

# comp_existmodels <- baggr_compare(bg5_p, bg5_f)

# test_that("baggr comparison method works for Full model", {

#   expect_is(comp_pl, "baggr_compare")
#   expect_is(comp_pr, "baggr_compare")
#   expect_is(comp_existmodels, "baggr_compare")

#   expect_is(testthat::capture_output(print(comp_pl)), "character")
#   expect_is(testthat::capture_output(print(comp_pr)), "character")
#   expect_is(testthat::capture_output(print(comp_existmodels)), "character")

#   expect_gt(length(comp_pl), 0)
#   expect_gt(length(comp_pr), 0)
#   expect_gt(length(comp_existmodels), 0)

#   expect_is(plot(comp_pl), "gg")
#   expect_is(plot(comp_pl, grid_models = TRUE), "gtable")
#   expect_is(plot(comp_pr), "gg")
#   expect_is(plot(comp_pr, grid_models = TRUE), "gtable")
#   expect_is(plot(comp_existmodels), "gg")
#   expect_is(plot(comp_existmodels, grid_models = TRUE), "gtable")

# })



# # Setting control pooling, control priors -----

# test_that("Prior specifications for baselines work", {

#   skip_on_cran()

#   bg1 <- expect_warning(baggr(df_binary, "logit", pooling = "none",
#                               pooling_control = "partial",
#                               iter = 150, chains = 2, refresh = 0,
#                               show_messages = F))
#   bg2 <- expect_warning(baggr(df_binary, "logit", pooling = "none",
#                               pooling_control = "partial", prior_control = normal(0, 5),
#                               iter = 150, chains = 2, refresh = 0,
#                               show_messages = F))
#   bg3 <- expect_warning(baggr(df_binary, "logit", pooling = "none",
#                               pooling_control = "partial", prior_control = normal(0, 5), prior_control_sd = uniform(0, 1),
#                               iter = 150, chains = 2, refresh = 0,
#                               show_messages = F))

#   expect_is(bg1, "baggr")
#   expect_is(bg2, "baggr")
#   expect_is(bg3, "baggr")

#   expect_error(baggr(df_binary, "logit", pooling = "none",
#                               pooling_control = "partial", prior_control = multinormal(c(0,0), diag(2)),
#                               iter = 150, chains = 2, refresh = 0,
#                               show_messages = F))
# })

