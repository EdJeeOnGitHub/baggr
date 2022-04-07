library(tidyverse) # data manipulation
library(tidybayes) # data manipulation... but Bayesian
library(data.table) # data manipulation... but fast
library(furrr) # like purrr but in parallel thanks to future
library(cmdstanr) # stan


################################################################################
##SBC Functions woooooooooooooo#################################################
################################################################################

#' Generate "Data" for simulations
#' 
#' Here data really means hyperparameters such as number of observations, 
#' number of covariates, and who is censored (any data that's not 'directly' 
#' modelled) 
#' 
#' #TODO: Create a meta gen_data function to vary hyper params 
meta_gen_data_function = function(N = 1000,
                                  nc = 3,
                                  nsc = 3,
                                  beta_mean = rep(0, nc),
                                  beta_sigma = rep(1, nc)){
    function(seed) {
        set.seed(seed + 1e6)
        N = N
        nc = nc
        nsc = nsc 
        beta_mean = beta_mean
        beta_sigma = beta_sigma

        censoring = rbernoulli(N, p = 0.5)

        return(list(
            N = N,
            nc = nc,
            nsc = nsc,
            beta_mean = beta_mean,
            beta_sigma = beta_sigma,
            censoring = censoring
        ))
    }


}

#' Generate parameters 
#' 
#' These are the parameters we want to estimate and correspond to the models
#' `parameters` block
gen_params <- function(seed,
                       data){
  set.seed(seed + 2e6)
    beta = rnorm(data$nc, data$beta_mean, data$beta_sigma)

    lshape = rnorm(1, 2, 0.5)
    shape = exp(lshape)
    beta_transformed = -1*beta/shape
  
  return(list(
      beta = beta,
      shape = shape,
      lshape = lshape,
      beta_transformed = beta_transformed
  ))
  
}
# N = 1000
# lshape = rnorm(N, 2, 0.5)
# shape = exp(lshape)
# beta = rnorm(N, 0, 1)
# bt = -1*beta/shape
# rweibull(N, shape, exp(bt)) %>% hist()


#' Generate Modelled Data
#'
#'
#' Using hyper parameters from `gen_data` and parameters from `gen_params`,
#' create our modelled data. That is, create left and interval censored 
#' survival data according to the generative model. 
gen_modeled_data <- function(seed,
                             data,
                             params) {
  set.seed(seed + 3e6)
    

    # X indep normals 0, 1 
    X = matrix(
        rnorm(data$N*data$nc),
        nrow = data$N,
        ncol = data$nc
    )
    
    scale = exp(X %*% params$beta_transformed)
    censoring = data$censoring
    N = data$N
    # Preallocate these to fill, should really just use apply or something
    #  but this matches stan file which makes debugging easier
    interval_left = matrix(data = NA, 
                           nrow = N,
                           ncol = 1)
    interval_right = matrix(data = NA, 
                            nrow = N,
                            ncol = 1)
   for (i in 1:N) {
       if (censoring[i] == 1) {
           interval_left[i, 1] = rweibull(1, params$shape, scale[i]) 
            # maybe should set this to -Inf to ensure likelihood blows up if 
            # we've incorrectly indexed censoring[i]
           interval_right[i, 1] = Inf 
       } 
       else {
           #! Not sure if this is correct
           # Just draw from two weibulls and take smallest as left censor 
           # largest as right censor
           interval = sort(rweibull(2, params$shape, scale[i]))
           interval_left[i, 1] = interval[1]
           interval_right[i, 1] = interval[2]
       } 
   } 



  return(list(
      interval_left = interval_left[, 1],
      interval_right = interval_right[, 1],
      X = X
  ))
}


#' Sampling wrapper function
#' 
#' N.B. cmdstanr syntax `modelobject$sample()`
sample_from_model <- function(seed,
                              data,
                              params,
                              modeled_data,
                              iters) {
    data_for_stan <- c(data, modeled_data)
    fit = surv_model$sample(
                  data = data_for_stan,
                  seed = seed,
                  chains = 1,
                  refresh = 0
    )
    return(fit)
}


#' Create a 'comparison df'
#' 
#' We want to compare our \theta_0 to \hat{\theta} draws and check if model 
#' can recover it's own parameters.
#' 
#' Therefore, we extract model draws and tidy them up. Then we extract 
#' generative prior draw and tidy up. Finally, we left_join
create_comp_df = function(model_fit, params){

    draw_df = gather_draws(model_fit, beta[j], shape) %>%
        to_broom_names() %>%
        mutate(j = replace_na(j, 1))

    param_df = bind_rows(
        enframe(params$beta, name = "j") %>% mutate(term = "beta"),
        enframe(params$shape)  %>% mutate(term = "shape") %>% mutate(j = 1)
    ) %>% rename(prior_draw = value) %>% select(-name)

    comp_df = left_join(
        draw_df, 
        param_df,
        by = c("term", "j")
    )

    if (nrow(comp_df) != nrow(draw_df)) {
        stop("`create_comp_df()` left_join failed.")
    }

    return(comp_df)

}


#' Create Rank Statistic DF
#' 
#' Calculate how often \hat{\theta} < \theta_0 across posterior draws
#' 
#' N.B. This is hardcoded to index \beta[j] at the moment and will have to 
#' be adjusted if we add e.g. \beta[j, t] or some multiple index covariance 
#' matrix
create_rank_stat = function(comp_df) {
    
    rank_df = comp_df %>%
        group_by(term, j) %>%
        summarise(rank_stat = mean(estimate < prior_draw))
    return(rank_df)

}

#' Extract N Effective Draws from cmdstanr fitted object
#' 
#'  Creates a table with each parameters number of effective draws.
#' 
#' I am not a smart man - this is not a good implementation.
#' 
#' For whatever reason, the R side summary functions don't seem to report N_eff
#' but the cmdstan side does. However, N_eff is only passed to us through 
#' `stdout` so I wrote this horrific code to extract the string as a table 
#' and then manipulate this string table into something resembling data. 
#' 
#' An alternative would be to calculate N_eff directly (1 /sum(\hat{\rho}_j)) 
#' or whatever or figure out where on earth this is reported (it must be 
#' somewhere?).
#' #TODO: Fix this bigly edward.
extract_N_eff_badly = function(model_fit){
    # This always prints to stdout even with `invisible` :(
    mod_summary = model_fit$cmdstan_summary()

    # convert to string table :) ... :(
    string_table = read.table(
        textConnection(mod_summary$stdout),
        sep = "\t",
        skip = 5
    ) %>% as_tibble()

    clean_N_eff = string_table %>% 
        # let's hope there's only ever 9 useless params reported at the fron
        slice(9:(nrow(string_table) - 4)) %>% 
        # split if we see a space followed by a digit, word, +, or -
        separate(V1, sep = "\\s(?=(\\w+|\\d+|[\\+\\-]))", into = paste0("col_", 1:10)) %>%
        select(col_1, col_8) %>%
        # clean up the columnns we've grabbed
        mutate(N_eff = as.numeric(col_8), term = str_remove_all(col_1, " ")) %>%
        select(-col_8, -col_1) %>%
        mutate(j = str_extract(term, "\\d+"), j = as.numeric(replace_na(j, 1))) %>%
        mutate(term = str_remove(term, "\\[.*$"))
    return(clean_N_eff)
}
#' Simulation Based Calibration Draw
#' 
#' Takes in a seed and spits out one simulated draw.
#' 
#'  SBC samples from prior, uses model to create data, fits model to data 
#' and compares posterior to sampled prior.
#' 
#' Option to `thin` or not because of the horrific N_eff extraction.
#' If something is breaking turning thin off will be a good idea.
#' 
#'  THINNING
#' Why are we thinning? Because markov chains tend to be auto correlated and 
#' so we want to thin so we're comparing independent quantiles so our rank 
#' stat is well behaved.
#' 
#' We do that by essentially taking N/N_eff every other draw.
#' 
#' Unfortunately NUTS is so good we actually get N_eff > N often because 
#' we find a great spot (thanks HMC) and then NUTS (No U-Turn Sampler)
#'  swings us in the other direction and we end up with negatively 
#' autocorrelated draws => N_eff > N.
#' 
#' Gelman and co/others keep running till N_eff < N and thin every other draw 
#' anyway but I have opted to skip this step. Possible #TODO
sbc_draw = function(seed = 1234, thin = TRUE, ...){
    gen_data = meta_gen_data_function(...)
    data = gen_data(seed)
    params = gen_params(seed,
                        data)
    modeled_data = gen_modeled_data(seed,
                                    data,
                                    params)
    fit = sample_from_model(seed,
                            data,
                            params,
                            modeled_data)
    comp_df = create_comp_df(fit, params)
    


    # thinning code
    if (thin == TRUE) {
        # Only extract parameters we care about, if we want to ignore a parameter
        # give it an underscore name or call it TRANSformed_variable.
        N_eff_df = extract_N_eff_badly(
            fit
            ) %>%
            filter(!str_detect(term, "trans|_"))


        comp_df = comp_df %>%
            left_join(
                N_eff_df %>%
                    mutate(N = fit$metadata()$iter_sampling, 
                        thin = ceiling(N/N_eff)),
                by = c("j", "term")
            ) %>%
            mutate(thin_draw = .iteration %% thin )

    } else {
        comp_df$thin_draw = 0
    }


    rank_stat = comp_df %>%
        filter(thin_draw == 0) %>%
        create_rank_stat()


  return(rank_stat)
}






################################################################################
###Simulation Time!!!!!!!!!!!###################################################
################################################################################


# Load our vanilla survival model
surv_model = cmdstan_model("inst/stan/for_ed.stan")
plan(multicore, workers = 8)

draws <- 1:200 %>%
  future_map_dfr(~sbc_draw(.x) %>% mutate(draw = .x),
                 .options = furrr_options(
                     seed = TRUE,
                     scheduling = 2L # This means we move jobs dynamically
                 ),
                 .progress = TRUE)  

draws %>% 
    ggplot(aes(sample = rank_stat,
           colour = term)) +
  stat_qq(distribution = stats::qunif,
          alpha = 1) +
  stat_qq_line(distribution = stats::qunif,
               colour = "black",linetype = "longdash") +
  facet_wrap(term~j ,
             scales = "free") +
  theme_bw() +
  guides(colour = "none") +
  labs(
      title = "Weibull Survival Model Simulation Based Calibration",
      subtitle = "Left and Interval Censoring", 
      caption = "Uniform CDF indicates well calibrated posterior coverage.",
      x = "Theoretical Quantile", 
      y = "Realised Quantile"
  )

ggsave(
    "data/plots/weibull-sbc.png",
    width = 8,
    height = 6
)


################ Scratchpad Stuff ##############################################
# Delete for future use #
######## Tests ##########  
# seed = 1
# seed = seed + 1
# sbc_draw(seed)
# # print(seed)
# # # seed = 1
# data <- gen_data(seed)
# params <- gen_params(seed,
#                     data)
# modelled_data <- gen_modeled_data(seed,
#                                 data,
#                                 params)
# dt = data.table(
#     il = modelled_data$interval_left,
#     ir = modelled_data$interval_right,
#     censored = data$censoring
# )
# dt[, c("x_1", "x_2") :=  as.data.table(modelled_data$X)]

# dt
# dt[censored == FALSE] %>%
#     ggplot(aes( 
#         x = ir
#     )) +
#     geom_histogram()
# params
# rweibull(100, params$shape, 100)
# data_for_stan <- c(data, modelled_data)

#     fit = surv_model$sample(
#                   data = data_for_stan,
#                   seed = seed,
#                   chains = 1
#     )



# ed = extract_N_eff_badly(fit)
# ed = fit$cmdstan_summary()
# tab = read.table(textConnection(ed$stdout), sep = "\t", skip = 5) %>% as_tibble()
# tab
# tab %>%
# tab


# # assuming 'fit' is from CmdStanR
# stanfit <- rstan::read_stan_csv(fit$output_files())
# launch_shinystan(stanfit)
# ed %>%
#     gather_draws(beta[j], shape[j])
# comp_df %>%
#     group_by(term, j) %>%
#     summarise(rank_stat = mean(estimate < prior_draw))

# comp_df = create_rank_df(ed, params)


# param_df 
# param_df

# param_df = map(params, ~t(tibble(paste0(names(.x), ))))


# tibble(c(1, 2, 3), t(params$beta))
# enframe(params$beta)
