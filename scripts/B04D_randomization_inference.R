# # # Emma Krasovich Southworth, emmars@stanford.edu, Jeff Wen
# # # Last Updated: April 19, 2024
# # # Description: run regression on burned structures and smoke speciation + plot response curves at different levels of smoke
#
# loadd(c(burned_struc_smoke_spec_df, parameter_categories), cache = drake_cache)

# ----------------------------------------------------------------
#  RUN PERMUTATION TEST TO TEST SHARP NULL
# ----------------------------------------------------------------

randomization_inference_permutation_test <- function(burned_struc_smoke_spec_df, parameter_categories) {
# permutation/ randomization inference ----
## sampling within state, across years
withinAcrossResampler <- function(data, within_condition, across_condition, interest_var){

  data_dt <- data %>% as.data.table()

  # iterate through within condition
  outer_dt <- lapply(data[[within_condition]] %>% unique(), function(within_unit){

    # get list of possible values of across condition
    across_condition_list <-  data_dt[get(within_condition)==within_unit][[across_condition]] %>%
      unique()

    across_condition_dt <- lapply(across_condition_list, function(across_unit){

      # number to sample needs to match how many obs we are trying to replace
      inner_dt <- data_dt[get(within_condition)==within_unit &
                            get(across_condition)==across_unit]

      num_samples <- nrow(inner_dt)

      # filter for within unit but not across unit, then sample without replacement
      # from set of possible options
      possible_new_values <- data_dt[get(within_condition)==within_unit &
                                       get(across_condition)!=across_unit
      ][[interest_var]]

      # get num_samples from possible_new_values
      new_values <- possible_new_values[sample.int(length(possible_new_values), 
                                                   size=num_samples)]

      inner_dt[[stringr::str_glue('{interest_var}_new')]] <- new_values

      return(inner_dt)
    }) %>% bind_rows()

    return(across_condition_dt)
  }) %>% bind_rows()

  return(outer_dt)
}


permTest <- function(i, outcome, rhs, data, treatment, interest, cluster_err, across_condition, sample_group=NULL){

  ## sample and create new data with shuffled treatment

  if(!is.null(sample_group)){
    new_data <- withinAcrossResampler(data, 
                                      within_condition=sample_group, 
                                      across_condition=across_condition,
                                      interest_var=treatment)

  } else{
    new_data <- data %>%
      mutate(!!paste0(treatment,"_new") := sample(get(c(treatment)), size=n())) %>%
      ungroup()
  }


  fixest_df <- lapply(outcome, function(x){
    fmla <- paste0(x,"~",rhs)
    fixe_est <- fixest::feols(as.formula(stringr::str_replace_all(fmla, treatment, paste0(treatment,"_new"))),
                              data=new_data)

    ## no clustering SEs here because we are just interested in the coef estimates
    data.frame("coef"=c(summary(fixe_est)$coefficients[[paste0(interest,"_new")]]),
               "species"=c(as.character(fixe_est$fml)[2]))
  }) %>%
    bind_rows()

  return(fixest_df)
}

TREATMENT_VAR <- "contrib_daily_structures_destroyed"
INTEREST_VAR <- "smokePM:contrib_daily_structures_destroyed"
LINE_COLOR <- "firebrick"
DEFAULT_COLOR <- "gray60"

# set up regression dataframe
reg_df <- burned_struc_smoke_spec_df %>%
  filter(contrib_daily_structures_destroyed > 0) %>%
  #filter(smokePM > 0) %>%
  mutate(monitor_month = paste0(site_id, "-", month))

# Run a regression using your sample
sample_mod = feols(c(TI, MN, PB, CU, NI, ZN, AS, CR, MG)
                   ~ smokePM*contrib_daily_structures_destroyed |
                     monitor_month + year,
                   reg_df, cluster = 'site_id')
# get coeffs
sample_coeffs <- coeftable(sample_mod) %>%
  rename(pval = 'Pr(>|t|)',
         se = 'Std. Error',
         species = 'lhs') %>%
  mutate(pval = round(pval, digits = 3)) %>%
  dplyr::select(-id)  %>%
  left_join(confint(
    sample_mod) %>%
      rename(species = 'lhs',
             CI25 = '2.5 %',
             CI975 = '97.5 %') %>%
      dplyr::select(-id))

# parellelize
future::plan(multisession)

start_time <- Sys.time()
# set seed
set.seed(123)
result_df <- lapply(1:1000, permTest,
                    outcome=c('TI', 'MN','PB', 'CU', 'NI', 'ZN', 'AS', 'CR', 'MG'),
                    rhs="smokePM:contrib_daily_structures_destroyed | monitor_month + year",
                    data=reg_df,
                    treatment=TREATMENT_VAR,
                    interest=INTEREST_VAR,
                    across_condition="year",
                    sample_group="state_name") %>%
  bind_rows()

perm_result_df <- result_df %>%
  left_join(sample_coeffs %>% 
              filter(coefficient == 'smokePM:contrib_daily_structures_destroyed') %>% 
              dplyr::select(species, obs_est = 'Estimate', pval))

rm(sample_mod)

end_time <- Sys.time()
future::plan(NULL)


return(perm_result_df)
}
