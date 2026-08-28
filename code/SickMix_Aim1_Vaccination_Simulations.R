################################################################################
#                                                                              #
#                     PROJECT: NSF IHBEM SickMix Aim 1                         #
#                     CODED BY: Jessica Ibiebele                               #
#                     E-MAIL: jibiebel@bu.edu                                  #
#                                                                              #
################################################################################

#------------------------------------------------------------------------------#
# Description
#------------------------------------------------------------------------------#

# This code script simulates vaccination-based public health interventions using
# behavior-dynamic and behavior-naive models.  

#------------------------------------------------------------------------------#
# Set up
#------------------------------------------------------------------------------#

# Set working directory
setwd('')

# Load packages
library(readxl)
library(tidyverse)
library(ggplot2)
library(deSolve)
library(dplyr)

# Load data 
sim_population <- read_xlsx("sim_population.xlsx")
contact_matrix_acute <- as.matrix(read_xlsx('contact_matrix_acute.xlsx'))
contact_matrix_late <- as.matrix(read_xlsx('contact_matrix_late.xlsx'))
contact_matrix_naive <- as.matrix(read_xlsx('contact_matrix_naive.xlsx'))
contact_matrix_home <- as.matrix(read_xlsx('contact_matrix_home.xlsx'))

# Data manipulation
## Add rownames to all matrices, if you haven't already 
rownames(contact_matrix_acute) <- c("<1", "1–4", "5–17", "18–29", "30–44", "45–64", "65+")
rownames(contact_matrix_late) <- c("<1", "1–4", "5–17", "18–29", "30–44", "45–64", "65+")
rownames(contact_matrix_naive) <- c("<1", "1–4", "5–17", "18–29", "30–44", "45–64", "65+")
rownames(contact_matrix_home) <- c("<1", "1–4", "5–17", "18–29", "30–44", "45–64", "65+")

# ******************* 50% of age group <1 is vaccinated ******************* #

# Assign initial fraction to enter into susceptible_vaccinated compartment. 
# Reduce beta by 50% and lambda by 50% for those who are vaccinated. 
# Vaccinated individuals can move through Ia_vaccinated, Il_vaccinated, and R_vaccinated.  

#------------------------------------------------------------------------------#
#                       SIR Model - Behavior-dynamic model 
#------------------------------------------------------------------------------#

SIR_dynamic_9 <- function (t, x, params) {
  # There are 7 compartments in each state to reflect the 7 age groups   
  S=x[1:7] # Susceptible 
  Ia=x[8:14] # Infectious (acute) (first 3 days of illness)
  Il=x[15:21] # Infectious (late) (next 7 days)
  R=x[22:28] # Recovered
  
  # Add vaccinated compartments 
  S_vax=x[29:35] # Susceptible_vaccinated 
  Ia_vax=x[36:42] # Infectious (acute) (first 3 days of illness)
  Il_vax=x[43:49] # Infectious (late) (next 7 days)
  R_vax=x[50:56] # Recovered
  
  beta  = params$beta # transmission rate (same for acute and late phases)
  beta_vax  = params$beta_vax # transmission rate (vaccinated) 
  gamma_a = params$gamma_a # acute infectious period = first 3 days
  gamma_l = params$gamma_l # late infectious period = next 7 days
  C_acute = params$C_acute # contact matrix for acute period
  C_late = params$C_late  # contact matrix for late period 
  
  # Population size, N
  N = S + Ia + Il + R + S_vax + Ia_vax + Il_vax + R_vax
  
  # Force of infection 
  ## Accounts for reduced transmissibility from vaccinated individuals and 
  ## reduced susceptibility of vaccinated individuals 
  infectious_contacts_unvax <- 
    t(C_acute) %*% Ia +
    t(C_late)  %*% Il
  
  infectious_contacts_vax <- 
    t(C_acute) %*% Ia_vax +
    t(C_late)  %*% Il_vax
  
  foi_sum_terms_unvax <- as.numeric(infectious_contacts_unvax) / N
  foi_sum_terms_vax   <- as.numeric(infectious_contacts_vax) / N
  
  lambda <- beta * foi_sum_terms_unvax + beta_vax * foi_sum_terms_vax 
  # This lambda is applied to all unvaccinated individuals and accounts for reduced
  # transmissibility of vaccinated individuals 
  lambda_vax <- 0.5*lambda 
  # This lambda is specific to vaccinated individuals and additionally accounts 
  # for reduced susceptibility due to vaccination 
  
  # Create a set of ODEs for the SIR model
  dxdt <- c(
    # Unvaccinated 
    - lambda * S, # S compartment
    lambda * S - gamma_a * Ia, # S -> Ia compartment
    gamma_a * Ia - gamma_l * Il, # Ia -> Il compartment
    gamma_l * Il, # Il -> R compartment
    
    # Vaccinated 
    ## Includes lambda_vax instead of lambda
    -lambda_vax * S_vax,
    lambda_vax * S_vax - gamma_a * Ia_vax,
    gamma_a * Ia_vax - gamma_l * Il_vax,
    gamma_l * Il_vax
  ) 
  
  return(list(dxdt))
}

# Set parameter values
params <- list(
  beta = 0.01,     # transmission rate
  beta_vax = 0.005, # vax transmission rate
  gamma_a = 1/3,     # acute infectious period = first 3 days 
  gamma_l = 1/7,    # late infectious period = next 7 days
  C_acute = contact_matrix_acute, # contact matrix for acute period
  C_late  = contact_matrix_late # contact matrix for late period  
)

# Make sure contact matrices are stored as numeric matrices 
params$C_acute <- as.matrix(contact_matrix_acute)
params$C_late  <- as.matrix(contact_matrix_late)

storage.mode(params$C_acute) <- "double"
storage.mode(params$C_late)  <- "double"

# Incorporate age structure
age_group_nums <- with(sim_population, setNames(age_group_total, age_group))

# Set initial conditions using population data
initial_infec_fraction <- 1e-5

## Assign the same proportion of infected individuals to each age group 

# Split population
S0_total <- age_group_nums * (1 - initial_infec_fraction)
Ia0_total <- age_group_nums * initial_infec_fraction

# 50% vaccinated
vax_coverage <- c(0.5, 0, 0, 0, 0, 0, 0) # only 50% of age group <1 is vaccinated
# Apply vaccination by age group
S0_vax <- vax_coverage * S0_total
S0     <- (1 - vax_coverage) * S0_total

Ia0_vax <- vax_coverage * Ia0_total
Ia0     <- (1 - vax_coverage) * Ia0_total

Il0      <- rep(0, 7)
Il0_vax  <- rep(0, 7)

R0       <- rep(0, 7)
R0_vax   <- rep(0, 7)

# Combine all compartments
initial_state <- c(
  S0, Ia0, Il0, R0,
  S0_vax, Ia0_vax, Il0_vax, R0_vax
)

# Check initial values 
sum(c(S0, Ia0, Il0, R0))
sum(sim_population$age_group_total)
sum(Ia0)

# Set time steps
times <- seq(0, 1000, by = 1)

# Run ODE (deterministic)
out_ODE_9 <- ode(
  y = initial_state,
  times = times,
  func = SIR_dynamic_9,
  parms = params
)

out_ODE_9 <- as.data.frame(out_ODE_9)

# Check output
head(out_ODE_9)
tail(out_ODE_9)

final_9 <- tail(out_ODE_9, 1)
R_unvax_9 <- as.numeric(final_9[23:29])
R_vax_9   <- as.numeric(final_9[51:57])
total_final_cases_9 <- sum(R_unvax_9 + R_vax_9)
print(total_final_cases_9)

# cases per 100,000
total_final_cases_9/(sum(initial_state))*100000

# ---- Extract per-age-group final cases from ode() output ----

# -----------------------------
# Age group labels and populations
# -----------------------------
age_groups <- c("<1", "1-4", "5-17", "18-29", "30-44", "45-64", "65+")
N_age <- age_group_nums  # vector of population per age group (same order as model)

get_age_group_cases <- function(out_ODE_9, age_groups) {
  final_row <- tail(out_ODE_9, 1)
  
  # +1 offset everywhere because column 1 of out_ODE is "time"
  # Unvaccinated
  S_unvax_final  <- as.numeric(final_row[1 + (1:7)])
  Ia_unvax_final <- as.numeric(final_row[1 + (8:14)])
  Il_unvax_final <- as.numeric(final_row[1 + (15:21)])
  R_unvax_final  <- as.numeric(final_row[1 + (22:28)])
  
  # Vaccinated
  S_vax_final  <- as.numeric(final_row[1 + (29:35)])
  Ia_vax_final <- as.numeric(final_row[1 + (36:42)])
  Il_vax_final <- as.numeric(final_row[1 + (43:49)])
  R_vax_final  <- as.numeric(final_row[1 + (50:56)])
  
  # Combine vax + unvax per age group
  S_total  <- S_unvax_final  + S_vax_final
  Ia_total <- Ia_unvax_final + Ia_vax_final
  Il_total <- Il_unvax_final + Il_vax_final
  R_total  <- R_unvax_final  + R_vax_final
  
  N <- S_total + Ia_total + Il_total + R_total
  
  # Total cases = everyone no longer susceptible (captures R plus any
  # still infectious at the end of the simulation)
  cases <- N - S_total
  
  data.frame(
    age_group = age_groups,
    N = N,
    S_final = S_total,
    Ia_final = Ia_total,
    Il_final = Il_total,
    R_final = R_total,
    cases = cases
  )
}

# ---- Apply to existing run ----
cases_by_age_9 <- get_age_group_cases(out_ODE_9, age_groups)
print(cases_by_age_9)


#------------------------------------------------------------------------------#
#                       SIR Model - Behavior-naive model 
#------------------------------------------------------------------------------#

SIR_naive_10 <- function (t, x, params) {
  # There are 7 compartments in each state to reflect the 7 age groups   
  S=x[1:7] # Susceptible 
  Ia=x[8:14] # Infectious (acute) (first 3 days of illness)
  Il=x[15:21] # Infectious (late) (next 7 days)
  R=x[22:28] # Recovered
  
  # Add vaccinated compartments 
  S_vax=x[29:35] # Susceptible_vaccinated 
  Ia_vax=x[36:42] # Infectious (acute) (first 3 days of illness)
  Il_vax=x[43:49] # Infectious (late) (next 7 days)
  R_vax=x[50:56] # Recovered
  
  beta  = params$beta # transmission rate (same for acute and late phases)
  beta_vax  = params$beta_vax # transmission rate (vaccinated) 
  gamma_a = params$gamma_a # acute infectious period = first 3 days
  gamma_l = params$gamma_l # late infectious period = next 7 days
  C_naive = params$C_naive # contact matrix for acute & late periods
  
  # Population size, N
  N = S + Ia + Il + R + S_vax + Ia_vax + Il_vax + R_vax
  
  # Force of infection 
  ## Accounts for reduced transmissibility from vaccinated individuals and 
  ## reduced susceptibility of vaccinated individuals 
  infectious_contacts_unvax <- 
    t(C_naive) %*% Ia +
    t(C_naive)  %*% Il
  
  infectious_contacts_vax <- 
    t(C_naive) %*% Ia_vax +
    t(C_naive)  %*% Il_vax
  
  foi_sum_terms_unvax <- as.numeric(infectious_contacts_unvax) / N
  foi_sum_terms_vax   <- as.numeric(infectious_contacts_vax) / N
  
  lambda <- beta * foi_sum_terms_unvax + beta_vax * foi_sum_terms_vax 
  # This lambda is applied to all unvaccinated individuals and accounts for reduced
  # transmissibility of vaccinated individuals 
  lambda_vax <- 0.5*lambda 
  # This lambda is specific to vaccinated individuals and additionally accounts 
  # for reduced susceptibility due to vaccination 
  
  # Create a set of ODEs for the SIR model
  dxdt <- c(
    # Unvaccinated 
    - lambda * S, # S compartment
    lambda * S - gamma_a * Ia, # S -> Ia compartment
    gamma_a * Ia - gamma_l * Il, # Ia -> Il compartment
    gamma_l * Il, # Il -> R compartment
    
    # Vaccinated 
    ## Includes lambda_vax instead of lambda
    -lambda_vax * S_vax,
    lambda_vax * S_vax - gamma_a * Ia_vax,
    gamma_a * Ia_vax - gamma_l * Il_vax,
    gamma_l * Il_vax
  ) 
  
  return(list(dxdt))
}

# Set parameter values
params <- list(
  beta = 0.01,     # transmission rate
  beta_vax = 0.005, # vax transmission rate
  gamma_a = 1/3,     # acute infectious period = first 3 days 
  gamma_l = 1/7,    # late infectious period = next 7 days
  C_naive = contact_matrix_naive # contact matrix for acute & late periods
)

# Make sure contact matrices are stored as numeric matrices 
params$C_naive <- as.matrix(contact_matrix_naive)

storage.mode(params$C_naive) <- "double"

# Incorporate age structure
age_group_nums <- with(sim_population, setNames(age_group_total, age_group))

# Set initial conditions using population data
initial_infec_fraction <- 1e-5

## Assign the same proportion of infected individuals to each age group 

# Split population
S0_total <- age_group_nums * (1 - initial_infec_fraction)
Ia0_total <- age_group_nums * initial_infec_fraction

# 50% vaccinated
vax_coverage <- c(0.5, 0, 0, 0, 0, 0, 0) # only 50% of age group <1 is vaccinated
# Apply vaccination by age group
S0_vax <- vax_coverage * S0_total
S0     <- (1 - vax_coverage) * S0_total

Ia0_vax <- vax_coverage * Ia0_total
Ia0     <- (1 - vax_coverage) * Ia0_total

Il0      <- rep(0, 7)
Il0_vax  <- rep(0, 7)

R0       <- rep(0, 7)
R0_vax   <- rep(0, 7)

# Combine all compartments
initial_state <- c(
  S0, Ia0, Il0, R0,
  S0_vax, Ia0_vax, Il0_vax, R0_vax
)

# Check initial values 
sum(c(S0, Ia0, Il0, R0))
sum(sim_population$age_group_total)
sum(Ia0)

# Set time steps
times <- seq(0, 1000, by = 1)

# Run ODE (deterministic)
out_ODE_10 <- ode(
  y = initial_state,
  times = times,
  func = SIR_naive_10,
  parms = params
)

out_ODE_10 <- as.data.frame(out_ODE_10)

# Check output
head(out_ODE_10)
tail(out_ODE_10)

final_10 <- tail(out_ODE_10, 1)
R_unvax_10 <- as.numeric(final_10[23:29])
R_vax_10   <- as.numeric(final_10[51:57])
total_final_cases_10 <- sum(R_unvax_10 + R_vax_10)
print(total_final_cases_10)

# cases per 100,000
total_final_cases_10/(sum(initial_state))*100000

# ---- Extract per-age-group final cases from an ode() output ----

get_age_group_cases <- function(out_ODE_10, age_groups) {
  final_row <- tail(out_ODE_10, 1)
  
  # +1 offset everywhere because column 1 of out_ODE is "time"
  # Unvaccinated
  S_unvax_final  <- as.numeric(final_row[1 + (1:7)])
  Ia_unvax_final <- as.numeric(final_row[1 + (8:14)])
  Il_unvax_final <- as.numeric(final_row[1 + (15:21)])
  R_unvax_final  <- as.numeric(final_row[1 + (22:28)])
  
  # Vaccinated
  S_vax_final  <- as.numeric(final_row[1 + (29:35)])
  Ia_vax_final <- as.numeric(final_row[1 + (36:42)])
  Il_vax_final <- as.numeric(final_row[1 + (43:49)])
  R_vax_final  <- as.numeric(final_row[1 + (50:56)])
  
  # Combine vax + unvax per age group
  S_total  <- S_unvax_final  + S_vax_final
  Ia_total <- Ia_unvax_final + Ia_vax_final
  Il_total <- Il_unvax_final + Il_vax_final
  R_total  <- R_unvax_final  + R_vax_final
  
  N <- S_total + Ia_total + Il_total + R_total
  
  # Total cases = everyone no longer susceptible (captures R plus any
  # still infectious at the end of the simulation)
  cases <- N - S_total
  
  data.frame(
    age_group = age_groups,
    N = N,
    S_final = S_total,
    Ia_final = Ia_total,
    Il_final = Il_total,
    R_final = R_total,
    cases = cases
  )
}

# ---- Apply to existing run ----
cases_by_age_10 <- get_age_group_cases(out_ODE_10, age_groups)
print(cases_by_age_10)


# ******************* 50% of age group 1-4 is vaccinated ******************* #

# Assign initial fraction to enter into susceptible_vaccinated compartment. 
# Reduce beta by 50% and lambda by 50% for those who are vaccinated. 
# Vaccinated individuals can move through Ia_vaccinated, Il_vaccinated, and R_vaccinated. 

#------------------------------------------------------------------------------#
#                       SIR Model - Behavior-dynamic model 
#------------------------------------------------------------------------------#

SIR_dynamic_11 <- function (t, x, params) {
  # There are 7 compartments in each state to reflect the 7 age groups   
  S=x[1:7] # Susceptible 
  Ia=x[8:14] # Infectious (acute) (first 3 days of illness)
  Il=x[15:21] # Infectious (late) (next 7 days)
  R=x[22:28] # Recovered
  
  # Add vaccinated compartments 
  S_vax=x[29:35] # Susceptible_vaccinated 
  Ia_vax=x[36:42] # Infectious (acute) (first 3 days of illness)
  Il_vax=x[43:49] # Infectious (late) (next 7 days)
  R_vax=x[50:56] # Recovered
  
  beta  = params$beta # transmission rate (same for acute and late phases)
  beta_vax  = params$beta_vax # transmission rate (vaccinated) 
  gamma_a = params$gamma_a # acute infectious period = first 3 days
  gamma_l = params$gamma_l # late infectious period = next 7 days
  C_acute = params$C_acute # contact matrix for acute period
  C_late = params$C_late  # contact matrix for late period 
  
  # Population size, N
  N = S + Ia + Il + R + S_vax + Ia_vax + Il_vax + R_vax
  
  # Force of infection 
  ## Accounts for reduced transmissibility from vaccinated individuals and 
  ## reduced susceptibility of vaccinated individuals 
  infectious_contacts_unvax <- 
    t(C_acute) %*% Ia +
    t(C_late)  %*% Il
  
  infectious_contacts_vax <- 
    t(C_acute) %*% Ia_vax +
    t(C_late)  %*% Il_vax
  
  foi_sum_terms_unvax <- as.numeric(infectious_contacts_unvax) / N
  foi_sum_terms_vax   <- as.numeric(infectious_contacts_vax) / N
  
  lambda <- beta * foi_sum_terms_unvax + beta_vax * foi_sum_terms_vax 
  # This lambda is applied to all unvaccinated individuals and accounts for reduced
  # transmissibility of vaccinated individuals 
  lambda_vax <- 0.5*lambda 
  # This lambda is specific to vaccinated individuals and additionally accounts 
  # for reduced susceptibility due to vaccination 
  
  # Create a set of ODEs for the SIR model
  dxdt <- c(
    # Unvaccinated 
    - lambda * S, # S compartment
    lambda * S - gamma_a * Ia, # S -> Ia compartment
    gamma_a * Ia - gamma_l * Il, # Ia -> Il compartment
    gamma_l * Il, # Il -> R compartment
    
    # Vaccinated 
    ## Includes lambda_vax instead of lambda
    -lambda_vax * S_vax,
    lambda_vax * S_vax - gamma_a * Ia_vax,
    gamma_a * Ia_vax - gamma_l * Il_vax,
    gamma_l * Il_vax
  ) 
  
  return(list(dxdt))
}

# Set parameter values
params <- list(
  beta = 0.01,     # transmission rate
  beta_vax = 0.005, # vax transmission rate
  gamma_a = 1/3,     # acute infectious period = first 3 days 
  gamma_l = 1/7,    # late infectious period = next 7 days
  C_acute = contact_matrix_acute, # contact matrix for acute period
  C_late  = contact_matrix_late # contact matrix for late period  
)

# Make sure contact matrices are stored as numeric matrices 
params$C_acute <- as.matrix(contact_matrix_acute)
params$C_late  <- as.matrix(contact_matrix_late)

storage.mode(params$C_acute) <- "double"
storage.mode(params$C_late)  <- "double"

# Incorporate age structure
age_group_nums <- with(sim_population, setNames(age_group_total, age_group))

# Set initial conditions using population data
initial_infec_fraction <- 1e-5

## Assign the same proportion of infected individuals to each age group 

# Split population
S0_total <- age_group_nums * (1 - initial_infec_fraction)
Ia0_total <- age_group_nums * initial_infec_fraction

# 50% vaccinated
vax_coverage <- c(0, 0.5, 0, 0, 0, 0, 0) # only 50% of age group 1-4 is vaccinated
# Apply vaccination by age group
S0_vax <- vax_coverage * S0_total
S0     <- (1 - vax_coverage) * S0_total

Ia0_vax <- vax_coverage * Ia0_total
Ia0     <- (1 - vax_coverage) * Ia0_total

Il0      <- rep(0, 7)
Il0_vax  <- rep(0, 7)

R0       <- rep(0, 7)
R0_vax   <- rep(0, 7)

# Combine all compartments
initial_state <- c(
  S0, Ia0, Il0, R0,
  S0_vax, Ia0_vax, Il0_vax, R0_vax
)

# Check initial values 
sum(c(S0, Ia0, Il0, R0))
sum(sim_population$age_group_total)
sum(Ia0)

# Set time steps
times <- seq(0, 1000, by = 1)

# Run ODE (deterministic)
out_ODE_11 <- ode(
  y = initial_state,
  times = times,
  func = SIR_dynamic_11,
  parms = params
)

out_ODE_11 <- as.data.frame(out_ODE_11)

# Check output
head(out_ODE_11)
tail(out_ODE_11)

final_11 <- tail(out_ODE_11, 1)
R_unvax_11 <- as.numeric(final_11[23:29])
R_vax_11   <- as.numeric(final_11[51:57])
total_final_cases_11 <- sum(R_unvax_11 + R_vax_11)
print(total_final_cases_11)

# cases per 100,000
total_final_cases_11/(sum(initial_state))*100000

# ---- Extract per-age-group final cases from an ode() output ----

get_age_group_cases <- function(out_ODE_11, age_groups) {
  final_row <- tail(out_ODE_11, 1)
  
  # +1 offset everywhere because column 1 of out_ODE is "time"
  # Unvaccinated
  S_unvax_final  <- as.numeric(final_row[1 + (1:7)])
  Ia_unvax_final <- as.numeric(final_row[1 + (8:14)])
  Il_unvax_final <- as.numeric(final_row[1 + (15:21)])
  R_unvax_final  <- as.numeric(final_row[1 + (22:28)])
  
  # Vaccinated
  S_vax_final  <- as.numeric(final_row[1 + (29:35)])
  Ia_vax_final <- as.numeric(final_row[1 + (36:42)])
  Il_vax_final <- as.numeric(final_row[1 + (43:49)])
  R_vax_final  <- as.numeric(final_row[1 + (50:56)])
  
  # Combine vax + unvax per age group
  S_total  <- S_unvax_final  + S_vax_final
  Ia_total <- Ia_unvax_final + Ia_vax_final
  Il_total <- Il_unvax_final + Il_vax_final
  R_total  <- R_unvax_final  + R_vax_final
  
  N <- S_total + Ia_total + Il_total + R_total
  
  # Total cases = everyone no longer susceptible (captures R plus any
  # still infectious at the end of the simulation)
  cases <- N - S_total
  
  data.frame(
    age_group = age_groups,
    N = N,
    S_final = S_total,
    Ia_final = Ia_total,
    Il_final = Il_total,
    R_final = R_total,
    cases = cases
  )
}

# ---- Apply to existing run ----
cases_by_age_11 <- get_age_group_cases(out_ODE_11, age_groups)
print(cases_by_age_11)


#------------------------------------------------------------------------------#
#                       SIR Model - Behavior-naive model 
#------------------------------------------------------------------------------#

SIR_naive_12 <- function (t, x, params) {
  # There are 7 compartments in each state to reflect the 7 age groups   
  S=x[1:7] # Susceptible 
  Ia=x[8:14] # Infectious (acute) (first 3 days of illness)
  Il=x[15:21] # Infectious (late) (next 7 days)
  R=x[22:28] # Recovered
  
  # Add vaccinated compartments 
  S_vax=x[29:35] # Susceptible_vaccinated 
  Ia_vax=x[36:42] # Infectious (acute) (first 3 days of illness)
  Il_vax=x[43:49] # Infectious (late) (next 7 days)
  R_vax=x[50:56] # Recovered
  
  beta  = params$beta # transmission rate (same for acute and late phases)
  beta_vax  = params$beta_vax # transmission rate (vaccinated) 
  gamma_a = params$gamma_a # acute infectious period = first 3 days
  gamma_l = params$gamma_l # late infectious period = next 7 days
  C_naive = params$C_naive # contact matrix for acute & late periods
  
  # Population size, N
  N = S + Ia + Il + R + S_vax + Ia_vax + Il_vax + R_vax
  
  # Force of infection 
  ## Accounts for reduced transmissibility from vaccinated individuals and 
  ## reduced susceptibility of vaccinated individuals 
  infectious_contacts_unvax <- 
    t(C_naive) %*% Ia +
    t(C_naive)  %*% Il
  
  infectious_contacts_vax <- 
    t(C_naive) %*% Ia_vax +
    t(C_naive)  %*% Il_vax
  
  foi_sum_terms_unvax <- as.numeric(infectious_contacts_unvax) / N
  foi_sum_terms_vax   <- as.numeric(infectious_contacts_vax) / N
  
  lambda <- beta * foi_sum_terms_unvax + beta_vax * foi_sum_terms_vax 
  # This lambda is applied to all unvaccinated individuals and accounts for reduced
  # transmissibility of vaccinated individuals 
  lambda_vax <- 0.5*lambda 
  # This lambda is specific to vaccinated individuals and additionally accounts 
  # for reduced susceptibility due to vaccination 
  
  # Create a set of ODEs for the SIR model
  dxdt <- c(
    # Unvaccinated 
    - lambda * S, # S compartment
    lambda * S - gamma_a * Ia, # S -> Ia compartment
    gamma_a * Ia - gamma_l * Il, # Ia -> Il compartment
    gamma_l * Il, # Il -> R compartment
    
    # Vaccinated 
    ## Includes lambda_vax instead of lambda
    -lambda_vax * S_vax,
    lambda_vax * S_vax - gamma_a * Ia_vax,
    gamma_a * Ia_vax - gamma_l * Il_vax,
    gamma_l * Il_vax
  ) 
  
  return(list(dxdt))
}

# Set parameter values
params <- list(
  beta = 0.01,     # transmission rate
  beta_vax = 0.005, # vax transmission rate
  gamma_a = 1/3,     # acute infectious period = first 3 days 
  gamma_l = 1/7,    # late infectious period = next 7 days
  C_naive = contact_matrix_naive # contact matrix for acute & late periods
)

# Make sure contact matrices are stored as numeric matrices 
params$C_naive <- as.matrix(contact_matrix_naive)

storage.mode(params$C_naive) <- "double"

# Incorporate age structure
age_group_nums <- with(sim_population, setNames(age_group_total, age_group))

# Set initial conditions using population data
initial_infec_fraction <- 1e-5

## Assign the same proportion of infected individuals to each age group 

# Split population
S0_total <- age_group_nums * (1 - initial_infec_fraction)
Ia0_total <- age_group_nums * initial_infec_fraction

# 50% vaccinated
vax_coverage <- c(0, 0.5, 0, 0, 0, 0, 0) # only 50% of age group 1-4 is vaccinated
# Apply vaccination by age group
S0_vax <- vax_coverage * S0_total
S0     <- (1 - vax_coverage) * S0_total

Ia0_vax <- vax_coverage * Ia0_total
Ia0     <- (1 - vax_coverage) * Ia0_total

Il0      <- rep(0, 7)
Il0_vax  <- rep(0, 7)

R0       <- rep(0, 7)
R0_vax   <- rep(0, 7)

# Combine all compartments
initial_state <- c(
  S0, Ia0, Il0, R0,
  S0_vax, Ia0_vax, Il0_vax, R0_vax
)

# Check initial values 
sum(c(S0, Ia0, Il0, R0))
sum(sim_population$age_group_total)
sum(Ia0)

# Set time steps
times <- seq(0, 1000, by = 1)

# Run ODE (deterministic)
out_ODE_12 <- ode(
  y = initial_state,
  times = times,
  func = SIR_naive_12,
  parms = params
)

out_ODE_12 <- as.data.frame(out_ODE_12)

# Check output
head(out_ODE_12)
tail(out_ODE_12)

final_12 <- tail(out_ODE_12, 1)
R_unvax_12 <- as.numeric(final_12[23:29])
R_vax_12   <- as.numeric(final_12[51:57])
total_final_cases_12 <- sum(R_unvax_12 + R_vax_12)
print(total_final_cases_12)

# cases per 100,000
total_final_cases_12/(sum(initial_state))*100000

# ---- Extract per-age-group final cases from an ode() output ----

get_age_group_cases <- function(out_ODE_12, age_groups) {
  final_row <- tail(out_ODE_12, 1)
  
  # +1 offset everywhere because column 1 of out_ODE is "time"
  # Unvaccinated
  S_unvax_final  <- as.numeric(final_row[1 + (1:7)])
  Ia_unvax_final <- as.numeric(final_row[1 + (8:14)])
  Il_unvax_final <- as.numeric(final_row[1 + (15:21)])
  R_unvax_final  <- as.numeric(final_row[1 + (22:28)])
  
  # Vaccinated
  S_vax_final  <- as.numeric(final_row[1 + (29:35)])
  Ia_vax_final <- as.numeric(final_row[1 + (36:42)])
  Il_vax_final <- as.numeric(final_row[1 + (43:49)])
  R_vax_final  <- as.numeric(final_row[1 + (50:56)])
  
  # Combine vax + unvax per age group
  S_total  <- S_unvax_final  + S_vax_final
  Ia_total <- Ia_unvax_final + Ia_vax_final
  Il_total <- Il_unvax_final + Il_vax_final
  R_total  <- R_unvax_final  + R_vax_final
  
  N <- S_total + Ia_total + Il_total + R_total
  
  # Total cases = everyone no longer susceptible (captures R plus any
  # still infectious at the end of the simulation)
  cases <- N - S_total
  
  data.frame(
    age_group = age_groups,
    N = N,
    S_final = S_total,
    Ia_final = Ia_total,
    Il_final = Il_total,
    R_final = R_total,
    cases = cases
  )
}

# ---- Apply to existing run ----
cases_by_age_12 <- get_age_group_cases(out_ODE_12, age_groups)
print(cases_by_age_12)


# ******************* 50% of age group 5-17 is vaccinated ******************* #

# Assign initial fraction to enter into susceptible_vaccinated compartment. 
# Reduce beta by 50% and lambda by 50% for those who are vaccinated. 
# Vaccinated individuals can move through Ia_vaccinated, Il_vaccinated, and R_vaccinated. 

#------------------------------------------------------------------------------#
#                       SIR Model - Behavior-dynamic model 
#------------------------------------------------------------------------------#

SIR_dynamic_13 <- function (t, x, params) {
  # There are 7 compartments in each state to reflect the 7 age groups   
  S=x[1:7] # Susceptible 
  Ia=x[8:14] # Infectious (acute) (first 3 days of illness)
  Il=x[15:21] # Infectious (late) (next 7 days)
  R=x[22:28] # Recovered
  
  # Add vaccinated compartments 
  S_vax=x[29:35] # Susceptible_vaccinated 
  Ia_vax=x[36:42] # Infectious (acute) (first 3 days of illness)
  Il_vax=x[43:49] # Infectious (late) (next 7 days)
  R_vax=x[50:56] # Recovered
  
  beta  = params$beta # transmission rate (same for acute and late phases)
  beta_vax  = params$beta_vax # transmission rate (vaccinated) 
  gamma_a = params$gamma_a # acute infectious period = first 3 days
  gamma_l = params$gamma_l # late infectious period = next 7 days
  C_acute = params$C_acute # contact matrix for acute period
  C_late = params$C_late  # contact matrix for late period 
  
  # Population size, N
  N = S + Ia + Il + R + S_vax + Ia_vax + Il_vax + R_vax
  
  # Force of infection 
  ## Accounts for reduced transmissibility from vaccinated individuals and 
  ## reduced susceptibility of vaccinated individuals 
  infectious_contacts_unvax <- 
    t(C_acute) %*% Ia +
    t(C_late)  %*% Il
  
  infectious_contacts_vax <- 
    t(C_acute) %*% Ia_vax +
    t(C_late)  %*% Il_vax
  
  foi_sum_terms_unvax <- as.numeric(infectious_contacts_unvax) / N
  foi_sum_terms_vax   <- as.numeric(infectious_contacts_vax) / N
  
  lambda <- beta * foi_sum_terms_unvax + beta_vax * foi_sum_terms_vax 
  # This lambda is applied to all unvaccinated individuals and accounts for reduced
  # transmissibility of vaccinated individuals 
  lambda_vax <- 0.5*lambda 
  # This lambda is specific to vaccinated individuals and additionally accounts 
  # for reduced susceptibility due to vaccination 
  
  # Create a set of ODEs for the SIR model
  dxdt <- c(
    # Unvaccinated 
    - lambda * S, # S compartment
    lambda * S - gamma_a * Ia, # S -> Ia compartment
    gamma_a * Ia - gamma_l * Il, # Ia -> Il compartment
    gamma_l * Il, # Il -> R compartment
    
    # Vaccinated 
    ## Includes lambda_vax instead of lambda
    -lambda_vax * S_vax,
    lambda_vax * S_vax - gamma_a * Ia_vax,
    gamma_a * Ia_vax - gamma_l * Il_vax,
    gamma_l * Il_vax
  ) 
  
  return(list(dxdt))
}

# Set parameter values
params <- list(
  beta = 0.01,     # transmission rate
  beta_vax = 0.005, # vax transmission rate
  gamma_a = 1/3,     # acute infectious period = first 3 days 
  gamma_l = 1/7,    # late infectious period = next 7 days
  C_acute = contact_matrix_acute, # contact matrix for acute period
  C_late  = contact_matrix_late # contact matrix for late period  
)

# Make sure contact matrices are stored as numeric matrices 
params$C_acute <- as.matrix(contact_matrix_acute)
params$C_late  <- as.matrix(contact_matrix_late)

storage.mode(params$C_acute) <- "double"
storage.mode(params$C_late)  <- "double"

# Incorporate age structure
age_group_nums <- with(sim_population, setNames(age_group_total, age_group))

# Set initial conditions using population data
initial_infec_fraction <- 1e-5

## Assign the same proportion of infected individuals to each age group 

# Split population
S0_total <- age_group_nums * (1 - initial_infec_fraction)
Ia0_total <- age_group_nums * initial_infec_fraction

# 50% vaccinated
vax_coverage <- c(0, 0, 0.5, 0, 0, 0, 0) # only 50% of age group 5-17 is vaccinated
# Apply vaccination by age group
S0_vax <- vax_coverage * S0_total
S0     <- (1 - vax_coverage) * S0_total

Ia0_vax <- vax_coverage * Ia0_total
Ia0     <- (1 - vax_coverage) * Ia0_total

Il0      <- rep(0, 7)
Il0_vax  <- rep(0, 7)

R0       <- rep(0, 7)
R0_vax   <- rep(0, 7)

# Combine all compartments
initial_state <- c(
  S0, Ia0, Il0, R0,
  S0_vax, Ia0_vax, Il0_vax, R0_vax
)

# Check initial values 
sum(c(S0, Ia0, Il0, R0))
sum(sim_population$age_group_total)
sum(Ia0)

# Set time steps
times <- seq(0, 1000, by = 1)

# Run ODE (deterministic)
out_ODE_13 <- ode(
  y = initial_state,
  times = times,
  func = SIR_dynamic_13,
  parms = params
)

out_ODE_13 <- as.data.frame(out_ODE_13)

# Check output
head(out_ODE_13)
tail(out_ODE_13)

final_13 <- tail(out_ODE_13, 1)
R_unvax_13 <- as.numeric(final_13[23:29])
R_vax_13   <- as.numeric(final_13[51:57])
total_final_cases_13 <- sum(R_unvax_13 + R_vax_13)
print(total_final_cases_13)

# cases per 100,000
total_final_cases_13/(sum(initial_state))*100000

# ---- Extract per-age-group final cases from an ode() output ----

get_age_group_cases <- function(out_ODE_13, age_groups) {
  final_row <- tail(out_ODE_13, 1)
  
  # +1 offset everywhere because column 1 of out_ODE is "time"
  # Unvaccinated
  S_unvax_final  <- as.numeric(final_row[1 + (1:7)])
  Ia_unvax_final <- as.numeric(final_row[1 + (8:14)])
  Il_unvax_final <- as.numeric(final_row[1 + (15:21)])
  R_unvax_final  <- as.numeric(final_row[1 + (22:28)])
  
  # Vaccinated
  S_vax_final  <- as.numeric(final_row[1 + (29:35)])
  Ia_vax_final <- as.numeric(final_row[1 + (36:42)])
  Il_vax_final <- as.numeric(final_row[1 + (43:49)])
  R_vax_final  <- as.numeric(final_row[1 + (50:56)])
  
  # Combine vax + unvax per age group
  S_total  <- S_unvax_final  + S_vax_final
  Ia_total <- Ia_unvax_final + Ia_vax_final
  Il_total <- Il_unvax_final + Il_vax_final
  R_total  <- R_unvax_final  + R_vax_final
  
  N <- S_total + Ia_total + Il_total + R_total
  
  # Total cases = everyone no longer susceptible (captures R plus any
  # still infectious at the end of the simulation)
  cases <- N - S_total
  
  data.frame(
    age_group = age_groups,
    N = N,
    S_final = S_total,
    Ia_final = Ia_total,
    Il_final = Il_total,
    R_final = R_total,
    cases = cases
  )
}

# ---- Apply to existing run ----
cases_by_age_13 <- get_age_group_cases(out_ODE_13, age_groups)
print(cases_by_age_13)


#------------------------------------------------------------------------------#
#                       SIR Model - Behavior-naive model 
#------------------------------------------------------------------------------#

SIR_naive_14 <- function (t, x, params) {
  # There are 7 compartments in each state to reflect the 7 age groups   
  S=x[1:7] # Susceptible 
  Ia=x[8:14] # Infectious (acute) (first 3 days of illness)
  Il=x[15:21] # Infectious (late) (next 7 days)
  R=x[22:28] # Recovered
  
  # Add vaccinated compartments 
  S_vax=x[29:35] # Susceptible_vaccinated 
  Ia_vax=x[36:42] # Infectious (acute) (first 3 days of illness)
  Il_vax=x[43:49] # Infectious (late) (next 7 days)
  R_vax=x[50:56] # Recovered
  
  beta  = params$beta # transmission rate (same for acute and late phases)
  beta_vax  = params$beta_vax # transmission rate (vaccinated) 
  gamma_a = params$gamma_a # acute infectious period = first 3 days
  gamma_l = params$gamma_l # late infectious period = next 7 days
  C_naive = params$C_naive # contact matrix for acute & late periods
  
  # Population size, N
  N = S + Ia + Il + R + S_vax + Ia_vax + Il_vax + R_vax
  
  # Force of infection 
  ## Accounts for reduced transmissibility from vaccinated individuals and 
  ## reduced susceptibility of vaccinated individuals 
  infectious_contacts_unvax <- 
    t(C_naive) %*% Ia +
    t(C_naive)  %*% Il
  
  infectious_contacts_vax <- 
    t(C_naive) %*% Ia_vax +
    t(C_naive)  %*% Il_vax
  
  foi_sum_terms_unvax <- as.numeric(infectious_contacts_unvax) / N
  foi_sum_terms_vax   <- as.numeric(infectious_contacts_vax) / N
  
  lambda <- beta * foi_sum_terms_unvax + beta_vax * foi_sum_terms_vax 
  # This lambda is applied to all unvaccinated individuals and accounts for reduced
  # transmissibility of vaccinated individuals 
  lambda_vax <- 0.5*lambda 
  # This lambda is specific to vaccinated individuals and additionally accounts 
  # for reduced susceptibility due to vaccination 
  
  # Create a set of ODEs for the SIR model
  dxdt <- c(
    # Unvaccinated 
    - lambda * S, # S compartment
    lambda * S - gamma_a * Ia, # S -> Ia compartment
    gamma_a * Ia - gamma_l * Il, # Ia -> Il compartment
    gamma_l * Il, # Il -> R compartment
    
    # Vaccinated 
    ## Includes lambda_vax instead of lambda
    -lambda_vax * S_vax,
    lambda_vax * S_vax - gamma_a * Ia_vax,
    gamma_a * Ia_vax - gamma_l * Il_vax,
    gamma_l * Il_vax
  ) 
  
  return(list(dxdt))
}

# Set parameter values
params <- list(
  beta = 0.01,     # transmission rate
  beta_vax = 0.005, # vax transmission rate
  gamma_a = 1/3,     # acute infectious period = first 3 days 
  gamma_l = 1/7,    # late infectious period = next 7 days
  C_naive = contact_matrix_naive # contact matrix for acute & late periods
)

# Make sure contact matrices are stored as numeric matrices 
params$C_naive <- as.matrix(contact_matrix_naive)

storage.mode(params$C_naive) <- "double"

# Incorporate age structure
age_group_nums <- with(sim_population, setNames(age_group_total, age_group))

# Set initial conditions using population data
initial_infec_fraction <- 1e-5

## Assign the same proportion of infected individuals to each age group 

# Split population
S0_total <- age_group_nums * (1 - initial_infec_fraction)
Ia0_total <- age_group_nums * initial_infec_fraction

# 50% vaccinated
vax_coverage <- c(0, 0, 0.5, 0, 0, 0, 0) # only 50% of age group 5-17 is vaccinated
# Apply vaccination by age group
S0_vax <- vax_coverage * S0_total
S0     <- (1 - vax_coverage) * S0_total

Ia0_vax <- vax_coverage * Ia0_total
Ia0     <- (1 - vax_coverage) * Ia0_total

Il0      <- rep(0, 7)
Il0_vax  <- rep(0, 7)

R0       <- rep(0, 7)
R0_vax   <- rep(0, 7)

# Combine all compartments
initial_state <- c(
  S0, Ia0, Il0, R0,
  S0_vax, Ia0_vax, Il0_vax, R0_vax
)

# Check initial values 
sum(c(S0, Ia0, Il0, R0))
sum(sim_population$age_group_total)
sum(Ia0)

# Set time steps
times <- seq(0, 1000, by = 1)

# Run ODE (deterministic)
out_ODE_14 <- ode(
  y = initial_state,
  times = times,
  func = SIR_naive_14,
  parms = params
)

out_ODE_14 <- as.data.frame(out_ODE_14)

final_14 <- tail(out_ODE_14, 1)
R_unvax_14 <- as.numeric(final_14[23:29])
R_vax_14   <- as.numeric(final_14[51:57])
total_final_cases_14 <- sum(R_unvax_14 + R_vax_14)
print(total_final_cases_14)

# cases per 100,000
total_final_cases_14/(sum(initial_state))*100000

# ---- Extract per-age-group final cases from an ode() output ----

get_age_group_cases <- function(out_ODE_14, age_groups) {
  final_row <- tail(out_ODE_14, 1)
  
  
  # +1 offset everywhere because column 1 of out_ODE is "time"
  # Unvaccinated
  S_unvax_final  <- as.numeric(final_row[1 + (1:7)])
  Ia_unvax_final <- as.numeric(final_row[1 + (8:14)])
  Il_unvax_final <- as.numeric(final_row[1 + (15:21)])
  R_unvax_final  <- as.numeric(final_row[1 + (22:28)])
  
  # Vaccinated
  S_vax_final  <- as.numeric(final_row[1 + (29:35)])
  Ia_vax_final <- as.numeric(final_row[1 + (36:42)])
  Il_vax_final <- as.numeric(final_row[1 + (43:49)])
  R_vax_final  <- as.numeric(final_row[1 + (50:56)])
  
  # Combine vax + unvax per age group
  S_total  <- S_unvax_final  + S_vax_final
  Ia_total <- Ia_unvax_final + Ia_vax_final
  Il_total <- Il_unvax_final + Il_vax_final
  R_total  <- R_unvax_final  + R_vax_final
  
  N <- S_total + Ia_total + Il_total + R_total
  
  # Total cases = everyone no longer susceptible (captures R plus any
  # still infectious at the end of the simulation)
  cases <- N - S_total
  
  data.frame(
    age_group = age_groups,
    N = N,
    S_final = S_total,
    Ia_final = Ia_total,
    Il_final = Il_total,
    R_final = R_total,
    cases = cases
  )
}

# ---- Apply to existing run ----
cases_by_age_14 <- get_age_group_cases(out_ODE_14, age_groups)
print(cases_by_age_14)


# ******************* 50% of age group 65+ is vaccinated ******************* #

# Assign initial fraction to enter into susceptible_vaccinated compartment. 
# Reduce beta by 50% and lambda by 50% for those who are vaccinated. 
# Vaccinated individuals can move through Ia_vaccinated, Il_vaccinated, and R_vaccinated. 

#------------------------------------------------------------------------------#
#                       SIR Model - Behavior-dynamic model 
#------------------------------------------------------------------------------#

SIR_dynamic_15 <- function (t, x, params) {
  # There are 7 compartments in each state to reflect the 7 age groups   
  S=x[1:7] # Susceptible 
  Ia=x[8:14] # Infectious (acute) (first 3 days of illness)
  Il=x[15:21] # Infectious (late) (next 7 days)
  R=x[22:28] # Recovered
  
  # Add vaccinated compartments 
  S_vax=x[29:35] # Susceptible_vaccinated 
  Ia_vax=x[36:42] # Infectious (acute) (first 3 days of illness)
  Il_vax=x[43:49] # Infectious (late) (next 7 days)
  R_vax=x[50:56] # Recovered
  
  beta  = params$beta # transmission rate (same for acute and late phases)
  beta_vax  = params$beta_vax # transmission rate (vaccinated) 
  gamma_a = params$gamma_a # acute infectious period = first 3 days
  gamma_l = params$gamma_l # late infectious period = next 7 days
  C_acute = params$C_acute # contact matrix for acute period
  C_late = params$C_late  # contact matrix for late period 
  
  # Population size, N
  N = S + Ia + Il + R + S_vax + Ia_vax + Il_vax + R_vax
  
  # Force of infection 
  ## Accounts for reduced transmissibility from vaccinated individuals and 
  ## reduced susceptibility of vaccinated individuals 
  infectious_contacts_unvax <- 
    t(C_acute) %*% Ia +
    t(C_late)  %*% Il
  
  infectious_contacts_vax <- 
    t(C_acute) %*% Ia_vax +
    t(C_late)  %*% Il_vax
  
  foi_sum_terms_unvax <- as.numeric(infectious_contacts_unvax) / N
  foi_sum_terms_vax   <- as.numeric(infectious_contacts_vax) / N
  
  lambda <- beta * foi_sum_terms_unvax + beta_vax * foi_sum_terms_vax 
  # This lambda is applied to all unvaccinated individuals and accounts for reduced
  # transmissibility of vaccinated individuals 
  lambda_vax <- 0.5*lambda 
  # This lambda is specific to vaccinated individuals and additionally accounts 
  # for reduced susceptibility due to vaccination 
  
  # Create a set of ODEs for the SIR model
  dxdt <- c(
    # Unvaccinated 
    - lambda * S, # S compartment
    lambda * S - gamma_a * Ia, # S -> Ia compartment
    gamma_a * Ia - gamma_l * Il, # Ia -> Il compartment
    gamma_l * Il, # Il -> R compartment
    
    # Vaccinated 
    ## Includes lambda_vax instead of lambda
    -lambda_vax * S_vax,
    lambda_vax * S_vax - gamma_a * Ia_vax,
    gamma_a * Ia_vax - gamma_l * Il_vax,
    gamma_l * Il_vax
  ) 
  
  return(list(dxdt))
}

# Set parameter values
params <- list(
  beta = 0.01,     # transmission rate
  beta_vax = 0.005, # vax transmission rate
  gamma_a = 1/3,     # acute infectious period = first 3 days 
  gamma_l = 1/7,    # late infectious period = next 7 days
  C_acute = contact_matrix_acute, # contact matrix for acute period
  C_late  = contact_matrix_late # contact matrix for late period  
)

# Make sure contact matrices are stored as numeric matrices 
params$C_acute <- as.matrix(contact_matrix_acute)
params$C_late  <- as.matrix(contact_matrix_late)

storage.mode(params$C_acute) <- "double"
storage.mode(params$C_late)  <- "double"

# Incorporate age structure
age_group_nums <- with(sim_population, setNames(age_group_total, age_group))

# Set initial conditions using population data
initial_infec_fraction <- 1e-5

## Assign the same proportion of infected individuals to each age group 

# Split population
S0_total <- age_group_nums * (1 - initial_infec_fraction)
Ia0_total <- age_group_nums * initial_infec_fraction

# 50% vaccinated
vax_coverage <- c(0, 0, 0, 0, 0, 0, 0.5) # only 50% of age group 65+ is vaccinated
# Apply vaccination by age group
S0_vax <- vax_coverage * S0_total
S0     <- (1 - vax_coverage) * S0_total

Ia0_vax <- vax_coverage * Ia0_total
Ia0     <- (1 - vax_coverage) * Ia0_total

Il0      <- rep(0, 7)
Il0_vax  <- rep(0, 7)

R0       <- rep(0, 7)
R0_vax   <- rep(0, 7)

# Combine all compartments
initial_state <- c(
  S0, Ia0, Il0, R0,
  S0_vax, Ia0_vax, Il0_vax, R0_vax
)

# Check initial values 
sum(c(S0, Ia0, Il0, R0))
sum(sim_population$age_group_total)
sum(Ia0)

# Set time steps
times <- seq(0, 1000, by = 1)

# Run ODE (deterministic)
out_ODE_15 <- ode(
  y = initial_state,
  times = times,
  func = SIR_dynamic_15,
  parms = params
)

out_ODE_15 <- as.data.frame(out_ODE_15)

# Check output
head(out_ODE_15)
tail(out_ODE_15)

final_15 <- tail(out_ODE_15, 1)
R_unvax_15 <- as.numeric(final_15[23:29])
R_vax_15   <- as.numeric(final_15[51:57])
total_final_cases_15 <- sum(R_unvax_15 + R_vax_15)
print(total_final_cases_15)

# cases per 100,000
total_final_cases_15/(sum(initial_state))*100000

# ---- Extract per-age-group final cases from an ode() output ----

get_age_group_cases <- function(out_ODE_15, age_groups) {
  final_row <- tail(out_ODE_15, 1)
  
  
  # +1 offset everywhere because column 1 of out_ODE is "time"
  # Unvaccinated
  S_unvax_final  <- as.numeric(final_row[1 + (1:7)])
  Ia_unvax_final <- as.numeric(final_row[1 + (8:14)])
  Il_unvax_final <- as.numeric(final_row[1 + (15:21)])
  R_unvax_final  <- as.numeric(final_row[1 + (22:28)])
  
  # Vaccinated
  S_vax_final  <- as.numeric(final_row[1 + (29:35)])
  Ia_vax_final <- as.numeric(final_row[1 + (36:42)])
  Il_vax_final <- as.numeric(final_row[1 + (43:49)])
  R_vax_final  <- as.numeric(final_row[1 + (50:56)])
  
  # Combine vax + unvax per age group
  S_total  <- S_unvax_final  + S_vax_final
  Ia_total <- Ia_unvax_final + Ia_vax_final
  Il_total <- Il_unvax_final + Il_vax_final
  R_total  <- R_unvax_final  + R_vax_final
  
  N <- S_total + Ia_total + Il_total + R_total
  
  # Total cases = everyone no longer susceptible (captures R plus any
  # still infectious at the end of the simulation)
  cases <- N - S_total
  
  data.frame(
    age_group = age_groups,
    N = N,
    S_final = S_total,
    Ia_final = Ia_total,
    Il_final = Il_total,
    R_final = R_total,
    cases = cases
  )
}

# ---- Apply to existing run ----
cases_by_age_15 <- get_age_group_cases(out_ODE_15, age_groups)
print(cases_by_age_15)


#------------------------------------------------------------------------------#
#                       SIR Model - Behavior-naive model 
#------------------------------------------------------------------------------#

SIR_naive_16 <- function (t, x, params) {
  # There are 7 compartments in each state to reflect the 7 age groups   
  S=x[1:7] # Susceptible 
  Ia=x[8:14] # Infectious (acute) (first 3 days of illness)
  Il=x[15:21] # Infectious (late) (next 7 days)
  R=x[22:28] # Recovered
  
  # Add vaccinated compartments 
  S_vax=x[29:35] # Susceptible_vaccinated 
  Ia_vax=x[36:42] # Infectious (acute) (first 3 days of illness)
  Il_vax=x[43:49] # Infectious (late) (next 7 days)
  R_vax=x[50:56] # Recovered
  
  beta  = params$beta # transmission rate (same for acute and late phases)
  beta_vax  = params$beta_vax # transmission rate (vaccinated) 
  gamma_a = params$gamma_a # acute infectious period = first 3 days
  gamma_l = params$gamma_l # late infectious period = next 7 days
  C_naive = params$C_naive # contact matrix for acute & late periods
  
  # Population size, N
  N = S + Ia + Il + R + S_vax + Ia_vax + Il_vax + R_vax
  
  # Force of infection 
  ## Accounts for reduced transmissibility from vaccinated individuals and 
  ## reduced susceptibility of vaccinated individuals 
  infectious_contacts_unvax <- 
    t(C_naive) %*% Ia +
    t(C_naive)  %*% Il
  
  infectious_contacts_vax <- 
    t(C_naive) %*% Ia_vax +
    t(C_naive)  %*% Il_vax
  
  foi_sum_terms_unvax <- as.numeric(infectious_contacts_unvax) / N
  foi_sum_terms_vax   <- as.numeric(infectious_contacts_vax) / N
  
  lambda <- beta * foi_sum_terms_unvax + beta_vax * foi_sum_terms_vax 
  # This lambda is applied to all unvaccinated individuals and accounts for reduced
  # transmissibility of vaccinated individuals 
  lambda_vax <- 0.5*lambda 
  # This lambda is specific to vaccinated individuals and additionally accounts 
  # for reduced susceptibility due to vaccination 
  
  # Create a set of ODEs for the SIR model
  dxdt <- c(
    # Unvaccinated 
    - lambda * S, # S compartment
    lambda * S - gamma_a * Ia, # S -> Ia compartment
    gamma_a * Ia - gamma_l * Il, # Ia -> Il compartment
    gamma_l * Il, # Il -> R compartment
    
    # Vaccinated 
    ## Includes lambda_vax instead of lambda
    -lambda_vax * S_vax,
    lambda_vax * S_vax - gamma_a * Ia_vax,
    gamma_a * Ia_vax - gamma_l * Il_vax,
    gamma_l * Il_vax
  ) 
  
  return(list(dxdt))
}

# Set parameter values
params <- list(
  beta = 0.01,     # transmission rate
  beta_vax = 0.005, # vax transmission rate
  gamma_a = 1/3,     # acute infectious period = first 3 days 
  gamma_l = 1/7,    # late infectious period = next 7 days
  C_naive = contact_matrix_naive # contact matrix for acute & late periods
)

# Make sure contact matrices are stored as numeric matrices 
params$C_naive <- as.matrix(contact_matrix_naive)

storage.mode(params$C_naive) <- "double"

# Incorporate age structure
age_group_nums <- with(sim_population, setNames(age_group_total, age_group))

# Set initial conditions using population data
initial_infec_fraction <- 1e-5

## Assign the same proportion of infected individuals to each age group 

# Split population
S0_total <- age_group_nums * (1 - initial_infec_fraction)
Ia0_total <- age_group_nums * initial_infec_fraction

# 50% vaccinated
vax_coverage <- c(0, 0, 0, 0, 0, 0, 0.5) # only 50% of age group 65+ is vaccinated
# Apply vaccination by age group
S0_vax <- vax_coverage * S0_total
S0     <- (1 - vax_coverage) * S0_total

Ia0_vax <- vax_coverage * Ia0_total
Ia0     <- (1 - vax_coverage) * Ia0_total

Il0      <- rep(0, 7)
Il0_vax  <- rep(0, 7)

R0       <- rep(0, 7)
R0_vax   <- rep(0, 7)

# Combine all compartments
initial_state <- c(
  S0, Ia0, Il0, R0,
  S0_vax, Ia0_vax, Il0_vax, R0_vax
)

# Check initial values 
sum(c(S0, Ia0, Il0, R0))
sum(sim_population$age_group_total)
sum(Ia0)

# Set time steps
times <- seq(0, 1000, by = 1)

# Run ODE (deterministic)
out_ODE_16 <- ode(
  y = initial_state,
  times = times,
  func = SIR_naive_16,
  parms = params
)

out_ODE_16 <- as.data.frame(out_ODE_16)

# Check output
head(out_ODE_16)
tail(out_ODE_16)

final_16 <- tail(out_ODE_16, 1)
R_unvax_16 <- as.numeric(final_16[23:29])
R_vax_16   <- as.numeric(final_16[51:57])
total_final_cases_16 <- sum(R_unvax_16 + R_vax_16)
print(total_final_cases_16)

# cases per 100,000
total_final_cases_16/(sum(initial_state))*100000

# ---- Extract per-age-group final cases from an ode() output ----

get_age_group_cases <- function(out_ODE_16, age_groups) {
  final_row <- tail(out_ODE_16, 1)
  
  
  # +1 offset everywhere because column 1 of out_ODE is "time"
  # Unvaccinated
  S_unvax_final  <- as.numeric(final_row[1 + (1:7)])
  Ia_unvax_final <- as.numeric(final_row[1 + (8:14)])
  Il_unvax_final <- as.numeric(final_row[1 + (15:21)])
  R_unvax_final  <- as.numeric(final_row[1 + (22:28)])
  
  # Vaccinated
  S_vax_final  <- as.numeric(final_row[1 + (29:35)])
  Ia_vax_final <- as.numeric(final_row[1 + (36:42)])
  Il_vax_final <- as.numeric(final_row[1 + (43:49)])
  R_vax_final  <- as.numeric(final_row[1 + (50:56)])
  
  # Combine vax + unvax per age group
  S_total  <- S_unvax_final  + S_vax_final
  Ia_total <- Ia_unvax_final + Ia_vax_final
  Il_total <- Il_unvax_final + Il_vax_final
  R_total  <- R_unvax_final  + R_vax_final
  
  N <- S_total + Ia_total + Il_total + R_total
  
  # Total cases = everyone no longer susceptible (captures R plus any
  # still infectious at the end of the simulation)
  cases <- N - S_total
  
  data.frame(
    age_group = age_groups,
    N = N,
    S_final = S_total,
    Ia_final = Ia_total,
    Il_final = Il_total,
    R_final = R_total,
    cases = cases
  )
}

# ---- Apply to existing run ----
cases_by_age_16 <- get_age_group_cases(out_ODE_16, age_groups)
print(cases_by_age_16)


# ******************* 50% of age group >1 is vaccinated ******************* #

# Assign initial fraction to enter into susceptible_vaccinated compartment. 
# Reduce beta by 50% and lambda by 50% for those who are vaccinated. 
# Vaccinated individuals can move through Ia_vaccinated, Il_vaccinated, and R_vaccinated. 

#------------------------------------------------------------------------------#
#                       SIR Model - Behavior-dynamic model 
#------------------------------------------------------------------------------#

SIR_dynamic_17 <- function (t, x, params) {
  # There are 7 compartments in each state to reflect the 7 age groups   
  S=x[1:7] # Susceptible 
  Ia=x[8:14] # Infectious (acute) (first 3 days of illness)
  Il=x[15:21] # Infectious (late) (next 7 days)
  R=x[22:28] # Recovered
  
  # Add vaccinated compartments 
  S_vax=x[29:35] # Susceptible_vaccinated 
  Ia_vax=x[36:42] # Infectious (acute) (first 3 days of illness)
  Il_vax=x[43:49] # Infectious (late) (next 7 days)
  R_vax=x[50:56] # Recovered
  
  beta  = params$beta # transmission rate (same for acute and late phases)
  beta_vax  = params$beta_vax # transmission rate (vaccinated) 
  gamma_a = params$gamma_a # acute infectious period = first 3 days
  gamma_l = params$gamma_l # late infectious period = next 7 days
  C_acute = params$C_acute # contact matrix for acute period
  C_late = params$C_late  # contact matrix for late period 
  
  # Population size, N
  N = S + Ia + Il + R + S_vax + Ia_vax + Il_vax + R_vax
  
  # Force of infection 
  ## Accounts for reduced transmissibility from vaccinated individuals and 
  ## reduced susceptibility of vaccinated individuals 
  infectious_contacts_unvax <- 
    t(C_acute) %*% Ia +
    t(C_late)  %*% Il
  
  infectious_contacts_vax <- 
    t(C_acute) %*% Ia_vax +
    t(C_late)  %*% Il_vax
  
  foi_sum_terms_unvax <- as.numeric(infectious_contacts_unvax) / N
  foi_sum_terms_vax   <- as.numeric(infectious_contacts_vax) / N
  
  lambda <- beta * foi_sum_terms_unvax + beta_vax * foi_sum_terms_vax 
  # This lambda is applied to all unvaccinated individuals and accounts for reduced
  # transmissibility of vaccinated individuals 
  lambda_vax <- 0.5*lambda 
  # This lambda is specific to vaccinated individuals and additionally accounts 
  # for reduced susceptibility due to vaccination 
  
  # Create a set of ODEs for the SIR model
  dxdt <- c(
    # Unvaccinated 
    - lambda * S, # S compartment
    lambda * S - gamma_a * Ia, # S -> Ia compartment
    gamma_a * Ia - gamma_l * Il, # Ia -> Il compartment
    gamma_l * Il, # Il -> R compartment
    
    # Vaccinated 
    ## Includes lambda_vax instead of lambda
    -lambda_vax * S_vax,
    lambda_vax * S_vax - gamma_a * Ia_vax,
    gamma_a * Ia_vax - gamma_l * Il_vax,
    gamma_l * Il_vax
  ) 
  
  return(list(dxdt))
}

# Set parameter values
params <- list(
  beta = 0.01,     # transmission rate
  beta_vax = 0.005, # vax transmission rate
  gamma_a = 1/3,     # acute infectious period = first 3 days 
  gamma_l = 1/7,    # late infectious period = next 7 days
  C_acute = contact_matrix_acute, # contact matrix for acute period
  C_late  = contact_matrix_late # contact matrix for late period  
)

# Make sure contact matrices are stored as numeric matrices 
params$C_acute <- as.matrix(contact_matrix_acute)
params$C_late  <- as.matrix(contact_matrix_late)

storage.mode(params$C_acute) <- "double"
storage.mode(params$C_late)  <- "double"

# Incorporate age structure
age_group_nums <- with(sim_population, setNames(age_group_total, age_group))

# Set initial conditions using population data
initial_infec_fraction <- 1e-5

## Assign the same proportion of infected individuals to each age group 

# Split population
S0_total <- age_group_nums * (1 - initial_infec_fraction)
Ia0_total <- age_group_nums * initial_infec_fraction

# 50% vaccinated
vax_coverage <- c(0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5) # 50% of everyone over the age of 1 is vaccinated
# Apply vaccination by age group
S0_vax <- vax_coverage * S0_total
S0     <- (1 - vax_coverage) * S0_total

Ia0_vax <- vax_coverage * Ia0_total
Ia0     <- (1 - vax_coverage) * Ia0_total

Il0      <- rep(0, 7)
Il0_vax  <- rep(0, 7)

R0       <- rep(0, 7)
R0_vax   <- rep(0, 7)

# Combine all compartments
initial_state <- c(
  S0, Ia0, Il0, R0,
  S0_vax, Ia0_vax, Il0_vax, R0_vax
)

# Check initial values 
sum(c(S0, Ia0, Il0, R0))
sum(sim_population$age_group_total)
sum(Ia0)

# Set time steps
times <- seq(0, 1500, by = 1)

# Run ODE (deterministic)
out_ODE_17 <- ode(
  y = initial_state,
  times = times,
  func = SIR_dynamic_17,
  parms = params
)

final_17 <- tail(out_ODE_17, 1)
R_unvax_17 <- as.numeric(final_17[23:29])
R_vax_17   <- as.numeric(final_17[51:57])
total_final_cases_17 <- sum(R_unvax_17 + R_vax_17)
print(total_final_cases_17)

# cases per 100,000
total_final_cases_17/(sum(initial_state))*100000

# ---- Extract per-age-group final cases from an ode() output ----

get_age_group_cases <- function(out_ODE_17, age_groups) {
  final_row <- tail(out_ODE_17, 1)
  
  
  # +1 offset everywhere because column 1 of out_ODE is "time"
  # Unvaccinated
  S_unvax_final  <- as.numeric(final_row[1 + (1:7)])
  Ia_unvax_final <- as.numeric(final_row[1 + (8:14)])
  Il_unvax_final <- as.numeric(final_row[1 + (15:21)])
  R_unvax_final  <- as.numeric(final_row[1 + (22:28)])
  
  # Vaccinated
  S_vax_final  <- as.numeric(final_row[1 + (29:35)])
  Ia_vax_final <- as.numeric(final_row[1 + (36:42)])
  Il_vax_final <- as.numeric(final_row[1 + (43:49)])
  R_vax_final  <- as.numeric(final_row[1 + (50:56)])
  
  # Combine vax + unvax per age group
  S_total  <- S_unvax_final  + S_vax_final
  Ia_total <- Ia_unvax_final + Ia_vax_final
  Il_total <- Il_unvax_final + Il_vax_final
  R_total  <- R_unvax_final  + R_vax_final
  
  N <- S_total + Ia_total + Il_total + R_total
  
  # Total cases = everyone no longer susceptible (captures R plus any
  # still infectious at the end of the simulation)
  cases <- N - S_total
  
  data.frame(
    age_group = age_groups,
    N = N,
    S_final = S_total,
    Ia_final = Ia_total,
    Il_final = Il_total,
    R_final = R_total,
    cases = cases
  )
}

# ---- Apply to existing run ----
cases_by_age_17 <- get_age_group_cases(out_ODE_17, age_groups)
print(cases_by_age_17)


#------------------------------------------------------------------------------#
#                       SIR Model - Behavior-naive model 
#------------------------------------------------------------------------------#

SIR_naive_18 <- function (t, x, params) {
  # There are 7 compartments in each state to reflect the 7 age groups   
  S=x[1:7] # Susceptible 
  Ia=x[8:14] # Infectious (acute) (first 3 days of illness)
  Il=x[15:21] # Infectious (late) (next 7 days)
  R=x[22:28] # Recovered
  
  # Add vaccinated compartments 
  S_vax=x[29:35] # Susceptible_vaccinated 
  Ia_vax=x[36:42] # Infectious (acute) (first 3 days of illness)
  Il_vax=x[43:49] # Infectious (late) (next 7 days)
  R_vax=x[50:56] # Recovered
  
  beta  = params$beta # transmission rate (same for acute and late phases)
  beta_vax  = params$beta_vax # transmission rate (vaccinated) 
  gamma_a = params$gamma_a # acute infectious period = first 3 days
  gamma_l = params$gamma_l # late infectious period = next 7 days
  C_naive = params$C_naive # contact matrix for acute & late periods
  
  # Population size, N
  N = S + Ia + Il + R + S_vax + Ia_vax + Il_vax + R_vax
  
  # Force of infection 
  ## Accounts for reduced transmissibility from vaccinated individuals and 
  ## reduced susceptibility of vaccinated individuals 
  infectious_contacts_unvax <- 
    t(C_naive) %*% Ia +
    t(C_naive)  %*% Il
  
  infectious_contacts_vax <- 
    t(C_naive) %*% Ia_vax +
    t(C_naive)  %*% Il_vax
  
  foi_sum_terms_unvax <- as.numeric(infectious_contacts_unvax) / N
  foi_sum_terms_vax   <- as.numeric(infectious_contacts_vax) / N
  
  lambda <- beta * foi_sum_terms_unvax + beta_vax * foi_sum_terms_vax 
  # This lambda is applied to all unvaccinated individuals and accounts for reduced
  # transmissibility of vaccinated individuals 
  lambda_vax <- 0.5*lambda 
  # This lambda is specific to vaccinated individuals and additionally accounts 
  # for reduced susceptibility due to vaccination 
  
  # Create a set of ODEs for the SIR model
  dxdt <- c(
    # Unvaccinated 
    - lambda * S, # S compartment
    lambda * S - gamma_a * Ia, # S -> Ia compartment
    gamma_a * Ia - gamma_l * Il, # Ia -> Il compartment
    gamma_l * Il, # Il -> R compartment
    
    # Vaccinated 
    ## Includes lambda_vax instead of lambda
    -lambda_vax * S_vax,
    lambda_vax * S_vax - gamma_a * Ia_vax,
    gamma_a * Ia_vax - gamma_l * Il_vax,
    gamma_l * Il_vax
  ) 
  
  return(list(dxdt))
}

# Set parameter values
params <- list(
  beta = 0.01,     # transmission rate
  beta_vax = 0.005, # vax transmission rate
  gamma_a = 1/3,     # acute infectious period = first 3 days 
  gamma_l = 1/7,    # late infectious period = next 7 days
  C_naive = contact_matrix_naive # contact matrix for acute & late periods
)

# Make sure contact matrices are stored as numeric matrices 
params$C_naive <- as.matrix(contact_matrix_naive)

storage.mode(params$C_naive) <- "double"

# Incorporate age structure
age_group_nums <- with(sim_population, setNames(age_group_total, age_group))

# Set initial conditions using population data
initial_infec_fraction <- 1e-5

## Assign the same proportion of infected individuals to each age group 

# Split population
S0_total <- age_group_nums * (1 - initial_infec_fraction)
Ia0_total <- age_group_nums * initial_infec_fraction

# 50% vaccinated
vax_coverage <- c(0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5) # 50% of everyone over the age of 1 is vaccinated
# Apply vaccination by age group
S0_vax <- vax_coverage * S0_total
S0     <- (1 - vax_coverage) * S0_total

Ia0_vax <- vax_coverage * Ia0_total
Ia0     <- (1 - vax_coverage) * Ia0_total

Il0      <- rep(0, 7)
Il0_vax  <- rep(0, 7)

R0       <- rep(0, 7)
R0_vax   <- rep(0, 7)

# Combine all compartments
initial_state <- c(
  S0, Ia0, Il0, R0,
  S0_vax, Ia0_vax, Il0_vax, R0_vax
)

# Check initial values 
sum(c(S0, Ia0, Il0, R0))
sum(sim_population$age_group_total)
sum(Ia0)

# Set time steps
times <- seq(0, 1000, by = 1)

# Run ODE (deterministic)
out_ODE_18 <- ode(
  y = initial_state,
  times = times,
  func = SIR_naive_18,
  parms = params
)

final_18 <- tail(out_ODE_18, 1)
R_unvax_18 <- as.numeric(final_18[23:29])
R_vax_18   <- as.numeric(final_18[51:57])
total_final_cases_18 <- sum(R_unvax_18 + R_vax_18)
print(total_final_cases_18)

# cases per 100,000
total_final_cases_18/(sum(initial_state))*100000

# ---- Extract per-age-group final cases from an ode() output ----

get_age_group_cases <- function(out_ODE_18, age_groups) {
  final_row <- tail(out_ODE_18, 1)
  
  
  # +1 offset everywhere because column 1 of out_ODE is "time"
  # Unvaccinated
  S_unvax_final  <- as.numeric(final_row[1 + (1:7)])
  Ia_unvax_final <- as.numeric(final_row[1 + (8:14)])
  Il_unvax_final <- as.numeric(final_row[1 + (15:21)])
  R_unvax_final  <- as.numeric(final_row[1 + (22:28)])
  
  # Vaccinated
  S_vax_final  <- as.numeric(final_row[1 + (29:35)])
  Ia_vax_final <- as.numeric(final_row[1 + (36:42)])
  Il_vax_final <- as.numeric(final_row[1 + (43:49)])
  R_vax_final  <- as.numeric(final_row[1 + (50:56)])
  
  # Combine vax + unvax per age group
  S_total  <- S_unvax_final  + S_vax_final
  Ia_total <- Ia_unvax_final + Ia_vax_final
  Il_total <- Il_unvax_final + Il_vax_final
  R_total  <- R_unvax_final  + R_vax_final
  
  N <- S_total + Ia_total + Il_total + R_total
  
  # Total cases = everyone no longer susceptible (captures R plus any
  # still infectious at the end of the simulation)
  cases <- N - S_total
  
  data.frame(
    age_group = age_groups,
    N = N,
    S_final = S_total,
    Ia_final = Ia_total,
    Il_final = Il_total,
    R_final = R_total,
    cases = cases
  )
}

# ---- Apply to existing run ----
cases_by_age_18 <- get_age_group_cases(out_ODE_18, age_groups)
print(cases_by_age_18)

#------------------------------------------------------------------------------#
#                              Visualize Results 
#------------------------------------------------------------------------------#

# Figure 7. Percent reduction in simulated cumulative cases within age groups targeted by vaccination
df_bar <- data.frame(
  intervention = rep(c("<1 Year", "1-4 Years", "5-17 Years", "65+ Years", ">1 Year"), each = 2),
  model = rep(c("Behavior-Dynamic", "Behavior-Naive"), 5),
  percent_reduction = c(24.5, 22.1,
                        35.2, 33.6,   
                        65.7, 28.3,   
                        27.5, 28.6,
                        85.5, 57.9)   
)

# Bar plot
ggplot(df_bar, aes(x = intervention, y = percent_reduction, fill = model)) +
  geom_col(
    position = position_dodge(width = 0.75),
    width = 0.65
  ) +
  geom_text(
    aes(label = percent_reduction),
    position = position_dodge(width = 0.75),
    vjust = -0.4,
    size = 4
  ) +
  scale_fill_manual(values = c("#969696", "#525252")) +
  labs(
    x = "Vaccinated Age Group",
    y = "Percent Reduction",
    fill = "Model"
  ) +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "top",
    axis.text.x = element_text(angle = 25, hjust = 1)
  ) +
  expand_limits(y = max(df_bar$percent_reduction) * 1.15)

# Export to high-resolution figure
ggsave("Fig7_vax_percent_reduction.png", width = 7, height = 5, dpi = 300)


# Supplementary Figure 12. Percent reduction in cumulative cases for age-targeted vaccination, by model type
df_bar <- data.frame(
  intervention = rep(c("<1 Year", "1-4 Years", "5-17 Years", "65+ Years", ">1 Year"), each = 2),
  model = rep(c("Behavior-Dynamic", "Behavior-Naive"), 5),
  percent_reduction = c(0.6, 0.5,
                        3.7, 4.3,   
                        66.4, 29.6,   
                        1.3, 2.2,
                        85.5, 58.0)   
)

# Bar plot
ggplot(df_bar, aes(x = intervention, y = percent_reduction, fill = model)) +
  geom_col(
    position = position_dodge(width = 0.75),
    width = 0.65
  ) +
  geom_text(
    aes(label = percent_reduction),
    position = position_dodge(width = 0.75),
    vjust = -0.4,
    size = 4
  ) +
  scale_fill_manual(values = c("#969696", "#525252")) +
  labs(
    x = "Vaccinated Age Group",
    y = "Percent Reduction",
    fill = "Model"
  ) +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "top",
    axis.text.x = element_text(angle = 25, hjust = 1)
  ) +
  expand_limits(y = max(df_bar$percent_reduction) * 1.15)

# Export to high-resolution figure
ggsave("Fig7_vax_percent_reduction.png", width = 7, height = 5, dpi = 300)




