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

# This code script simulates isolation-based public health interventions using
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
library(ggrepel)

# Load data 
sim_population <- read_xlsx("sim_population.xlsx")
contact_matrix_acute <- as.matrix(read_xlsx('cases_index_matrix.xlsx'))
contact_matrix_late <- as.matrix(read_xlsx('cases_1week_matrix.xlsx'))
contact_matrix_naive <- as.matrix(read_xlsx('cases_2week_matrix.xlsx'))
contact_matrix_home <- as.matrix(read_xlsx('cases_index_home_matrix.xlsx'))

# Data manipulation
## Add rownames to all matrices, if you haven't already 
rownames(contact_matrix_acute) <- c("<1", "1–4", "5–17", "18–29", "30–44", "45–64", "65+")
rownames(contact_matrix_late) <- c("<1", "1–4", "5–17", "18–29", "30–44", "45–64", "65+")
rownames(contact_matrix_naive) <- c("<1", "1–4", "5–17", "18–29", "30–44", "45–64", "65+")
rownames(contact_matrix_home) <- c("<1", "1–4", "5–17", "18–29", "30–44", "45–64", "65+")


# ********************* Age 0-4 years isolation (first 3 days) ********************* #
# NOTE: this is simulating isolation of ill individuals aged 0-4 years during acute infectious phase.

# Create contact matrices with decreased contacts only for ages 0-4.  
# The dynamic_acute contact matrix will be the SickMix respondent day 0 matrix, 
# except the row for the <1 and 1-4 groups will be from the SickMix case at home day 0 matrix.
# The dynamic_late contact matrix will be the original SickMix respondent week 1 matrix.  
# The naive_acute matrix will be the the naive matrix, except the row for the 
# <1 and 1-4 groups will be from the SickMix case at home day 0 matrix. 
# The naive_late matrix will be the original naive matrix (respondent week 2 matrix).  

C_naive_acute_0to4 <- contact_matrix_naive
C_naive_late_0to4 <- contact_matrix_naive
C_dynamic_acute_0to4 <- contact_matrix_acute
C_dynamic_late_0to4 <- contact_matrix_late

row_name_1 <- "<1"
row_name_2 <- "1–4"

C_naive_acute_0to4[row_name_1, ] <- contact_matrix_home[row_name_1, ]
C_naive_acute_0to4[row_name_2, ] <- contact_matrix_home[row_name_2, ]
C_naive_late_0to4 <- C_naive_late_0to4
C_dynamic_acute_0to4[row_name_1, ] <- contact_matrix_home[row_name_1, ]
C_dynamic_acute_0to4[row_name_2, ] <- contact_matrix_home[row_name_2, ]
C_dynamic_late_0to4 <- C_dynamic_late_0to4

#------------------------------------------------------------------------------#
#                   Model 3: SIR Model - Behavior-dynamic  
#------------------------------------------------------------------------------#

# SIR model 
SIR_dynamic_3 <- function (t, x, params) {
  # There are 7 compartments in each state to reflect the 7 age groups   
  S=x[1:7] # Susceptible 
  Ia=x[8:14] # Infectious (acute) (first 3 days of illness)
  Il=x[15:21] # Infectious (late) (next 7 days)
  R=x[22:28] # Recovered
  
  beta  = params$beta # transmission rate (same for acute and late phases)
  gamma_a = params$gamma_a # acute infectious period = first 3 days
  gamma_l = params$gamma_l # late infectious period = next 7 days
  C_dynamic_acute_0to4 = params$C_dynamic_acute_0to4 # contact matrix for acute period
  C_dynamic_late_0to4 = params$C_dynamic_late_0to4 # contact matrix for late period
  
  # Population size, N
  N = S + Ia + Il + R 
  
  # Total infectious contacts received by each recipient age group (j) 
  ## transpose contact matrices because rows are infectors (i), cols are contacts (j)
  infectious_contacts_to_j <- 
    t(C_dynamic_acute_0to4) %*% Ia + 
    t(C_dynamic_late_0to4) %*% Il 
  
  # Convert to per-capita exposure rate for contacts (divide by N_j)
  foi_sum_terms <- as.numeric(infectious_contacts_to_j) / N
  
  # Force of infection (lambda) 
  lambda <- beta * foi_sum_terms
  
  # foi by age group
  foi_vec <- beta*foi_sum_terms 
  output <- data.frame(age = c("<1", "1–4", "5–17", "18–29", "30–44", "45–64", "65+")
  )
  output$N <- N
  output$foi <- foi_vec # daily FOI
  
  # Create a set of ODEs for the SIR model
  dxdt <- c(
    - lambda * S, # S compartment
    lambda * S - gamma_a * Ia, # S -> Ia compartment
    gamma_a * Ia - gamma_l * Il, # Ia -> Il compartment
    gamma_l * Il # Il -> R compartment
  ) 
  
  return(list(dxdt))
}

# Set parameter values
params <- list(
  beta = 0.01,     # transmission rate
  gamma_a = 1/3,     # acute infectious period = first 3 days 
  gamma_l = 1/7,    # late infectious period = next 7 days
  C_dynamic_acute_0to4 = C_dynamic_acute_0to4, # contact matrix for acute period
  C_dynamic_late_0to4 = C_dynamic_late_0to4 # contact matrix for late period
)

# Make sure contact matrices are stored as numeric matrices 
params$C_dynamic_acute_0to4 <- as.matrix(C_dynamic_acute_0to4)
params$C_dynamic_late_0to4 <- as.matrix(C_dynamic_late_0to4)

storage.mode(params$C_dynamic_acute_0to4) <- "double"
storage.mode(params$C_dynamic_late_0to4) <- "double"

# Incorporate age structure
age_group_nums <- with(sim_population, setNames(age_group_total, age_group))

# Set initial conditions using population data
initial_infec_fraction <- 1e-5
## Assign the same proportion of infected individuals to each age group 
S0  <- age_group_nums * (1 - initial_infec_fraction)
Ia0 <- age_group_nums * initial_infec_fraction
Il0 <- rep(0, 7)
R0  <- rep(0, 7)

initial_state <- c(S0, Ia0, Il0, R0)

# Check initial values 
sum(c(S0, Ia0, Il0, R0))
sum(sim_population$age_group_total)
sum(Ia0)

# Set time steps
times <- seq(0, 1000, by = 1)

# Run ODE
out_ODE_3 <- ode(
  y = initial_state,
  times = times,
  func = SIR_dynamic_3,
  parms = params
)

out_ODE_3 <- as.data.frame(out_ODE_3)

# Check output
head(out_ODE_3)
tail(out_ODE_3)

final_3 <- tail(out_ODE_3, 1)
R_final_3 <- as.numeric(final_3[23:29])
total_final_cases_3 <- sum(R_final_3)
print(total_final_cases_3)

# ---- Extract per-age-group final cases from ode() output ----

# Age group labels and populations
age_groups <- c("<1", "1–4", "5–17", "18–29", "30–44", "45–64", "65+")
N_age <- age_group_nums  # vector of population per age group (same order as model)

get_age_group_cases <- function(out_ODE_3, age_groups) {
  final_row <- tail(out_ODE_3, 1)
  
  # +1 offset everywhere because column 1 of out_ODE is "time"
  S_final  <- as.numeric(final_row[1 + (1:7)])
  Ia_final <- as.numeric(final_row[1 + (8:14)])
  Il_final <- as.numeric(final_row[1 + (15:21)])
  R_final  <- as.numeric(final_row[1 + (22:28)])
  
  N <- S_final + Ia_final + Il_final + R_final
  
  # Total cases = everyone no longer susceptible (captures R plus any
  # still infectious at the end of the simulation)
  cases <- N - S_final
  
  data.frame(
    age_group = age_groups,
    N = N,
    S_final = S_final,
    Ia_final = Ia_final,
    Il_final = Il_final,
    R_final = R_final,
    cases = cases
  )
}

# ---- Apply to existing run ----
cases_by_age_3 <- get_age_group_cases(out_ODE_3, age_groups)
print(cases_by_age_3)


#------------------------------------------------------------------------------#
#                   Model 4: SIR Model - Behavior-naive   
#------------------------------------------------------------------------------#

# SIR model 
SIR_naive_4 <- function (t, x, params) {
  # There are 7 compartments in each state to reflect the 7 age groups   
  S=x[1:7] # Susceptible 
  Ia=x[8:14] # Infectious (acute) (first 3 days of illness)
  Il=x[15:21] # Infectious (late) (next 7 days)
  R=x[22:28] # Recovered
  
  beta  = params$beta # transmission rate (same for acute and late phases)
  gamma_a = params$gamma_a # acute infectious period = first 3 days
  gamma_l = params$gamma_l # late infectious period = next 7 days
  C_naive_acute_0to4 = params$C_naive_acute_0to4 # contact matrix for acute period
  C_naive_late_0to4 = params$C_naive_late_0to4 # contact matrix for late period
  
  # Population size, N
  N = S + Ia + Il + R 
  
  # Total infectious contacts received by each recipient age group (j) 
  ## transpose contact matrices because rows are infectors (i), cols are contacts (j)
  infectious_contacts_to_j <- 
    t(C_naive_acute_0to4) %*% Ia + 
    t(C_naive_late_0to4) %*% Il 
  
  # Convert to per-capita exposure rate for contacts (divide by N_j)
  foi_sum_terms <- as.numeric(infectious_contacts_to_j) / N
  
  # Force of infection (lambda) 
  lambda <- beta * foi_sum_terms
  
  # foi by age group
  foi_vec <- beta*foi_sum_terms 
  output <- data.frame(age = c("<1", "1–4", "5–17", "18–29", "30–44", "45–64", "65+")
  )
  output$N <- N
  output$foi <- foi_vec # daily FOI
  
  # Create a set of ODEs for the SIR model
  dxdt <- c(
    - lambda * S, # S compartment
    lambda * S - gamma_a * Ia, # S -> Ia compartment
    gamma_a * Ia - gamma_l * Il, # Ia -> Il compartment
    gamma_l * Il # Il -> R compartment
  ) 
  
  return(list(dxdt))
}

# Set parameter values
params <- list(
  beta = 0.01,     # transmission rate
  gamma_a = 1/3,     # acute infectious period = first 3 days 
  gamma_l = 1/7,    # late infectious period = next 7 days
  C_naive_acute_0to4 = C_naive_acute_0to4, # contact matrix for acute period
  C_naive_late_0to4 = C_naive_late_0to4 # contact matrix for late period
)

# Make sure contact matrices are stored as numeric matrices 
params$C_naive_acute_0to4 <- as.matrix(C_naive_acute_0to4)
params$C_naive_late_0to4 <- as.matrix(C_naive_late_0to4)

storage.mode(params$C_naive_acute_0to4) <- "double"
storage.mode(params$C_naive_late_0to4) <- "double"

# Incorporate age structure
age_group_nums <- with(sim_population, setNames(age_group_total, age_group))

# Set initial conditions using population data
initial_infec_fraction <- 1e-5
## Assign the same proportion of infected individuals to each age group 
S0  <- age_group_nums * (1 - initial_infec_fraction)
Ia0 <- age_group_nums * initial_infec_fraction
Il0 <- rep(0, 7)
R0  <- rep(0, 7)

initial_state <- c(S0, Ia0, Il0, R0)

# Check initial values 
sum(c(S0, Ia0, Il0, R0))
sum(sim_population$age_group_total)
sum(Ia0)

# Set time steps
times <- seq(0, 1000, by = 1)

# Run ODE
out_ODE_4 <- ode(
  y = initial_state,
  times = times,
  func = SIR_naive_4,
  parms = params
)

out_ODE_4 <- as.data.frame(out_ODE_4)

# Check output
head(out_ODE_4)
tail(out_ODE_4)

final_4 <- tail(out_ODE_4, 1)
R_final_4 <- as.numeric(final_4[23:29])
total_final_cases_4 <- sum(R_final_4)
print(total_final_cases_4)

# ---- Extract per-age-group final cases from ode() output ----

get_age_group_cases <- function(out_ODE_4, age_groups) {
  final_row <- tail(out_ODE_4, 1)
  
  # +1 offset everywhere because column 1 of out_ODE is "time"
  S_final  <- as.numeric(final_row[1 + (1:7)])
  Ia_final <- as.numeric(final_row[1 + (8:14)])
  Il_final <- as.numeric(final_row[1 + (15:21)])
  R_final  <- as.numeric(final_row[1 + (22:28)])
  
  N <- S_final + Ia_final + Il_final + R_final
  
  # Total cases = everyone no longer susceptible (captures R plus any
  # still infectious at the end of the simulation)
  cases <- N - S_final
  
  data.frame(
    age_group = age_groups,
    N = N,
    S_final = S_final,
    Ia_final = Ia_final,
    Il_final = Il_final,
    R_final = R_final,
    cases = cases
  )
}

# ---- Apply to existing run ----
cases_by_age_4 <- get_age_group_cases(out_ODE_4, age_groups)
print(cases_by_age_4)


# ********************* Age 5-17 years isolation (first 3 days) only ********************* #
# NOTE: this is simulating isolation of ill individuals aged 5-17 years during acute infectious phase.

# Create contact matrices with decreased contacts only for ages 5-17.  
# The dynamic_acute contact matrix will be the SickMix respondent day 0 matrix, 
# except the row for the 5-17 group will be from the SickMix case at home day 0 matrix.
# The dynamic_late contact matrix will be the original SickMix respondent week 1 matrix.  
# The naive_acute matrix will be the the naive matrix, except the row for the 
# 5-17 group will be from the SickMix case at home day 0 matrix. 
# The naive_late matrix will be the regular naive matrix (week 2 matrix).  

C_naive_acute_5to17 <- contact_matrix_naive
C_naive_late_5to17 <- contact_matrix_naive
C_dynamic_acute_5to17 <- contact_matrix_acute
C_dynamic_late_5to17 <- contact_matrix_late

row_name <- "5–17"

C_naive_acute_5to17[row_name, ] <- contact_matrix_home[row_name, ]
C_naive_late_5to17 <- C_naive_late_5to17
C_dynamic_acute_5to17[row_name, ] <- contact_matrix_home[row_name, ]
C_dynamic_late_5to17 <- C_dynamic_late_5to17

#------------------------------------------------------------------------------#
#                   Model 5: SIR Model - Behavior-dynamic   
#------------------------------------------------------------------------------#

# SIR model 
SIR_dynamic_5 <- function (t, x, params) {
  # There are 7 compartments in each state to reflect the 7 age groups   
  S=x[1:7] # Susceptible 
  Ia=x[8:14] # Infectious (acute) (first 3 days of illness)
  Il=x[15:21] # Infectious (late) (next 7 days)
  R=x[22:28] # Recovered
  
  beta  = params$beta # transmission rate (same for acute and late phases)
  gamma_a = params$gamma_a # acute infectious period = first 3 days
  gamma_l = params$gamma_l # late infectious period = next 7 days
  C_dynamic_acute_5to17 = params$C_dynamic_acute_5to17 # contact matrix for acute period
  C_dynamic_late_5to17 = params$C_dynamic_late_5to17 # contact matrix for late period
  
  # Population size, N
  N = S + Ia + Il + R 
  
  # Total infectious contacts received by each recipient age group (j) 
  ## transpose contact matrices because rows are infectors (i), cols are contacts (j)
  infectious_contacts_to_j <- 
    t(C_dynamic_acute_5to17) %*% Ia + 
    t(C_dynamic_late_5to17) %*% Il 
  
  # Convert to per-capita exposure rate for contacts (divide by N_j)
  foi_sum_terms <- as.numeric(infectious_contacts_to_j) / N
  
  # Force of infection (lambda) 
  lambda <- beta * foi_sum_terms
  
  # foi by age group
  foi_vec <- beta*foi_sum_terms 
  output <- data.frame(age = c("<1", "1–4", "5–17", "18–29", "30–44", "45–64", "65+")
  )
  output$N <- N
  output$foi <- foi_vec # daily FOI
  
  # Create a set of ODEs for the SIR model
  dxdt <- c(
    - lambda * S, # S compartment
    lambda * S - gamma_a * Ia, # S -> Ia compartment
    gamma_a * Ia - gamma_l * Il, # Ia -> Il compartment
    gamma_l * Il # Il -> R compartment
  ) 
  
  return(list(dxdt))
}

# Set parameter values
params <- list(
  beta = 0.01,     # transmission rate
  gamma_a = 1/3,     # acute infectious period = first 3 days 
  gamma_l = 1/7,    # late infectious period = next 7 days
  C_dynamic_acute_5to17 = C_dynamic_acute_5to17, # contact matrix for acute period
  C_dynamic_late_5to17 = C_dynamic_late_5to17 # contact matrix for late period
)

# Make sure contact matrices are stored as numeric matrices 
params$C_dynamic_acute_5to17 <- as.matrix(C_dynamic_acute_5to17)
params$C_dynamic_late_5to17 <- as.matrix(C_dynamic_late_5to17)

storage.mode(params$C_dynamic_acute_5to17) <- "double"
storage.mode(params$C_dynamic_late_5to17) <- "double"

# Incorporate age structure
age_group_nums <- with(sim_population, setNames(age_group_total, age_group))

# Set initial conditions using population data
initial_infec_fraction <- 1e-5
## Assign the same proportion of infected individuals to each age group 
S0  <- age_group_nums * (1 - initial_infec_fraction)
Ia0 <- age_group_nums * initial_infec_fraction
Il0 <- rep(0, 7)
R0  <- rep(0, 7)

initial_state <- c(S0, Ia0, Il0, R0)

# Check initial values 
sum(c(S0, Ia0, Il0, R0))
sum(sim_population$age_group_total)
sum(Ia0)

# Set time steps
times <- seq(0, 1000, by = 1)

# Run ODE
out_ODE_5 <- ode(
  y = initial_state,
  times = times,
  func = SIR_dynamic_5,
  parms = params
)

out_ODE_5 <- as.data.frame(out_ODE_5)

# Check output
head(out_ODE_5)
tail(out_ODE_5)

final_5 <- tail(out_ODE_5, 1)
R_final_5 <- as.numeric(final_5[23:29])
total_final_cases_5 <- sum(R_final_5)
print(total_final_cases_5)

# ---- Extract per-age-group final cases from ode() output ----

get_age_group_cases <- function(out_ODE_5, age_groups) {
  final_row <- tail(out_ODE_5, 1)
  
  # +1 offset everywhere because column 1 of out_ODE is "time"
  S_final  <- as.numeric(final_row[1 + (1:7)])
  Ia_final <- as.numeric(final_row[1 + (8:14)])
  Il_final <- as.numeric(final_row[1 + (15:21)])
  R_final  <- as.numeric(final_row[1 + (22:28)])
  
  N <- S_final + Ia_final + Il_final + R_final
  
  # Total cases = everyone no longer susceptible (captures R plus any
  # still infectious at the end of the simulation)
  cases <- N - S_final
  
  data.frame(
    age_group = age_groups,
    N = N,
    S_final = S_final,
    Ia_final = Ia_final,
    Il_final = Il_final,
    R_final = R_final,
    cases = cases
  )
}

# ---- Apply to existing run ----
cases_by_age_5 <- get_age_group_cases(out_ODE_5, age_groups)
print(cases_by_age_5)

sum(cases_by_age_5$cases)
total_final_cases_5

#------------------------------------------------------------------------------#
#                   Model 6: SIR Model - Behavior-naive   
#------------------------------------------------------------------------------#

# SIR model 
SIR_naive_6 <- function (t, x, params) {
  # There are 7 compartments in each state to reflect the 7 age groups   
  S=x[1:7] # Susceptible 
  Ia=x[8:14] # Infectious (acute) (first 3 days of illness)
  Il=x[15:21] # Infectious (late) (next 7 days)
  R=x[22:28] # Recovered
  
  beta  = params$beta # transmission rate (same for acute and late phases)
  gamma_a = params$gamma_a # acute infectious period = first 3 days
  gamma_l = params$gamma_l # late infectious period = next 7 days
  C_naive_acute_5to17 = params$C_naive_acute_5to17 # contact matrix for acute periods
  C_naive_late_5to17 = params$C_naive_late_5to17 # contact matrix for late period
  
  # Population size, N
  N = S + Ia + Il + R 
  
  # Total infectious contacts received by each recipient age group (j) 
  ## transpose contact matrices because rows are infectors (i), cols are contacts (j)
  infectious_contacts_to_j <- 
    t(C_naive_acute_5to17) %*% Ia + 
    t(C_naive_late_5to17) %*% Il 
  
  # Convert to per-capita exposure rate for contacts (divide by N_j)
  foi_sum_terms <- as.numeric(infectious_contacts_to_j) / N
  
  # Force of infection (lambda) 
  lambda <- beta * foi_sum_terms
  
  # foi by age group
  foi_vec <- beta*foi_sum_terms 
  output <- data.frame(age = c("<1", "1–4", "5–17", "18–29", "30–44", "45–64", "65+")
  )
  output$N <- N
  output$foi <- foi_vec # daily FOI
  
  # Create a set of ODEs for the SIR model
  dxdt <- c(
    - lambda * S, # S compartment
    lambda * S - gamma_a * Ia, # S -> Ia compartment
    gamma_a * Ia - gamma_l * Il, # Ia -> Il compartment
    gamma_l * Il # Il -> R compartment
  ) 
  
  return(list(dxdt))
}

# Set parameter values
params <- list(
  beta = 0.01,     # transmission rate
  gamma_a = 1/3,     # acute infectious period = first 3 days 
  gamma_l = 1/7,    # late infectious period = next 7 days
  C_naive_acute_5to17 = C_naive_acute_5to17, # contact matrix for acute period
  C_naive_late_5to17 = C_naive_late_5to17 # contact matrix for late period 
)

# Make sure contact matrices are stored as numeric matrices 
params$C_naive_acute_5to17 <- as.matrix(C_naive_acute_5to17)
params$C_naive_late_5to17 <- as.matrix(C_naive_late_5to17)

storage.mode(params$C_naive_acute_5to17) <- "double"
storage.mode(params$C_naive_late_5to17) <- "double"

# Incorporate age structure
age_group_nums <- with(sim_population, setNames(age_group_total, age_group))

# Set initial conditions using population data
initial_infec_fraction <- 1e-5
## Assign the same proportion of infected individuals to each age group 
S0  <- age_group_nums * (1 - initial_infec_fraction)
Ia0 <- age_group_nums * initial_infec_fraction
Il0 <- rep(0, 7)
R0  <- rep(0, 7)

initial_state <- c(S0, Ia0, Il0, R0)

# Check initial values 
sum(c(S0, Ia0, Il0, R0))
sum(sim_population$age_group_total)
sum(Ia0)

# Set time steps
times <- seq(0, 1000, by = 1)

# Run ODE 
out_ODE_6 <- ode(
  y = initial_state,
  times = times,
  func = SIR_naive_6,
  parms = params
)

out_ODE_6 <- as.data.frame(out_ODE_6)

# Check output
head(out_ODE_6)
tail(out_ODE_6)

final_6 <- tail(out_ODE_6, 1)
R_final_6 <- as.numeric(final_6[23:29])
total_final_cases_6 <- sum(R_final_6)
print(total_final_cases_6)

# ---- Extract per-age-group final cases from ode() output ----

get_age_group_cases <- function(out_ODE_6, age_groups) {
  final_row <- tail(out_ODE_6, 1)
  
  # +1 offset everywhere because column 1 of out_ODE is "time"
  S_final  <- as.numeric(final_row[1 + (1:7)])
  Ia_final <- as.numeric(final_row[1 + (8:14)])
  Il_final <- as.numeric(final_row[1 + (15:21)])
  R_final  <- as.numeric(final_row[1 + (22:28)])
  
  N <- S_final + Ia_final + Il_final + R_final
  
  # Total cases = everyone no longer susceptible (captures R plus any
  # still infectious at the end of the simulation)
  cases <- N - S_final
  
  data.frame(
    age_group = age_groups,
    N = N,
    S_final = S_final,
    Ia_final = Ia_final,
    Il_final = Il_final,
    R_final = R_final,
    cases = cases
  )
}

# ---- Apply to existing run ----
cases_by_age_6 <- get_age_group_cases(out_ODE_6, age_groups)
print(cases_by_age_6)

sum(cases_by_age_6$cases)
total_final_cases_6

# ***************** 65+ isolation - 3 days only ***************** #

# Create contact matrices with decreased contacts only for age 65+.  
# The dynamic acute contact matrix will be the SickMix respondent day 0 matrix, 
# except the row for the 65+ group will be from the SickMix case at home day 0 matrix.
# The dynamic late contact matrix will be the SickMix respondent week 1 matrix.  
# The naive acute matrix will be the naive matrix except the row for the 65+ group 
# will be from the SickMix case at home day 0 matrix. 
# The naive late matrix will be the original naive matrix.  

C_naive_acute_over65 <- contact_matrix_naive
C_naive_late_over65 <- contact_matrix_naive
C_dynamic_acute_over65 <- contact_matrix_acute
C_dynamic_late_over65 <- contact_matrix_late

row_name <- "65+"

C_naive_acute_over65[row_name, ] <- contact_matrix_home[row_name, ]
C_naive_late_over65 <- C_naive_late_over65
C_dynamic_acute_over65[row_name, ] <- contact_matrix_home[row_name, ]
C_dynamic_late_over65 <- C_dynamic_late_over65

#------------------------------------------------------------------------------#
#                   Model 7: SIR Model - Behavior-dynamic   
#------------------------------------------------------------------------------#

# SIR model 
SIR_dynamic_7 <- function (t, x, params) {
  # There are 7 compartments in each state to reflect the 7 age groups   
  S=x[1:7] # Susceptible 
  Ia=x[8:14] # Infectious (acute) (first 3 days of illness)
  Il=x[15:21] # Infectious (late) (next 7 days)
  R=x[22:28] # Recovered
  
  beta  = params$beta # transmission rate (same for acute and late phases)
  gamma_a = params$gamma_a # acute infectious period = first 3 days
  gamma_l = params$gamma_l # late infectious period = next 7 days
  C_dynamic_acute_over65 = params$C_dynamic_acute_over65 # contact matrix for acute period
  C_dynamic_late_over65 = params$C_dynamic_late_over65 # contact matrix for late period
  
  # Population size, N
  N = S + Ia + Il + R 
  
  # Total infectious contacts received by each recipient age group (j) 
  ## transpose contact matrices because rows are infectors (i), cols are contacts (j)
  infectious_contacts_to_j <- 
    t(C_dynamic_acute_over65) %*% Ia + 
    t(C_dynamic_late_over65) %*% Il 
  
  # Convert to per-capita exposure rate for contacts (divide by N_j)
  foi_sum_terms <- as.numeric(infectious_contacts_to_j) / N
  
  # Force of infection (lambda) 
  lambda <- beta * foi_sum_terms
  
  # foi by age group
  foi_vec <- beta*foi_sum_terms 
  output <- data.frame(age = c("<1", "1–4", "5–17", "18–29", "30–44", "45–64", "65+")
  )
  output$N <- N
  output$foi <- foi_vec # daily FOI
  
  # Create a set of ODEs for the SIR model
  dxdt <- c(
    - lambda * S, # S compartment
    lambda * S - gamma_a * Ia, # S -> Ia compartment
    gamma_a * Ia - gamma_l * Il, # Ia -> Il compartment
    gamma_l * Il # Il -> R compartment
  ) 
  
  return(list(dxdt))
}

# Set parameter values
params <- list(
  beta = 0.01,     # transmission rate
  gamma_a = 1/3,     # acute infectious period = first 3 days 
  gamma_l = 1/7,    # late infectious period = next 7 days
  C_dynamic_acute_over65 = C_dynamic_acute_over65, # contact matrix for acute period
  C_dynamic_late_over65 = C_dynamic_late_over65 # contact matrix for late period
)

# Make sure contact matrices are stored as numeric matrices 
params$C_dynamic_acute_over65 <- as.matrix(C_dynamic_acute_over65)
params$C_dynamic_late_over65 <- as.matrix(C_dynamic_late_over65)

storage.mode(params$C_dynamic_acute_over65) <- "double"
storage.mode(params$C_dynamic_late_over65) <- "double"

# Incorporate age structure
age_group_nums <- with(sim_population, setNames(age_group_total, age_group))

# Set initial conditions using population data
initial_infec_fraction <- 1e-5
## Assign the same proportion of infected individuals to each age group 
S0  <- age_group_nums * (1 - initial_infec_fraction)
Ia0 <- age_group_nums * initial_infec_fraction
Il0 <- rep(0, 7)
R0  <- rep(0, 7)

initial_state <- c(S0, Ia0, Il0, R0)

# Check initial values 
sum(c(S0, Ia0, Il0, R0))
sum(sim_population$age_group_total)
sum(Ia0)

# Set time steps
times <- seq(0, 1000, by = 1)

# Run ODE
out_ODE_7 <- ode(
  y = initial_state,
  times = times,
  func = SIR_dynamic_7,
  parms = params
)

out_ODE_7 <- as.data.frame(out_ODE_7)

# Check output
head(out_ODE_7)
tail(out_ODE_7)

final_7 <- tail(out_ODE_7, 1)
R_final_7 <- as.numeric(final_7[23:29])
total_final_cases_7 <- sum(R_final_7)
print(total_final_cases_7)

# ---- Extract per-age-group final cases from ode() output ----

get_age_group_cases <- function(out_ODE_7, age_groups) {
  final_row <- tail(out_ODE_7, 1)
  
  # +1 offset everywhere because column 1 of out_ODE is "time"
  S_final  <- as.numeric(final_row[1 + (1:7)])
  Ia_final <- as.numeric(final_row[1 + (8:14)])
  Il_final <- as.numeric(final_row[1 + (15:21)])
  R_final  <- as.numeric(final_row[1 + (22:28)])
  
  N <- S_final + Ia_final + Il_final + R_final
  
  # Total cases = everyone no longer susceptible (captures R plus any
  # still infectious at the end of the simulation)
  cases <- N - S_final
  
  data.frame(
    age_group = age_groups,
    N = N,
    S_final = S_final,
    Ia_final = Ia_final,
    Il_final = Il_final,
    R_final = R_final,
    cases = cases
  )
}

# ---- Apply to existing run ----
cases_by_age_7 <- get_age_group_cases(out_ODE_7, age_groups)
print(cases_by_age_7)

sum(cases_by_age_7$cases)
total_final_cases_7

#------------------------------------------------------------------------------#
#                   Model 8: SIR Model - Behavior-naive   
#------------------------------------------------------------------------------#

# SIR model 
SIR_naive_8 <- function (t, x, params) {
  # There are 7 compartments in each state to reflect the 7 age groups   
  S=x[1:7] # Susceptible 
  Ia=x[8:14] # Infectious (acute) (first 3 days of illness)
  Il=x[15:21] # Infectious (late) (next 7 days)
  R=x[22:28] # Recovered
  
  beta  = params$beta # transmission rate (same for acute and late phases)
  gamma_a = params$gamma_a # acute infectious period = first 3 days
  gamma_l = params$gamma_l # late infectious period = next 7 days
  C_naive_acute_over65 = params$C_naive_acute_over65 # contact matrix for acute period
  C_naive_late_over65 = params$C_naive_late_over65 # contact matrix for late period
  
  # Population size, N
  N = S + Ia + Il + R 
  
  # Total infectious contacts received by each recipient age group (j) 
  ## transpose contact matrices because rows are infectors (i), cols are contacts (j)
  infectious_contacts_to_j <- 
    t(C_naive_acute_over65) %*% Ia + 
    t(C_naive_late_over65) %*% Il 
  
  # Convert to per-capita exposure rate for contacts (divide by N_j)
  foi_sum_terms <- as.numeric(infectious_contacts_to_j) / N
  
  # Force of infection (lambda) 
  lambda <- beta * foi_sum_terms
  
  # foi by age group
  foi_vec <- beta*foi_sum_terms 
  output <- data.frame(age = c("<1", "1–4", "5–17", "18–29", "30–44", "45–64", "65+")
  )
  output$N <- N
  output$foi <- foi_vec # daily FOI
  
  # Create a set of ODEs for the SIR model
  dxdt <- c(
    - lambda * S, # S compartment
    lambda * S - gamma_a * Ia, # S -> Ia compartment
    gamma_a * Ia - gamma_l * Il, # Ia -> Il compartment
    gamma_l * Il # Il -> R compartment
  ) 
  
  return(list(dxdt))
}

# Set parameter values
params <- list(
  beta = 0.01,     # transmission rate
  gamma_a = 1/3,     # acute infectious period = first 3 days 
  gamma_l = 1/7,    # late infectious period = next 7 days
  C_naive_acute_over65 = C_naive_acute_over65, # contact matrix for acute period
  C_naive_late_over65 = C_naive_late_over65 # contact matrix for late period
)

# Make sure contact matrices are stored as numeric matrices 
params$C_naive_acute_over65 <- as.matrix(C_naive_acute_over65)
params$C_naive_late_over65 <- as.matrix(C_naive_late_over65)

storage.mode(params$C_naive_acute_over65) <- "double"
storage.mode(params$C_naive_late_over65) <- "double"

# Incorporate age structure
age_group_nums <- with(sim_population, setNames(age_group_total, age_group))

# Set initial conditions using population data
initial_infec_fraction <- 1e-5
## Assign the same proportion of infected individuals to each age group 
S0  <- age_group_nums * (1 - initial_infec_fraction)
Ia0 <- age_group_nums * initial_infec_fraction
Il0 <- rep(0, 7)
R0  <- rep(0, 7)

initial_state <- c(S0, Ia0, Il0, R0)

# Check initial values 
sum(c(S0, Ia0, Il0, R0))
sum(sim_population$age_group_total)
sum(Ia0)

# Set time steps
times <- seq(0, 1000, by = 1)

# Run ODE
out_ODE_8 <- ode(
  y = initial_state,
  times = times,
  func = SIR_naive_8,
  parms = params
)

out_ODE_8 <- as.data.frame(out_ODE_8)

# Check output
head(out_ODE_8)
tail(out_ODE_8)

final_8 <- tail(out_ODE_8, 1)
R_final_8 <- as.numeric(final_8[23:29])
total_final_cases_8 <- sum(R_final_8)
print(total_final_cases_8)

# ---- Extract per-age-group final cases from ode() output ----

get_age_group_cases <- function(out_ODE_8, age_groups) {
  final_row <- tail(out_ODE_8, 1)
  
  # +1 offset everywhere because column 1 of out_ODE is "time"
  S_final  <- as.numeric(final_row[1 + (1:7)])
  Ia_final <- as.numeric(final_row[1 + (8:14)])
  Il_final <- as.numeric(final_row[1 + (15:21)])
  R_final  <- as.numeric(final_row[1 + (22:28)])
  
  N <- S_final + Ia_final + Il_final + R_final
  
  # Total cases = everyone no longer susceptible (captures R plus any
  # still infectious at the end of the simulation)
  cases <- N - S_final
  
  data.frame(
    age_group = age_groups,
    N = N,
    S_final = S_final,
    Ia_final = Ia_final,
    Il_final = Il_final,
    R_final = R_final,
    cases = cases
  )
}

# ---- Apply to existing run ----
cases_by_age_8 <- get_age_group_cases(out_ODE_8, age_groups)
print(cases_by_age_8)

sum(cases_by_age_8$cases)
total_final_cases_8

#------------------------------------------------------------------------------#
#                             Visualize Results   
#------------------------------------------------------------------------------#

# Figure 6. Percent reduction in cases among target age groups due to isolation interventions, by model type
df_bar <- data.frame(
  intervention = rep(c("0-4 Years", "5-17 Years", "65+ Years"), each = 2),
  model = rep(c("Behavior-Dynamic", "Behavior-Naive"), 3),
  percent_reduction = c(8.7, 12.2,   
                        9.4, 11.1,   
                        1.3, 3.2)   
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
  scale_fill_manual(values = c("Behavior-Dynamic" = "#969696", "Behavior-Naive" = "#525252")) +
  labs(
    x = "Intervention Target Age Group",
    y = "Percent Reduction",
    fill = "Model"
  ) +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "top",
    axis.text.x = element_text(angle = 25, hjust = 1)
  ) +
  expand_limits(y = max(df_bar$percent_reduction) * 1.15)

ggsave("Figure6_isolation_interventions.png", width = 7, height = 5, dpi = 300)


# Supplementary Figure 11. Percent reduction in cumulative cases for age-targeted 
# isolation interventions, by model type
df_bar <- data.frame(
  intervention = rep(c("0-4 Years", "5-17 Years", "65+ Years"), each = 2),
  model = rep(c("Behavior-Dynamic", "Behavior-Naive"), 3),
  percent_reduction = c(1.4, 2.3,   
                        12.6, 16.5,   
                        0.3, 0.7)   
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
  scale_fill_manual(values = c("Behavior-Dynamic" = "#969696", "Behavior-Naive" = "#525252")) +
  labs(
    x = "Intervention Target Age Group",
    y = "Percent Reduction",
    fill = "Model"
  ) +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "top",
    axis.text.x = element_text(angle = 25, hjust = 1)
  ) +
  expand_limits(y = max(df_bar$percent_reduction) * 1.15)

ggsave("SFig11_isolation_interventions.png", width = 7, height = 5, dpi = 300)

