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

# This code script creates disease-agnostic transmission models that incorporate 
# age-structured contact matrices. 

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
pop_age_dist <- read_xlsx('United_States_country_level_age_distribution_85.xlsx')

#------------------------------------------------------------------------------#
#  Data Manipulation
#------------------------------------------------------------------------------#

# Calculate U.S. population age distributions using data from Mistry et al. 
# (https://github.com/mobs-lab/mixing-patterns/tree/main/data/age_distributions)
pop_age_dist$Age <- as.numeric(pop_age_dist$Age)

age_dist <- pop_age_dist %>%
  mutate(
    age_group = case_when(
      Age < 1  ~ "<1",
      Age >= 1  & Age <= 4  ~ "1–4",
      Age >= 5  & Age <= 17 ~ "5–17",
      Age >= 18 & Age <= 29 ~ "18–29",
      Age >= 30 & Age <= 44 ~ "30–44",
      Age >= 45 & Age <= 64 ~ "45–64",
      Age >= 65              ~ "65+",
      TRUE ~ NA_character_
    )
  ) %>%
  group_by(age_group) %>%
  summarise(
    people = sum(as.numeric(Count), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    proportion = people / sum(people),
    age_group = factor(
      age_group,
      levels = c("<1", "1–4", "5–17", "18–29", "30–44", "45–64", "65+")
    )
  ) %>%
  arrange(age_group)
# Check
sum(age_dist$proportion)
## proportions sum to 1 

# Create a simulated population data frame with the actual U.S. age distribution
## 2021 census = 331449281
sim_population <- age_dist %>%
  mutate(
    age_group_total = round(proportion * 331449281)
  ) %>%
  mutate(age_group_total = age_group_total +
           (331449281 - sum(age_group_total)) * (age_group_total == max(age_group_total))
  ) %>%
  dplyr::select(-2) # drop second column which is num individuals from Mistry data frame

# Save 
write_xlsx(sim_population, "sim_population.xlsx")

#------------------------------------------------------------------------------#
#  Load contact matrices
#------------------------------------------------------------------------------#
# NOTE: THE CONTACT MATRICES ARE CONSTRUCED SUCH THAT 
# ROW = RESPONDENT/INFECTED AND COLUMN = CONTACT/SUSCEPTIBLE

contact_matrix_acute <- read_xlsx("contact_matrix_acute.xlsx")
contact_matrix_late <- read_xlsx("contact_matrix_late.xlsx")
contact_matrix_naive <- read_xlsx("contact_matrix_naive.xlsx")
contact_matrix_home <- read_xlsx("contact_matrix_home.xlsx")

#------------------------------------------------------------------------------#
#                     SIR Model - Behavior-dynamic model 
#------------------------------------------------------------------------------#

SIR_behavior_dynamic <- function (t, x, params) {
  # There are 7 compartments in each state to reflect the 7 age groups   
  S=x[1:7] # Susceptible 
  Ia=x[8:14] # Infectious (acute) (first 3 days of illness)
  Il=x[15:21] # Infectious (late) (next 7 days)
  R=x[22:28] # Recovered
  
  beta  = params$beta # transmission rate (same for acute and late phases)
  gamma_a = params$gamma_a # acute infectious period = first 3 days
  gamma_l = params$gamma_l # late infectious period = next 7 days
  C_acute = params$C_acute # contact matrix for acute period
  C_late = params$C_late  # contact matrix for late period 
  
  # Population size, N
  N = S + Ia + Il + R 
  
  # Total infectious contacts received by each recipient age group (j) 
  ## transpose contact matrices because rows are infectors (i), cols are contacts (j)
  infectious_contacts_to_j <- 
    t(C_acute) %*% Ia + 
    t(C_late) %*% Il 
  
  # Convert to per-capita exposure rate for contacts (divide by N_j)
  foi_sum_terms <- as.numeric(infectious_contacts_to_j) / N
  
  # Force of infection (lambda) 
  lambda <- beta * foi_sum_terms
  
  # FOI by age group
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
out_ODE <- ode(
  y = initial_state,
  times = times,
  func = SIR_behavior_dynamic,
  parms = params
)

out_ODE <- as.data.frame(out_ODE)

# Check output
head(out_ODE)
tail(out_ODE)

final_1 <- tail(out_ODE, 1)
R_final_1 <- as.numeric(final_1[23:29])
total_final_cases_1 <- sum(R_final_1)
print(total_final_cases_1)

# cases per 100,000
total_final_cases_1/(sum(initial_state))*100000

# ==================
### Summary Statistics ###
# ==================

# ---- Extract per-age-group final cases from ode() output ----

# Age group labels and populations
age_groups <- c("<1", "1–4", "5–17", "18–29", "30–44", "45–64", "65+")
N_age <- age_group_nums  # vector of population per age group (same order as model)

get_age_group_cases <- function(out_ODE, age_groups) {
  final_row <- tail(out_ODE, 1)
  
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
cases_by_age_1 <- get_age_group_cases(out_ODE, age_groups)
print(cases_by_age_1)

# ==================
### Plot Results ###
# ==================

#----- Time series plot for infected individuals per 100,000 by age group -----#

## Supplementary Figure 10. Age-stratified time series plots for infected individuals 
## per 100,000 population, by model type

## Extract Ia and Il matrices
Ia_mat <- as.matrix(out_ODE[, 9:15])
Il_mat <- as.matrix(out_ODE[, 16:22])
total_I <- rowSums(out_ODE[, 9:22]) # Ia + Il 

# Age group populations
pop_vec <- age_group_nums

# Convert infections to per 100k
I_per100k <- sweep(Ia_mat + Il_mat, 2, pop_vec, "/") * 100000

colnames(I_per100k) <- names(age_group_nums)

I_by_age <- as.data.frame(I_per100k)
I_by_age$time <- out_ODE$time

# Convert to long format
df_age <- I_by_age %>%
  pivot_longer(
    cols = -time,
    names_to = "age_group",
    values_to = "infected_per100k"
  )

# Compute peak infections 
peak_df <- df_age %>%
  group_by(age_group) %>%
  filter(infected_per100k == max(infected_per100k)) %>%
  slice(1) %>%
  ungroup()

df_age <- df_age %>%
  mutate(age_group = factor(age_group, 
                            levels = c("<1", "1–4", "5–17", "18–29", "30–44", "45–64", "65+")))

# Sort peaks by magnitude
peak_df <- peak_df %>%
  arrange(infected_per100k)

# Plot
p <- ggplot(df_age, aes(x = time, y = infected_per100k, color = age_group)) +
  
  geom_line(linewidth = 1.2) +
  
  scale_color_manual(
    values = c(
      "<1" = "#E69F00" , 
      "1–4" = "#56B4E9", 
      "5–17" = "#009E73", 
      "18–29" = "#C4A000",
      "30–44" = "#0072B2", 
      "45–64" = "#D55E00", 
      "65+" = "#CC79A7"
    ),
    labels = c(
      "<1" = "<1 Year",
      "1–4" = "1–4 Years",
      "5–17" = "5–17 Years",
      "18–29" = "18–29 Years",
      "30–44" = "30–44 Years",
      "45–64" = "45–64 Years",
      "65+" = "65+ Years"
    )
  ) +
  
  labs(
    title = NULL,
    x = "Days",
    y = "Active Infections per 100,000",
    color = "Age Group"
  ) +
  
  scale_y_continuous(
    breaks = c(2500, 5000, 7500, 10000, 12500),
    expand = expansion(mult = c(0.02, 0.08))
    ) +  
  
  coord_cartesian(xlim = c(0, 300)) +
  
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "right",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"
  )
  )

print(p)

# Export to high-resolution figure
ggsave("SFig10_timeseries_age_dynamic.png", width = 7, height = 5, dpi = 300)

# ====================================================
### Next generation matrix and dominant eigenvalue ###
# ====================================================

# NOTE: C_acute and C_late are per-capita contact rates,
# so no additional division by population size is required
# Contact matrices were defined as per-capita rates; therefore, 
# population size normalization is implicit in the next-generation matrix.

# Extract parameters
beta    <- params$beta
gamma_a <- params$gamma_a
gamma_l <- params$gamma_l

C_acute <- params$C_acute
C_late  <- params$C_late

N_age <- age_group_nums

K_dynamic <- beta * (
  sweep(C_acute, 2, N_age, "/") / gamma_a +
    sweep(C_late,  2, N_age, "/") / gamma_l
)

# Build next generation matrix
K_dynamic <- beta * (
  t(C_acute) / gamma_a +
    t(C_late)  / gamma_l
)

# Check dimensions 
dim(K_dynamic)
# should be 7 x 7

# Compute the dominant eigenvalue (=R0) 
eig_dynamic <- eigen(K_dynamic)

R0_dynamic <- max(Re(eig_dynamic$values))
R0_dynamic


#------------------------------------------------------------------------------#
#                        SIIR Model - Behavior-naive model
#------------------------------------------------------------------------------#

SIR_behavior_naive <- function (t, x, params) {
  # There are 7 compartments in each state to reflect the 7 age groups   
  S=x[1:7] # Susceptible 
  Ia=x[8:14] # Infectious (acute) (first 3 days of illness)
  Il=x[15:21] # Infectious (late) (next 7 days)
  R=x[22:28] # Recovered
  
  beta  = params$beta # transmission rate (same for acute and late phases)
  gamma_a = params$gamma_a # acute infectious period = first 3 days
  gamma_l = params$gamma_l # late infectious period = next 7 days
  C_naive = params$C_naive # contact matrix for both periods
  
  # Population size, N
  N = S + Ia + Il + R 
  
  # Total infectious contacts received by each recipient age group (j) 
  ## transpose contact matrices because rows are infectors (i), cols are contacts (j)
  infectious_contacts_to_j <- 
    t(C_naive) %*% Ia + 
    t(C_naive) %*% Il 
  
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
  C_naive = contact_matrix_naive # contact matrix for both periods
)

# Make sure contact matrices are stored as numeric matrices 
params$C_naive <- as.matrix(contact_matrix_naive)

storage.mode(params$C_naive) <- "double"

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
out_ODE_2 <- ode(
  y = initial_state,
  times = times,
  func = SIR_behavior_naive,
  parms = params
)

out_ODE_2 <- as.data.frame(out_ODE_2)

# Check output
head(out_ODE_2)
tail(out_ODE_2)

final_2 <- tail(out_ODE_2, 1)
R_final_2 <- as.numeric(final_2[23:29])
total_final_cases_2 <- sum(R_final_2)
print(total_final_cases_2)

# cases per 100,000
total_final_cases_2/(sum(initial_state))*100000

# ==================
### Summary Statistics ###
# ==================

# ---- Extract per-age-group final cases from an ode() output ----

get_age_group_cases <- function(out_ODE_2, age_groups) {
  final_row <- tail(out_ODE_2, 1)
  
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
cases_by_age_2 <- get_age_group_cases(out_ODE_2, age_groups)
print(cases_by_age_2)


# ==================
### Plot Results ###
# ==================

#----- Time series plot for infected individuals per 100,000 by age group -----#
## Supplementary Figure 10. Age-stratified time series plots for infected individuals 
## per 100,000 population, by model type

## Extract Ia and Il matrices
Ia_mat_2 <- as.matrix(out_ODE_2[, 9:15])
Il_mat_2 <- as.matrix(out_ODE_2[, 16:22])
total_I_2 <- rowSums(out_ODE_2[, 9:22]) # Ia + Il 

# Age group populations
pop_vec <- age_group_nums

# Convert infections to per 100k
I_per100k_2 <- sweep(Ia_mat_2 + Il_mat_2, 2, pop_vec, "/") * 100000

colnames(I_per100k_2) <- names(age_group_nums)

I_by_age_2 <- as.data.frame(I_per100k_2)
I_by_age_2$time <- out_ODE_2$time

# Convert to long format
df_age_2 <- I_by_age_2 %>%
  pivot_longer(
    cols = -time,
    names_to = "age_group",
    values_to = "infected_per100k"
  )

df_age_2 <- df_age_2 %>%
  mutate(age_group = factor(age_group, 
                            levels = c("<1", "1–4", "5–17", "18–29", "30–44", "45–64", "65+")))

# Compute peak infections 
peak_df_2 <- df_age_2 %>%
  group_by(age_group) %>%
  filter(infected_per100k == max(infected_per100k)) %>%
  slice(1) %>%
  ungroup()

# Sort peaks by magnitude
peak_df_2 <- peak_df_2 %>%
  arrange(infected_per100k)

# Plot
p2 <- ggplot(df_age_2, aes(x = time, y = infected_per100k, color = age_group)) +
  
  geom_line(linewidth = 1.2) +
  
  coord_cartesian(xlim = c(0, 300)) + 
  
  scale_color_manual(
    values = c(
      "<1" = "#E69F00" , 
      "1–4" = "#56B4E9", 
      "5–17" = "#009E73", 
      "18–29" = "#C4A000",
      "30–44" = "#0072B2", 
      "45–64" = "#D55E00", 
      "65+" = "#CC79A7"
    ),
    labels = c(
      "<1" = "<1 Year",
      "1–4" = "1–4 Years",
      "5–17" = "5–17 Years",
      "18–29" = "18–29 Years",
      "30–44" = "30–44 Years",
      "45–64" = "45–64 Years",
      "65+" = "65+ Years"
    )
  ) +
  
  labs(
    title = NULL,
    x = "Days",
    y = "Active Infections per 100,000",
    color = "Age Group"
  ) +
  
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "right",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"
  )
  )

print(p2)

# Export to high-resolution figure
ggsave("SFig10_timeseries_age_naive.png", width = 7, height = 5, dpi = 300)

# ====================================================
### Next generation matrix and dominant eigenvalue ###
# ====================================================

# Extract parameters
beta    <- params$beta
gamma_a <- params$gamma_a
gamma_l <- params$gamma_l

C_naive <- params$C_naive

N_age <- age_group_nums

# Build next generation matrix
K_naive <- beta * (
  t(C_naive) / gamma_a +
    t(C_naive)  / gamma_l
)

# Check dimensions 
dim(K_naive)
# should be 7 x 7

# Compute the dominant eigenvalue (=R0) 
eig_naive <- eigen(K_naive)

R0_naive <- max(Re(eig_naive$values))
R0_naive

# =================================================
### RO Ratio (behavior dynamic / behavior naive) ###
# =================================================
R0_dynamic / R0_naive


#--- Model Comparison: Time series plot for infected individuals per 100,000 ---#

## Figure 5. Time series plot for the total number of infected 
## individuals each day (per 100,000 population) for the behavior-dynamic 
## and behavior-naïve model simulations 

N <- 331449281

# Combine datasets
df1 <- data.frame(
  time = out_ODE$time,
  infected = total_I,
  model = "Behavior-Dynamic"
)

df2 <- data.frame(
  time = out_ODE_2$time,
  infected = total_I_2,
  model = "Behavior-Naive"
)

df <- bind_rows(df1, df2)

# Convert infections to rate per 100,000
df$infected_per100k <- (df$infected / N) * 100000

# Find peak infection point for each model
peak_df <- df %>%
  group_by(model) %>%
  filter(infected_per100k == max(infected_per100k)) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(peak_label = paste0("Peak: ", round(infected_per100k,1), "\nDay ", time))

# Plot
ggplot(df, aes(x = time, y = infected_per100k, color = model, linetype = model)) +
  
  geom_line(linewidth = 1.4) +
  
  scale_color_manual(values = c(
    "Behavior-Dynamic" = "#5B84B1",
    "Behavior-Naive" = "#8CB369"
  )) +
  
  scale_linetype_manual(values = c(
    "Behavior-Dynamic" = "solid",
    "Behavior-Naive" = "dashed"
  )) +
  
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.02))) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.2))) + # extra space for labels
  
  labs(
    title = "Infected Individuals Over Time",
    x = "Days",
    y = "Infections per 100,000 Population",
    color = "Model",
    linetype = "Model"
  ) +
  
  theme_minimal(base_size = 14) +
  
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "bottom",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black")
  ) +
  
  coord_cartesian(xlim = c(0, 300), clip = "off")

# Export to high-resolution figure
ggsave("Fig5_TimeSeries_Model_Comparison.png", width = 7, height = 5, dpi = 300)





