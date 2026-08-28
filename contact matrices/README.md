# **Contact Matrices**

This folder contains contact matrices for mathematical modeling and visualization of age-structured contact patterns. 

Contact matrices are based on data from the SickMix study. 

Matrices corresponding to acute gastroenteritis (AGE) and acute respiratory infection (ARI) cases represent self-reported contacts across three survey time points, 
weighted using MarketScan benchmarks according to age, sex, and diagnosis (AGE/ARI). 

Matrices corresponding to household members are based on self-reported contacts of household members of acute gastroenteritis and acute respiratory infection cases across two survey time points. 

Each cell of the contact matrices $M_{ij}$ represents the weighted mean number of contacts with individuals in age group j (column) reported by a participant in age group i (row). 

$$
M_{ij} = \frac{\sum_{r=1}^{N_i} C_{ij}^r}{N_j}
$$

- $i$ is the age group of the respondent
- $j$ is the age group of the contact
- $r$ is the respondent in age group i
- $C_{ij}^r$ is the number of contacts reported by respondent r in age group i with a contact in age group j
- $N_{j}$ is the number of contacts in age group j

The folder case_matrices contains these matrices for cases, and the folder HHM_matrices contains these matrices for household members. 
