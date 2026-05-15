library(tidyverse)
library(ggplot2)
library(R2jags)

# Define parameters
J <- 5  # Number of dose levels
N <- 9 # max patients
target_tox_rate <- 0.3
beta0_prior_mean <- 0
beta0_prior_sd <- 1
beta1_prior_rate <- 1
cohort_size <- 3
max_titrations <- 3
c <- 0.5  # Early termination threshold for high toxicity at the lowest dose


# Define initial skeleton of toxicity probabilities 
skeleton <-  c(0.02, 0.12, 0.3, 0.5, 0.65)

#True toxicity probs
scenario1 <- c(0.1, 0.3, 0.55, 0.65, 0.75) #true mtd at dose level 2
scenario2 <- c(0.05, 0.15, 0.3, 0.5, 0.8) #true mtd at dose level 3
scenario3 <- c(0.02, 0.05, 0.1, 0.3, 0.5) #true mtd at dose level 4
scenario4 <- c(0.05, 0.3, 0.2, 0.4, 0.6) #true mtd at dose level 2f but not monotonically increasing


simulate_data <- function(n_patients, dose_level, toxicity_probs) {
  tox_prob <- toxicity_probs[dose_level]  # Probability of toxicity at this dose
  
  # Generate toxicities for each patient in the cohort
  toxicity <- rbinom(n_patients, 1, tox_prob)  # Toxicity for each patient (1 = toxic, 0 = no toxic)
  
  # Create a dataframe with patient ID, dose level, and the simulated toxicity
  data <- data.frame(
    patient_id = 1:n_patients,
    dose_level = rep(dose_level, n_patients),
    toxicity = toxicity
  )
  
  return(data)  
}




# Define the function to sample from the posterior distribution
gibbs <- function(data, beta0_prior_mean, beta0_prior_sd, beta1_prior_rate) {
  
  # Define the model in JAGS syntax
  model_code <- "
    model {
      # Prior distributions for beta0 and beta1
      beta0 ~ dnorm(beta0_prior_mean, 1 / (beta0_prior_sd^2))
      beta1 ~ dexp(beta1_prior_rate)
      
      # Likelihood
      for (i in 1:N) {
        toxicity[i] ~ dbern(p[i])
        logit(p[i]) <- beta0 + beta1 * dose_level[i]
      }
    }
  "
  
  # Define the data
  jags_data <- list(
    N = length(data$toxicity),
    toxicity = data$toxicity,
    dose_level = data$dose_level,
    beta0_prior_mean = beta0_prior_mean,
    beta0_prior_sd = beta0_prior_sd,
    beta1_prior_rate = beta1_prior_rate
  )
  
  # Define the parameters to monitor
  parameters <- c("beta0", "beta1")
  
  # Fit the model
  model_fit <- jags(data = jags_data, 
                    parameters.to.save = parameters, 
                    model.file = textConnection(model_code), 
                    n.chains = 4, 
                    n.iter = 1000,
                    n.burnin = 200,
                    n.thin = 2)
  
  beta0_samples <- model_fit$BUGSoutput$sims.list$beta0
  beta1_samples <- model_fit$BUGSoutput$sims.list$beta1
  
  # Return the posterior samples
  return(list(beta0_samples = beta0_samples, beta1_samples = beta1_samples))
}








# MTD Calculation


estimate_mtd <- function(beta0_samples, beta1_samples, target_tox_rate, J) {
  # Matrix to store the estimated toxicity probabilities for each dose level
  estimated_tox_probs <- matrix(0, nrow = length(beta0_samples), ncol = J)
  
  # Calculate estimated toxicity probabilities for each dose level using the posterior samples
  for (i in 1:J) {
    estimated_tox_probs[, i] <- 1 / (1 + exp(-(beta0_samples + beta1_samples * (i))))
  }
  
  # Calculate the mean toxicity probability for each dose level
  mean_tox_probs <- apply(estimated_tox_probs, 2, mean)
  print("Mean Toxicity Probabilities by Dose Level:")
  print(mean_tox_probs) 
  
  # Find the dose level closest to the target toxicity rate
  closest_dose <- which.min(abs(mean_tox_probs - target_tox_rate))
  
  return(closest_dose)
}






#IP-CRM simulation
simulate_trial_ip_crm <- function(J, N, cohort_size, max_titrations, scenario, target_tox_rate, c, skeleton) {
  # Initialize data frame to store results from all patients
  all_data <- data.frame(patient_id = integer(), dose_level = integer(), toxicity = integer(), escalation_count = integer())
  dose_level <- 1  # Start with the lowest dose level
  cohort <- 1
  trial_terminated <- FALSE
  final_mtd <- 0
  toxic_patients <- 0  # Initialize the number of toxic patients
  
  while (cohort <= ceiling(N / cohort_size) && !trial_terminated) {
    print(paste("Cohort", cohort, "starting at dose level", dose_level))
    
    # Simulate the data for the current cohort
    data <- simulate_data(cohort_size, dose_level, scenario)
    
    # Count toxic patients for this cohort
    toxic_patients <- toxic_patients + sum(data$toxicity)
    
    # Initialize the escalation count for each patient (set to 0 initially)
    data$escalation_count <- 0
    
    # Add the data for the cohort to the overall trial data
    all_data <- rbind(all_data, data)
    
    # Check if the first dose level is too toxic (i.e., toxicity > c) for any patient
    if (dose_level == 1 && mean(data$toxicity) > c) {
      print("Trial terminated early due to high toxicity at the first dose level.")
      trial_terminated <- TRUE
      break
    }
    
    # Implement intra-patient dose escalation: check toxicity and escalate doses for patients
    for (i in 1:nrow(data)) {
      while (data$toxicity[i] == 0 && data$escalation_count[i] < max_titrations && data$dose_level[i] < J) {
        new_dose <- min(data$dose_level[i] + 1, J)
        print(paste("Patient", data$patient_id[i], "escalated from dose", data$dose_level[i], "to dose", new_dose))
        
        # Simulate toxicity at the new dose level
        data$toxicity[i] <- rbinom(1, 1, scenario[new_dose])  
        
        # Update the patient's current dose and increment escalation count
        data$dose_level[i] <- new_dose  
        data$escalation_count[i] <- data$escalation_count[i] + 1
        
        # Update all_data with the new treatment cycle results
        all_data <- rbind(all_data, data[i, ,drop = FALSE])
        
        # Count toxic patients for this escalation
        toxic_patients <- toxic_patients + sum(data$toxicity[i] == 1)
      }
    }
    
    # Re-estimate the MTD using logistic regression (based on all data so far)
    gibbs_results <- gibbs(all_data, beta0_prior_mean, beta0_prior_sd, beta1_prior_rate)
    final_mtd <- estimate_mtd(gibbs_results$beta0_samples, gibbs_results$beta1_samples, target_tox_rate, J)
    print(paste("Updated MTD estimate: Dose level", final_mtd))
    
    # Update the starting dose for the next cohort based on the new MTD
    if (final_mtd > dose_level) {
      dose_level <- min(final_mtd, J)
    } else if (final_mtd < dose_level) {
      dose_level <- max(final_mtd, 1)
    } else {
      dose_level <- final_mtd
    }
    
    cohort <- cohort + 1
  }
  
  # Calculate the safety index: proportion of patients with toxic doses
  safety_index <- toxic_patients / nrow(all_data)
  return(list(all_data = all_data, final_mtd = final_mtd, safety_index = safety_index))
}





# Traditional CRM Simulation
simulate_trial_traditional_crm <- function(J, N, cohort_size, scenario, target_tox_rate, skeleton) {
  # Initialize data frame to store results from all patients
  all_data <- data.frame(patient_id = integer(), dose_level = integer(), toxicity = integer())
  dose_level <- 1  # Start with the lowest dose level
  cohort <- 1
  trial_terminated <- FALSE
  toxic_patients <- 0  # Initialize the number of toxic patients
  
  while (cohort <= ceiling(N / cohort_size) && !trial_terminated) {
    print(paste("Cohort", cohort, "starting at dose level", dose_level))
    
    # Simulate the data for the current cohort
    data <- simulate_data(cohort_size, dose_level, scenario)
    
    # Count toxic patients for this cohort
    toxic_patients <- toxic_patients + sum(data$toxicity)
    
    # Add the data for the cohort to the overall trial data
    all_data <- rbind(all_data, data)
    
    # Re-estimate the MTD using logistic regression (based on all data so far)
    gibbs_results <- gibbs(all_data, beta0_prior_mean, beta0_prior_sd, beta1_prior_rate)
    final_mtd <- estimate_mtd(gibbs_results$beta0_samples, gibbs_results$beta1_samples, target_tox_rate, J)
    print(paste("Updated MTD estimate: Dose level", final_mtd))
    
    # Update the starting dose for the next cohort based on the new MTD
    if (final_mtd > dose_level) {
      dose_level <- min(final_mtd, J)
    } else if (final_mtd < dose_level) {
      dose_level <- max(final_mtd, 1)
    }
    
    cohort <- cohort + 1
  }
  
  # Calculate the safety index: proportion of patients with toxic doses
  safety_index <- toxic_patients / nrow(all_data)
  return(list(all_data = all_data, final_mtd = final_mtd, safety_index = safety_index))
}





# 3+3 method simulation
simulate_trial_3_plus_3 <- function(J, N, cohort_size, scenario, target_tox_rate) {
  # Initialize data frame to store results from all patients
  all_data <- data.frame(patient_id = integer(), dose_level = integer(), toxicity = integer())
  dose_level <- 1  # Start with the lowest dose level
  cohort <- 1
  trial_terminated <- FALSE
  toxic_patients <- 0  # Initialize the number of toxic patients
  
  while (cohort <= ceiling(N / cohort_size) && !trial_terminated) {
    print(paste("Cohort", cohort, "starting at dose level", dose_level))
    
    # Simulate the data for the current cohort (3 patients at a time)
    data <- simulate_data(cohort_size, dose_level, scenario)
    
    # Count toxic patients for this cohort
    toxic_patients <- toxic_patients + sum(data$toxicity)
    
    # Add the data for the cohort to the overall trial data
    all_data <- rbind(all_data, data)
    
    # Apply 3+3 rule to determine the next dose level
    num_toxicities <- sum(data$toxicity)
    if (num_toxicities == 0) {
      dose_level <- min(dose_level + 1, J)
    } else if (num_toxicities == 1) {
      print("1 patient with toxicity, stay at the same dose level")
    } else {
      dose_level <- max(dose_level - 1, 1)
    }
    
    # If the dose level is 1 and too toxic, terminate the trial early
    if (dose_level == 1 && mean(data$toxicity) > target_tox_rate) {
      print("Trial terminated early due to high toxicity at the first dose level.")
      trial_terminated <- TRUE
    }
    
    cohort <- cohort + 1
  }
  
  # Calculate the safety index: proportion of patients with toxic doses
  safety_index <- toxic_patients / nrow(all_data)
  return(list(all_data = all_data, final_dose_level = dose_level, safety_index = safety_index))
}



# Initialize a list to store all results
all_results <- data.frame(
  method = character(), 
  final_mtd = numeric(), 
  safety_index = numeric(),
  scenario = character(),
  stringsAsFactors = FALSE
)

# Run the trials for each scenario
num_trials <- 2000
for (scenario_index in 1:4) {
  # Select the scenario based on the scenario index
  scenario <- list(scenario1, scenario2, scenario3, scenario4)[[scenario_index]]
  
  # Initialize lists to store results for each method
  ip_crm_results <- vector("list", num_trials)
  ip_crm_safety_indices <- numeric(num_trials)  
  
  traditional_crm_results <- vector("list", num_trials)
  traditional_crm_safety_indices <- numeric(num_trials) 
  
  results_3_plus_3 <- vector("list", num_trials)
  safety_indices_3_plus_3 <- numeric(num_trials)  
  
  # Run IP-CRM simulation
  for (i in 1:num_trials) {
    ip_crm_results[[i]] <- simulate_trial_ip_crm(J, N, cohort_size, max_titrations, scenario, target_tox_rate, c, skeleton)
    ip_crm_safety_indices[i] <- ip_crm_results[[i]]$safety_index
  }
  
  # Run Traditional CRM simulation
  for (i in 1:num_trials) {
    traditional_crm_results[[i]] <- simulate_trial_traditional_crm(J, N, cohort_size, scenario, target_tox_rate, skeleton)
    traditional_crm_safety_indices[i] <- traditional_crm_results[[i]]$safety_index
  }
  
  # Run 3+3 simulation
  for (i in 1:num_trials) {
    results_3_plus_3[[i]] <- simulate_trial_3_plus_3(J, N, cohort_size, scenario, target_tox_rate)
    safety_indices_3_plus_3[i] <- results_3_plus_3[[i]]$safety_index
  }
  
  # Extract final MTDs for each method
  ip_crm_mtds <- sapply(ip_crm_results, function(x) x$final_mtd)
  traditional_crm_mtds <- sapply(traditional_crm_results, function(x) x$final_mtd)
  mtds_3_plus_3 <- sapply(results_3_plus_3, function(x) x$final_dose_level)
  
  # Combine results into a data frame for this scenario
  scenario_df <- data.frame(
    method = rep(c("IP-CRM", "Traditional CRM", "3+3"), each = num_trials),
    safety_index = c(ip_crm_safety_indices, traditional_crm_safety_indices, safety_indices_3_plus_3),
    final_mtd = c(ip_crm_mtds, traditional_crm_mtds, mtds_3_plus_3),
    scenario = rep(paste("Scenario", scenario_index), num_trials * 3)
  )
  
  # Append the current scenario's results to the main results dataframe
  all_results <- rbind(all_results, scenario_df)
}



#Visualization

true_mtd_values <- c(2, 3, 4, 2)

 
ggplot(all_results, aes(x = method, y = final_mtd, fill = method)) +
  geom_boxplot() +
  facet_wrap(~ scenario) +  # Create separate plots for each scenario
  geom_hline(data = data.frame(scenario = paste("Scenario", 1:4), true_mtd = true_mtd_values),
             aes(yintercept = true_mtd), color = "red", linetype = "dashed", size = 1) +
  labs(
    title = "Final MTD Estimates by Method and Scenario",
    x = "Method",
    y = "Final MTD",
    fill = "Method"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Calculate the selection percentages for each dose level and method
selection_percentages <- all_results %>%
  group_by(scenario, method, final_mtd) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(scenario, method) %>%
  mutate(percentage = (count/2000) * 100) %>%
  ungroup() %>%
  # Complete all combinations of scenario, method, and final_mtd
  complete(scenario, method, final_mtd, fill = list(count = 0, percentage = 0)) %>%
  ungroup()

# Create the bar plot with consistent bar widths and show all dose levels
ggplot(selection_percentages, aes(x = factor(final_mtd), y = percentage, fill = method)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.8) +  # Set width for consistent bar size
  facet_wrap(~ scenario, scales = "free_x") +  
  labs(
    title = "Selection Percentage by Dose Level and Method",
    x = "Dose Level",
    y = "Selection Percentage (%)",
    fill = "Method"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(hjust = 1))


# Calculate the mean safety index for each scenario and method
safety_indices <- all_results %>% 
  group_by(method, scenario) %>%
  summarise(safety_index = mean(safety_index))
  
