library(ggplot2)
library(truncnorm)

# Set seed for reproducibility
set.seed(40)

# Distribution of Prior Knowledge
dPriorKnowledge <- function(n_individuals, distribution_type = "truncnorm", params = list()) {
  switch(distribution_type,
         "truncnorm" = {
           # Truncated normal distribution between between 0 and 1
           # with specified mean and sd
           mean_val <- if(is.null(params$mean)) 0.5 else params$mean
           sd_val <- if(is.null(params$sd)) 0.15 else params$sd
           truncnorm::rtruncnorm(n_individuals, a = 0, b = 1, mean = mean_val, sd = sd_val)
         }
  )
}



# Distribution of Individual First Estimates
dIndividualFirstEstimate <- function(prior_knowledge, true_value) {
  # The higher prior knowledge, the smaller the distribution of individual estimates
  # For maximum prior knowledge = 1, sd_of_log = 0.05
  # For minimum prior knowledge = 0, sd_of_log = 1
  sd_of_log = -0.95 * prior_knowledge + 1
  mean_of_log = log(true_value) - sd_of_log^2/2
 
  # Generate a sample from this distribution
  estimate <- rlnorm(1, meanlog = mean_of_log, sdlog = sd_of_log)
  
  return(estimate)
}

# Distribution of Group First Estimates
dGroupFirstEstimate <- function(prior_knowledge_vector, true_value) {
  # Generate first estimates for all individuals
  first_estimates <- sapply(prior_knowledge_vector, dIndividualFirstEstimate, true_value = true_value)
  return(first_estimates)
}

# Function to determine weight of advice based on perceived advice quality
getWeightOnAdvice <- function(advice_quality) {
  weight <- switch(advice_quality,
         "low" = 0.32,
         "neutral" = 0.37,
         "high" = 0.48,
         stop("Invalid advice quality. Use 'low', 'neutral', or 'high'."))

  return(weight)
}


# Psi function for integrating social information
psi <- function(first_estimate, social_info, advice_quality) {
  weight_on_advice <- getWeightOnAdvice(advice_quality)
  self_weight <- 1 - weight_on_advice
  expected_second_estimate <- weight_on_advice * social_info + self_weight * first_estimate
  return(expected_second_estimate)
  # The expected value for the second estimate based on weight on advice and self weight
  # Interpreted as a persons tendency to a second estimate value, it can be seen
  # as the central tendency of her second estimate distribution
}

# Distribution of Individual Second Estimates
dIndividualSecondEstimate <- function(first_estimate, social_info, advice_quality) {
  # Returns the second estimate after applying the Psi function and random noise
  expected_second_estimate <- psi(first_estimate, social_info, advice_quality)
  second_estimate <- rnorm(1, mean = expected_second_estimate, sd = 1)
  return(second_estimate)
}

# Distribution of Group Second Estimates
dGroupSecondEstimates <- function(first_estimate_vector, social_info, advice_quality) {
  # Generate second estimates for all individuals
  second_estimates <- mapply(dIndividualSecondEstimate,
                             first_estimate = first_estimate_vector,
                             advice_quality = advice_quality,
                             MoreArgs = list(social_info = social_info))
  return(second_estimates)
}

# Function to determine WOC with social influence
determineWOCSocialInfluence <- function(first_estimate_vector, second_estimate_vector, true_value) {
  # Calculate mean estimates
  mean_first <- mean(first_estimate_vector)
  mean_second <- mean(second_estimate_vector)
  
  # Calculate individual means (average of individual estimates)
  mean_individual_first <- mean(first_estimate_vector)
  mean_individual_second <- mean(second_estimate_vector)
  
  # Calculate individual average distances (average error of individuals)
  mean_individual_distance_first <- mean(abs(first_estimate_vector - true_value))
  mean_individual_distance_second <- mean(abs(second_estimate_vector - true_value))
  
  # Calculate group distances (error of the aggregated estimate)
  group_distance_first <- abs(mean_first - true_value)
  group_distance_second <- abs(mean_second - true_value)
  
  # Calculate WOC benefit in absolute units (how many units closer the group is vs avg individual)
  woc_benefit_first <- mean_individual_distance_first - group_distance_first
  woc_benefit_second <- mean_individual_distance_second - group_distance_second
  
  # Change in WOC benefit
  woc_benefit_change <- woc_benefit_second - woc_benefit_first
  
  # Percentage WOC benefit (what % of individual error is eliminated by aggregation)
  woc_percent_benefit_first <- (woc_benefit_first / mean_individual_distance_first) * 100
  woc_percent_benefit_second <- (woc_benefit_second / mean_individual_distance_second) * 100
  woc_percent_benefit_change <- woc_percent_benefit_second - woc_percent_benefit_first
  
  # Also calculate the squared error metrics for completeness
  error_first <- (mean_first - true_value)^2
  error_second <- (mean_second - true_value)^2
  individual_error_first <- mean((first_estimate_vector - true_value)^2)
  individual_error_second <- mean((second_estimate_vector - true_value)^2)
  
  return(list(
    # Group estimates
    group_mean_first = mean_first,
    group_mean_second = mean_second,
    group_distance_first = group_distance_first,
    group_distance_second = group_distance_second,
    
    # Individual averages
    individual_mean_distance_first = mean_individual_distance_first,
    individual_mean_distance_second = mean_individual_distance_second,
    
    # WOC benefits in original units
    woc_benefit_first = woc_benefit_first,
    woc_benefit_second = woc_benefit_second,
    woc_benefit_change = woc_benefit_change,
    
    # WOC benefits as percentages
    woc_percent_benefit_first = woc_percent_benefit_first,
    woc_percent_benefit_second = woc_percent_benefit_second,
    woc_percent_benefit_change = woc_percent_benefit_change,
    
    # Original squared error metrics
    error_first = error_first,
    error_second = error_second,
    individual_error_first = individual_error_first,
    individual_error_second = individual_error_second,
    improvement = error_first - error_second,
    individual_improvement = individual_error_first - individual_error_second
  ))
}

# Function to run a simulation
runSimulation <- function(n_individuals, true_value, advice_quality, knowledge_distribution = "truncnorm", 
                          knowledge_params = list(), social_info_type = "mean", 
                          manipulated_social_info = NULL, n_trials = 1) {
  
  results <- list()
  
  for(i in 1:n_trials) {
    # Generate prior knowledge
    prior_knowledge <- dPriorKnowledge(n_individuals, knowledge_distribution, knowledge_params)
    
    # Generate first estimates
    first_estimates <- dGroupFirstEstimate(prior_knowledge, true_value)
    
    # Calculate personalized social information and second estimates
    if(!is.null(manipulated_social_info)) {
      # If social_info is manipulated, use the same value for everyone
      social_info <- manipulated_social_info
      second_estimates <- dGroupSecondEstimates(first_estimates, social_info, advice_quality)
    } else {
      # Create vectors to store personalized social info and second estimates
      social_info_vector <- numeric(n_individuals)
      second_estimates <- numeric(n_individuals)
      
      # Calculate personalized social info for each individual (excluding their own estimate)
      for(j in 1:n_individuals) {
        # Get all estimates except this person's
        other_estimates <- first_estimates[-j]
        
        
        # Calculate social information based on specified type
        switch(social_info_type,
               "mean" = {
                 social_info_vector[j] <- mean(other_estimates)
               }
        )
        
        # Calculate second estimate for this person using their personalized social info
        second_estimates[j] <- dIndividualSecondEstimate(
          first_estimate = first_estimates[j],
          social_info = social_info_vector[j],
          advice_quality = advice_quality
        )
      }
      
      # For compatibility with the rest of the code, store the average social info
      social_info <- mean(social_info_vector)
    }
    
    # Determine WOC effects
    woc_effects <- determineWOCSocialInfluence(first_estimates, second_estimates, true_value)
    
    # Store results for this trial
    results[[i]] <- list(
      prior_knowledge = prior_knowledge,
      advice_quality = advice_quality,
      first_estimates = first_estimates,
      social_info = social_info,
      second_estimates = second_estimates,
      woc_effects = woc_effects
    )
  }
  
  if(n_trials == 1) {
    return(results[[1]])
  } else {
    aggregate_results <- list(
      # Advice Quality
      advice_quality = advice_quality,
      
      # Group estimates
      mean_group_first = mean(sapply(results, function(x) x$woc_effects$group_mean_first)),
      mean_group_second = mean(sapply(results, function(x) x$woc_effects$group_mean_second)),
      mean_group_distance_first = mean(sapply(results, function(x) x$woc_effects$group_distance_first)),
      mean_group_distance_second = mean(sapply(results, function(x) x$woc_effects$group_distance_second)),
      
      # Individual averages
      mean_individual_distance_first = mean(sapply(results, function(x) x$woc_effects$individual_mean_distance_first)),
      mean_individual_distance_second = mean(sapply(results, function(x) x$woc_effects$individual_mean_distance_second)),
      
      # WOC benefits in original units
      mean_woc_benefit_first = mean(sapply(results, function(x) x$woc_effects$woc_benefit_first)),
      mean_woc_benefit_second = mean(sapply(results, function(x) x$woc_effects$woc_benefit_second)),
      mean_woc_benefit_change = mean(sapply(results, function(x) x$woc_effects$woc_benefit_change)),
      
      # WOC benefits as percentages
      mean_woc_percent_benefit_first = mean(sapply(results, function(x) x$woc_effects$woc_percent_benefit_first)),
      mean_woc_percent_benefit_second = mean(sapply(results, function(x) x$woc_effects$woc_percent_benefit_second)),
      mean_woc_percent_benefit_change = mean(sapply(results, function(x) x$woc_effects$woc_percent_benefit_change)),
      
      # Proportion of trials showing improvement
      proportion_group_improved = mean(sapply(results, function(x) x$woc_effects$group_distance_first > x$woc_effects$group_distance_second)),
      proportion_individual_improved = mean(sapply(results, function(x) x$woc_effects$individual_mean_distance_first > x$woc_effects$individual_mean_distance_second)),
      proportion_woc_strengthened = mean(sapply(results, function(x) x$woc_effects$woc_benefit_change > 0))
    )
    
    return(list(
      trial_results = results,
      aggregate_results = aggregate_results
    ))
  }
}

# Simulate group of 10 individuals with high prior knowledge (0.9)
## Advice quality = "low"

lowAQ_group_results_high <- runSimulation(
  n_individuals = 10,              
  true_value = 100, 
  advice_quality = "low",
  social_info_type = "mean",
  knowledge_distribution = "truncnorm", 
  knowledge_params = list(
    mean = 0.9,                   
    sd = 0.15                      
  ),
  n_trials = 100
)
## Advice quality = "neutral"

neutralAQ_group_results_high <- runSimulation(
  n_individuals = 10,              
  true_value = 100, 
  advice_quality = "neutral",
  social_info_type = "mean",
  knowledge_distribution = "truncnorm", 
  knowledge_params = list(
    mean = 0.9,                   
    sd = 0.15                       
  ),
  n_trials = 100
)

## Advice quality = "high"

highAQ_group_results_high <- runSimulation(
  n_individuals = 10,              
  true_value = 100, 
  advice_quality = "high",
  social_info_type = "mean",
  knowledge_distribution = "truncnorm", 
  knowledge_params = list(
    mean = 0.9,                   
    sd = 0.15                     
  ),
  n_trials = 100
)

print(lowAQ_group_results_high$aggregate_results)
print(neutralAQ_group_results_high$aggregate_results)
print(highAQ_group_results_high$aggregate_results)



# Simulate group of 10 individuals with low prior knowledge (0.1)
## Advice quality = "low"

lowAQ_group_results_low <- runSimulation(
  n_individuals = 10,              
  true_value = 100, 
  advice_quality = "low",
  social_info_type = "mean",
  knowledge_distribution = "truncnorm", 
  knowledge_params = list(
    mean = 0.1,                   
    sd = 0.15                      
  ),
  n_trials = 100
)

## Advice quality = "neutral"

neutralAQ_group_results_low <- runSimulation(
  n_individuals = 10,              
  true_value = 100, 
  advice_quality = "neutral",
  social_info_type = "mean",
  knowledge_distribution = "truncnorm", 
  knowledge_params = list(
    mean = 0.1,                   
    sd = 0.15                       
  ),
  n_trials = 100
)

## Advice quality = "high"

highAQ_group_results_low <- runSimulation(
  n_individuals = 10,              
  true_value = 100, 
  advice_quality = "high",
  social_info_type = "mean",
  knowledge_distribution = "truncnorm", 
  knowledge_params = list(
    mean = 0.1,                   
    sd = 0.15                     
  ),
  n_trials = 100
)

print(lowAQ_group_results_low$aggregate_results)
print(neutralAQ_group_results_low$aggregate_results)
print(highAQ_group_results_low$aggregate_results)



# Print WOC-benefit_changes
print(lowAQ_group_results_high$aggregate_results$mean_woc_percent_benefit_change)
print(neutralAQ_group_results_high$aggregate_results$mean_woc_percent_benefit_change)
print(highAQ_group_results_high$aggregate_results$mean_woc_percent_benefit_change)

print(highAQ_group_results_low$aggregate_results$mean_woc_percent_benefit_change)
print(neutralAQ_group_results_low$aggregate_results$mean_woc_percent_benefit_change)
print(highAQ_group_results_low$aggregate_results$mean_woc_percent_benefit_change)

#Als Tabelle
woc_percent_benefit_change_table <- data.frame(
  Advice_Quality = c("Low", "Neutral", "High"),  # Die Zeilen sind die Advice Quality Kategorien
  `High Prior Knowledge` = c(
    lowAQ_group_results_high$aggregate_results$mean_woc_percent_benefit_change,
    neutralAQ_group_results_high$aggregate_results$mean_woc_percent_benefit_change,
    highAQ_group_results_high$aggregate_results$mean_woc_percent_benefit_change
  ),
  `Low Prior Knowledge` = c(
    lowAQ_group_results_low$aggregate_results$mean_woc_percent_benefit_change,
    neutralAQ_group_results_low$aggregate_results$mean_woc_percent_benefit_change,
    highAQ_group_results_low$aggregate_results$mean_woc_percent_benefit_change
  )
)

# Ausgabe der Tabelle mit prozentualer Veränderung der WOC
print(woc_percent_benefit_change_table)





# Effektstärke von prior knowledge (low vs. high) auf den Mittelwert der ersten Schätzung
## Funktion zur Berechnung von Cohen's d
calculate_cohens_d <- function(group1, group2) {
  # Entferne NA-Werte
  group1 <- na.omit(group1)
  group2 <- na.omit(group2)
  
  mean_group1 <- mean(group1)
  mean_group2 <- mean(group2)
  sd_group1 <- sd(group1)
  sd_group2 <- sd(group2)
  n1 <- length(group1)
  n2 <- length(group2)
  
  # Gemeinsame Standardabweichung
  sp <- sqrt(((n1 - 1) * sd_group1^2 + (n2 - 1) * sd_group2^2) / (n1 + n2 - 2))
  
  # Cohen's d
  d <- (mean_group1 - mean_group2) / sp
  return(d)
}

## Extrahiere first_estimates für high und low Prior Knowledge
first_estimates_high <- unlist(lapply(lowAQ_group_results_high$trial_results, function(x) x$first_estimates))
first_estimates_low <- unlist(lapply(lowAQ_group_results_low$trial_results, function(x) x$first_estimates))

## Berechne Mittelwerte ohne NA-Werte
mean_first_estimates_high <- mean(first_estimates_high, na.rm = TRUE)
mean_first_estimates_low <- mean(first_estimates_low, na.rm = TRUE)

## Berechne Cohen's d für die Mittelwerte der first_estimates (low vs high Prior Knowledge)
cohens_d_mean <- calculate_cohens_d(first_estimates_high, first_estimates_low)

## Tabelle mit der berechneten Effektstärke (Cohen's d) für die Mittelwerte der first estimates
effect_size_table <- data.frame(
  Measure = "Mean First Estimate",
  Cohen_d = cohens_d_mean
)

## Ausgabe der Tabelle
print(effect_size_table)





# Berechnung der Effektstärken für Advice Quality (low vs high) bei mittleren Prior Knowledge (0.5) 

## Funktion zur Berechnung von Cohen's d
calculate_cohens_d <- function(group1, group2) {
  # Entferne NA-Werte
  group1 <- na.omit(group1)
  group2 <- na.omit(group2)
  
  mean_group1 <- mean(group1)
  mean_group2 <- mean(group2)
  sd_group1 <- sd(group1)
  sd_group2 <- sd(group2)
  n1 <- length(group1)
  n2 <- length(group2)
  
  # Gemeinsame Standardabweichung
  sp <- sqrt(((n1 - 1) * sd_group1^2 + (n2 - 1) * sd_group2^2) / (n1 + n2 - 2))
  
  # Cohen's d
  d <- (mean_group1 - mean_group2) / sp
  return(d)
}

## Simuliere eine Gruppe mit mittlerem Prior Knowledge (0.5)
## Advice quality = "low"
lowAQ_group_results_medium <- runSimulation(
  n_individuals = 10,              
  true_value = 100, 
  advice_quality = "low",
  social_info_type = "mean",
  knowledge_distribution = "truncnorm", 
  knowledge_params = list(
    mean = 0.5,                   
    sd = 0.15                      
  ),
  n_trials = 100
)

## Advice quality = "neutral"
neutralAQ_group_results_medium <- runSimulation(
  n_individuals = 10,              
  true_value = 100, 
  advice_quality = "neutral",
  social_info_type = "mean",
  knowledge_distribution = "truncnorm", 
  knowledge_params = list(
    mean = 0.5,                   
    sd = 0.15                      
  ),
  n_trials = 100
)

## Advice quality = "high"
highAQ_group_results_medium <- runSimulation(
  n_individuals = 10,              
  true_value = 100, 
  advice_quality = "high",
  social_info_type = "mean",
  knowledge_distribution = "truncnorm", 
  knowledge_params = list(
    mean = 0.5,                   
    sd = 0.15                      
  ),
  n_trials = 100
)

## Extrahiere second_estimates für low und high Advice Quality bei mittlerem Prior Knowledge
second_estimates_low_medium <- unlist(lapply(lowAQ_group_results_medium$trial_results, function(x) x$second_estimates))
second_estimates_high_medium <- unlist(lapply(highAQ_group_results_medium$trial_results, function(x) x$second_estimates))

## Berechnung der Effektstärke (Cohen's d) für low vs high Advice Quality bei mittlerem Prior Knowledge
cohens_d_advice_quality_medium <- calculate_cohens_d(second_estimates_low_medium, second_estimates_high_medium)

## Tabelle mit den berechneten Effektstärken
effect_size_table_advice_quality_medium <- data.frame(
  Advice_Quality_Comparison = "Low vs High (Medium Prior Knowledge)",
  Cohen_d = cohens_d_advice_quality_medium
)

## Ausgabe der Tabelle mit der berechneten Effektstärke
print(effect_size_table_advice_quality_medium)






## Plot der Verteilung der ersten Schätzung in der Gruppe abhängig von Prior Knowledge in der Gruppe

## Extrahiere die first_estimates für jede Gruppe (high und low Prior Knowledge)
first_estimates_high_pknw <- unlist(lapply(lowAQ_group_results_high$trial_results, function(x) x$first_estimates))
first_estimates_low_pknw <- unlist(lapply(lowAQ_group_results_low$trial_results, function(x) x$first_estimates))

## Erstelle DataFrame für den Plot
data_for_first_estimates_plot <- data.frame(
  Estimate = c(first_estimates_high_pknw, first_estimates_low_pknw),
  Prior_Knowledge = rep(c("High Prior Knowledge", "Low Prior Knowledge"), 
                        each = length(first_estimates_high_pknw))
)

## Plot erstellen
ggplot(data_for_first_estimates_plot, aes(x = Estimate, fill = Prior_Knowledge)) +
  geom_density(alpha = 0.6) +
  labs(title = "Verteilung der ersten Schätzungen in der Gruppe abhängig von Prior Knowledge",
       x = "Erste Schätzungen",
       y = "Dichte") +
  theme_minimal() +
  theme(legend.position = "top") +
  xlim(0, 200)






# Plot der Verteilung der zweiten Schätzung in der Gruppe abhängig von Prior Knowledge in der Gruppe und der gegebenen advice quality

## Extrahiere die second_estimates für jede Gruppe (high und low Prior Knowledge und verschiedene Advice Quality)
second_estimates_high_pknw_low_advice <- unlist(lapply(lowAQ_group_results_high$trial_results, function(x) x$second_estimates))
second_estimates_high_pknw_neutral_advice <- unlist(lapply(neutralAQ_group_results_high$trial_results, function(x) x$second_estimates))
second_estimates_high_pknw_high_advice <- unlist(lapply(highAQ_group_results_high$trial_results, function(x) x$second_estimates))

second_estimates_low_pknw_low_advice <- unlist(lapply(lowAQ_group_results_low$trial_results, function(x) x$second_estimates))
second_estimates_low_pknw_neutral_advice <- unlist(lapply(neutralAQ_group_results_low$trial_results, function(x) x$second_estimates))
second_estimates_low_pknw_high_advice <- unlist(lapply(highAQ_group_results_low$trial_results, function(x) x$second_estimates))

## Erstelle DataFrame für den Plot
data_for_second_estimates_plot <- data.frame(
  Estimate = c(
    second_estimates_high_pknw_low_advice, second_estimates_high_pknw_neutral_advice, second_estimates_high_pknw_high_advice,
    second_estimates_low_pknw_low_advice, second_estimates_low_pknw_neutral_advice, second_estimates_low_pknw_high_advice
  ),
  Prior_Knowledge = rep(c("High Prior Knowledge", "Low Prior Knowledge"), 
                        each = length(second_estimates_high_pknw_low_advice) + length(second_estimates_high_pknw_neutral_advice) + length(second_estimates_high_pknw_high_advice)),
  Advice_Quality = rep(c("Low", "Neutral", "High"), 
                       each = length(second_estimates_high_pknw_low_advice), times = 2)  # 2 Gruppen (high, low Prior Knowledge)
)

## Plot erstellen
ggplot(data_for_second_estimates_plot, aes(x = Estimate, color = Advice_Quality, linetype = Advice_Quality)) +
  geom_density(linewidth = 1.2) +  # Linien statt Füllung
  facet_wrap(~ Prior_Knowledge) +  # Facet für verschiedene Prior Knowledge Levels
  labs(title = "Verteilung der zweiten Schätzungen in der Gruppe abhängig von Prior Knowledge und advice quality",
       x = "Zweite Schätzungen",
       y = "Dichte") +
  theme_minimal() +
  theme(legend.position = "top") +
  xlim(0, 200)
  scale_color_manual(values = c("Low" = "red", "Neutral" = "blue", "High" = "green")) +  # Farben nach Advice Quality
  scale_linetype_manual(values = c("Low" = "solid", "Neutral" = "dashed", "High" = "dotted"))  # Verschiedene Linientypen für Advice Quality


