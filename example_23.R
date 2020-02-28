# From https://github.com/richelbilderbeek/pirouette_article/issues/58 :
#
# Show the true and twin errors for DD trees of different likelihoods
# using a BD tree prior,
suppressMessages(library(pirouette))
suppressMessages(library(ggplot2))

################################################################################
# Constants
################################################################################
is_testing <- is_on_travis()
example_no <- 23

# The total number of DD trees to simulate
n_all_trees <- 12 * 10

# 'data' is a list, of which each element has
#  * 'phylogeny': a reconstructed phylogeny
#  * 'log_likelihood': the log likelihood of that tree
data <- list()

################################################################################
# Creates phylogenies of a known log-likelihood
################################################################################
for (i in seq(1, n_all_trees)) {
  print(paste(i, "/", n_all_trees))
  # Create a list of trees

  speciation_rate <- 0.8 # lambda
  extinction_rate <- 0.1 # mu
  carrying_capacity <- 40 # clade-level
  crown_age <- 10
  dd_parameters <- c(speciation_rate, extinction_rate, carrying_capacity)
  ddmodel <- 1 # linear dependence in speciation rate with parameter K
  set.seed(i)
  dd_sim_result <- DDD::dd_sim(pars = dd_parameters, age  = crown_age, ddmodel = ddmodel)
  phylogeny <- dd_sim_result$tes # Only extant species

  max_num_species <- 2 * carrying_capacity
  conditioning <- 1 # crown age and non-extinction of the phylogeny
  likelihood_branching_times_or_phylogeny <- 0 # branching times
  verbose <- 0 # be silent
  stem_age_or_crown_age <- 2 # crown age
  model_settings <- c(
    max_num_species,
    ddmodel,
    conditioning,
    likelihood_branching_times_or_phylogeny,
    verbose,
    stem_age_or_crown_age
  )
  log_likelihood <- DDD::dd_loglik(
    pars1 = dd_parameters,
    pars2 = model_settings,
    brts = dd_sim_result$brts,
    missnumspec = 0 # Number of missing/unsampled species
  )
  data[[i]] <- list()
  data[[i]]$log_likelihood <- log_likelihood
  data[[i]]$phylogeny <- phylogeny
  testit::assert("phylogeny" %in% names(data[[i]]))
  testit::assert("log_likelihood" %in% names(data[[i]]))
  testit::assert(class(data[[i]]$phylogeny) == "phylo")
  testit::assert(class(data[[i]]$log_likelihood) == "numeric")
}

# Show the distribution of log likelihoods
ggplot2::ggplot(
  data.frame(log_likelihood = sapply(data,'[[', 1)),
  aes(x = log_likelihood)
) + geom_density() + ggsave("likelihoods.png")

# Sort 'data' by log-likelihood
sorted_data <- data[order(sapply(data,'[[',1))]

# Check if really sorted on log-likelihood
lowest <- sorted_data[[1]]$log_likelihood
for (i in seq_along(sorted_data)) {
  testit::assert(lowest <= sorted_data[[i]]$log_likelihood)
  lowest <- sorted_data[[i]]$log_likelihood
}


# The indices of the phylogenies to use:
# 1, quarter, middle, third-quarter, last
indices <- c(
  1,
  round(1 * length(data) / 4),
  round(length(data) / 2),
  round(3 * length(data) / 4),
  length(data)
)

# Show the distribution of log likelihoods
ggplot2::ggplot(
  data.frame(log_likelihood = sapply(data,'[[', 1)),
  aes(x = log_likelihood)
) +
  geom_vline(
    data = data.frame(index = as.factor(seq(1, 5))),
    aes(xintercept = sapply(sorted_data[indices],'[[', 1), colour = index),
    size = 2
  ) +
  geom_density() + ggsave("likelihoods.png")

for (i in seq_along(indices)) {

  # First RNG seed must be 314
  rng_seed <- 314 + i - 1
  testit::assert(rng_seed >= 314)
  testit::assert(rng_seed <= 318)
  print(rng_seed)

  folder_name <- file.path(paste0("example_", example_no, "_", rng_seed))

  set.seed(rng_seed)
  phylogeny <- data[[ indices[i] ]]$phylogeny

  pir_params <- create_std_pir_params(folder_name = folder_name)

  if (is_testing) {
    pir_params <- shorten_pir_params(pir_params)
  }

  errors <- pir_run(
    phylogeny,
    pir_params = pir_params
  )

  utils::write.csv(
    x = errors,
    file = file.path(folder_name, "errors.csv"),
    row.names = FALSE
  )

  pir_plot(errors) +
    ggsave(file.path(folder_name, "errors.png"), width = 7, height = 7)

  pir_to_pics(
    phylogeny = phylogeny,
    pir_params = pir_params,
    folder = folder_name
  )

  pir_to_tables(
    pir_params = pir_params,
    folder = folder_name
  )
}
