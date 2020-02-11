# From https://github.com/richelbilderbeek/pirouette_article/issues/58 :
#
# Show the true and twin errors for DD trees of different likelihoods
# using a BD tree prior,
suppressMessages(library(pirouette))
suppressMessages(library(ggplot2))

root_folder <- getwd()
example_no <- 23

testit::assert(is_beast2_installed())

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

  example_folder <- file.path(root_folder, paste0("example_", example_no, "_", rng_seed))
  dir.create(example_folder, showWarnings = FALSE, recursive = TRUE)
  setwd(example_folder)

  set.seed(rng_seed)
  phylogeny <- data[[ indices[i] ]]$phylogeny

  alignment_params <- create_alignment_params(
    sim_tral_fun = get_sim_tral_with_std_nsm_fun(
      mutation_rate = pirouette::create_standard_mutation_rate(
        phylogeny = phylogeny
      )
    ),
    root_sequence = create_blocked_dna(length = 1000),
    rng_seed = rng_seed,
    fasta_filename = "true_alignment.fas"
  )

  # No need for candidate models here
  experiment <- create_gen_experiment()
  experiment$beast2_options$input_filename <- "true_alignment_gen.xml"
  experiment$beast2_options$output_state_filename <- "true_alignment_gen.xml.state"
  experiment$inference_model$mcmc$tracelog$filename <- "true_alignment_gen.log"
  experiment$inference_model$mcmc$treelog$filename <- "true_alignment_gen.trees"
  experiment$inference_model$mcmc$screenlog$filename <- "true_alignment_gen.csv"
  experiment$errors_filename <- "true_errors_gen.csv"
  experiments <- list(experiment)

  # Set the RNG seed
  for (i in seq_along(experiments)) {
    experiments[[i]]$beast2_options$rng_seed <- rng_seed
  }

  # Shorter on Travis
  if (is_on_travis()) {
    for (i in seq_along(experiments)) {
      experiments[[i]]$inference_model$mcmc$chain_length <- 3000
      experiments[[i]]$inference_model$mcmc$store_every <- 1000
    }
  }

  twinning_params <- create_twinning_params(
    rng_seed_twin_tree = rng_seed,
    sim_twin_tree_fun = get_sim_bd_twin_tree_fun(),
    rng_seed_twin_alignment = rng_seed,
    sim_twal_fun = get_sim_twal_with_std_nsm_fun(
      mutation_rate = pirouette::create_standard_mutation_rate(
        phylogeny
      )
    ),
    twin_tree_filename = "twin_tree.newick",
    twin_alignment_filename = "twin_alignment.fas",
    twin_evidence_filename = "twin_evidence.csv"
  )

  error_measure_params <- pirouette::create_error_measure_params(
    error_fun = pirouette::get_nltt_error_fun()
  )

  pir_params <- create_pir_params(
    alignment_params = alignment_params,
    experiments = experiments,
    twinning_params = twinning_params,
    error_measure_params = error_measure_params
  )

  rm_pir_param_files(pir_params)

  errors <- pir_run(
    phylogeny,
    pir_params = pir_params
  )

  utils::write.csv(
    x = errors,
    file = file.path(example_folder, "errors.csv"),
    row.names = FALSE
  )

  pir_plot(errors) +
    ggsave(file.path(example_folder, "errors.png"), width = 7, height = 7)

  pir_to_pics(
    phylogeny = phylogeny,
    pir_params = pir_params,
    folder = example_folder
  )

  pir_to_tables(
    pir_params = pir_params,
    folder = example_folder
  )
}
