# From https://github.com/richelbilderbeek/pirouette_article/issues/58 :
#
# A pirouette example
# that shows the true and twin errors for DD trees of different likelihoods,
# using a BD tree prior
#
# It does so by simulating 'n_trees' DD trees from the same parameters.
# From each of those trees, the likelihood is known.
# Afterwards, the simulated DD trees are sorted by that likelihood
#
# From those sorted trees, it takes the top 'n_replicates' trees to
# get a pirouette error distribution.
# Likewise for the bottom 'n_replicates' trees.
# Likewise for the middle 'n_replicates' trees.

library(pirouette)
library(beautier)
library(beastier)
library(testthat)
library(ggplot2)

# Constants
is_testing <- is_on_travis()
example_no <- 23
# The total number of DD trees to simulate
n_trees <- 1000
n_replicates <- 5

if (is_testing) {
  n_replicates <- 2
  n_trees <- n_replicates * 5
}

################################################################################
# Create list of phylogenies and likelihoods
# 'data' is a list, of which each element has
#  * 'phylogeny': a reconstructed phylogeny
#  * 'log_likelihood': the log likelihood of that tree
################################################################################
data <- list()
# Creates phylogenies of a known log-likelihood
for (i in seq(1, n_trees)) {
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

################################################################################
# Sort 'data', create 'sorted_data'
################################################################################
sorted_data <- data[order(sapply(data,'[[',1))]
# Check if really sorted on log-likelihood
lowest <- sorted_data[[1]]$log_likelihood
for (i in seq_along(sorted_data)) {
  testit::assert(lowest <= sorted_data[[i]]$log_likelihood)
  lowest <- sorted_data[[i]]$log_likelihood
}

################################################################################
# Get the phylogenies
################################################################################
phylogenies <- list()
indices <- c(
  seq(1, n_replicates),
  seq(floor((n_trees / 2) - (n_replicates / 2)) + 1, length.out = n_replicates),
  seq(n_trees - n_replicates + 1, length.out = n_replicates)
)
phylogenies <- sapply(data[indices], "[[", 1)
expect_equal(length(phylogenies), 3 * n_replicates)

################################################################################
# Create pirouette parameter sets
################################################################################
pir_paramses <- create_std_pir_paramses(n = 3 * n_replicates)
expect_equal(length(pir_paramses), length(phylogenies))
if (is_testing) {
  pir_paramses <- shorten_pir_paramses(pir_paramses)
}

# Do the runs
pir_outs <- pir_runs(
  phylogenies = phylogenies,
  pir_paramses = pir_paramses
)

# Save
for (i in seq_along(pir_outs)) {
  pir_save(
    phylogeny = phylogenies[[i]],
    pir_params = pir_paramses[[i]],
    pir_out = pir_outs[[i]],
    folder_name = dirname(pir_paramses[[i]]$alignment_params$fasta_filename)
  )
}

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
