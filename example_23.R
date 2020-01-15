# Code of example 1
#
# Works under Windows
#
#
#

suppressMessages(library(pirouette))
suppressMessages(library(ggplot2))


#if (!is_on_travis()) {
#  setwd("~/GitHubs/pirouette_example_1")
#}

root_folder <- getwd()
example_no <- 1
rng_seed <- 314
example_folder <- file.path(root_folder, paste0("example_", example_no, "_", rng_seed))
dir.create(example_folder, showWarnings = FALSE, recursive = TRUE)
setwd(example_folder)
set.seed(rng_seed)
testit::assert(is_beast2_installed())
phylogeny <- create_yule_tree(n_taxa = 6, crown_age = 10)

alignment_params <- create_alignment_params(
  sim_tral_fun = get_sim_tral_with_std_nsm_fun(
    mutation_rate = 0.1,
    site_model = beautier::create_jc69_site_model()
  ),
  root_sequence = create_blocked_dna(length = 1000),
  rng_seed = rng_seed,
  fasta_filename = "true_alignment.fas"
)

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

# Testing on local machine
#if (!is_on_travis()) {
#  for (i in seq_along(experiments)) {
#    experiments[[i]]$inference_model$mcmc <- create_mcmc(chain_length = 20000, store_every = 1000)
#  }
#}

pir_params <- create_pir_params(
  alignment_params = alignment_params,
  experiments = experiments
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

