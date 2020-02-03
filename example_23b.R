
print("Start")
# beast?
if (beastier::is_beast2_installed() == FALSE) {
  beastier::install_beast2()
}
if (mauricer::is_beast2_ns_pkg_installed() == FALSE) {
  mauricer::install_beast2_pkg("NS")
}
stopifnot(beastier::is_beast2_installed())
stopifnot(mauricer::is_beast2_ns_pkg_installed())

print("parsetting")
# I suggest to do simulations with DDD_sim with crown age 5 and crown age 15,
# and extinction rates 0 and 0.4. Speciation rate is 0.8 and K is 40. And then 10-100 trees for each parameter setting.
parses <- data.frame(
  lambda = 0.8, # matches
  mu = 0.1, # matches
  kk = 40, # matches
  crown_age = 10, # matches
  cond = 1,
  ddmodel = 1,
  n_0 = 2
)

# simulate trees
sim_data <- vector("list", length(parses))
for (i in seq_along(parses)){
  pars <- parses
  loglik <- rep(NA, n_replicates)
  for (seed in seq_len(n_replicates)) {
    set.seed(seed)
    sim_data[[i]][[seed]] <- DDD::dd_sim(
      pars = c(pars$lambda, pars$mu, pars$kk), # matches
      age = pars$crown_age, # matches
      ddmodel = pars$ddmodel # mathes
    )
    loglik[seed] <- DDD::dd_loglik(
      pars1 = c(pars$lambda, pars$mu, pars$kk),
      pars2 = c(
        2 * pars$kk, # matches
        pars$ddmodel, # matches
        pars$cond, # matches
        0, # matches
        0, # ?
        pars$n_0 # matches
      ),
      brts = sim_data[[i]][[seed]]$brts, # matches
      missnumspec = 0 # matches
    )
    sim_data[[i]][[seed]]$loglik <- loglik[seed]
    sim_data[[i]][[seed]]$gamma <- phytools::gammatest(
      phytools::ltt(sim_data[[i]][[seed]]$tes, plot = FALSE, gamma = FALSE)
    )$gamma
  }
}

# create pir_params
pir_paramseses <- vector("list", length(parses))
for (i in seq_along(parses)){
  for (seed in seq_len(n_replicates)) {
    phylogeny1 <- sim_data[[i]][[1]]$tes
    pir_paramseses[[i]][[seed]] <- pirouette::create_test_pir_params(
      alignment_params = pirouette::create_alignment_params(
        sim_tral_fun = pirouette::get_sim_tral_with_std_nsm_fun(
          mutation_rate = pirouette::create_standard_mutation_rate(
            phylogeny = phylogeny1
          ),
          site_model = beautier::create_jc69_site_model()
        ),
        root_sequence = pirouette::create_blocked_dna(length = 1e3),
      ),
      twinning_params = pirouette::create_twinning_params(
        rng_seed_twin_tree = seed,
        rng_seed_twin_alignment = seed,
        sim_twin_tree_fun = pirouette::get_sim_bd_twin_tree_fun(),
        sim_twal_fun = pirouette::get_sim_twal_with_std_nsm_fun(
          mutation_rate = pirouette::create_standard_mutation_rate(
            phylogeny = phylogeny1
          ),
          site_model = beautier::create_jc69_site_model()
        )
      ),
      experiments = pirouette::create_all_experiments(),
      error_measure_params = pirouette::create_error_measure_params(
        error_fun = pirouette::get_gamma_error_fun()
      )
    )
  }
}

# pir run!
pir_outs <- vector("list", length(parses))
for (i in seq_along(parses)){
  phylogenies <- lapply(sim_data[[i]], function(x) x$tes)
  pir_paramses <- pir_paramseses[[i]]
  pir_outs[[i]] <- pirouette::pir_runs(
    phylogenies = phylogenies,
    pir_paramses = pir_paramses
  )
}
pir_outs
