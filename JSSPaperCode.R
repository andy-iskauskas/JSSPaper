#### Loading Package ####
## Installation from CRAN
# install.packages('hmer')
## For newest version, install from github
# require(devtools)
# devtools::install_github("andy-iskauskas/hmer")

### Load in the required and suggested packages
# deSolve is only required to run ode_results
library(hmer)
library(lhs)
library(deSolve)
library(ggplot2)
set.seed(12)

################# HELPER FUNCTIONS ###################

# Simple code for generating the output of the SIR deterministic model
ode_results <- function(parms, end_time = 25) {
  des = function(time, state, parms) {
    with(as.list(c(state, parms)), {
      dS <- aSR*R-aSI*I*S/(S+I+R)
      dI <- aSI*I*S/(S+I+R)-aIR*I
      dR <- aIR*I-aSR*R
      return(list(c(dS, dI, dR)))
    })
  }
  yini = c(S = 950, I = 50, R = 0)
  times = seq(0, end_time, by = 1)
  out = ode(yini, times, des, parms)
  return(out)
}

# Collates the results to give the outputs for emulation
get_res <- function(inputs) {
  ode_out <- data.frame(t(apply(inputs, 1, function(x) ode_results(x)[11, c('S', 'I', 'R')])))
  return(cbind(inputs, ode_out) |> setNames(c(names(inputs), c('nS', 'nI', 'nR'))))
}

# Very simple wrapper for lhs::randomLHS that scales the inputs to our ranges
get_lhs_design <- function(npoints, ranges) {
  initial_design <- lhs::optimumLHS(npoints, length(ranges))
  return(setNames(data.frame(t(apply(initial_design, 1, function(x) x*purrr::map_dbl(ranges, diff)+purrr::map_dbl(ranges, ~.[[1]])))), names(ranges)))
}

## Usage of the above functions ##
# initial_points <- get_lhs_design(90, ranges)
# initial_results <- get_res(initial_points)

########### EMULATION SET UP ##############

# Define the training and validation runs
training1 <- SIRSample$training
valid1 <- SIRSample$validation
# Set up the ranges
ranges <- list(aSI = c(0.1, 0.8), aIR = c(0, 0.5), aSR = c(0, 0.05))
targets <- SIREmulators$targets
# Just a little plot function
plot(do.call('rbind.data.frame', SIRSample)[,names(ranges)], pch = c(rep(16, nrow(SIRSample$training)), rep(4, nrow(SIRSample$validation))))

# Train the emulators using the emulator_from_data defaults
ems_wave1 <- emulator_from_data(SIRSample$training, names(targets), ranges)
# Examine the emulator structure
ems_wave1$nI

# Demonstrating the behaviour of emulators at training points
ems_wave1$nS$get_exp(SIRSample$training[c(7, 19), ])
SIRSample$training[c(7, 19), 'nS']
ems_wave1$nS$get_cov(SIRSample$training[c(7, 19), ])

# Diagnostics on the emulators
validation_diagnostics(ems_wave1, targets, SIRSample$validation)

############### PLOTTING ###############

# Plot of emulator expectations
emulator_plot(ems_wave1)

# Demonstration of the difference between an untrained and trained emulator
plot(ems_wave1$nR)
plot(ems_wave1$nR$o_em)

# nth implausibility plot
emulator_plot(ems_wave1$nS, plot_type = 'imp', targets = targets, ppd = 40)
emulator_plot(ems_wave1, plot_type = "nimp", targets = targets, ppd = 40, cb = TRUE)

# Plot lattice: minimum implausibility and optical depth
plot_lattice(ems_wave1, targets, ppd = 35)

# Diagnostic for 'optimal' cutoff
space_removed(ems_wave1, targets, ppd = 20) +
  geom_vline(xintercept = 3, lty = 2)

# Same idea, but w.r.t. emulator variance
# space_removed(ems_wave1, targets, ppd = 20, modified = 'var', u_mod = seq(0.5, 1.5, by = 0.1))

############# PROPOSING POINTS ####################
# Simple generate_new_runs call
proposal1 <- generate_new_runs(ems_wave1, 500, targets)

# Plotting the proposal
plot_wrap(proposal1, ranges)

############# MULTIWAVE PLOTTING ################

# Plot the outputs across waves
wave_values(SIRMultiWaveData, targets, l_wid = 0.8)

# Examining the emulator uncertainties
sapply(SIRMultiWaveEmulators[[3]], function(x) x$u_sigma)

# Simulator plot: checking that we can hit all targets with some model runs
simulator_plot(SIRMultiWaveData, targets, barcol = 'black')

sessionInfo()
