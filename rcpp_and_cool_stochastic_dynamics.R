library(Rcpp)
library(ggplot2)
library(data.table)

sourceCpp("noise/rcpp_stochastic_dynamics/gillespie_algorithm.cpp")

load('noise/rcpp_stochastic_dynamics/data.RData')

# Oscillations ----
## Rock paper scissors (unstable internal equilibrium with oscillations - infinite size: real part of Jacobian leading eigenvalue at internal equilibrium > 0; complex part != 0) ----
rock_paper_scissors <- gillespie_simulation(states = c(100,  # rock
                                                       100,  # paper
                                                       100), # scissors
                                            max_simulation_time = 10000,
                                            number_reagents = 3,
                                            number_reactions = 3,
                                            omega = 300,
                                            reaction_rates = c(1, 1, 1),
                                            # columns: rock, paper, scissors
                                            stoichiometry_reagents = matrix(c(1, 1, 0,
                                                                              0, 1, 1,
                                                                              1, 0, 1), 
                                                                            ncol = 3, 
                                                                            byrow = T),
                                            stoichiometry_products = matrix(c(2, 0, 0,
                                                                              0, 2, 0,
                                                                              0, 0, 2), 
                                                                            ncol = 3,
                                                                            byrow = T))

ggplot(as.data.table(rock_paper_scissors)) +
  geom_step(aes(time, species_0), color = "red") +
  geom_step(aes(time, species_1), color = "blue") +
  geom_step(aes(time, species_2), color = "green") +
  ylab("Number") +
  ggtitle("Rock-paper-scissors")

ggsave("noise/rcpp_stochastic_dynamics/figures/rps.png",
       height = 5,
       width = 7,
       bg = "white")

## Lotka-Voltera predator-prey (stable limit cycles / neutral stability - infinite size: real part of Jacobian leading eigenvalue at internal equilibrium = 0; complex part != 0) ----
lotka_volterra_predator_prey <- gillespie_simulation(states = c(1000,  # prey
                                                                1000), # predators
                                                     max_simulation_time = 500,
                                                     number_reagents = 2,
                                                     number_reactions = 4,
                                                     omega = 1000,
                                                     reaction_rates = c(0.5,  # prey reproduction
                                                                        0.2,  # prey capture
                                                                        0.3,  # predator reproduction
                                                                        0.2), # predator death
                                                     # columns: prey, predator
                                                     stoichiometry_reagents = matrix(c(1, 0,  # prey reproduction
                                                                                       1, 1,  # prey capture
                                                                                       1, 1,  # predator reproduction
                                                                                       0, 1), # predator death
                                                                                     ncol = 2, 
                                                                                     byrow = T) ,
                                                     stoichiometry_products = matrix(c(2, 0,  # prey reproduction
                                                                                       0, 1,  # prey capture
                                                                                       0, 2,  # predator reproduction
                                                                                       0, 0), # predator death
                                                                                     ncol = 2,
                                                                                     byrow = T))

ggplot(as.data.table(lotka_volterra_predator_prey)) +
  geom_step(aes(time, species_0), color = "blue") +
  geom_step(aes(time, species_1), color = "red") +
  ylab("Number") +
  ggtitle("Lotka-Volterra predator-prey")

ggsave("noise/rcpp_stochastic_dynamics/figures/predator_prey.png",
       height = 5,
       width = 7,
       bg = "white")

# Multi-stability ----
## Ternary interactions collective decision-making model (infinite size: real part of Jacobian leading eigenvalue at internal equilibrium > 0; complex part = 0) ----
data_ternary <- data.table()
for(replicate in 1:100) {
  print(replicate)
  ternary_model <- gillespie_simulation(states = c(100,  # prefer A
                                                   100), # prefer B
                                        max_simulation_time = 50,
                                        number_reagents = 2,
                                        number_reactions = 4,
                                        omega = 200,
                                        reaction_rates = c(1,     # conversion: 2A + B -> 3A
                                                           1,     # conversion: A + 2B -> 3B
                                                           0.01,  # noise: A -> B
                                                           0.01), # noise: B -> A
                                        # columns: prefer A, prefer B
                                        stoichiometry_reagents = matrix(c(2, 1,  # conversion: 2A + B -> 3A
                                                                          1, 2,  # conversion: A + 2B -> 3B
                                                                          1, 0,  # noise: A -> B
                                                                          0, 1), # noise: B -> A
                                                                        ncol = 2, 
                                                                        byrow = T),
                                        stoichiometry_products = matrix(c(3, 0,  # conversion: 2A + B -> 3A
                                                                          0, 3,  # conversion: A + 2B -> 3B
                                                                          0, 1,  # noise: A -> B
                                                                          1, 0), # noise: B -> A
                                                                        ncol = 2, 
                                                                        byrow = T))
  
  data_ternary <- rbind(data_ternary,
                        as.data.table(ternary_model)[, ":="(replicate = replicate)])
}

ggplot(data_ternary) +
  geom_step(aes(time, species_0, group = replicate, color = replicate)) +
  ylab("Number prefering option X") +
  ggtitle("Multistability of ternary model (species 1)") 

ggsave("noise/rcpp_stochastic_dynamics/figures/voter.png",
       height = 5,
       width = 7,
       bg = "white")


## Cross-inhibition collective decision-making model (infinite size: real part of Jacobian leading eigenvalue at internal equilibrium > 0 for small sigma; complex part = 0) ----
bifurcation_anlaysis <- data.table()
for(sigma in seq(0, 0.4, 0.01)) {
  print(sigma)
  for(replicate in 1:100) {
    print(replicate)
    cross_inhibition <- gillespie_simulation(states = c(100,  # prefer A
                                                        100,  # prefer B
                                                        0),   # undecided U 
                                             max_simulation_time = 200,
                                             number_reagents = 3,
                                             number_reactions = 8,
                                             omega = 200,
                                             reaction_rates = c(1,       # inhibition: A + B -> A + U
                                                                1,       # inhibition: B + A -> B + U
                                                                1,       # conversion: A + U -> A + A
                                                                1,       # conversion: B + U -> B + B
                                                                sigma,   # noise: A -> B
                                                                sigma,   # noise: B -> A
                                                                0.01,    # noise: A -> U
                                                                0.01),   # noise: B -> U
                                             # columns: prefer A, prefer B, undecided U 
                                             stoichiometry_reagents = matrix(c(1, 1, 0,  # inhibition: A + B -> A + U
                                                                               1, 1, 0,  # inhibition: B + A -> B + U
                                                                               1, 0, 1,  # conversion: A + U -> A + A
                                                                               0, 1, 1,  # conversion: B + U -> B + B
                                                                               1, 0, 0,  # noise: A -> B
                                                                               0, 1, 0,  # noise: B -> A
                                                                               1, 0, 0,  # noise: A -> U
                                                                               0, 1, 0), # noise: B -> U
                                                                             ncol = 3, 
                                                                             byrow = T),
                                             stoichiometry_products = matrix(c(1, 0, 1,  # inhibition: A + B -> A + U
                                                                               0, 1, 1,  # inhibition: B + A -> B + U
                                                                               2, 0, 0,  # conversion: A + U -> A + A
                                                                               0, 2, 0,  # conversion: B + U -> B + B
                                                                               0, 1, 0,  # noise: A -> B
                                                                               1, 0, 0,  # noise: B -> A
                                                                               0, 0, 1,  # noise: A -> U
                                                                               0, 0, 1), # noise: B -> U
                                                                             ncol = 3, 
                                                                             byrow = T))
    
    bifurcation_anlaysis <- rbind(bifurcation_anlaysis,
                                  as.data.table(cross_inhibition)[time > 25 # steady state only, nice and cleean
                                  ][, ":="(diff = species_0 - species_1,
                                           duration = c(waiting_time[-1], NA))
                                  ][-(.N) # eliminate last state of which we do not know the duration
                                  ][, ":="(max_duration = sum(duration))
                                  ][, .(prob = sum(duration) / unique(max_duration),
                                        replicate = replicate,
                                        sigma = sigma),
                                    by = diff])
  }
}

bifurcation_anlaysis_s <- bifurcation_anlaysis[, .(probability = mean(prob)),
                                               by = list(sigma, diff)]

ggplot(bifurcation_anlaysis_s) +
  geom_tile(aes(sigma, diff, fill = log(probability))) +
  geom_hline(aes(yintercept = 0), linewidth = 2, lty = "dashed") +
  scale_fill_viridis_c(option = "plasma") +
  #ggtitle("Bifurcation diagram") +
  ylab("Preferences (from 200 prefering Y to 200 prefering X)") +
  xlab("Switch rate")

ggsave("noise/rcpp_stochastic_dynamics/figures/biff.png",
       height = 5,
       width = 7,
       bg = "white")
  
## Genetic switch ----
a <- 0.65
b <- 0.65
omega <- 500

# Infinte size approximation 
# This is the solution of deterministic system of differential equations:
# dx/dt = -a * x + y
# dy/dt = x^2 / (1 + x^2) - b * y
x_isocline <- data.table(x = seq(0, 2.5, l = 1000),
                         y = a * seq(0, 2.5, l = 1000))
y_isocline <- data.table(x = seq(0, 2.5, l = 1000),
                         y = seq(0, 2.5, l = 1000)^2 / ((seq(0, 2.5, l = 1000)^2 + 1)*b))

id <- 1
gene_switch_data <- data.table()
for(x_initial in seq(0.2, 2.5, l = 5)) {
  for(y_initial in seq(0.2, 1.5, l = 5)) {
    print(x_initial)
    print(y_initial)
    gene_switch <- data.table()
    gene_switch <- gillespie_simulation(states = c(x_initial * omega,  # protein X
                                                   y_initial * omega,  # mRNA Y
                                                   500,  # promoter P
                                                   0),   # promoter and protein dimer PX2
                                        max_simulation_time = 20,
                                        number_reagents = 4,
                                        number_reactions = 6,
                                        omega = omega,
                                        reaction_rates = c(a,  # X -> 0
                                                           b,  # Y -> 0
                                                           1,    # Y -> X + Y
                                                           # functional response type III (outputs equilbrium proportion of P that are binded with X2 so that it can promote synthesis of Y)
                                                           10,    # P + 2X -> PX2
                                                           10,    # PX2 -> P + 2X
                                                           # effect of PX2 
                                                           1    # PX2 -> PX2 + Y
                                        ),  
                                        # columns: protein X, mRNA Y, promoter P, promoter and protein dimer PX2
                                        stoichiometry_reagents = matrix(c(1, 0, 0, 0, # X -> 0
                                                                          0, 1, 0, 0, # Y -> 0
                                                                          0, 1, 0, 0, # Y -> X + Y
                                                                          # functional response type III (outputs equilbrium proportion of P that are binded with X2 so that it can promote synthesis of Y)
                                                                          2, 0, 1, 0, # P + 2X -> PX2
                                                                          0, 0, 0, 1, # PX2 -> P + 2X
                                                                          # effect of PX2 
                                                                          0, 0, 0, 1  # PX2 -> PX2 + Y
                                        ),
                                        ncol = 4, 
                                        byrow = T),
                                        stoichiometry_products = matrix(c(0, 0, 0, 0, # X -> 0
                                                                          0, 0, 0, 0, # Y -> 0
                                                                          1, 1, 0, 0, # Y -> X + Y
                                                                          # functional response type III (outputs equilbrium proportion of P that are binded with X2 so that it can promote synthesis of Y)
                                                                          0, 0, 0, 1, # P + 2X -> PX2
                                                                          2, 0, 1, 0, # PX2 -> P + 2X
                                                                          # effect of PX2 
                                                                          0, 1, 0, 1  # PX2 -> PX2 + Y
                                        ), 
                                        ncol = 4, 
                                        byrow = T))
    
    gene_switch_data <- rbind(gene_switch_data,
                              as.data.table(gene_switch)[, ":="(x_initial = x_initial,
                                                                y_initial = y_initial,
                                                                id = id)]) 
    id <- id + 1
  }
}

plot <- ggplot(gene_switch_data) +
  geom_point(aes(species_0,
                species_1, 
                color = time),
             size = 0.1) +
  # protein isocline
  geom_line(data = x_isocline * omega,
            aes(x, y),
            color = "red", linewidth = 1) +
  # mRNA isocline
  geom_line(data = y_isocline * omega,
            aes(x, y),
            color = "blue", linewidth = 1) +
  # stable equilibrium (SWITCH OFF)
  geom_point(aes(0, 0), 
             size = 3) +
  # unstable equilibrium 
  geom_point(aes((1 - sqrt(1 - (4 * a^2 * b^2))) / (2 * a * b) * omega, 
                 (1 - sqrt(1 - (4 * a^2 * b^2))) / (2 * b) * omega), 
             size = 3, 
             fill = "white",
             shape = 21) +
  # stable equilibrium (SWITCH ON)
  geom_point(aes((1 + sqrt(1 - (4 * a^2 * b^2))) / (2 * a * b) * omega, 
                 (1 + sqrt(1 - (4 * a^2 * b^2))) / (2 * b) * omega),
             size = 3) +
  scale_color_viridis_c() +
  ggtitle("Phase plot") +
  xlab("Protein number") +
  ylab("mRNA number") +
  facet_wrap(~id)

ggsave(plot,
       "noise/rcpp_stochastic_dynamics/figures/gene_switch.png",
       height = 30,
       width = 30,
       bg = "white")
