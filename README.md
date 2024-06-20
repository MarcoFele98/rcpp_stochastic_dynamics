# Background
R is very slow and C++ is very fast. Despite the speed gains, a cumbersone C++ aspect is running simulation replicates or exploring a parameter space, becuase you have to pass and read external text files and/or rely on bash scripting. Rcpp gives me the speed of C++, and I can swiftly use R to interface with my simulation to run replicates. Rcpp also has r-like functions that easen writing code, but it is a bit trickier to do fancy C++ things like meta-programming. For info and rcpp installation see https://teuder.github.io/rcpp4everyone_en/. 

Here I implement the Gillespie algorithm in rcpp, and use it to demonstrate oscillations and mulitstability. The reason I did this is to force me to think about the microscopic reactions needed to obtain interesting dyamical outcomes. I still do not know whether it possible to design micro-level interactions starting from macro-level dynamical desiderata, so in this document I only report observations. Examples of oscillations are 1) non-transitive interactions with a rock-paper-scissors game, 2) the Lotka-Volterra predator-prey system. Examples of multistability are 1) a voter model that includes second order interactions for collective decision-making, 2) a honey-bee-inspired collective decision-making model based on cross-inhibition, 3) a gene switch in which protein transcription is governed by an autocatalytic sytem regulated by protein, mRNA, and promoter number. 

## Rock-paper-scissors
An simple way to obtain oscillations is through non-transitive interactions. A transitive interactions is one where A wins on B, B wins on C, and hence A wins on C. Rock-paper-scissors is non-transitive so the last interaction is inverted. The microscopic kinetics are A + B ⟶ 2B, B + C ⟶ 2C, C + A ⟶ 2A. Note that reaction rates are not shown in all of the examples provided. The dynamics is somewhat disappointing (in my opinion) becuase unstable. 

## Lotka-volterra predator-prey
The microscopic kinetics are prey (X) clonal reproduction X ⟶ 2X, prey death by predator (Y) X + Y ⟶ Y, predator reproduction by eating prey X + Y ⟶ 2Y, and predator death Y ⟶ $\varnothing$. The oscillations are neutrally stable, meaning the system has infinte attraction cycles instead of just one. This is caused by eigenvalues hanging out at the boundary between stablity and instability. I think that this means that in a stochastic system the oscillation cycle will randomly drift, increasing the probability that predators or preys will eventually go extinct.

## Second order interactions collective decision-making model
There are two options between which the population has to chose. The first microscopic reaction is two individuals that prefer X convert one that prefers Y (could be seen as the power of a social majority) 2X + Y ⟶ 3X. The order of this reaction is two, since the stoichiometry of the greater number of reagents involved is two. It is the order of this reaction that gives rise to mulitstability. Another reaction is that X spontaneusly switches opionion X ⟶ Y. The full set of reactions includes ones were X and Y are inveretd (symmetric reactions). The population reaches consenus, where everyone either prefers option X or option Y. 

## Cross-inhibition collective decision-making model
Instead of increasing the order of the reactions like in the previous model, multistability can be achived by adding reagents. In this model (from a cool paper https://www.nature.com/articles/s42005-023-01345-3) there are individuals that prefer X, individuals that prefer Y, and undecided individuals U. Reactions are conversion X + U ⟶ X, inhibition X + Y ⟶ X + U, and some random switching X ⟶ Y and X ⟶ U (again, the full set of reactions are symmetric). I show the bifurcation diagram where on the x axis there is the rate of random switching X ⟶ Y, on the y axis the population state, and the color rappresents the probability the population is found in that state. Multistability is only found for low values of random switching.

## Gene switch
Here I present an autocatalytic gene switch, where protein transcription can be turned on or off depending on a external input of proteins or mRNA. Multistability is an awsome way to have a fixed genetic code respond in a flexible way to the environemnt, virtually increasing the number of phenotypes expressed. In this example, multistability is achived through a type III functional response. There are four molecules involved in this switch: a protein X, a mRNA Y, a promoter P, and a promoter PX2 where two X proteins are binded to its allosteric regulation site, enhancing mRNA synthesis. The reactions are protein degradation X ⟶ $\varnothing$, mRNA degradation Y ⟶ $\varnothing$, and protein transcription by the mRNA Y ⟶ Y + X. The type III functional response is given by the association and dissociation of the proteins to the promoter P + 2X ⟶ PX2 and PX2 ⟶ P + 2X. mRNA synthesis by the promoter is given by PX2 ⟶ PX2 + Y, which closes the autocatalytic loop. This system can be reduced from four to only two dimensions (X and Y) by assuming steady state (rate of change = 0) of the promoter in both its forms (P and PX2). This is sensible if the dynamics of the type III functional response are much faster than the others, i.e., there is a separation of timescales. With this approximation, mRNA synthesis will be proportional to PX2, and $$PX2 = frac{x^2}{k_{-}/k{+}+x^2}$$. Nevertheless all the reactions still occur from a microscopic point of view. I compare the stochastic dynamics in the phase plot to the system of differential equations derived from steady steate approximation and infinite size assumption, which is used to find the mRNA and protein nullclines and the three system equilibria. Interestingly, the stochastic dynamics are comparable to the deterministic ones even if association and dissociation rates are comparable to the rates of the other reactions. 
