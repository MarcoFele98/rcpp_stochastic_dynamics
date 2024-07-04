# Background
R is very slow and C++ is very fast. Despite the speed gains, a cumbersone aspect of C++ is running simulation replicates or exploring a parameter space, becuase you have to pass and read external text files and/or rely on bash scripting. Rcpp gives me the speed of C++, and I can swiftly use R to interface with my simulation to run replicates. Rcpp also has r-like functions that easen writing code, but it is a bit trickier to do fancy C++ things like meta-programming. For info and rcpp installation see https://teuder.github.io/rcpp4everyone_en/. 

Here I implement a general version of the Gillespie algorithm in rcpp, and use it to demonstrate oscillations and multistability. The reason I did this is to force me to think about the microscopic reactions needed to obtain interesting dyamical outcomes. I still do not know whether it possible to design micro-level interactions starting from macro-level dynamical desiderata. I implement the Gillespie algorithm to simulate five classic models from theoretical ecology, collective decision-making, and systems biology.

## Rock-paper-scissors
An simple way to obtain oscillations is through non-transitive interactions. A transitive interactions is one where A wins on B, B wins on C, and hence A wins on C. Rock-paper-scissors is non-transitive so the last interaction is inverted. The microscopic reactions are A + B ⟶ 2B, B + C ⟶ 2C, C + A ⟶ 2A. Note that reaction rates are not shown in all of the examples provided. The dynamics is unstable. 
![rps](https://github.com/MarcoFele98/rcpp_stochastic_dynamics/assets/122376407/eeb13548-6d58-40bc-be4b-0a61f5d6e90b)

## Lotka-volterra predator-prey
The microscopic reactions are prey clonal reproduction X ⟶ 2X, prey death by predator X + Y ⟶ Y, predator reproduction by eating prey X + Y ⟶ 2Y, and predator death Y ⟶ $\varnothing$. The oscillations are neutrally stable, meaning the system has infinte attraction cycles. This is caused by eigenvalues hanging out at the boundary between stability and instability. I think that this means that in a stochastic system the oscillation cycle will randomly drift in magnitude, increasing the probability that predators or preys will eventually go extinct.
![predator_prey](https://github.com/MarcoFele98/rcpp_stochastic_dynamics/assets/122376407/921df9de-b5a4-4315-bb02-e787eeddfa8a)

## Second order interactions collective decision-making model
Assume that a population of individuals has to choose between two available options. The first microscopic reaction is that two individuals that prefer X convert one that prefers Y (could be seen as social conformity) 2X + Y ⟶ 3X. The order of this reaction is two, since the stoichiometry of the highest number of reagents involved is two. It is the order of this reaction that gives rise to mulitstability. Another reaction is that X spontaneusly switches opionion X ⟶ Y. The full set of reactions includes ones were X and Y are inverted (symmetric reactions). The population reaches consenus, where everyone either prefers option X or option Y. 
![voter](https://github.com/MarcoFele98/rcpp_stochastic_dynamics/assets/122376407/b97c34e7-cbc9-4a34-beca-623971c46017)

## Cross-inhibition collective decision-making model
Instead of increasing the order of the reactions like in the previous model, multistability can be achived by adding reagents. In this model (from a cool paper https://www.nature.com/articles/s42005-023-01345-3) there are individuals that prefer X, individuals that prefer Y, and undecided individuals U. Reactions are conversion X + U ⟶ X, inhibition X + Y ⟶ X + U, and some random switching X ⟶ Y and X ⟶ U (again, the full set of reactions are symmetric). I show the bifurcation diagram where on the x axis there is the rate of random switching X ⟶ Y, on the y axis the population state, and the color rappresents the probability the population is found in that state. Multistability is only found for low values of random switching.
![biff](https://github.com/MarcoFele98/rcpp_stochastic_dynamics/assets/122376407/778db95e-34fe-42b9-80d8-be45b335e569)

## Gene switch
Here I present an autocatalytic gene switch, where protein expression can be turned on or off depending on a external input of proteins or mRNA. Multistability is an awsome way to have a fixed genetic code respond in a flexible way to the environemnt, virtually increasing the number of phenotypes expressed. In this example, multistability is achived through a type III functional response. There are four molecules involved in this switch: a protein X, a mRNA Y, a promoter P, and a promoter PX2 where two X proteins are binded to the allosteric regulation site of P, enhancing mRNA transcription. The reactions are protein degradation X ⟶ $\varnothing$, mRNA degradation Y ⟶ $\varnothing$, and protein translation by the mRNA Y ⟶ Y + X. The type III functional response is given by the association and dissociation of the proteins to the promoter P + 2X ⟶ PX2 (rate $k_{+}$) and PX2 ⟶ P + 2X (rate $k_{-}$). mRNA transcription by the promoter is given by PX2 ⟶ PX2 + Y, which closes the autocatalytic loop. I compare the stochastic dynamics in the phase plot to the system of differential equations derived from steady steate approximation and infinite size assumption, which is used to find the mRNA and protein nullclines and the three system equilibria. This system can be reduced from four to only two dimensions (X and Y) by assuming steady state (rate of change = 0) of the promoter in both its forms (P and PX2). This is sensible if the dynamics of the type III functional response are much faster than the others, i.e., there is a separation of timescales. With this approximation, we calculate the equilibrium number of PX2 as $\frac{x^2}{k_{-}/k_{+}+x^2}$. Interestingly, the stochastic dynamics are comparable to the deterministic predictions even if association and dissociation rates are comparable to the rates of the other reactions. 
![gene_switch](https://github.com/MarcoFele98/rcpp_stochastic_dynamics/assets/122376407/4e3b52cf-d597-414c-b31b-fbf26c2021e8)

