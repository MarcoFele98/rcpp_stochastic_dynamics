# Background
R is very slow and C++ is very fast. Despite the speed gains, a cumbersone C++ aspect is running simulation replicates or exploring a parameter space, becuase you have to pass and read external text files and/or rely on bash scripting. Rcpp gives me the speed of C++, and I can swiftly use R to interface with my simulation to run replicates. Rcpp also has r-like functions that easen writing code, but it is a bit trickier to do fancy C++ things like meta-programming. For info and rcpp installation see https://teuder.github.io/rcpp4everyone_en/. 

Here I implement the Gillespie algorithm in rcpp, and use it to demonstrate oscillations and mulitstability. The reason I did this is to force me to think about the microscopic interactions needed to obtain interesting dyamical outcomes. I still do not know whether it possible to design micro-level interactions starting from macro-level dynamical desiderata. Examples of oscillations are 1) non-transitive interactions with a rock-paper-scissors game, 2) the Lotka-Volterra predator-prey system. Examples of multistability are 1) a voter model that includes higher order interactions for collective decision-making, 2) a honey-bee-inspired collective decision-making model based on cross-inhibition, 3) a gene switch in which protein transcription is governed by an autocatalytic sytem regulated by protein, mRNA, and promoter number. 

## Rock-paper-scissors
An simple way to obtain oscillations is through non-transitive interactions. A transitive interactions is one where A wins on B, B wins on C, and hence A wins on C. Rock-paper-scissors is non-transitive and the last interaction is inverted. The microscopic kinetics are A + B ⟶ 2B, B + C ⟶ 2C, C + A ⟶ 2A. The dynamics is somewhat disappointing (in my opinion) becuase unstable. 

## Lotka-volterra predator-prey
The microscopic kinetics are X ⟶ 2X, X + Y ⟶ Y, X + Y ⟶ 2Y, Y ⟶ $$\varnothing$$
