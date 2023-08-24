# Chaotic-APSO
This is a CAPSO algorithm specifically adjusted to the optimization problem concerning the electromagnetic cloaking of a perfectly electric conducting or dielectric core. The considered cloaking model can be provided upon request; see below.

This code file is accompanied by many comments providing clarity regarding the complex optimization problem it has provided feasible solutions for. The code can be adjusted to work for any objective function, and the test function of a sphere is provided, alongside with basic initialization and boundary checking functions. However, the focus of this experiment is not related to benchmark testing. If in search for an APSO demo, you could check ref. [2].

In case you are interested in testing the code alongside the model or seeing the results of the said experiment, please contact us:

Nikolaos Tsitsas, Aristotle University of Thessaloniki: ntsitsas@csd.auth.gr
Alkmini Michaloglou, Aristotle University of Thessaloniki: malkmini@csd.auth.gr
The CAPSO algorithm is a proposed variant of the APSO algorithm, and it utilizes chaotic maps for parameter tuning, and useful randomness. Reference:

[1] Gandomi, A.H.; Yun, G.J.; Yang, X.S.; Talatahari, S. Chaos-enhanced accelerated particle swarm optimization. Communications in Nonlinear Science and Numerical Simulation 2013,18, 327–340.

The code is built upon XS Yang's APSO demo available at (given the github upload date):

[2] XS Yang (2021). Accelerated Particle Swarm Optimization (APSO) (https://www.mathworks.com/matlabcentral/fileexchange/74766-accelerated-particle-swarm-optimization-apso), MATLAB Central File Exchange. Retrieved July 21, 2021.

APSO code and description can also be found in (and the several other editions of this book):

[3] Xin-She Yang, Nature-Inspired Metaheuristic Algorithms, Second Edition, Luniver Press, (2010). www.luniver.com

2022 and 2023 Changes: This has been edited and many commented out unused variables and redundant comments have been sorted out.
