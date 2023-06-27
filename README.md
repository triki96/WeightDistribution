# WeightDistribution
Estimators for the weight distribution of linear codes. 
The repository consists of three folders: LDPC, RandomCodes and Utils.


## Utils
  Contains support code for performing weight estimation tests.
  
## RandomCodes
The folder contains two executable files, *estimatorsComparison.sage* and *progressiveComparison.sage*.

* *estimatorsComparison.sage* allows to choose a weight range $[w_{min},w_{max}]$ within which we want to estimate the weight distribution for the code under analysis. For each $w \in [w_{min},w_{max}]$ it is possible to set the parameter *max_attempt*, i.e. the maximum number of calls we want to make to ISD. The code will then make *max_attempt* attempts at ISD for every $w$ in the choosen set.

* *progressiveComparison.sage* allows to choose a fixed weight ($w$) and a maximum number of calls to ISD (*max_attempt*). The code then progressively shows how the estimates of $N_C(w)$ improve as the number of calls to ISD increases.

## LDPC
The LDPC folder contains the code relating to the comparison between the new estimators and the old one. Tests are performed on LDPC codes taken from https://www.inference.org.uk/mackay/CodesFiles.html. The folder contains one executable file, *estimatorsComparison.sage*, which allows to choose a LDPC code among the available ones, together with Stern's parameters $w,p,l$ and the number $t$ of iterations of the test.
The code will then make $t$ attempts at ISD for that specific code, with the choosen parameters set.
