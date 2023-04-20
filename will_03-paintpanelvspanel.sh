#!/bin/sh

## So, in ancientvsancient.cp we will define the parameters for chromopainter, the phase and recombination files, and then importantly, the donoridfile and popidfile
## The donoridfile tells the code who the donor and recipient populations are
## The popidfile tells which individuals are in which populations (and which are to be used in the analysis)
## These are exactly as we'd use for GLOBETROTTER

bash will_paint_withinpanel.sh full `pwd`/ancientvsancient will_ancientvsancient.cp

## I think you:
## Run this
## Do the jobs it gives you
## Run it again to combine the output
## Clean up if desired
