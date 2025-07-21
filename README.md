# SLiM-code
The script options here are the stepping stones for my overall project to be able to simulate oncogenesis, the first goal was to recreate a case study patient with clonal hematopoiesis
For beginners to SLiM, I'd highly recommend visiting the messerlab.org website and downloading the SLiM manual. 
SLiM download is also avaliable on the messerlab.org website 

Key: 

NM: Non-mutant

M: Mutant 

NWF: non wright-fischer

WF: Wright-fischer 

*** Tip for python when doing a newick conversion 
in order for python to be able to run a conversion of the tree file make sure the following add-ons are downloaded:
- homebrew: https://brew.sh/
- pip: https://pip.pypa.io/en/stable/installation/  (or conda) 
- pyslim: https://tskit.dev/pyslim/docs/latest/installation.html
- tskit: https://tskit.dev/tskit/docs/stable/installation.html
  
All are open source on the web

*** You cannot add mutations to a individual >0 year old in a non WF simulation, this can be negated by specifically adding a child to the simulation at the year you wish to add the mutation (check case study code), however there is potential for this to interfer with the accuracy of the final phylogeny because of unclear familial orgin. 

*** For clarity in my final phylogenies many simulation scripts show mutation rate at a rate of 0, this is intentional and does not effect any new mutations added to the population in later generations 

happy simulating :) 
