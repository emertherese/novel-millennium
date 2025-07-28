This repository contains the R code to replicate all analyses and figures from our paper:

> Cunningham, E. T., T. L. Staples, I. R. Butler, M. Lepore, H. Markham Summers, G. Roff, and J. M. Pandolfi. 2024. Taxonomic novelty emerges more frequently and independently of functional novelty in historical coral communities. Oikos:e10912.

#### Data

This project considers historical corals from reef matrix cores, collected across the east of so-called Australia. Thank you very much to Drs Hannah Markham Summers, George Roff, Mauro Lepore, and Ian Butler for their field work and processing.

These cores constitute a larger *Queensland Quaternary Core Collection* funded by the Australian Research Council and National Environmental Research Program. 

#### Acknowledging Country

The corals we studied lived long lives in Sea Country. Over millennia, they have built incredible underwater places and felt changes to their homes. We pay our respect to corals, their Sea Country, and the Traditional Custodians across the Reef who are working towards a healthier future on Country.

We acknowledge Mandingalby Yidinji and Gungandji Sea Country (Frankland Islands), Bwgcolman Sea Country (Palm Islands), Woppaburra Sea Country (Keppel Islands), and Butchulla Sea Country (Hervey Bay).

#### Workflow

Open `novel-millennium.Rproj` and run scripts in sequence from `1-data-preparation.R`. Each script has also been designed to run out of sequence by calling on previous code, such that you could simply open and run `3-figures.R` and create figure outputs.

#### Graphical abstract

![Graphical abstract for this paper including key findings and coral imagery.](/abstract/graphical-abstract.png)

#### Any questions?

Please contact me ([e.cunningham@uqconnect.edu.au](mailto:e.cunningham@uqconnect.edu.au?subject=Novel millennium enquiry)) if you would like to learn more! I am happy to chat about anything, including the raw reef matrix core data that were processed into `sorted-matrices.RDS`.
