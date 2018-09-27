# NC-cells

Introduction 

This code simulates the migration of neural crest cells. There are two phenotypes of cells, namely leaders and followers. The main difference between these two cell phenotypes is the mechanism driving motility and invasion: the leaders follow a cell-induced gradient of a chemoattractant, while the trailing cells follow the leaders. We use a two-dimensional individual-based model to represent cells and a continuum reaction-diffusion model for the dynamics of the chemoattractant. We assume two guidance mechansisms from leaders to followers: by contacts or via creation of tunnels. 

We provide two versions of the code: NC-model, where the interactions between the cells in NC-model are only by contact; NC-model-tunnelling where interactions are by contact and tunnelling mechanisms


The guidelines how to install Aboria can be found https://martinjrobins.github.io/Aboria/aboria/installation_and_getting_started.html.
