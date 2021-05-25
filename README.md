# R_sim_scripts

###
# by M. Chas Sumner (mcnu)
 Used for analyses in Sumner, Torrisi, Brickner, Brickner (2021, eLife)
###

# TOC

## License (MIT license)

## multimsd.R

 Function that takes multiple merged MSD files from smoothMSD.R to plot MSD for different conditions within the same figure (primary method of MSD plotting in the paper)

## Particle_Tracker.Rmd

 Markdown with two functions to segment and track single particles within images across multiple .tif files as recognized by the .nd output from Metamorph acquisition software

## peri_dist_trfr.Rmd

 Markdown similar to Particle_Tracker but with added functionality to detect distance from the nuclear boundary, and can utilize either a separate color channel for the membrane or the same color as the tracked particle soluble fluorophore is concentrated in the nucleus.

## SimClustFunc.R

 Function used to perform pairwise analysis of two tracked particles and produce the data for figures 6 and 7

## SimLin_Repos.R

 Function used to analyze and plot repositioning data in figure 5. Takes considerably longer trajectories and requires a separate look up table containing the frame where first contact between the chromatin particle and nuclear boundary occur

## Simulated_MSD.R

 Initially ad MSD calculator for only simulated output from YGRW, now standard MSD function used in our article

## smoothMSD.R

 General plot function to visualize MSD for a single condition then merge all the individual MSD files together for multiMSD.R input

## Trajectory_Analysis.Rmd

 Initially the markdown used to analyze simulation output, this markdown grew to source and run all functions used in the paper and can perform all displacement data analysis. In addition, this markdown contains heuristics and tests used to compare simulations to their respective in vivo particles and  conditions




