# Plot_SSPSE
Additional plotting functions to complement sspse R package

- rtp.simple() creates a simplified version of the Reingold-Tilford Plot that places recruitment chains horizontally with seeds all at the same level.
- post.plot() creates a ggplot version of posterior distribution for N.
- visplots() creates ggplot versions of visibility diagnostic plots. The following types of plot can be specified: 
  - "meanvis" mean of visibility distribution
  - "sdvis" standard deviation of visibility distribution
  - "visdist" visibility distribution
  - "visdistoverlay" visibility distribution overlaid with observed (augmented) network sizes
  - "visall" visibility vs network size by person
  - "pvdist" visibility distributions for random sample of 20
