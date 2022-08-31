library(sspse)
library(ggplot2)
library(igraph)
library(ggsci)
library(ggpubr)
library(tidyr)

# Simplified version of Reingold-Tilford Plot that places recruitment chains horizontally with seeds all at the same level
# Has the option to color points based on a provided column in the data
# Cre: Katherine McLaughlin

## rdf = name of rds data frame
## char = name of column in data to color points by
## title = name of plot

## EXAMPLE
# library(RDS)
# data(fauxmadrona)
# rtp.simple(fauxmadrona, title="Example Plot", char="seed")

rtp.simple <- function (rdf, 
                        char = NULL, 
                        title = "Recruitment Plot") 
{
  rdf$seed <- get.seed.id(rdf)
  rdf$wave <- get.wave(rdf)
  rdf <- rdf[order(rdf$wave), ]
  tp.el <- matrix(c(get.rid(rdf), get.id(rdf)), ncol = 2)
  g <- graph_from_edgelist(tp.el)
  lyt <- layout.reingold.tilford(g)
  lyt <- lyt[-which.max(lyt[,2]), ]
  seed.names <- unique(rdf$seed)
  n.seeds <- length(seed.names)
  for (i in 2:n.seeds) {
    sid.p <- seed.names[i - 1]
    sid.c <- seed.names[i]
    cids.p <- which(rdf$seed == sid.p)
    cids.c <- which(rdf$seed == sid.c)
    max.x.p <- max(lyt[cids.p, 1])
    min.x.c <- min(lyt[cids.c, 1])
    to.add <- 0.5 - (min.x.c - max.x.p)
    sid.o <- seed.names[i:n.seeds]
    cids.o <- which(rdf$seed %in% sid.o)
    lyt[cids.o, 1] <- lyt[cids.o, 1] + to.add
  }
  
  s.mins <- c()
  s.maxes <- c()
  seeds <- unique(get.seed.id(rdf))
  for (i in 2:n.seeds) {
    ids <- which(rdf$seed == seeds[i])
    lyt[ids, 1] <- lyt[ids, 1] + (i - 1) * 2
    idsp <- which(rdf$seed == seeds[i - 1])
  }
  for (i in 1:n.seeds) {
    ids <- which(rdf$seed == seeds[i])
    s.mins[i] <- min(lyt[ids, 1])
    s.maxes[i] <- max(lyt[ids, 1])
  }
  
  n.seg <- dim(tp.el)[1] - n.seeds
  x0 <- rep(0, n.seg)
  y0 <- rep(0, n.seg)
  x1 <- rep(0, n.seg)
  y1 <- rep(0, n.seg)
  sd.id <- which(tp.el[,1] == "seed")
  tp.el.new <- tp.el[-sd.id, ]
  for (i in 1:n.seg) {
    x0[i] <- lyt[which(get.id(rdf) == tp.el.new[i, 1]), 1]
    y0[i] <- lyt[which(get.id(rdf) == tp.el.new[i, 1]), 2]
    x1[i] <- lyt[which(get.id(rdf) == tp.el.new[i, 2]), 1]
    y1[i] <- lyt[which(get.id(rdf) == tp.el.new[i, 2]), 2]
  }
  if (!is.null(char)) {
    charcol <- which(names(rdf) == char)
    if (length(charcol) == 0) {
      stop("Error: covariate column name does not exist")
    }
    n.cat <- length(unique(rdf[, charcol]))
    if (n.cat <= 10) {
      col.num <- pal_d3("category10", alpha=0.8)(n.cat)
    } else {
      col.num <- pal_d3("category20", alpha=0.8)(n.cat)
    }
    col.id <- col.num[as.numeric(factor(rdf[, charcol]))]
  } else {
    col.id <- rep("#2780E3", nrow(rdf))
  }
  n.wav <- max(rdf$wave)
  
  par(mar = c(0, 0.5, 1.5, 0.5))
  
  plot(lyt, axes = FALSE, main = title, xlab = "", ylab = "", 
       xlim = c(min(lyt[, 1]), max(lyt[, 1])+5)-5,
       ylim = c(-1.5, n.wav + 1), col="white", xpd=TRUE)
  
  segments(x0 = x0, y0 = y0, x1 = x1, y1 = y1, col = "grey50", 
           lwd = 1.5)
  points(lyt, col="black", bg = col.id, pch = 21)
  
  par(mar = c(5.1, 4.1, 4.1, 2.1))
}

# ggplot version of posterior distribution for N
# Cre: Katherine McLaughlin, Laura Gamble

## sspse.obj = name of fitted sspse object
## xmax = upper bound for x (on plot)
## band.post = controls smoothing
## ptitle = plot title

## EXAMPLE
# library(sspse)
# data(fauxmadrona)
# fit <- posteriorsize(fauxmadrona, median.prior.size=1000, visibility=FALSE)
# post.plot(fit, ptitle="Posterior Distribution for N")

post.plot <- function(sspse.obj, xmax=NULL, band.post = 1, ptitle=NULL){
  n = sspse.obj$n
  
  prior <- exp(sspse.obj$lpriorm - max(sspse.obj$lpriorm))
  denom <- sum(prior[((0:99) + (1:100))/200*length(prior)]*length(prior)/100)
  prior.n <- prior/denom
  
  sumry <- summary(sspse.obj, HPD.level = 0.9)
  df <- data.frame("x" =  n:(length(prior) + n - 1), 
                   "Prior" = prior.n)
  df2 <- data.frame("Median" = sumry$Median[2], 
                    "Mean" = sumry$Mean[2],
                    "5P" = sumry$'5%'[2],
                    "95P" = sumry$'95%'[2])
  dens <- data.frame("x" = density(sspse.obj$sample[,"N"], 
                                   adjust = band.post, from = n)$x, 
                     "y" = density(sspse.obj$sample[,"N"], 
                                   adjust = band.post, from = n)$y)
  
  if(is.null(xmax)) {
    x.lim <- sort(sspse.obj$sample[,"N"])[.975*10000]
  } else {
    x.lim <- xmax
  }
  y.lim <- max(dens$y, df$Prior)
  
  labs <- c("Posterior", "Prior", "Median", "Mean")
  labs.c <- c("Posterior" = "black", "Prior" = "black", "Median" = "darkred", 
              "Mean" = "green4")
  labs.l <- c("Posterior" = 1, "Prior" = 2, 
              "Median" = 1, "Mean" = 2)
  labs.txt <- data.frame("x" = c(n, df2$X5P, df2$Median,
                                 df2$Mean, df2$X95P), 
                         "y" = rep(-.05*y.lim, 5),
                         "label" = round(c(n, df2$X5P, df2$Median,
                                           df2$Mean, df2$X95P)))
  if (is.null(ptitle)) {
    ptitle <- "Posterior for Population Size"
  }
  
  p.plot <- ggplot(df, aes(x = x, y = Prior)) + 
    geom_line(aes(color = 'Prior', lty = 'Prior')) +
    geom_line(data = dens, aes(x = x, y = y, color = 'Posterior', 
                               lty = 'Posterior')) +
    geom_segment(aes(x = n, y = 0,
                     xend = n, yend = dens$y[1],
                     color = 'Posterior', 
                     lty = 'Posterior')) +
    geom_area(data = subset(dens, (x >= df2$X5P & x <= df2$X95P)), 
              aes(x = x, y = y, fill = '90% Interval')) + 
    geom_segment(aes(x = n, y = 0, xend = x.lim, yend = 0)) + 
    geom_segment(data = df2, aes(x = Median, y = 0, 
                                 xend = Median, yend = 1.1*y.lim,
                                 color = 'Median', lty = 'Median'), size = 1) +
    geom_segment(data = df2, aes(x = Mean, y = 0, 
                                 xend = Mean, yend = 1.1*y.lim,
                                 color = 'Mean', lty = 'Mean'), size = 1) + 
    xlim(c(0, x.lim)) +
    scale_colour_manual(name = "",
                        breaks = names(labs.c),
                        values = labs.c) +
    scale_linetype_manual(name = "",
                          breaks = names(labs.l),
                          values = labs.l) +
    scale_fill_manual(name = "", labels = "90% Interval", 
                      values = alpha("dodgerblue4", 0.5)) +
    scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) +
    geom_text(data = labs.txt, aes(x = x, y = y, label = label), size = 2,
              angle = 45, color = c("black", "dodgerblue4", "darkred", 
                                    "green4", "dodgerblue4")) +
    xlab("Population Size") + ylab("Density") +
    ggtitle(ptitle) + theme_bw()
  #theme(plot.title = element_text(hjust = 0.5))
  
  suppressWarnings(print(p.plot))
}


# ggplot version of visibility diagnostic plots
# Cre: Katherine McLaughlin

## sspse.obj = name of fitted sspse object
## type = vector containing one or more of "meanvis", "sdvis", "visdist", "visdistoverlay", "visall", "pvdist"

## EXAMPLE
# library(sspse)
# data(fauxmadrona)
# fitv <- posteriorsize(fauxmadrona, median.prior.size=1000, visibility=TRUE)
# visplots(fitv, type=c("visdistoverlay", "visall"))

visplots <- function(fit, type=c("visdistoverlay", "visall")) {
  
if (fit$visibility == FALSE) {
  stop("Error: visibility plots can only be made with a visibility sspse object.")
}

  if (length(type) < length(match(type, c("meanvis", "sdvis", "visdist", "visdistoverlay", "visall", "pvdist")))) {
    stop("Error: invalid plot type specified")
  }
  
  plist <- list()
  
  out <- data.frame(fit$sample)
  
  # Visibility mean
  if ("meanvis" %in% type) {
    meanvis <- ggplot(out, aes(x=mu)) + 
      geom_density(fill="black", alpha=0.5) +
      theme_bw() +
      labs(x="Mean visibility", title="Posterior for mean visibility")
    plist <- c(plist, list(meanvis))
  }
  
  # Visibility standard deviation
  if ("sdvis" %in% type) {
    sdvis <- ggplot(out, aes(x=sigma)) + 
      geom_density(fill="black", alpha=0.5) +
      theme_bw() +
      labs(x="s.d. visibility", title="Posterior for s.d. of the visibility")
    plist <- c(plist, list(sdvis))
  }
  
  # Visibility distribution (plain)
  df.visdist <- data.frame(Visibility=seq_along(fit$predictive.visibility),
                           Probability=fit$predictive.visibility)
  
  if ("visdist" %in% type) {
    visdist <- ggplot(df.visdist, aes(x=Visibility, y=Probability)) +
      geom_bar(stat="identity") +
      theme_bw() +
      labs(title="Visibility Distribution")
    plist <- c(plist, list(visdist))
  }
  
  
  # Visibility distribution with network sizes
  network.size <- as.numeric(fit$data[[attr(fit$data,"network.size.variable")]])
  
  # Augment the reported network size by the number of recruits and the recruiter (if not seed).
  # If network size is less than #recruits+recruiter, replace
  nr <- RDS::get.number.of.recruits(fit$data)
  is.seed <- (RDS::get.rid(fit$data)=="seed")
  network.size <- pmax(network.size,nr+!is.seed)
  
  
  remvalues <- is.na(network.size)
  if (sum(remvalues) > 0) {print("There is at least one missing network size.")}
  
  Kmax <- max(seq_along(fit$predictive.visibility))
  ns.prob <- tabulate(network.size[!remvalues],nbins=Kmax)
  ns.prob1 <- ns.prob/length(network.size)
  
  df.visdist$ns.prob <- ns.prob1 
  df.visdist2 <- df.visdist %>% gather("Probability", "ns.prob", key=Method, value=Probability)
  df.visdist2$Method <- ifelse(df.visdist2$Method=="Probability", "Visibility", "Network Size")
  
  n.extreme <- length(network.size) - sum(ns.prob)
  per.extreme <- round(n.extreme/length(network.size), digits=3)*100
  
  if ("visdistoverlay" %in% type) {
    visdistoverlay <- ggplot(df.visdist2, aes(x=Visibility, y=Probability, fill=Method)) +
      geom_bar(stat="identity", position="identity", width=rep(c(0.9,0.4),each=Kmax)) +
      theme_bw() +
      labs(title="Visibility Distribution with Network Sizes Overlaid") +
      scale_fill_manual(values=c("red", "gray35")) +
      annotate(geom = 'text', label = paste("Observed network size > K:\n", per.extreme, "%", sep=""), x = Kmax, y = max(df.visdist2$Probability), color="red", size=4, hjust=1, vjust=1)
    plist <- c(plist, list(visdistoverlay))
  }
  
  # Visibility vs network size by person
  df.visall <- data.frame(Visibility=fit$visibilities, Network=network.size)
  qs <- apply(fit$vsample,2,stats::quantile,probs=c(0.25,0.75))
  df.visall$q1 <- qs[1,]
  df.visall$q3 <- qs[2,]
  
  if ("visall" %in% type) {
    visall <- ggplot(df.visall, aes(x=Network, y=Visibility)) +
      theme_bw() +
      geom_errorbar(aes(ymin=q1, ymax=q3), width=2, alpha=0.4) +
      geom_count(color="black", fill="red", shape=21) +
      labs(title="Estimated visibilities for each respondent", x="Network Size (Augmented)")
    plist <- c(plist, list(visall))
  }
  
  # Visibility dist for random participants
  if ("pvdist" %in% type) {
    vsamp <- data.frame(fit$vsample)
    vsamp2 <- vsamp %>% gather(sample(names(vsamp), 20), key="ID", value="Visibility")
    
    pvdist <- ggplot(vsamp2, aes(x=Visibility, group=ID)) + 
      geom_density(fill="black", alpha=0.08, bw=0.5) +
      theme_bw() +
      labs(title="Visibility distributions for 20 randomly selected participants", y="Density")
    plist <- c(plist, list(pvdist))
  }
  
  nc <- ifelse(length(plist) >= 2, 2, 1)
  nr <- ceiling(length(plist)/2)
  
  do.call(ggarrange, c(plist, ncol=nc, nrow=nr))
  
}
