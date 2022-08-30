# Simplified version of Reingold-Tilford Plot that places recruitment chains horizontally with seeds all at the same level
# Has the option to color points based on a provided column in the data

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
