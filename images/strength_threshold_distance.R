library(grooveFinder)
library(bulletxtrctr)
library(x3ptools)
library(tidyverse)

data.x3p <- read_x3p("../bulletQuality/manual_ident/phoenix/Gun 1-A9/B1/L4.x3p")
sggh <- safely(get_grooves_hough)
length_between <- function(data.x3p, qu, crosscut){
  grooves <-sggh(x3p_to_df(data.x3p), qu)
  if(is.null(grooves$error)){
    distance <- grooves$result$right.groove.fit(crosscut) - grooves$result$left.groove.fit(crosscut)
  }
  else{
    distance <- NA
  }
  return(distance)
  
}


qu <- c(0.999, 0.995, 0.99, 0.95, 0.9, 0.85, 0.8, 0.75, 0.7)
dist <- c()
yep <- function(data.x3p, qu){
  dist <- c()
  ccdata <- x3p_crosscut_optimize(data.x3p)
  for(i in 1:length(qu)){
    dist[i] <- length_between(data.x3p, qu[i], crosscut)
  }
  return(dist)
}
yep(data.x3p, qu)
