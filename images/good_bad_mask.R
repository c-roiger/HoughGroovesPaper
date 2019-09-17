library(grooveFinder)
library(bulletxtrctr)
library(tidyverse)
library(x3ptools)

# Pre-prep the bullets 
b1 <- read_bullet(urllist = hamby252demo[[1]])
b2 <- read_bullet(urllist = hamby252demo[[2]])

b1$bullet <- 1
b2$bullet <- 2
b1$land <- 1:6
b2$land <- 1:6
bullets <- rbind(b1, b2)


bullets$x3p[[1]]$header.info$incrementY
bullets$x3p[[1]]$header.info$incrementX
summary(as.vector(bullets$x3p[[1]]$surface.matrix))

# get hough grooves

bullets <- bullets %>% mutate(
  x3p = x3p %>% purrr::map(.f = rotate_x3p, 90),
  x3p = x3p %>% purrr::map(.f = x3p_m_to_mum),
  ccdata = x3p %>% purrr::map(.f = x3p_to_df)
)

# Good example
grooves <- get_grooves_hough(bullets$ccdata[[1]], 0.999)
a <- get_mask_hough(bullets$x3p[[1]], grooves)

image_x3p(a)

# Bad example

grooves <- get_grooves_hough(bullets$ccdata[[6]], 0.9)
a <- get_mask_hough(bullets$x3p[[6]], grooves)
image_x3p(a)
