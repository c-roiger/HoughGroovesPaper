# Load Libraries
library(tidyverse)
library(x3ptools)
library(imager)

source("/Users/charlotteroiger/Documents/GitHub/bulletQuality/charlotte_code/rho_to_ab.R")

# Load in data 

files <- list.files( "/Users/charlotteroiger/Documents/CSAFE_bullet_quality/data/Houston_barf/", pattern = "*.x3p", full.names = T, recursive = T)

houston.bar1 <- data.frame(file.names = files) %>%
  mutate(x3p = purrr::map(as.character(file.names), .f = function(file.names){
    x3ptools::read_x3p(file.names)
  }))

bullets <- houston.bar1 %>% mutate(
  ccdata = x3p %>% purrr::map(.f = x3p_to_df),
  x3p = x3p %>% purrr::map(.f = y_flip_x3p)
)




strong.b1.l1 <- as.cimg(bullets$x3p[[1]]$surface.matrix)
strong.b1.l2 <- as.cimg(bullets$x3p[[2]]$surface.matrix)
strong.b1.l3 <- as.cimg(bullets$x3p[[3]]$surface.matrix)
strong.b1.l4 <- as.cimg(bullets$x3p[[4]]$surface.matrix)
strong.b1.l5 <- as.cimg(bullets$x3p[[5]]$surface.matrix)
strong.b1.l6 <- as.cimg(bullets$x3p[[6]]$surface.matrix)
strong.b2.l1 <- as.cimg(bullets$x3p[[7]]$surface.matrix)
strong.b2.l2 <- as.cimg(bullets$x3p[[8]]$surface.matrix)
strong.b2.l3 <- as.cimg(bullets$x3p[[9]]$surface.matrix)
strong.b2.l4 <- as.cimg(bullets$x3p[[10]]$surface.matrix)
strong.b2.l5 <- as.cimg(bullets$x3p[[11]]$surface.matrix)
strong.b2.l6 <- as.cimg(bullets$x3p[[12]]$surface.matrix)

### Create x3p image 
x3p_image(bullets$x3p[[1]], multiply = 2, zoom = 0.5, file = "/Users/charlotteroiger/Documents/GitHub/HoughGroovesPaper/images/Houston_BarrelF_Bullet1_Land1_Scan.png")


### Create cimage of scan
png(file = "images/Houston_BarrelF_Bullet1_cimg.png", width = 600, height = 350)
plot(as.cimg(bullets$x3p[[1]]$surface.matrix))
dev.off()

### post edge detection image
dx <- imgradient(strong.b1.l1, "x")
dy <- imgradient(strong.b1.l1, "y")
grad.mag <- sqrt(dx^2+dy^2)
strong.b1.l1 <- grad.mag > quantile(grad.mag, .99, na.rm = TRUE)

png(file = "images/Houston_BarrelF_Bullet1_Strong_Edge.png", width = 600, height = 350)
plot(strong.b1.l1)
dev.off()

### Create image with hough lines
df.strong.b1.l1 <- hough_line(strong.b1.l1, data.frame = TRUE)

png("images/Houston_BarrelF_Bullet1_Hough_Bin100.png")
plot(strong.b1.l1)
with(subset(df.strong.b1.l1,score > quantile(score, .999) & (theta < pi/4)) ,nfline(theta,rho,col="red"))
dev.off()

### Find Hough lines closes to the middle two thirds

segments <- rho_to_ab(df = df.strong.b1.l1)
segments <- segments %>%
  mutate(
         xaverage = ifelse(theta == 0, xintercept, ((0-yintercept)/slope + (height(strong.b1.l1) - yintercept)/slope)/2),
         pixset.intercept = ifelse(theta == 0, xintercept, (height(strong.b1.l1) - yintercept)/slope))

good_vertical_segs <- segments %>%
  filter(score > quantile(score, 0.9975) & theta < pi/4) %>%
  extract2("xaverage")


lthird <- width(strong.b1.l1)/6
uthird <- 5*width(strong.b1.l1)/6

closelthird <- good_vertical_segs[which.min(abs(good_vertical_segs - lthird))]
closeuthird <- good_vertical_segs[which.min(abs(good_vertical_segs - uthird))]


bestfit <- segments %>%
  filter(xaverage %in% c(closelthird, closeuthird))

png("images/Houston_BarrelF_Bullet1_BestFit.png")
plot(strong.b1.l1)
with(bestfit, nfline(theta, rho, col = "red"))
abline(v = lthird, col = "green")
abline(v = uthird, col = "green")
dev.off()


