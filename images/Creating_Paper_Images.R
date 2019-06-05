# Load Libraries
library(tidyverse)
library(x3ptools)
library(imager)
library(bulletxtrctr)
library(fixedpoints)

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

cimg <- as.cimg(bullets$x3p[[1]]$surface.matrix)
### Create cimage of scan
png(file = "images/Houston_BarrelF_Bullet1_cimg.png", width = 800, height = 550)
plot(as.cimg(bullets$x3p[[1]]$surface.matrix))
dev.off()

### post edge detection image
dx <- imgradient(strong.b1.l1, "x")
dy <- imgradient(strong.b1.l1, "y")
grad.mag <- sqrt(dx^2+dy^2)
strong.b1.l1 <- grad.mag > quantile(grad.mag, .99, na.rm = TRUE)

png(file = "images/Houston_BarrelF_Bullet1_Strong_Edge.png", width = 800, height = 550)
plot(strong.b1.l1)
dev.off()

# With Canny Edge
## Hysteresis
t2 <- quantile(grad.mag, .99, na.rm = TRUE)
t1 <- quantile(grad.mag, .8, na.rm = TRUE)

weak.b1.l1 <- grad.mag %inr% c(t1,t2)

overlap <- grow(strong.b1.l1, 3) & weak.b1.l1
strong.b1.l1.new <- strong.b1.l1 | overlap

expand.strong <- function(ws){
  overlap <- grow(ws$strong, 3) & ws$weak
  ws$strong[overlap] <- TRUE
  ws$weak[overlap <- FALSE]
  ws 
}

hystFP <- fp(expand.strong)
a <- Sys.time()
out <- list(strong = strong.b1.l1, weak = weak.b1.l1) %>% hystFP
b <- Sys.time(); b-a

png(file = "images/Houston_BarrelF_Bullet1_Canny_Edge.png", width = 800, height = 550)
plot(out$strong)
dev.off()




### Create image with hough lines
df.strong.b1.l1 <- hough_line(strong.b1.l1, data.frame = TRUE)

png("images/Houston_BarrelF_Bullet1_Hough_Bin10.png", width = 800, height = 550)
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

png("images/Houston_BarrelF_Bullet1_middle_twothirds.png", width = 800, height = 550)
plot(strong.b1.l1)
abline(v = lthird, col = "green", lwd = 3)
abline(v = uthird, col = "green", lwd = 3)
dev.off()


closelthird <- good_vertical_segs[which.min(abs(good_vertical_segs - lthird))]
closeuthird <- good_vertical_segs[which.min(abs(good_vertical_segs - uthird))]


bestfit <- segments %>%
  filter(xaverage %in% c(closelthird, closeuthird))

png("images/Houston_BarrelF_Bullet1_BestFit.png", width = 800, height = 550)
plot(strong.b1.l1)
with(bestfit, nfline(theta, rho, col = "red", lwd = 3))
abline(v = lthird, col = "green", lwd = 3)
abline(v = uthird, col = "green", lwd = 3)
dev.off()

# Houstong Barrel F Bullet1 BestFit 
#########################################################


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
  ccdata = x3p %>% purrr::map(.f = x3p_to_df)
)

####### Land 3 images for hamby 252

land3 <- bullets$ccdata[[3]]
land3 <- df_to_x3p(land3)

land3.cimg <- as.cimg(land3$surface.matrix)

dx <- imgradient(land3.cimg, "x")
dy <- imgradient(land3.cimg, "y")
grad.mag <- sqrt(dx^2 + dy^2)

strong.hamby.l3 <- grad.mag > quantile(grad.mag, .99, na.rm = TRUE )

png(file = "images/Hamby252_Bullet1_Land3_Strong_edge.png", width = 800, height = 550)
plot(strong.hamby.l3)
dev.off()

weak.hamby.l3 <- grad.mag %inr% c(t1,t2)

# Canny Edge With Hamby
a <- Sys.time()
out <- list(strong = strong.hamby.l3, weak = weak.hamby.l3) %>% hystFP
b <- Sys.time(); b-a

png(file = "images/Hamby252_Bullet1_Land3_Canny_Edge.png", width = 800, height = 550)
plot(out$strong)
dev.off()



df.strong.l3 <- hough_line(strong.hamby.l3, data.frame = TRUE)

df.strong.l3 <- df.strong.l3 %>%
  mutate(theta = ifelse(theta <= pi, theta, theta - 2*pi)) %>%
  filter(score > quantile(score, .999),
         theta > (-pi/4),
         theta < (pi/4))

segments <- rho_to_ab(df = df.strong.l3)

segments <- segments %>%
  mutate(pixset.intercept = ifelse(theta==0, xintercept, (height(strong.hamby.l3) - yintercept)/slope),
         xaverage = ifelse(theta==0, xintercept, ((0-yintercept)/slope + (height(strong.hamby.l3) - yintercept)/slope)/2))

good_vertical_segs <- segments %>%
  extract2("xaverage")

lthird <- width(strong.hamby.l3)/6
uthird <- 5*width(strong.hamby.l3)/6

# Find hough line index where
closelthird <- good_vertical_segs[which.min(abs(good_vertical_segs - lthird))]
closeuthird <- good_vertical_segs[which.min(abs(good_vertical_segs - uthird))]

png("images/Hamby_252_Bullet1_Land3_BestFit.png", width = 800, height = 550)
plot(strong.hamby.l3)
with(subset(segments, xaverage %in% c(closelthird, closeuthird)), nfline(theta, rho, col = "red", lwd = 3))
abline(v = lthird, col = "green", lwd = 3)
abline(v = uthird, col = "green", lwd = 3)
dev.off()
