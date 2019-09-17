library(grooveFinder)
library(bulletxtrctr)
library(x3ptools)
library(tidyverse)

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


land.x3p <- bullets$x3p[[1]]
pix_to_micron <- function(x) {
  (x-1) * x3p_get_scale(land.x3p)
}

cimg <- as.cimg(land.x3p$surface.matrix)
# }

# Create image gradient
dx <- imgradient(cimg, "x")
dy <- imgradient(cimg, "y")

grad.mag <- sqrt(dx^2 + dy^2)

strong <- grad.mag > quantile(grad.mag, .99, na.rm = TRUE)
# create the hough transform
hough.df <- hough_line(strong, data.frame = TRUE, shift = FALSE) # we want to get values with respect to (0,0) not (1,1)

# Subset based on score and angle rotation
hough.df <- hough.df %>%
  dplyr::mutate(theta = ifelse(theta <= pi, theta, theta - 2 * pi)) %>%
  dplyr::filter(
    score > quantile(score, 0.999),
    theta > (-pi / 16), # identify only vertical(ish) lines
    theta < (pi / 16),
    (rho < abs(width(strong)) * 1 / 6 | rho > width(strong) * 5 / 6) # at either end of the LEA
  )


# get x and y intercepts  (in pixel dimensions)
segments <- grooveFinder:::rho_to_ab(df = hough.df)

# browser()
segments <- segments %>%
  dplyr::mutate(
    pixset.intercept = ifelse(theta == 0, xintercept, (height(strong) - yintercept) / slope),
    xaverage = ifelse(theta == 0, xintercept,
                      ((0 - yintercept) / slope + (height(strong) - yintercept) / slope) / 2
    ),
    xbottom = ifelse(theta == 0, xintercept, (height(strong) - yintercept) / slope),
    xtop = xintercept
  )

plot(strong)
abline(yint.left, slope.left, col = "green")
abline(yint.right, slope.right, col = "green")
abline(v = lthird, col = "red", lwd = 3)
abline(v = uthird, col = "red", lwd = 3)

# Find the middle 2/3rds
lthird <- width(strong) / 6
uthird <- 5 * width(strong) / 6

# Find best bottom and top index of groove for both sides
top.left <- segments$xintercept[which.min(abs(segments$xintercept - lthird))]
bottom.left <- segments$pixset.intercept[which.min(abs(segments$pixset.intercept - lthird))]


top.right <- segments$xintercept[which.min(abs(segments$xintercept - uthird))]
bottom.right <- segments$pixset.intercept[which.min(abs(segments$pixset.intercept - uthird))]

# Calculate equation of line for each side
slope.left <- -height(strong)/(top.left - bottom.left)
yint.left <- -(slope.left*top.left)

slope.right <- -height(strong)/(top.right - bottom.right)
yint.right <- -(slope.right*top.right)

strong.raster <- imager:::as.data.frame.cimg(strong)
ggplot() +
  geom_raster(data = strong.raster, aes(x = x, y = -y, fill = value)) +
  scale_fill_manual(values = c("black", "white", "grey"))+
  geom_vline(xintercept = lthird, col = "green") +
  geom_vline(xintercept = uthird, col = "green") +
  geom_abline(aes(slope = -slope, intercept = -yintercept), data = segments, col = "blue") + 
  # geom_abline(slope = -slope.left, intercept = -yint.left, col = "red") + 
  # geom_abline(slope = -slope.right, intercept = -yint.right, col = "red")+
  # geom_point(aes(x = top.left, y = 0), col = "red")+
  # geom_point(aes(x = top.right, y = 0), col = "red")+
  # geom_point(aes(x = bottom.left, y = -501), col = "red")+
  # geom_point(aes(x = bottom.right, y = -501), col = "red")+
  coord_fixed() +
  guides(fill = FALSE)
ggsave(filename = "images/Hamby252_Bullet1_Land1_hough_lines_process_step_one.png", width = 4, height = 3, unit = "in")


