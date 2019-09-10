rho_to_ab <- function(rho = NULL, theta = NULL, df = NULL) {
  if (is.null(df)) {
    df <- data.frame(rho = rho, theta = theta)
  }
  stopifnot(c("rho", "theta") %in% names(df))
  
  df <- df %>%
    mutate(
      yintercept = ifelse(theta == 0, NA, rho/sin(theta)),
      slope = -cos(theta)/sin(theta),
      xintercept = rho/cos(theta)) # cos(theta) == 0 for theta = ± pi/2 but realistically we don't get there
  df
}
