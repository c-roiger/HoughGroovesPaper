geom_nfline <- function(mapping = NULL, data = NULL,
                        ...,
                        theta,
                        rho,
                        na.rm = FALSE,
                        show.legend = NA) {
  
  # If nothing set, default to y = 0, x = 1
  if (is.null(mapping) && missing(theta) && missing(rho)) {
    rho <- 1
    theta <- 0
  }
  
  # Act like an annotation
  if (!missing(rho) || !missing(theta)) {
    
    # Warn if supplied mapping and/or data is going to be overwritten
    if (!is.null(mapping)) {
      warn_overwritten_args("geom_abline()", "mapping", c("rho", "theta"))
    }
    if (!is.null(data)) {
      warn_overwritten_args("geom_abline()", "data", c("rho", "theta"))
    }
    
    if (missing(rho)) rho <- 1
    if (missing(theta)) theta <- 0
    n_thetas <- max(length(theta), length(theta))
    
    data <- new_data_frame(list(
      theta = theta,
      rho = rho
    ), n = n_thetas)
    mapping <- aes(theta = theta, rho = rho)
    show.legend <- FALSE
  }
  
  layer(
    data = data,
    mapping = mapping,
    stat = StatIdentity,
    geom = GeomNfline,
    position = PositionIdentity,
    show.legend = show.legend,
    inherit.aes = FALSE,
    params = list(
      na.rm = na.rm,
      ...
    )
  )
}


GeomNfline <- ggproto("GeomNfline", Geom,
                      draw_panel = function(data, panel_params, coord) {
                        ranges <- coord$backtransform_range(panel_params)
                        
                        data$slope <- -cos(theta)/sin(theta)
                        data$intercept <- ifelse(theta == 0, NA, rho/sin(theta))
                        
                        data$x    <- ranges$x[1]
                        data$xend <- ranges$x[2]
                        data$y    <- ranges$x[1] * data$slope + data$intercept
                        data$yend <- ranges$x[2] * data$slope + data$intercept
                        
                        GeomSegment$draw_panel(unique(data), panel_params, coord)
                      },
                      
                      default_aes = aes(colour = "black", size = 0.5, linetype = 1, alpha = NA),
                      required_aes = c("theta", "rho"),
                      
                      draw_key = draw_key_abline
)

warn_overwritten_args <- function(fun_name, overwritten_arg, provided_args, plural_join = " and/or ") {
  overwritten_arg_text <- paste0("`", overwritten_arg, "`")
  
  n_provided_args <- length(provided_args)
  if (n_provided_args == 1) {
    provided_arg_text <- paste0("`", provided_args, "`")
    verb <- "was"
  } else if (n_provided_args == 2) {
    provided_arg_text <- paste0("`", provided_args, "`", collapse = plural_join)
    verb <- "were"
  } else {
    provided_arg_text <- paste0(
      paste0("`", provided_args[-n_provided_args], "`", collapse = ", "),
      ",", plural_join,
      "`", provided_args[n_provided_args], "`"
    )
    verb <- "were"
  }
  
  warning(
    sprintf(
      "%s: Ignoring %s because %s %s provided.",
      fun_name,
      overwritten_arg_text,
      provided_arg_text,
      verb
    ),
    call. = FALSE
  )
}


