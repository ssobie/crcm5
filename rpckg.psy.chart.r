library(ggplot2)
library(weights)

enthalpy_intersect <- function(enthalpy, alt = 0) {
  #' Intersection between enthalpy line and 100 \% RH line
  #'
  #' Calculates the intersection of RH 100 \% line and the enthalpy line.
  #' @param enthalpy Vector of enthalpies [kJ/kg].
  #' @param alt Vector of altitudes [m]. Defaults to 0 m (sea level).
  #' @return Returns the dry-bulb temperature at which the enthalpy and
  #'  100 \% RH line intersects [degC].
  #' @keywords internal
  #' @export
  #' @examples
  #' enthalpy_intersect(enthalpy = 50)
  #' @author Christoffer Rasmussen

  temp <- seq(-75, 75, 0.1)

  df.test <- data.frame(
    temp = temp,
    enthalpy.diff =
      abs(-0.557341 * (temp - 0.994036 * enthalpy) / (temp + 1385.6) -
      sat_hum_ratio(temp, alt))
  )

  return(df.test$temp[which.min(df.test$enthalpy.diff)])
}

slope_rel_hum <- function(temp.1, temp.2, rel.hum) {
  #' Slope of relative humidity line
  #'
  #' Calculate slope of relative humidity line.
  #' @param temp.1 Temperature of first point [degC].
  #' @param temp.2 Temperature of second point [degC].
  #' @param rel.hum Relative humidity line to find slope of [\%].
  #' @return Returns the slope angle [deg].
  #' @keywords internal
  #' @export
  #' @examples
  #' slope_rel_hum(15, 20, 80)
  #' @author Christoffer Rasmussen

  x1 <- temp.1
  x2 <- temp.2
  y1 <- hum_ratio(rel.hum, temp.1)
  y2 <- hum_ratio(rel.hum, temp.2)

  angle <- atan((y2 - y1) * 1000 / (x2 - x1))

  return(angle * 180 / pi)
}

hum_ratio <- function(rel.hum, temp.air, alt = 0) {
  #' Humidity ratio (W)
  #'
  #' Calculates the humidity ratio from relative humidity and temperature.
  #' Eq. 25 & 12 - ASHRAE Fundamentals Handbook 2002, Psychrometrics.
  #' @param rel.hum Vector of relative humidities [\%].
  #' @param temp.air Vector of air temperatures [degC].
  #' @param alt Vector of altitudes [m]. Defaults to 0 m (sea level).
  #' @return Returns a vector of humidity ratios [kg/kg].
  #' @export
  #' @examples
  #' hum_ratio(rel.hum = 60, temp.air = 25, alt = 0)
  #' @author Christoffer Rasmussen

  # Humidity ratio (eq. 25 and 12, solved for W)
  rel.hum <- rel.hum / 100
  w.s <- sat_hum_ratio(temp.air, alt)
  p <- bar_press(alt)
  p.ws <- sat_w_press(temp.air)

  return(-rel.hum * w.s * (p - p.ws) / (rel.hum * p.ws - p))
}

hum_ratio_enthalpy <- function(enthalpy, temp.db, alt = 0) {
  #' Humidity ratio for enthalpy lines
  #'
  #' Calculates the humidity ratio based on the enthalpy.
  #' Eq. 32 - ASHRAE Fundamentals Handbook 2002, Psychrometrics.
  #' @param enthalpy Vector of enthalpies [kJ/kg].
  #' @param temp.db Vector of dry bulb temperatures [degC].
  #' @return Returns a vector of humidity ratios [kg/kg].
  #' @keywords internal
  #' @export
  #' @examples
  #' hum_ratio_enthalpy(enthalpy = 30, temp.db = 25)
  #' @author Christoffer Rasmussen

  # Humidity ratio (eq. 32, solved for W)
  w <- -0.557341 * (temp.db - 0.994036 * enthalpy) / (temp.db + 1385.6)

  cut.profile <- sat_hum_ratio(temp.db, alt)

  w <- ifelse(w > cut.profile, NA, ifelse(w < 0, NA, w))

  return(w)

}

enthalpy <- function(temp.db, hum.ratio) {
  #' Enthalpy (h)
  #'
  #' Calculates the enthalpy of a given moist air state.
  #' Eq. 32 - ASHRAE Fundamentals Handbook 2002, Psychrometrics.
  #' @param temp.db Vector of dry bulb temperatures [degC].
  #' @param hum.ratio Vector of humidity ratios [kg/kg].
  #' @return Returns a vector of enthalpies [kJ/kg].
  #' @export
  #' @examples
  #' enthalpy(25, 0.010)
  #'
  #' # Calculating cooling demand
  #' density <- 1.2 # density of air (kg/m3)
  #' (enthalpy(35, 0.020) - enthalpy(15, 0.010)) * density
  #' @author Christoffer Rasmussen

  # Enthalpy (eq. 32)
  return(1.006 * temp.db + hum.ratio * (2501 + 1.805 * temp.db))
}

my_round <- function(x, n, type) {
  #' Round to nearest number
  #'
  #' Round up or down to nearest user specified number.
  #' @param x Vector of numbers to round.
  #' @param n Nearest number to which the x-value should be rounded.
  #' @param type A string (either 'ceiling' or 'floor') to indicate if the
  #'  number should be rounded up or down.
  #' @return Returns a vector of rounded numbers.
  #' @keywords internal
  #' @export
  #' @examples
  #' my_round(24.9, 5, "floor")
  #' my_round(c(24.9, 20.1, 25), 5, "ceiling")
  #' @author Christoffer Rasmussen

  if (type == "ceiling") {
    x.rounded <- ceiling(x / n) * n
  } else if (type == "floor") {
    x.rounded <- floor(x / n) * n
  }

  return(x.rounded)
}

wetbulb_intersect <- function(temp.wb, alt = 0) {
  #' Wet-bulb intersection with dry-bulb temperature axis
  #'
  #' Calculates the drybulb temperature at which the wetbulb temperature line
  #' intersects in the psychrometric chart.
  #' Eq. 32 - ASHRAE Fundamentals Handbook 2002, Psychrometrics.
  #' @param temp.wb Vector of wet-bulb temperatures [degC].
  #' @param alt Vector of altitudes [m]. Defaults to 0 m (sea level).
  #' @return Returns a vector of dry-bulb temperatures at which the wet-bulb
  #'   temperature line intersects [degC].
  #' @keywords internal
  #' @export
  #' @examples
  #' wetbulb_intersect(15)
  #' wetbulb_intersect(15, alt = 1000)
  #' @author Christoffer Rasmussen

  # Calculate satuated humidity for wet-bulb temperature (kg/kg)
  w.s <- sat_hum_ratio(temp.wb, alt)

  # Return dry-bulb temperature for given wet-bulb temperature (Eq. 35)
  return(-2.36581 * ((temp.wb - 1050.84) * w.s - 0.422689 * temp.wb))
}

sat_w_press <- function(temp.air) {
  #' Satuated water vapour pressure (p_ws)
  #'
  #' Calculates the satuated water vapour pressure over water and ice.
  #' Eq.5 & 6 - ASHRAE Fundamentals Handbook 2002, Psychrometrics.
  #' @param temp.air Vector of air temperatures [degC].
  #' @return Returns a vector of satuated water pressure over water and
  #'  ice [kPa].
  #' @export
  #' @examples
  #' sat_w_press(temp.air = 25)
  #' @author Christoffer Rasmussen

  temp.air <- temp.air + 273.15

  # Coefficients for calculating the satuated pressure over water and ice
  c <- c(-5674.5359, 6.3925247, -0.009677843, 6.22157E-07, 2.07478E-09,
         -9.48402E-13, 4.1635019, -5800.2206, 1.3914993, -0.048640239,
         4.17648E-05, -1.44521E-08, 6.5459673)

  # Calculate satuated pressure over ice (Eq. 5) and water (Eq. 6)
  sat.w.press <- ifelse(temp.air < 273.15,
                 exp(c[1] / temp.air +
                     c[2] +
                     c[3] * temp.air +
                     c[4] * temp.air ** 2 +
                     c[5] * temp.air ** 3 +
                     c[6] * temp.air ** 4 +
                     c[7] * log(temp.air)),
                 exp(c[8] / temp.air +
                     c[9] +
                     c[10] * temp.air +
                     c[11] * temp.air ** 2 +
                     c[12] * temp.air ** 3 +
                     c[13] * log(temp.air)))

  # Return pressure in kPa
  return(sat.w.press / 1000)
}

sat_hum_ratio <- function(temp.air, alt = 0) {
  #' Humidity ratio at satuation (W_s)
  #'
  #' Eq. 23 - ASHRAE Fundamentals Handbook 2002, Psychrometrics.
  #' @param temp.air Vector of air temperatures [degC].
  #' @param alt Vector of altitudes [m]. Defaults to 0 m (sea level).
  #' @return Returns a vector of water vapour content at satuatuation [kg/kg].
  #' @export
  #' @examples
  #' sat_hum_ratio(temp.air = 25, alt = 0)
  #' @author Christoffer Rasmussen

  # Humidity ratio i kg/kg (eq. 23)
  return(0.62198 * sat_w_press(temp.air) / (bar_press(alt) -
         sat_w_press(temp.air)))
}

bar_press <- function(alt = 0) {
  #' Barometric/total pressure (p)
  #'
  #' Calculates barometric pressure based on altitude. This is also called total
  #'  pressure.
  #' Eq. 3 - ASHRAE Fundamentals Handbook 2002, Psychrometrics.
  #' @param alt Vector of altitudes [m]. Defaults to 0 m (sea level).
  #' @return Returns a vector of barometric pressure [kPa].
  #' @export
  #' @examples
  #' bar_press()
  #'
  #' bar_press(alt = 1000)
  #' @author Christoffer Rasmussen

  # Barometric pressure (eq. 3)
  return(101.325 * (1 - 2.25577e-05 * alt) ** 5.2559)
}

par_w_press <- function(hum.ratio, alt = 0) {
  #' Partial water vapour pressure (p_w)
  #'
  #' Calculates partial water vapour pressure.
  #' Eq. 36 - ASHRAE Fundamentals Handbook 2002, Psychrometrics.
  #' @param hum.ratio Vector of humidity ratios [kg/kg].
  #' @param alt Vector of altitudes [m]. Defaults to 0 m (sea level).
  #' @return Returns a vector of partial water vapour pressure [kPa].
  #' @export
  #' @examples
  #' par_w_press(hum.ratio = 0.010)
  #'
  #' hum.ratio <- rep(c(0.005, 0.010), 2)
  #' alt <- rep(c(0, 500), each = 2)
  #' par_w_press(hum.ratio, alt)
  #' @author Christoffer Rasmussen

  # Partial water vapor pressure (eq. 36)
  return(bar_press(alt) * hum.ratio / (0.62198 + hum.ratio))
}

dewpoint <- function(hum.ratio, alt = 0) {
  #' Dew-point temperature (t_d)
  #'
  #' Calculates the dew-point temperature based on humidity ratios.
  #' Eq. 37 & 38 - ASHRAE Fundamentals Handbook 2002, Psychrometrics.
  #' @param hum.ratio Vector of humidity ratios [kg/kg].
  #' @param alt Vector of altitudes [m]. Defaults to 0 m (sea level).
  #' @return Returns a vector of dew-point temperatures [degC].
  #' @export
  #' @examples
  #' dewpoint(hum.ratio = 0.010)
  #'
  #' hum.ratio = seq(0, 0.050, by = 0.001)
  #' plot(dewpoint(hum.ratio = hum.ratio), hum.ratio,
  #'      xlab = "Dew-point [degC]",
  #'      ylab = "Humidity ratio [kg/kg]")
  #' @author Christoffer Rasmussen

  par.w.press <- par_w_press(hum.ratio, alt)

  # Dewpoint temperature (eq. 37)
  dewpoint <- 6.54 + 14.526 * log(par.w.press) +
              0.7389 * log(par.w.press) ** 2 +
              0.09486 * log(par.w.press) ** 3 +
              0.4569 * par.w.press ** 0.1984

  dewpoint <- ifelse(dewpoint < 0,
                     6.09 + 12.608 * log(par.w.press) +
                     0.4959 * log(par.w.press) ** 2,
                     ifelse(dewpoint > 93, NA, dewpoint))

  return(dewpoint)
}

hum_ratio_rel_hum <- function(rel.hum, temp.air, hum.ratio.max, alt = 0) {
  #' Humidity ratio (W)
  #'
  #' Calculates the humidity ratio from relative humidity and temperature.
  #' Eq. 25 & 12 - ASHRAE Fundamentals Handbook 2002, Psychrometrics.
  #' @param rel.hum Vector of relative humidities [\%].
  #' @param temp.air Vector of air temperatures [degC].
  #' @param hum.ratio.max Max humidity ratio on psychrometric chart [kg/kg].
  #' @param alt Vector of altitudes [m]. Defaults to 0 m (sea level).
  #' @return Returns a vector of humidity ratios [kg/kg].
  #' @keywords internal
  #' @export
  #' @examples
  #' hum_ratio_rel_hum(rel.hum = 60, temp.air = 25, hum.ratio.max = 0.02)
  #' @author Christoffer Rasmussen

  # Humidity ratio (eq. 25 and 12, solved for W)
  rel.hum <- rel.hum / 100
  W_s <- sat_hum_ratio(temp.air, alt)
  p <- bar_press(alt)
  p_ws <- sat_w_press(temp.air)

  W <- (-rel.hum) * W_s * (p - p_ws) / (rel.hum * p_ws - p)

  W <- ifelse(W > hum.ratio.max, NA, W)

  return(W)

}

theme_parms <- function() {
  #' ggplot2 parameters
  #'
  #' Parameters for climateeng themes.
  #' @return Returns a list of ggplot2 theme parameters.
  #' @keywords internal
  #' @export
  #' @examples
  #' theme_parms()
  #' @author Christoffer Rasmussen

  # Check for package -------------------------------------------------------

  if (!requireNamespace("RColorBrewer", quietly = TRUE)) {
  stop("RColorBrewer is needed for this function to work. Please install it.",
       call. = FALSE)
  }

  # Set theme parameters ----------------------------------------------------

  palette <- RColorBrewer::brewer.pal("Greys", n = 9)

  l <- list(
    color.background = palette[1],
    color.grid.major = palette[4],
    color.axis.text = palette[6],
    color.axis.title = palette[7],
    color.title = palette[9],
    size.line = 0.25,
    text.size.axis = 12,
    text.size.axis.title = 12,
    text.size.plot.title = 12
  )

  return(l)
}


theme_climateeng_psy <- function(asp = NULL) {
  #' Custom ggplot2
  #'
  #' Custom theme for climateeng psychrometric chart.
  #' @param asp Aspect ratio of plot. Defaults to NULL.
  #' @return Returns a ggplot2 theme.
  #' @keywords internal
  #' @export
  #' @examples
  #' theme_climateeng_psy()
  #' @author Christoffer Rasmussen

  # Check for package -------------------------------------------------------

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is needed for this function to work. Please install it.",
         call. = FALSE)
  } else if (!requireNamespace("grid", quietly = TRUE)) {
    stop("grid is needed for this function to work. Please install it.",
         call. = FALSE)
  }

  # Define theme -----------------------------------------------------------

  # Generate the colors for the chart procedurally with RColorBrewer
  theme <- theme_parms()

  # Begin construction of chart
  ggplot2::theme_bw(base_size = 12) +

    # Set the entire chart region to a light gray color
    ggplot2::theme(panel.background = ggplot2::element_rect(
                     fill = theme$color.background,
                     color = theme$color.background)) +
    ggplot2::theme(plot.background = ggplot2::element_rect(
                     fill = theme$color.background,
                     color = theme$color.background)) +
    ggplot2::theme(panel.border = ggplot2::element_rect(
                     color = theme$color.background)) +

    # Format the grid
    ggplot2::theme(panel.grid = ggplot2::element_blank()) +
    ggplot2::theme(axis.ticks = ggplot2::element_line(NA)) +

    # Format the legend, but hide by default
    ggplot2::theme(legend.position = "none") +
    ggplot2::theme(legend.background = ggplot2::element_rect(
                     fill = theme$color.background)) +
    ggplot2::theme(legend.text = ggplot2::element_text(
                     size = 12,
                     color = theme$color.axis.title)) +

    # Set title and axis labels, and format these and tick marks
    ggplot2::theme(plot.title = ggplot2::element_text(
                     color = theme$color.title,
                     size = 12,
                     vjust = 1.25)) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(
                     size = 12,
                     color = theme$color.axis.text)) +
    ggplot2::theme(axis.text.y = ggplot2::element_text(
                     size = 12,
                     color = theme$color.axis.text)) +
    ggplot2::theme(axis.title.x = ggplot2::element_text(
                     size = 12,
                     color = theme$color.axis.title,
                     vjust = 0)) +
    ggplot2::theme(axis.title.y = ggplot2::element_text(
                     size = 12,
                     color = theme$color.axis.title,
                     vjust = 1.25)) +

    ggplot2::theme(aspect.ratio = asp) +

    # Plot margins
    ggplot2::theme(plot.margin = grid::unit(c(0.35, 0.2, 0.3, 0.35), "cm"))

}

psychro_chart <- function(temp.db = NULL, hum.ratio = NULL,
                          temp.min = -15, temp.max = 30,
                          humidity.max = 0.020, alt = 0,
                          mollier = FALSE, alpha = 0.25,
                          disable.warnings = TRUE) {
  #' Psychrometric chart
  #'
  #' Plot psychrometric chart with or with data.
  #' @param temp.db Vector of dry-bulb temperatures [degC]. Defaults to NULL.
  #' @param hum.ratio Vector of humidity ratios [kg/kg]. Defaults to NULL.
  #' @param temp.min Minimum value of temperature axis [degC]. Defaluts to
  #'  -15 degC.
  #' @param temp.max Maximum value of temperature axis [degC]. Defaluts to
  #'  30 degC.
  #' @param humidity.max Maximum value of humidity axis [kg/kg]. Defaults to
  #'  0.020 kg/kg.
  #' @param alt Vector of altitudes [m]. Defaults to 0 m (sea level).
  #' @param mollier Boolean operator. Plot Mollier chart (tx-chart) or tx-chart.
  #'  Defaults to FALSE.
  #' @param alpha Transparancy of the data points. Defaults to 0.25.
  #' @param disable.warnings Some ignorable warnings appear when plotting the
  #'  chart. To see these and eventually others disable.warnings should be
  #'  FALSE. Defaults to TRUE.
  #' @return Plots the psychrometric chart. If the error
  #' "TEXT_SHOW_BACKTRACE environmental variable.
  #' Error in grid.Call(L_textBounds, as.graphicsAnnot(x$label), x$x, x$y,  :
  #' polygon edge not found"
  #' shows up. Try run psychrometric_chart() again.
  #' @export
  #' @examples
  #' psychrometric_chart()
  #'
  #' temp.db <- c(20, 23, 18, 20, 10, 27)
  #' hum.ratio <- c(0.005, 0.009, 0.004, 0.011, 0.004, 0.01)
  #' psychrometric_chart(temp.db, hum.ratio, alpha = 1)
  #'
  #' psychrometric_chart(temp.db, hum.ratio, temp.min=5, mollier = TRUE,
  #'  alpha = 1)
  #' @author Christoffer Rasmussen

  # Check for package -------------------------------------------------------

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is needed for this function to work. Please install it.",
         call. = FALSE)
  } else if (!requireNamespace("weights", quietly = TRUE)) {
    stop("weights is needed for this function to work. Please install it.",
         call. = FALSE)
  }

  # Plot parameters ---------------------------------------------------------

  # Constants
  FONT.SCALE <- 2.834646  # Translate annotation font size to theme font size.
  LINE.MULTIPLIER <- 3  # Line size of boarders compares to rest of lines.
  N <- 2000  # Number og line segments for stat_function lines.

  # Round temperature and humidity limits
  temp.min <- floor(temp.min/5)*5
  temp.max <- ceiling(temp.max/5)*5
  humidity.max <- ceiling(humidity.max/0.005)*0.005

  # Temperature which RH annotation should be centered between
  temp.1 <- 20
  temp.2 <- 25

  # Temperature which RH annotation should be centered between
  temp.3 <- 10
  temp.4 <- 15

  # Wetbulb line to annotate
  wb <- 5

  # Aspect ratio
  if (mollier == F) {
    asp <- (humidity.max / 0.005) / ((temp.max - temp.min) / 5)
    deg <- 0
    k <- 1
  } else {
    asp <- ((temp.max - temp.min) / 5) / (humidity.max / 0.005)
    deg <- 90
    k <- -1
  }
  # Axis to filp. Depends on chart type.
  axis.to.flip <- "y"

  # Plot parameters
  theme <- theme_parms()


  # Create base plot --------------------------------------------------------


  p <- ggplot2::ggplot(data.frame(x = c(temp.min, temp.max)), ggplot2::aes(x)) +

    # Add theme
    theme_climateeng_psy(asp) +

    # Axis breaks
    ggplot2::scale_x_continuous(
      breaks = seq(temp.min, temp.max, 5)) +
    ggplot2::scale_y_continuous(
      breaks = seq(0, humidity.max, 0.005),
      labels = weights::rd(seq(0.0, humidity.max, 0.005), digits = 3, add = F),
      position = "right") +
    # Axis titles
    ggplot2::ylab(expression("Humidity ratio ("*kg[m]*"/"*kg[da]*")")) +
    ggplot2::xlab(expression("Dry-bulb temperature ("*degree*C*")"))


  # Trim chart --------------------------------------------------------------


  x_add <- (temp.max - temp.min) * 0.003
  y_add <- (humidity.max) * 0.003

  p <- p + ggplot2::coord_cartesian(
             xlim = c(temp.min - x_add, temp.max + x_add),
             ylim = c(-y_add, humidity.max + y_add),
             expand = F)


  # Converts tx-chart to xt-chart -------------------------------------------


  if (mollier == T) {
    p <- p + ggplot2::coord_flip(
               ylim = c(-y_add, humidity.max + y_add),
               xlim = c(temp.min - x_add, temp.max + x_add),
               expand = F)
    axis.to.flip <- "x"
  }


  # Relative humidity lines -------------------------------------------------


  for (i in seq(0, 100, 10)) {
    if (i %in% c(0, 100)) {
      # Border curve
      p <- p + ggplot2::stat_function(
                 fun = hum_ratio_rel_hum,
                 args = list(alt = alt,
                             rel.hum = i,
                             hum.ratio.max = humidity.max),
                 n = N,
                 size = theme$size.line * LINE.MULTIPLIER,
                 geom = "line",
                 col = theme$color.grid.major)
    } else {
      # Other RH curves
      p <- p + ggplot2::stat_function(
                 fun = hum_ratio_rel_hum,
                 args = list(alt = alt,
                             rel.hum = i,
                             hum.ratio.max = humidity.max),
                 n = N,
                 size = theme$size.line,
                 geom = "line",
                 col = theme$color.grid.major)
    }
  }


  # Absolute humidity lines -------------------------------------------------


  for (i in seq(0.005, humidity.max, 0.005)) {
    if (i == humidity.max) {
      # Border curve
      p <- p + ggplot2::geom_segment(
                 x = temp.max,
                 xend = dewpoint(hum.ratio = i, alt = alt),
                 y = i,
                 yend = i,
                 col = theme$color.grid.major,
                 size = theme$size.line * LINE.MULTIPLIER)
    } else {
      # Other RH curves
      p <- p + ggplot2::geom_segment(
                 x = temp.max,
                 xend = dewpoint(hum.ratio = i, alt = alt),
                 y = i,
                 yend = i,
                 size = theme$size.line,
                 col = theme$color.grid.major)
    }
  }


  # Dry bulb temperature lines ----------------------------------------------


  for (i in seq(temp.min, temp.max, 5)) {
    if (i %in% c(temp.min, temp.max)) {
      # Border curve
      p <- p + ggplot2::geom_segment(
                x = i,
                xend = i,
                y = 0,
                yend = min(sat_hum_ratio(i, alt), humidity.max),
                col = theme$color.grid.major,
                size = theme$size.line * LINE.MULTIPLIER)
    } else {
      # Other RH curves
      p <- p + ggplot2::geom_segment(
                x = i,
                xend = i,
                y = 0,
                yend = min(sat_hum_ratio(i, alt), humidity.max),
                size = theme$size.line,
                col = theme$color.grid.major)
    }
  }


  # Wet bulb temperature lines ----------------------------------------------


  temp_wb_min = temp.min
  while (wetbulb_intersect(temp_wb_min, alt) > temp.min) {
    temp_wb_min <- temp_wb_min - 5
  }

##  for (i in seq(temp_wb_min, temp.max - 5, 5)) {
##    p <- p + ggplot2::geom_segment(x = i,
##                          xend = wetbulb_intersect(i, alt),
##                          y = sat_hum_ratio(i, alt),
##                          yend = 0,
##                          size = theme$size.line,
##                          col = theme$color.grid.major,
##                          lty = "dotted")
##  }


  # Enthalpy lines ----------------------------------------------------------


  start <- my_round(enthalpy(temp.min, 0), 10, "ceiling")
  end <- my_round(enthalpy(temp.max, humidity.max), 10, "floor")

  # Add lines
##  for (i in seq(start, end, 5)) {
##    p <- p + ggplot2::stat_function(
##               fun = hum_ratio_enthalpy,
##               args = list(enthalpy = i,
##                           alt = alt),
##               n = N,
##               size = theme$size.line,
##               geom = "line",
##               col = theme$color.grid.major
##               )
##  }


  # Data points -------------------------------------------------------------


  # Add data to plot if present
  if (!is.null(temp.db) | !is.null(hum.ratio)) {

    # Create data frame
    df.points <- data.frame(
      temp = temp.db,
      hum = hum.ratio)

    # Add data points
    p <- p +
      ggplot2::scale_color_gradient(low = "#3F5151", high = "#9B110E") +
      ggplot2::geom_point(
        data = df.points,
        ggplot2::aes(x = temp, y = hum),
        size = 1.25,
        alpha = alpha)
  }


# Annotation of relative humidity lines -----------------------------------

  temp <- mean(c(temp.1, temp.2))

  if (mollier == F) {
    seq <- seq(10, 90, 20)
  } else if (mollier == T) {
    seq <- seq(20, 100, 20)
  }


  for (i in seq) {

    # Define label
    if (i == 10) {
      label <- "Relative humidity: 10 %\n"
    } else if (i == 100) {
      label <- "Relative humidity: 100 %\n"
    } else {
      label <- paste0(i, " %\n")
    }

    # Add label
    p <- p +
      ggplot2::annotate("text",
        x = temp,
        y = hum_ratio(i, temp),
        size = theme$text.size.axis / FONT.SCALE,
        label =  label,
        angle = deg + k * slope_rel_hum(temp.1, temp.2, i),
        col = theme$color.axis.text)
  }

  # Annotation of wetbulb lines -----------------------------------

  temp <- mean(c(temp.3, temp.4))

  # if (mollier == F) {
  #   seq <- seq(10, 90, 20)
  # } else if (mollier == T) {
  #   seq <- seq(20, 100, 20)
  # }


  # for (i in seq) {
  #
  #   # Define label
  #   if (i == 10) {
  #     label <- "Relative humidity: 10 %\n"
  #   } else if (i == 100) {
  #     label <- "Relative humidity: 100 %\n"
  #   } else {
  #     label <- paste0(i, " %\n")
  #   }

    # Add label
    label <- as.character(expression("Wet-bulb temperature: 5 "*degree*C*""))
    p <- p +
      ggplot2::annotate("text",
                        x = temp,
                        y = sat_hum_ratio(wb)/wetbulb_intersect(wb) * (temp-wb)-0.00015,
                        parse = TRUE,
                        size = theme$text.size.axis / FONT.SCALE,
                        label =  label,
                        angle = deg + k * atan(-(sat_hum_ratio(wb)*1000) / (wetbulb_intersect(wb)-wb)) * 180/pi,
                        col = theme$color.axis.text)
  # }




# Annotation of enthalpy lines --------------------------------------------


  start <- my_round(enthalpy(temp.min, 0), 10, "ceiling")
  end <- my_round(enthalpy(temp.max, humidity.max), 10, "floor")

  for (i in seq(start, end, 10)) {

    # Define label
    if (enthalpy_intersect(i, alt) < temp.min + 2 |
        enthalpy_intersect(i, alt) > dewpoint(humidity.max, alt) - 0.5) {
      label <- ""
    } else if (mollier == F) {
      if (i == 40) {
        label <- paste0(i, "\n\n")##paste0("Enthalpy (kJ/kg)") ##: ", i, " kJ/kg\n\n")
      } else {
        label <- paste0(i, "\n\n")
      }
    } else if (mollier == T) {
      if (i == 40) {
        label <- paste0(i, "\n\n")##"Enthalpy (kJ/kg)" ##paste0("\n\nEnthalpy: ", i, " kJ/kg")
      } else {
        label <- paste0("\n\n", i)
      }
    }

    intersect <- enthalpy_intersect(i, alt)

    # Add label
##    p <- p +
##      ggplot2::annotate("text",
##        x = intersect,
##        y = sat_hum_ratio(intersect, alt),
##        size = theme$text.size.axis / FONT.SCALE,
##        label =  label,
##        angle = deg + k * slope_rel_hum(intersect - 0.5,
##                                        intersect + 0.5,
##                                        100),
##        col = theme$color.axis.text)
  }

  # Plot --------------------------------------------------------------------

  if (disable.warnings == TRUE) {
    p
  } else {
    p
  }
}

