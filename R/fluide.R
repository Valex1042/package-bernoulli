#' @title Physical properties of standard fluids as a function of temperature
#'
#' @description
#' Returns the main physical properties of water, air, and oil
#' in fluid mechanics, with correction for the specified temperature.
#'
#'@details
#' This function provides the properties of three basic fluids:
#' \itemize{
#'   \item \code{name} : Full name of the fluid
#'   \item \code{rho} : Density in kg/m3
#'   \item \code{mu} : Dynamic viscosity in Pa.s
#'   \item \code{nu} : Kinematic viscosity in m2/s (nu/?)
#'  \item \code{gamma} : Specific weight in N/m3 (rho ? g)
#'   \item \code{T} : Actual temperature used in ?C
#'}
#'
#' @param fluid Character specifying the desired fluid. Options:
#'   \describe{
#'     \item{"water"}{Pure water}
#'    \item{"air"}{Dry air}
#'     \item{"oil"}{SAE 30 oil}
#'   }
#'   Default: "water"
#'
#' @param T Fluid temperature in degrees Celsius (?C).
#'   Recommended ranges:
#'   \itemize{
#'     \item Water: 0 to 100?C
#'     \item Air: -20 to 100?C
#'     \item Oil: 0 to 100?C
#'   }
#'   Default: 20
#'
#' @param g Gravitational acceleration in m/s2. Default: 9.81
#'
#' @return
#' A list of class "fluid_properties" containing:
#' \describe{
#'   \item{name}{Full name}
#'   \item{rho}{Density (kg/m3)}
#'   \item{mu}{Dynamic viscosity (Pa.s)}
#'   \item{nu}{Kinematic viscosity (m2/s)}
#'   \item{gamma}{Specific weight (N/m3)}
#'   \item{T}{Temperature (?C)}
#'   \item{fluid_type}{Fluid type}
#'}
#'
#'@export
#
#'@examples
#'Properties of water at different temperatures
#'fluid_properties("water", 20)
#'fluid_properties("water", 50)
#'fluid_properties("water", 100)
#'
#'Properties of air
#'fluid_properties("air", 0)
#'fluid_properties("air", 25)
#'
#'Properties of oil
#'fluid_properties("oil", 20)
#'fluid_properties("oil", 40)
#'
#' Use in a calculation
#' props <- fluid_properties("water", 80)
#' cat("Density of water at 80?C:", props$rho, "kg/m3\n")
#' cat("Viscosity:", props$mu, "Pa.s\n")
#'
#' @seealso
#' \code{\link{bernoulli_standard}}
#'
#' @references
#' White, F. M. (2011). \emph{Fluid Mechanics} (7th ed.). McGraw-Hill.
#'
#' @family fluid_properties_functions
#
#'@keywords fluid properties density viscosity temperature
#'
fluid_properties <- function(fluid = "water", T = 20, g = 9.81) {

  # Validation des entrees
  if (!is.numeric(T) || length(T) != 1) {
    stop("La temperature doit etre une valeur numerique unique")
  }

  if (!is.numeric(g) || g <= 0) {
    stop("Lacceleration gravitationnelle doit etre positive")
  }

  # Base de donnees des fluides avec correction de temperature
  fluids_db <- list(

    water = function(T) {
      # Proprietes de leau pure
      if (T < 0 || T > 100) {
        warning("Temperature hors plage recommandee pour leau (0-100?C)")
      }

      # Densite de leau (kg/m3) - approximation polynomiale
      # Valeur maximale a 4?C : 1000 kg/m3
      if (T <= 4) {
        # 0-4?C : densite augmente avec la temperature
        rho <- 999.87 + 0.032 * T - 0.0045 * T^2
      } else if (T <= 20) {
        # 4-20?C : densite diminue legerement
        rho <- 1000 - 0.2 * (T - 4) - 0.006 * (T - 4)^2
      } else {
        # 20-100?C : decroissance plus rapide
        rho <- 998.2 - 0.2 * (T - 20) - 0.006 * (T - 20)^2
      }

      # Viscosite dynamique (Pa.s) - equation simplifiee
      mu <- 0.00179 / (1 + 0.0337 * T + 0.000221 * T^2)

      list(
        name = "Eau pure",
        rho = max(rho, 958.4),  # Minimum a 100?C
        mu = mu,
        T_ref = T
      )
    },

    air = function(T) {
      # Proprietes de lair sec a pression atmospherique
      if (T < -20 || T > 100) {
        warning("Temperature hors plage recommandee pour lair (-20 a 100?C)")
      }

      # Densite (kg/m3) - loi des gaz parfaits
      # P = 101325 Pa, R = 287.05 J/(kg.K)
      T_k <- T + 273.15
      rho <- 101325 / (287.05 * T_k)

      # Viscosite dynamique (Pa.s) - formule simplifiee
      # Approximation lineaire entre -20?C et 100?C
      mu <- 1.72e-5 * (1 + 0.0025 * T)

      list(
        name = "Air sec",
        rho = rho,
        mu = mu,
        T_ref = T
      )
    },

    oil = function(T) {
      # Proprietes de lhuile SAE 30
      if (T < 0 || T > 100) {
        warning("Temperature hors plage recommandee pour lhuile (0-100?C)")
      }

      # Densite a 20?C
      rho_20 <- 875

      # Coefficient dexpansion thermique
      alpha <- 0.0007

      # Densite corrigee
      rho <- rho_20 * (1 - alpha * (T - 20))

      # Viscosite dynamique - variation exponentielle
      mu_20 <- 0.29
      mu <- mu_20 * exp(-0.035 * (T - 20))

      list(
        name = "Huile SAE 30",
        rho = max(rho, 800),  # Limite inferieure
        mu = max(mu, 0.01),   # Limite inferieure
        T_ref = T
      )
    }
  )

  # Verifier si le fluide est supporte
  if (!fluid %in% names(fluids_db)) {
    supported <- paste(names(fluids_db), collapse = ", ")
    stop("Fluide non supporte. Choisissez parmi: ", supported)
  }

  # Obtenir les proprietes a la temperature specifiee
  props <- fluids_db[[fluid]](T)

  # Calculer les proprietes derivees
  nu <- props$mu / props$rho  # Viscosite cinematique
  gamma <- props$rho * g      # Poids volumique

  # Retourner la liste complete
  result <- list(
    name = props$name,
    rho = round(props$rho, 2),
    mu = signif(props$mu, 3),
    nu = signif(nu, 4),
    gamma = round(gamma, 1),
    T = T,
    fluid_type = fluid,
    g = g
  )

  class(result) <- "fluid_properties"
  return(result)
}

#' @title Display fluid properties
#' @description Print method for fluid_properties objects
#' @param x fluid_properties object
#' @param ... Additional arguments
#' @export
print.fluid_properties <- function(x, ...) {
  cat("Proprietes du fluide :", x$name, "\n")
  cat("=", rep("=", nchar(x$name) + 20), "\n", sep = "")
  cat(sprintf("Temperature         : %6.1f ?C\n", x$T))
  cat(sprintf("Masse volumique (rho) : %6.2f kg/m3\n", x$rho))
  cat(sprintf("Viscosite dyn. (nu)  : %6.2e Pa.s\n", x$mu))
  cat(sprintf("Viscosite cin. (?)  : %6.2e m2/s\n", x$nu))
  cat(sprintf("Poids volumique  : %6.1f N/m3\n", x$gamma))
  cat("=", rep("=", nchar(x$name) + 20), "\n", sep = "")
}
fluid_properties("water")
