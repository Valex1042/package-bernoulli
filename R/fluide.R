#' @title Propriétés physiques des fluides standards en fonction de la température
#'
#' @description
#' Retourne les propriétés physiques principales de l'eau, l'air et l'huile
#' en mécanique des fluides, avec correction pour la température spécifiée.
#'
#' @details
#' Cette fonction fournit les propriétés des trois fluides de base :
#' \itemize{
#'   \item \code{name} : Nom complet du fluide
#'   \item \code{rho} : Masse volumique en kg/m³
#'   \item \code{mu} : Viscosité dynamique en Pa·s
#'   \item \code{nu} : Viscosité cinématique en m²/s (μ/ρ)
#'   \item \code{gamma} : Poids volumique en N/m³ (ρ × g)
#'   \item \code{T} : Température effectivement utilisée en °C
#' }
#'
#' @param fluid Caractère spécifiant le fluide souhaité. Options :
#'   \describe{
#'     \item{"water"}{Eau pure}
#'     \item{"air"}{Air sec}
#'     \item{"oil"}{Huile SAE 30}
#'   }
#'   Par défaut : "water"
#'
#' @param T Température du fluide en degrés Celsius (°C).
#'   Plages recommandées :
#'   \itemize{
#'     \item Eau : 0 à 100°C
#'     \item Air : -20 à 100°C
#'     \item Huile : 0 à 100°C
#'   }
#'   Par défaut : 20
#'
#' @param g Accélération gravitationnelle en m/s². Par défaut : 9.81
#'
#' @return
#' Une liste de classe "fluid_properties" contenant :
#' \describe{
#'   \item{name}{Nom complet}
#'   \item{rho}{Masse volumique (kg/m³)}
#'   \item{mu}{Viscosité dynamique (Pa·s)}
#'   \item{nu}{Viscosité cinématique (m²/s)}
#'   \item{gamma}{Poids volumique (N/m³)}
#'   \item{T}{Température (°C)}
#'   \item{fluid_type}{Type de fluide}
#' }
#'
#' @export
#'
#' @examples
#' # Propriétés de l'eau à différentes températures
#' fluid_properties("water", 20)
#' fluid_properties("water", 50)
#' fluid_properties("water", 100)
#'
#' # Propriétés de l'air
#' fluid_properties("air", 0)
#' fluid_properties("air", 25)
#'
#' # Propriétés de l'huile
#' fluid_properties("oil", 20)
#' fluid_properties("oil", 40)
#'
#' # Utilisation dans un calcul
#' props <- fluid_properties("water", 80)
#' cat("Densité de l'eau à 80°C:", props$rho, "kg/m³\n")
#' cat("Viscosité:", props$mu, "Pa·s\n")
#'
#' @seealso
#' \code{\link{bernoulli_standard}}
#'
#' @references
#' White, F. M. (2011). \emph{Fluid Mechanics} (7th ed.). McGraw-Hill.
#'
#' @family fluid_properties_functions
#'
#' @keywords fluid properties density viscosity temperature
fluid_properties <- function(fluid = "water", T = 20, g = 9.81) {

  # Validation des entrées
  if (!is.numeric(T) || length(T) != 1) {
    stop("La température doit être une valeur numérique unique")
  }

  if (!is.numeric(g) || g <= 0) {
    stop("L'accélération gravitationnelle doit être positive")
  }

  # Base de données des fluides avec correction de température
  fluids_db <- list(

    water = function(T) {
      # Propriétés de l'eau pure
      if (T < 0 || T > 100) {
        warning("Température hors plage recommandée pour l'eau (0-100°C)")
      }

      # Densité de l'eau (kg/m³) - approximation polynomiale
      # Valeur maximale à 4°C : 1000 kg/m³
      if (T <= 4) {
        # 0-4°C : densité augmente avec la température
        rho <- 999.87 + 0.032 * T - 0.0045 * T^2
      } else if (T <= 20) {
        # 4-20°C : densité diminue légèrement
        rho <- 1000 - 0.2 * (T - 4) - 0.006 * (T - 4)^2
      } else {
        # 20-100°C : décroissance plus rapide
        rho <- 998.2 - 0.2 * (T - 20) - 0.006 * (T - 20)^2
      }

      # Viscosité dynamique (Pa·s) - équation simplifiée
      mu <- 0.00179 / (1 + 0.0337 * T + 0.000221 * T^2)

      list(
        name = "Eau pure",
        rho = max(rho, 958.4),  # Minimum à 100°C
        mu = mu,
        T_ref = T
      )
    },

    air = function(T) {
      # Propriétés de l'air sec à pression atmosphérique
      if (T < -20 || T > 100) {
        warning("Température hors plage recommandée pour l'air (-20 à 100°C)")
      }

      # Densité (kg/m³) - loi des gaz parfaits
      # P = 101325 Pa, R = 287.05 J/(kg·K)
      T_k <- T + 273.15
      rho <- 101325 / (287.05 * T_k)

      # Viscosité dynamique (Pa·s) - formule simplifiée
      # Approximation linéaire entre -20°C et 100°C
      mu <- 1.72e-5 * (1 + 0.0025 * T)

      list(
        name = "Air sec",
        rho = rho,
        mu = mu,
        T_ref = T
      )
    },

    oil = function(T) {
      # Propriétés de l'huile SAE 30
      if (T < 0 || T > 100) {
        warning("Température hors plage recommandée pour l'huile (0-100°C)")
      }

      # Densité à 20°C
      rho_20 <- 875

      # Coefficient d'expansion thermique
      alpha <- 0.0007

      # Densité corrigée
      rho <- rho_20 * (1 - alpha * (T - 20))

      # Viscosité dynamique - variation exponentielle
      mu_20 <- 0.29
      mu <- mu_20 * exp(-0.035 * (T - 20))

      list(
        name = "Huile SAE 30",
        rho = max(rho, 800),  # Limite inférieure
        mu = max(mu, 0.01),   # Limite inférieure
        T_ref = T
      )
    }
  )

  # Vérifier si le fluide est supporté
  if (!fluid %in% names(fluids_db)) {
    supported <- paste(names(fluids_db), collapse = ", ")
    stop("Fluide non supporté. Choisissez parmi: ", supported)
  }

  # Obtenir les propriétés à la température spécifiée
  props <- fluids_db[[fluid]](T)

  # Calculer les propriétés dérivées
  nu <- props$mu / props$rho  # Viscosité cinématique
  gamma <- props$rho * g      # Poids volumique

  # Retourner la liste complète
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

#' @title Affichage des propriétés d'un fluide
#' @description Méthode print pour les objets fluid_properties
#' @param x Objet fluid_properties
#' @param ... Arguments supplémentaires
#' @export
print.fluid_properties <- function(x, ...) {
  cat("Propriétés du fluide :", x$name, "\n")
  cat("=", rep("=", nchar(x$name) + 20), "\n", sep = "")
  cat(sprintf("Température         : %6.1f °C\n", x$T))
  cat(sprintf("Masse volumique (ρ) : %6.2f kg/m³\n", x$rho))
  cat(sprintf("Viscosité dyn. (μ)  : %6.2e Pa·s\n", x$mu))
  cat(sprintf("Viscosité cin. (ν)  : %6.2e m²/s\n", x$nu))
  cat(sprintf("Poids volumique (γ) : %6.1f N/m³\n", x$gamma))
  cat("=", rep("=", nchar(x$name) + 20), "\n", sep = "")
}
fluid_properties("water")
