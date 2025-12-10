#' Standard Bernoulli Equation
#'
#' Solves the Bernoulli equation for perfect incompressible fluids.
#' The Bernoulli equation expresses the conservation of mechanical energy
#' for a steady, incompressible, and non-viscous flow.
#'
#' @param P1 Pressure at point 1 (Pascals, Pa). Can be NULL if its the variable to solve for.
#' @param v1 Velocity at point 1 (meters per second, m/s). Can be NULL if its the variable to solve for.
#' @param h1 Height at point 1 (meters, m). Can be NULL if its the variable to solve for.
#' @param P2 Pressure at point 2 (Pascals, Pa). Can be NULL if its the variable to solve for.
#' @param v2 Velocity at point 2 (meters per second, m/s). Can be NULL if its the variable to solve for.
#' @param h2 Height at point 2 (meters, m). Can be NULL if its the variable to solve for.
#' @param rho Fluid density (kilograms per cubic meter, kg/m3). Must be positive.
#' @param g Gravitational acceleration (meters per second squared, m/s2). Default value: 9.81.
#' @param solve_for Variable to solve for among: "P1", "P2", "v1", "v2", "h1", "h2".
#'
#' @return Numerical value of the unknown variable in appropriate units (Pa, m/s, or m).
#'
#' @details
#' The Bernoulli equation for a perfect incompressible fluid is written:
#' \deqn{P_1 + \frac{1}{2} \rho v_1^2 + \rho g h_1 = P_2 + \frac{1}{2} \rho v_2^2 + \rho g h_2}
#'
#' The function solves this equation for any unknown variable,
#' provided that the other 5 variables are specified.
#'
#'@section Warnings:
#'- The function assumes an incompressible and non-viscous fluid
#'- Head losses are not taken into account
#'- For velocity solutions, a negative term under the square root indicates
#'  physically impossible parameters
#'
#' @export
#'
#'@examples
#' Example 1: Calculation of P2 (downstream pressure)
#' bernoulli_standard(
#'   P1 = 101325,    # atmospheric pressure
#'   v1 = 2,         # initial velocity 2 m/s
#'   h1 = 5,         # initial height 5 m
#'   v2 = 5,         # final velocity 5 m/s
#'   h2 = 3,         # final height 3 m
#'   rho = 1000,     # water
#'   solve_for = "P2"
#' )
#'
#' # Example 2: Calculation of v2 (downstream velocity)
#' bernoulli_standard(
#'   P1 = 200000,
#'   v1 = 1,
#'   h1 = 0,
#'   P2 = 100000,
#'   h2 = 0,
#'   rho = 1000,
#'   solve_for = "v2"
#' )
#'
#' # Example 3: Calculation of h1 (initial height)
#' bernoulli_standard(
#'   P1 = 150000,
#'   v1 = 3,
#'   P2 = 101325,
#'   v2 = 8,
#'   h2 = 10,
#'   rho = 1000,
#'   solve_for = "h1"
#' )
#'
#' @references
#' Bernoulli, D. (1738). Hydrodynamica. Strasbourg.
#'
#' White, F. M. (2011). Fluid Mechanics (7th ed.). McGraw-Hill.
#'
#' @seealso
#' \code{\link{bernoulli_terms}} to calculate individual terms of the equation.

bernoulli_standard <- function(P1 = NULL, v1 = NULL, h1 = NULL,
                               P2 = NULL, v2 = NULL, h2 = NULL,
                               rho, g = 9.81, solve_for = "P2") {

  # Validation des entrees
  params <- list(P1 = P1, v1 = v1, h1 = h1, P2 = P2, v2 = v2, h2 = h2)
  null_count <- sum(sapply(params, is.null))

  if (null_count != 1) {
    stop("Exactement un parametre doit etre NULL (variable a resoudre)")
  }

  if (is.null(rho) || rho <= 0) {
    stop("La masse volumique doit etre positive")
  }

  # Equation de Bernoulli: P1 + ?rhov12 + rhogh1 = P2 + ?rhov22 + rhogh2
  bernoulli_eq <- function(P1, v1, h1, P2, v2, h2, rho, g) {
    P1 + 0.5 * rho * v1^2 + rho * g * h1 - (P2 + 0.5 * rho * v2^2 + rho * g * h2)
  }

  # Resolution selon la variable inconnue
  switch(solve_for,
         "P1" = {
           result <- P2 + 0.5 * rho * (v2^2 - v1^2) + rho * g * (h2 - h1)
         },
         "P2" = {
           result <- P1 + 0.5 * rho * (v1^2 - v2^2) + rho * g * (h1 - h2)
         },
         "v1" = {
           term <- (P2 - P1 + 0.5 * rho * v2^2 + rho * g * (h2 - h1)) / (0.5 * rho)
           if (term < 0) stop("Solution imaginaire - verifiez les parametres")
           result <- sqrt(term)
         },
         "v2" = {
           term <- (P1 - P2 + 0.5 * rho * v1^2 + rho * g * (h1 - h2)) / (0.5 * rho)
           if (term < 0) stop("Solution imaginaire - verifiez les parametres")
           result <- sqrt(term)
         },
         "h1" = {
           result <- h2 + (P2 - P1 + 0.5 * rho * (v2^2 - v1^2)) / (rho * g)
         },
         "h2" = {
           result <- h1 + (P1 - P2 + 0.5 * rho * (v1^2 - v2^2)) / (rho * g)
         },
         stop("Variable a resoudre non reconnue")
  )

  return(result)
}

#' Calculation of Bernoulli Terms
#'
#' Calculates and returns separately each term of the Bernoulli equation
#' as well as the total energy per unit volume.
#'
#' @param P Static pressure (Pascals, Pa).
#' @param v Flow velocity (meters per second, m/s).
#' @param h Height (elevation) relative to a reference (meters, m).
#' @param rho Fluid density (kilograms per cubic meter, kg/m3).
#' @param g Gravitational acceleration (meters per second squared, m/s2).
#'   Default value: 9.81.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{pressure_term}: Pressure term (Pa)
#'   \item \code{velocity_term}: Kinetic energy term (0.5 * rho * v2, in Pa)
#'   \item \code{elevation_term}: Potential energy term (rho * g * h, in Pa)
#'   \item \code{total_energy}: Total energy per unit volume (Pa)
#' }
#'
#' @details
#' The total energy per unit volume according to Bernoulli is the sum of three terms:
#' \deqn{E_{tot} = P + \frac{1}{2} \rho v^2 + \rho g h}
#'
#' This function allows analyzing the relative contribution of each term.
#'
#' @export
#'
#' @examples
#' Example 1: Water flow
#' termes <- bernoulli_terms(
#'   P = 101325,   # atmospheric pressure
#'   v = 5,        # velocity 5 m/s
#'   h = 10,       # height 10 m
#'   rho = 1000    # water
#' )
#'
#' Display results
#' str(termes)
#'
#' Percentage of each term
#' total <- termes$total_energy
#' pourcent_pression <- termes$pressure_term / total * 100
#' pourcent_cinetique <- termes$velocity_term / total * 100
#' pourcent_potentiel <- termes$elevation_term / total * 100
#'
#' cat(sprintf("Pressure: %.1f%%\n", pourcent_pression))
#' cat(sprintf("Kinetic: %.1f%%\n", pourcent_cinetique))
#' cat(sprintf("Potential: %.1f%%\n", pourcent_potentiel))
#'
#' Example 2: Analysis for different points of a flow
#' point1 <- bernoulli_terms(P = 200000, v = 2, h = 5, rho = 1000)
#' point2 <- bernoulli_terms(P = 150000, v = 8, h = 2, rho = 1000)
#'
#' Check conservation (without losses)
#' diff_energie <- point1$total_energy - point2$total_energy
#' cat(sprintf("Energy difference: %.2f Pa\n", diff_energie))
#'
#' @seealso
#'\code{\link[bernoulli]{bernoulli_standard}} for the complete Bernoulli equation.


bernoulli_terms <- function(P, v, h, rho, g = 9.81) {
  list(
    pressure_term = P,
    velocity_term = 0.5 * rho * v^2,
    elevation_term = rho * g * h,
    total_energy = P + 0.5 * rho * v^2 + rho * g * h
  )
}

# fonction des test du programme de base
test_bernoulli <- function() {
  cat("=== TESTS COMPLETS DE LEQUATION DE BERNOULLI DE BASE ===\n\n")

  # Test 1: Cas simple
  cat("Test 1 - Cas simple (P2):\n")
  res1 <- bernoulli_standard(P1 = 100000, v1 = 1, h1 = 0,
                             v2 = 2, h2 = 0, rho = 1000, solve_for = "P2")
  cat("  P2 =", round(res1, 2), "Pa\n")

  # Test 2: Avec hauteur
  cat("\nTest 2 - Avec difference de hauteur (v2):\n")
  res2 <- bernoulli_standard(P1 = 100000, v1 = 0, h1 = 10,
                             P2 = 100000, h2 = 0, rho = 1000, solve_for = "v2")
  cat("  v2 =", round(res2, 2), "m/s\n")

  # Test 3: Detail des termes
  cat("\nTest 3 - Detail des termes:\n")
  terms <- bernoulli_terms(P = 101325, v = 5, h = 10, rho = 1000)
  cat("  Pression:", round(terms$pressure_term/1000, 2), "kPa\n")
  cat("  Cinetique:", round(terms$velocity_term/1000, 2), "kPa\n")
  cat("  Hauteur:", round(terms$elevation_term/1000, 2), "kPa\n")
  cat("  Total:", round(terms$total_energy/1000, 2), "kPa\n")

  # Test 4: Conservation denergie
  cat("\nTest 4 - Conservation denergie:\n")
  P2_test <- bernoulli_standard(P1 = 150000, v1 = 2, h1 = 5,
                                v2 = 4, h2 = 3, rho = 1000, solve_for = "P2")
  E1 <- 150000 + 0.5*1000*2^2 + 1000*9.81*5
  E2 <- P2_test + 0.5*1000*4^2 + 1000*9.81*3
  cat("  Energie point 1:", round(E1, 2), "J/m3\n")
  cat("  Energie point 2:", round(E2, 2), "J/m3\n")
  cat("  Difference:", round(abs(E1-E2), 2), "J/m3\n")
  cat("  Test reussi:", abs(E1-E2) < 0.001, "\n")
}

# Executer tous les tests
test_bernoulli()
