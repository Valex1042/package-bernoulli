#' Flow velocity via Torricellis theorem
#'
#' Calculates the flow velocity of a fluid through an orifice
#' according to Torricellis theorem.
#'
#' @param h Height of the fluid column above the orifice (meters, m).
#' @param g Gravitational acceleration (meters per second squared, m/s2).
#' Default value: 9.81 (Earth).
#' @param P_atm Atmospheric pressure (Pascals, Pa). Default value: 101325.
#' This parameter is included for compatibility but does not affect the calculation
#' for a reservoir open to the atmosphere.
#'
#' @return Flow velocity in meters per second (m/s).
#'
#' @details
#' Torricellis theorem gives the flow velocity of a fluid
#' through a small orifice under the effect of gravity:
#' \deqn{v = \sqrt{2 g h}}
#'
#' This formula assumes:
#' \itemize{
#' \item Incompressible fluid
#' \item No friction losses
#' \item Orifice small compared to the reservoir
#' \item Identical pressure at the surface and at the outlet (open reservoir)
#' }
#'
#' @export
#'
#' @examples
#' Velocity for a height of 10 m on Earth
#' velocity_torricelli(h = 10)
#'
#' With different gravity (e.g., Mars)
#' velocity_torricelli(h = 10, g = 3.71)
#'
#' Velocity for different heights
#' hauteurs <- c(1, 5, 10, 20)
#' vitesses <- velocity_torricelli(h = hauteurs)
#' data.frame(Hauteur = hauteurs, Vitesse = round(vitesses, 2))
#'
#' @references
#' Torricelli, E. (1644). De motu aquarum.
#'
#' Munson, B. R., Young, D. F., & Okiishi, T. H. (2006).
#' Fundamentals of Fluid Mechanics (6th ed.). Wiley.
#'
#' @seealso
#' \code{\link{bernoulli_standard}} for the complete Bernoulli equation.
#' \code{\link{flow_rate_bernoulli}} for flow calculation.
velocity_torricelli <- function(h, g = 9.81, P_atm = 101325) {
  if (h < 0) stop("La hauteur ne peut pas etre negative")
  if (!is.numeric(h)) stop("La hauteur doit etre numerique")
  if (!is.numeric(g) || g <= 0) stop("Lacceleration gravitationnelle doit etre positive")
  sqrt(2 * g * h)
}

#' Pressure difference calculation in a Venturi
#'
#' Calculates the pressure difference between two sections of a Venturi tube
#' using the Bernoulli equation.
#'
#' @param rho Fluid density (kilograms per cubic meter, kg/m3).
#' @param v1 Velocity at the upstream (wide) section (meters per second, m/s).
#' @param v2 Velocity at the constricted section (meters per second, m/s).
#'
#' @return Pressure difference delta_P = P1 - P2 in Pascals (Pa).
#' A positive value indicates that P1 > P2 (pressure drop in the constricted section).
#'
#' @details
#' The Venturi effect describes the pressure reduction when a fluid
#' flows through a converging section:
#' \deqn{\Delta P = \frac{1}{2} \rho (v_2^2 - v_1^2)}
#'
#' Typical applications:
#' \itemize{
#' \item Venturi flow meters
#' \item Sprayers
#' \item Suction systems
#' }
#'
#' Assumptions:
#' \itemize{
#' \item Steady flow
#' \item Incompressible fluid
#' \item No head losses
#' \item No height variation
#' }
#'
#' @export
#'
#' @examples
#' Pressure difference for water
#' pressure_venturi(rho = 1000, v1 = 3, v2 = 5)
#'
#'  For air at different velocities
#' pressure_venturi(rho = 1.225, v1 = 10, v2 = 20)
#'
#' Analysis for a real Venturi
#' vitesses <- seq(1, 10, by = 1)
#' delta_P <- pressure_venturi(rho = 1000, v1 = 2, v2 = vitesses)
#' plot(vitesses, delta_P / 1000, type = "b",
#' xlab = "Vitesse v2 (m/s)", ylab = "delta_P (kPa)",
#' main = "Effet Venturi: delta_P vs v2")
#'
#' @seealso
#' \code{\link{velocity_pitot}} for the Pitot tube.
#' \code{\link{bernoulli_standard}} for the general Bernoulli equation.
#'
#' @references
#' Venturi, G. B. (1797). Recherches experimentales sur le principe de la
#' communication laterale du mouvement dans les fluides.
pressure_venturi <- function(rho, v1, v2) {
  0.5 * rho * (v2^2 - v1^2)
}

#' Velocity measured by a Pitot tube
#'
#' Calculates fluid velocity from total and static pressures
#' measured by a Pitot tube.
#'
#' @param P_total Total pressure (stagnation) measured by the Pitot tube
#' (Pascals, Pa).
#' @param P_static Static pressure measured on the wall (Pascals, Pa).
#' @param rho Fluid density (kilograms per cubic meter, kg/m3).
#'
#' @return Fluid velocity in meters per second (m/s).
#'
#' @details
#' The Pitot tube measures velocity by comparing the total pressure
#' (static pressure + dynamic pressure) and the static pressure:
#' \deqn{v = \sqrt{\frac{2(P_{total} - P_{static})}{\rho}}}
#'
#' Applications:
#' \itemize{
#' \item Aeronautical anemometry
#' \item Velocity measurement in hydraulics
#' \item Industrial instrumentation
#' }
#'
#' Assumptions:
#' \itemize{
#' \item Incompressible flow (Mach < 0.3)
#' \item Correct alignment with the flow
#' \item No local disturbances
#' }
#'
#' @export
#'
#' @examples
#' Measurement for an aircraft in flight
#' velocity_pitot(P_total = 102000, P_static = 101325, rho = 1.225)
#'
#'  For a water flow
#' velocity_pitot(P_total = 150000, P_static = 101325, rho = 1000)
#'
#' Sensitivity for different pressure differentials
#' delta_P <- seq(100, 1000, by = 100)
#' vitesses_air <- velocity_pitot(P_total = 101325 + delta_P,
#' P_static = 101325,
#' rho = 1.225)
#' vitesses_eau <- velocity_pitot(P_total = 101325 + delta_P,
#' P_static = 101325,
#' rho = 1000)
#'
#' plot(delta_P, vitesses_air, type = "b", col = "blue",
#' xlab = "delta_P(Pa)", ylab = "Vitesse (m/s)",
#' main = "Sensibilite du tube de Pitot")
#' lines(delta_P, vitesses_eau, type = "b", col = "red")
#' legend("topleft", legend = c("Air", "Eau"), col = c("blue", "red"), lty = 1)
#'
#' @references
#' Pitot, H. (1732). Description dune machine pour mesurer la vitesse
#' des eaux courantes et le sillage des vaisseaux.
#'
#' Anderson, J. D. (2010). Fundamentals of Aerodynamics (5th ed.). McGraw-Hill.
velocity_pitot <- function(P_total, P_static, rho) {
  sqrt(2 * (P_total - P_static) / rho)
}

#' Flow rate calculated by the Bernoulli equation
#'
#' Calculates the volumetric flow rate from Bernoulli parameters
#' between two sections of a pipe.
#'
#' @param A Cross-sectional area (square meters, m2).
#' @param P1 Pressure at section 1 (Pascals, Pa).
#' @param P2 Pressure at section 2 (Pascals, Pa).
#' @param rho Fluid density (kilograms per cubic meter, kg/m3).
#' @param h1 Height (elevation) of section 1 (meters, m).
#' Default value: 0.
#' @param h2 Height (elevation) of section 2 (meters, m).
#' Default value: 0.
#' @param g Gravitational acceleration (meters per second squared, m/s2).
#' Default value: 9.81.
#'
#' @return Volumetric flow rate in cubic meters per second (m3/s).
#'
#' @details
#' This function solves the Bernoulli equation to calculate the average
#' velocity, then the flow rate:
#' \deqn{Q = A \times \sqrt{2\left(\frac{P_1 - P_2}{\rho} + g(h_1 - h_2)\right)}}
#'
#' Applications:
#' \itemize{
#' \item Flow calculation in pipes
#' \item Pump sizing
#' \item Hydraulic network analysis
#' }
#'
#' Limitations:
#' \itemize{
#' \item Does not account for head losses
#' \item Uniform velocity across the section (no velocity profile)
#' \item Incompressible fluid
#' }
#'
#' @export
#'
#' @examples
#' Horizontal pipe with no height difference
#' flow_rate_bernoulli(A = 0.01, P1 = 200000, P2 = 100000, rho = 1000)
#'
#' Pipe with height difference
#' flow_rate_bernoulli(A = 0.005, P1 = 101325, P2 = 101325,
#' rho = 1000, h1 = 10, h2 = 0)
#'
#' Complete case with pressure and height
#' flow_rate_bernoulli(A = 0.02, P1 = 150000, P2 = 101325,
#' rho = 1000, h1 = 5, h2 = 2)
#'
#' Influence of cross-sectional area
#' aires <- seq(0.001, 0.01, by = 0.001)
#' debits <- sapply(aires, function(A) {
#' flow_rate_bernoulli(A = A, P1 = 200000, P2 = 101325, rho = 1000)
#' })
#'
#' plot(aires, debits, type = "b",
#' xlab = "Aire (m2)", ylab = "Debit (m3/s)",
#' main = "Debit en fonction de laire de section")
#'
#' @seealso
#' \code{\link{bernoulli_standard}} to solve the Bernoulli equation.
#' \code{\link{velocity_torricelli}} for the special case of reservoirs.

flow_rate_bernoulli <- function(A, P1, P2, rho, h1 = 0, h2 = 0, g = 9.81) {
  delta_P <- P1 - P2
  delta_h <- h1 - h2
  v <- sqrt(2 * (delta_P/rho + g * delta_h))
  A * v
}

test_application <- function(){
  cat("TEST DU THEOREME DE TORRECELLI\n")
  cat("valexxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n")
  test1 <- velocity_torricelli(10)
  cat("la vitesse de torrecelli est la suivant:", round(test1, 2), "\n")
  cat("effect venturi testtttttttttttttttttttttt\n")
  delta_p1 <- pressure_venturi(rho = 1000, v1 = 3, v2 = 5)
  cat("calul attendu: deltat_p1 <- 0.5* 1000*(5^2 - 3^2)")
  cat("la valeur de delta p1 est:", round(delta_p1, 2))

  #TEST PITOT
  message("\n---teste tube de pitot ----------------------\n")
  message("\n--Avion en vol\n")
  message("\n--Avion P_total=102 kPa P_statique=101.325 kPa air\n")
  message("\n--Vitesse theorique: sqrt(2*675/1.225) = 33.20 m/s\n")
  vitesse_pitot <- velocity_pitot(P_total = 102000, P_static = 101325, rho = 1.225)
  message("\n--Vitesse calculee :", round(vitesse_pitot, 2), "m/s\n")

  cat("\n TEST DU DEBIT PAR BERNOULLI\n")
  cat("---------------------------\n")
  cat("------------ sans difference de hauteur ---------------\n")
  Q1 <- flow_rate_bernoulli(A = 0.01, P1 = 150000, P2 = 100000, rho = 1000)
  cat("Vitesse: sqrt(2*50000/1000) = 10 m/s\n")
  cat("Debit theorique: 0.01 * 10 = 0.1 m3/s\n")
  cat("Debit calcule :", round(Q1, 4), "m3/s\n")
  cat("------------ avec difference de hauteur ---------------\n")
  Q2 <- flow_rate_bernoulli(A = 0.005, P1 = 101325, P2 = 101325,
                            rho = 1000, h1 = 10, h2 = 0)
  cat("\nReservoir, A=0.005 m2, deltaP=0, deltah=10 m, eau\n")
  cat("Vitesse: sqrt(2*9.81*10) = 14.007 m/s\n")
  cat("Debit theorique: 0.005 * 14.007 = 0.07004 m3/s\n")
  cat("Debit calcule :", round(Q2, 4), "m3/s\n")
  cat("------------ cas complet ---------------\n")
  Q3 <- flow_rate_bernoulli(A = 0.02, P1 = 120000, P2 = 100000,
                            rho = 1000, h1 = 5, h2 = 2)
  cat("\nCas complet, A=0.02 m2, deltaP=20 kPa, deltah=3 m, eau\n")
  cat("Vitesse: sqrt(2*(20000/1000 + 9.81*3)) = sqrt(40 + 58.86) = 9.94 m/s\n")
  cat("Debit theorique: 0.02 * 9.94 = 0.1988 m3/s\n")
  cat("Debit calcule :", round(Q3, 4), "m3/s\n")
}
# test global
test_application()
