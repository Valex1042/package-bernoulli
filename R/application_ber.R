#' Vitesse d'écoulement par le théorème de Torricelli
#'
#' Calcule la vitesse d'écoulement d'un fluide à travers un orifice
#' selon le théorème de Torricelli.
#'
#' @param h Hauteur de la colonne de fluide au-dessus de l'orifice (mètres, m).
#' @param g Accélération gravitationnelle (mètres par seconde carrée, m/s²).
#'   Valeur par défaut: 9.81 (Terre).
#' @param P_atm Pression atmosphérique (Pascals, Pa). Valeur par défaut: 101325.
#'   Ce paramètre est inclus pour compatibilité mais n'affecte pas le calcul
#'   pour un réservoir ouvert à l'atmosphère.
#'
#' @return Vitesse d'écoulement en mètres par seconde (m/s).
#'
#' @details
#' Le théorème de Torricelli donne la vitesse d'écoulement d'un fluide
#' par un petit orifice sous l'effet de la gravité :
#' \deqn{v = \sqrt{2 g h}}
#'
#' Cette formule suppose:
#' \itemize{
#'   \item Fluide incompressible
#'   \item Pas de pertes par frottement
#'   \item Orifice petit par rapport au réservoir
#'   \item Pression identique à la surface et à la sortie (réservoir ouvert)
#' }
#'
#' @export
#'
#' @examples
#' # Vitesse pour une hauteur de 10 m sur Terre
#' velocity_torricelli(h = 10)
#'
#' # Avec une gravité différente (ex: Mars)
#' velocity_torricelli(h = 10, g = 3.71)
#'
#' # Vitesse pour différentes hauteurs
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
#' \code{\link{bernoulli_standard}} pour l'équation de Bernoulli complète.
#' \code{\link{flow_rate_bernoulli}} pour le calcul du débit.
velocity_torricelli <- function(h, g = 9.81, P_atm = 101325) {
  if (h < 0) stop("La hauteur ne peut pas être négative")
  sqrt(2 * g * h)
}

#' Calcul de la différence de pression dans un venturi
#'
#' Calcule la différence de pression entre deux sections d'un tube
#' de Venturi en utilisant l'équation de Bernoulli.
#'
#' @param rho Masse volumique du fluide (kilogrammes par mètre cube, kg/m³).
#' @param v1 Vitesse à la section amont (large) (mètres par seconde, m/s).
#' @param v2 Vitesse à la section rétrécie (mètres par seconde, m/s).
#'
#' @return Différence de pression ΔP = P1 - P2 en Pascals (Pa).
#'   Une valeur positive indique que P1 > P2 (dépression dans la section rétrécie).
#'
#' @details
#' L'effet Venturi décrit la réduction de pression lorsqu'un fluide
#' s'écoule à travers une section convergente :
#' \deqn{\Delta P = \frac{1}{2} \rho (v_2^2 - v_1^2)}
#'
#' Applications typiques:
#' \itemize{
#'   \item Débitmètres Venturi
#'   \item Pulvérisateurs
#'   \item Systèmes d'aspiration
#' }
#'
#' Hypothèses:
#' \itemize{
#'   \item Écoulement permanent
#'   \item Fluide incompressible
#'   \item Pas de pertes de charge
#'   \item Pas de variation de hauteur
#' }
#'
#' @export
#'
#' @examples
#' # Différence de pression pour l'eau
#' pressure_venturi(rho = 1000, v1 = 3, v2 = 5)
#'
#' # Pour l'air à différentes vitesses
#' pressure_venturi(rho = 1.225, v1 = 10, v2 = 20)
#'
#' # Analyse pour un Venturi réel
#' vitesses <- seq(1, 10, by = 1)
#' delta_P <- pressure_venturi(rho = 1000, v1 = 2, v2 = vitesses)
#' plot(vitesses, delta_P / 1000, type = "b",
#'      xlab = "Vitesse v2 (m/s)", ylab = "ΔP (kPa)",
#'      main = "Effet Venturi: ΔP vs v2")
#'
#' @seealso
#' \code{\link{velocity_pitot}} pour le tube de Pitot.
#' \code{\link{bernoulli_standard}} pour l'équation de Bernoulli générale.
#'
#' @references
#' Venturi, G. B. (1797). Recherches expérimentales sur le principe de la
#' communication latérale du mouvement dans les fluides.
pressure_venturi <- function(rho, v1, v2) {
  0.5 * rho * (v2^2 - v1^2)
}

#' Vitesse mesurée par un tube de Pitot
#'
#' Calcule la vitesse d'un fluide à partir des pressions totale et statique
#' mesurées par un tube de Pitot.
#'
#' @param P_total Pression totale (stagnation) mesurée par le tube de Pitot
#'   (Pascals, Pa).
#' @param P_static Pression statique mesurée sur la paroi (Pascals, Pa).
#' @param rho Masse volumique du fluide (kilogrammes par mètre cube, kg/m³).
#'
#' @return Vitesse du fluide en mètres par seconde (m/s).
#'
#' @details
#' Le tube de Pitot mesure la vitesse en comparant la pression totale
#' (pression statique + pression dynamique) et la pression statique:
#' \deqn{v = \sqrt{\frac{2(P_{total} - P_{static})}{\rho}}}
#'
#' Applications:
#' \itemize{
#'   \item Anémométrie aéronautique
#'   \item Mesure de vitesse en hydraulique
#'   \item Instrumentation industrielle
#' }
#'
#' Hypothèses:
#' \itemize{
#'   \item Écoulement incompressible (Mach < 0.3)
#'   \item Alignement correct avec l'écoulement
#'   \item Pas de perturbations locales
#' }
#'
#' @export
#'
#' @examples
#' # Mesure pour un avion en vol
#' velocity_pitot(P_total = 102000, P_static = 101325, rho = 1.225)
#'
#' # Pour un écoulement d'eau
#' velocity_pitot(P_total = 150000, P_static = 101325, rho = 1000)
#'
#' # Sensibilité pour différentes pressions différentielles
#' delta_P <- seq(100, 1000, by = 100)
#' vitesses_air <- velocity_pitot(P_total = 101325 + delta_P,
#'                               P_static = 101325,
#'                               rho = 1.225)
#' vitesses_eau <- velocity_pitot(P_total = 101325 + delta_P,
#'                               P_static = 101325,
#'                               rho = 1000)
#'
#' plot(delta_P, vitesses_air, type = "b", col = "blue",
#'      xlab = "ΔP (Pa)", ylab = "Vitesse (m/s)",
#'      main = "Sensibilité du tube de Pitot")
#' lines(delta_P, vitesses_eau, type = "b", col = "red")
#' legend("topleft", legend = c("Air", "Eau"), col = c("blue", "red"), lty = 1)
#'
#' @references
#' Pitot, H. (1732). Description d'une machine pour mesurer la vitesse
#' des eaux courantes et le sillage des vaisseaux.
#'
#' Anderson, J. D. (2010). Fundamentals of Aerodynamics (5th ed.). McGraw-Hill.
velocity_pitot <- function(P_total, P_static, rho) {
  sqrt(2 * (P_total - P_static) / rho)
}

#' Débit calculé par l'équation de Bernoulli
#'
#' Calcule le débit volumique à partir des paramètres de Bernoulli
#' entre deux sections d'une conduite.
#'
#' @param A Aire de la section transversale (mètres carrés, m²).
#' @param P1 Pression à la section 1 (Pascals, Pa).
#' @param P2 Pression à la section 2 (Pascals, Pa).
#' @param rho Masse volumique du fluide (kilogrammes par mètre cube, kg/m³).
#' @param h1 Hauteur (élévation) de la section 1 (mètres, m).
#'   Valeur par défaut: 0.
#' @param h2 Hauteur (élévation) de la section 2 (mètres, m).
#'   Valeur par défaut: 0.
#' @param g Accélération gravitationnelle (mètres par seconde carrée, m/s²).
#'   Valeur par défaut: 9.81.
#'
#' @return Débit volumique en mètres cubes par seconde (m³/s).
#'
#' @details
#' Cette fonction résout l'équation de Bernoulli pour calculer la vitesse
#' moyenne, puis le débit:
#' \deqn{Q = A \times \sqrt{2\left(\frac{P_1 - P_2}{\rho} + g(h_1 - h_2)\right)}}
#'
#' Applications:
#' \itemize{
#'   \item Calcul de débit dans les conduites
#'   \item Dimensionnement de pompes
#'   \item Analyse de réseaux hydrauliques
#' }
#'
#' Limitations:
#' \itemize{
#'   \item Ne tient pas compte des pertes de charge
#'   \item Vitesse uniforme sur la section (pas de profil de vitesse)
#'   \item Fluide incompressible
#' }
#'
#' @export
#'
#' @examples
#' # Conduite horizontale sans différence de hauteur
#' flow_rate_bernoulli(A = 0.01, P1 = 200000, P2 = 100000, rho = 1000)
#'
#' # Conduite avec différence de hauteur
#' flow_rate_bernoulli(A = 0.005, P1 = 101325, P2 = 101325,
#'                    rho = 1000, h1 = 10, h2 = 0)
#'
#' # Cas complet avec pression et hauteur
#' flow_rate_bernoulli(A = 0.02, P1 = 150000, P2 = 101325,
#'                    rho = 1000, h1 = 5, h2 = 2)
#'
#' # Influence de l'aire de section
#' aires <- seq(0.001, 0.01, by = 0.001)
#' debits <- sapply(aires, function(A) {
#'   flow_rate_bernoulli(A = A, P1 = 200000, P2 = 101325, rho = 1000)
#' })
#'
#' plot(aires, debits, type = "b",
#'      xlab = "Aire (m²)", ylab = "Débit (m³/s)",
#'      main = "Débit en fonction de l'aire de section")
#'
#' @seealso
#' \code{\link{bernoulli_standard}} pour résoudre l'équation de Bernoulli.
#' \code{\link{velocity_torricelli}} pour le cas particulier des réservoirs.
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

  # TEST PITOT
  message("\n---teste tube de pitot ----------------------\n")
  message("\n--Avion en vol\n")
  message("\n--Avion P_total=102 kPa P_statique=101.325 kPa air\n")
  message("\n--Vitesse théorique: sqrt(2*675/1.225) = 33.20 m/s\n")
  vitesse_pitot <- velocity_pitot(P_total = 102000, P_static = 101325, rho = 1.225)
  message("\n--Vitesse calculée :", round(vitesse_pitot, 2), "m/s\n")


  cat("\n TEST DU DÉBIT PAR BERNOULLI\n")
  cat("---------------------------\n")
  cat("------------ sans difference de hauteur ---------------\n")
  Q1 <- flow_rate_bernoulli(A = 0.01, P1 = 150000, P2 = 100000, rho = 1000)
  cat("Vitesse: sqrt(2*50000/1000) = 10 m/s\n")
  cat("Débit théorique: 0.01 * 10 = 0.1 m³/s\n")
  cat("Débit calculé  :", round(Q1, 4), "m³/s\n")
  cat("------------ avec difference de hauteur ---------------\n")
  Q2 <- flow_rate_bernoulli(A = 0.005, P1 = 101325, P2 = 101325,
                            rho = 1000, h1 = 10, h2 = 0)
  cat("\nRéservoir, A=0.005 m², deltaP=0, deltah=10 m, eau\n")
  cat("Vitesse: sqrt(2*9.81*10) = 14.007 m/s\n")
  cat("Débit théorique: 0.005 * 14.007 = 0.07004 m³/s\n")
  cat("Débit calculé  :", round(Q2, 4), "m³/s\n")
  cat("------------ cas complet ---------------\n")
  Q3 <- flow_rate_bernoulli(A = 0.02, P1 = 120000, P2 = 100000,
                            rho = 1000, h1 = 5, h2 = 2)
  cat("\nCas complet, A=0.02 m², deltaP=20 kPa, deltah=3 m, eau\n")
  cat("Vitesse: sqrt(2*(20000/1000 + 9.81*3)) = sqrt(40 + 58.86) = 9.94 m/s\n")
  cat("Débit théorique: 0.02 * 9.94 = 0.1988 m³/s\n")
  cat("Débit calculé  :", round(Q3, 4), "m³/s\n")
}
#  teste gobale
test_application()
