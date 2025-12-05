#' bernoulli: Applications de l'Équation de Bernoulli en Mécanique des Fluides
#'
#' @description
#' Le package bernoulli implémente l'équation de Bernoulli et ses principales
#' applications en mécanique des fluides pour les fluides incompressibles.
#' Il fournit des outils pour calculer des vitesses, pressions, débits,
#' et analyser divers phénomènes fluidiques.
#'
#' @details
#' Le package bernoulli fournit des fonctions pour résoudre l'équation de Bernoulli
#' et ses applications classiques en mécanique des fluides.
#'
#' @section Équation de Bernoulli:
#' L'équation fondamentale implémentée est:
#' \deqn{P_1 + \frac{1}{2}\rho v_1^2 + \rho g h_1 = P_2 + \frac{1}{2}\rho v_2^2 + \rho g h_2}
#' où:
#' \itemize{
#'   \item \eqn{P} est la pression (Pa)
#'   \item \eqn{\rho} est la masse volumique (kg/m³)
#'   \item \eqn{v} est la vitesse (m/s)
#'   \item \eqn{g} est l'accélération gravitationnelle (9.81 m/s²)
#'   \item \eqn{h} est la hauteur (m)
#' }
#'
#' @section Fonctions principales:
#' \describe{
#'   \item{\code{\link{bernoulli_standard}}}{Résout l'équation de Bernoulli pour n'importe quelle variable}
#'   \item{\code{\link{velocity_torricelli}}}{Théorème de Torricelli pour les réservoirs}
#'   \item{\code{\link{pressure_venturi}}}{Effet Venturi pour les différences de pression}
#'   \item{\code{\link{velocity_pitot}}}{Tube de Pitot pour mesurer la vitesse}
#'   \item{\code{\link{flow_rate_bernoulli}}}{Calcule le débit volumique}
#'   \item{\code{\link{bernoulli_terms}}}{Calcule les termes individuels de Bernoulli}
#'   \item{\code{\link{plot_bernoulli_terms}}}{Visualise les contributions des termes}
#' }
#'
#' @section Fonctions auxiliaires:
#' \describe{
#'   \item{\code{\link{fluid_properties}}}{Propriétés des fluides communs}
#'   \item{\code{\link{validate_bernoulli}}}{Valide les hypothèses de Bernoulli}
#'   \item{\code{\link{plot_pressure_velocity}}}{Relation pression-vitesse}
#' }
#'
#' @section Applications typiques:
#' \itemize{
#'   \item Calcul de vitesse d'écoulement dans les réservoirs
#'   \item Mesure de débit dans les conduites
#'   \item Analyse de systèmes de pompage
#'   \item Dimensionnement de canalisations
#'   \item Études d'écoulements incompressibles
#'   \item Analyse de systèmes Venturi et Pitot
#'   \item Calculs hydrostatiques et hydrodynamiques
#' }
#'
#' @section Hypothèses et limitations:
#' Le package suppose:
#' \itemize{
#'   \item Fluides incompressibles (ρ constant)
#'   \item Écoulements permanents
#'   \item Pas de viscosité (fluides parfaits)
#'   \item Pas de transfert de chaleur
#'   \item Écoulement le long d'une ligne de courant
#'   \item Pas de pertes de charge par frottement
#' }
#'
#' @note
#' Pour les applications réelles avec pertes de charge, il est recommandé
#' d'utiliser un coefficient de correction ou de considérer les pertes
#' séparément.
#'
#' @author Djiogap \email{djiogapvalex1043@gmail.com}
#'
#' @references
#' Bernoulli, D. (1738). Hydrodynamica.
#'
#' White, F. M. (2011). Fluid Mechanics (7th ed.). McGraw-Hill.
#'
#' Munson, B. R., Young, D. F., & Okiishi, T. H. (2006).
#' Fundamentals of Fluid Mechanics (6th ed.). Wiley.
#'
#' Çengel, Y. A., & Cimbala, J. M. (2014). Fluid Mechanics:
#' Fundamentals and Applications (3rd ed.). McGraw-Hill.
#'
#' @seealso
#' Pour des analyses plus avancées:
#' \itemize{
#'   \item \pkg{fluids} : Calculs de mécanique des fluides avancés
#'   \item \pkg{hydraulics} : Fonctions pour l'ingénierie hydraulique
#'   \item \pkg{pracma} : Mathématiques pratiques pour l'ingénierie
#'   \item \pkg{units} : Gestion des unités physiques
#' }
#'
#' Sites web utiles:
#' \itemize{
#'   \item \url{https://en.wikipedia.org/wiki/Bernoulli%27s_principle}
#'   \item \url{https://www.engineeringtoolbox.com/}
#' }
#'
#' @examples
#' # Charger le package
#' library(bernoulli)
#'
#' # Exemple 1: Vitesse de sortie d'un réservoir (Torricelli)
#' v_torricelli <- velocity_torricelli(h = 10)
#' cat("Vitesse de Torricelli pour h = 10 m:", round(v_torricelli, 2), "m/s\n")
#' cat("Soit", round(v_torricelli * 3.6, 1), "km/h\n\n")
#'
#' # Exemple 2: Calcul de débit dans une conduite
#' Q <- flow_rate_bernoulli(
#'   A = 0.01,      # Aire de section: 0.01 m² (diamètre ~11.3 cm)
#'   P1 = 200000,   # Pression amont: 200 kPa
#'   P2 = 100000,   # Pression aval: 100 kPa
#'   rho = 1000     # Eau: 1000 kg/m³
#' )
#' cat("Débit dans la conduite:", round(Q, 4), "m³/s\n")
#' cat("Soit", round(Q * 1000, 1), "L/s\n\n")
#'
#' # Exemple 3: Effet Venturi
#' delta_P <- pressure_venturi(rho = 1000, v1 = 2, v2 = 8)
#' cat("Différence de pression dans un Venturi:", round(delta_P/1000, 1), "kPa\n\n")
#'
#' # Exemple 4: Tube de Pitot (mesure de vitesse d'avion)
#' vitesse_avion <- velocity_pitot(
#'   P_total = 102000,
#'   P_static = 101325,
#'   rho = 1.225
#' )
#' cat("Vitesse de l'avion:", round(vitesse_avion, 1), "m/s\n")
#' cat("Soit", round(vitesse_avion * 3.6, 0), "km/h\n\n")
#'
#' # Exemple 5: Analyse des termes de Bernoulli
#' termes <- bernoulli_terms(P = 150000, v = 5, h = 15, rho = 1000)
#' cat("Analyse des termes de Bernoulli:\n")
#' cat("  Terme de pression:", round(termes$pressure_term/1000, 1), "kPa\n")
#' cat("  Terme cinétique:", round(termes$velocity_term/1000, 1), "kPa\n")
#' cat("  Terme de hauteur:", round(termes$elevation_term/1000, 1), "kPa\n")
#' cat("  Énergie totale:", round(termes$total_energy/1000, 1), "kPa\n")
#'
#' # Exemple 6: Résolution de l'équation complète
#' P2 <- bernoulli_standard(
#'   P1 = 101325, v1 = 1, h1 = 10,
#'   v2 = 8, h2 = 2, rho = 1000,
#'   solve_for = "P2"
#' )
#' cat("\nPression au point 2:", round(P2/1000, 1), "kPa\n")
#'
#' # Pour plus d'exemples, voir les pages d'aide spécifiques
#' # ?velocity_torricelli
#' # ?bernoulli_standard
#' # ?plot_bernoulli_terms
#'
#' @keywords package
# ou pour des fonctions internes
#' @keywords internal

"_PACKAGE"
