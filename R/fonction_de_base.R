#' Équation de Bernoulli Standard
#'
#' Résout l'équation de Bernoulli pour les fluides parfaits incompressibles.
#' L'équation de Bernoulli exprime la conservation de l'énergie mécanique
#' pour un écoulement permanent, incompressible et non visqueux.
#'
#' @param P1 Pression au point 1 (Pascals, Pa). Peut être NULL si c'est la variable à résoudre.
#' @param v1 Vitesse au point 1 (mètres par seconde, m/s). Peut être NULL si c'est la variable à résoudre.
#' @param h1 Hauteur au point 1 (mètres, m). Peut être NULL si c'est la variable à résoudre.
#' @param P2 Pression au point 2 (Pascals, Pa). Peut être NULL si c'est la variable à résoudre.
#' @param v2 Vitesse au point 2 (mètres par seconde, m/s). Peut être NULL si c'est la variable à résoudre.
#' @param h2 Hauteur au point 2 (mètres, m). Peut être NULL si c'est la variable à résoudre.
#' @param rho Masse volumique du fluide (kilogrammes par mètre cube, kg/m³). Doit être positive.
#' @param g Accélération gravitationnelle (mètres par seconde carrée, m/s²). Valeur par défaut: 9.81.
#' @param solve_for Variable à résoudre parmi: "P1", "P2", "v1", "v2", "h1", "h2".
#'
#' @return Valeur numérique de la variable inconnue dans les unités appropriées (Pa, m/s, ou m).
#'
#' @details
#' L'équation de Bernoulli pour un fluide parfait incompressible s'écrit:
#' \deqn{P_1 + \frac{1}{2} \rho v_1^2 + \rho g h_1 = P_2 + \frac{1}{2} \rho v_2^2 + \rho g h_2}
#'
#' La fonction résout cette équation pour n'importe quelle variable inconnue,
#' à condition que les 5 autres variables soient spécifiées.
#'
#' @section Avertissements:
#' - La fonction suppose un fluide incompressible et non visqueux
#' - Les pertes de charge ne sont pas prises en compte
#' - Pour les solutions de vitesse, un terme négatif sous la racine indique
#'   des paramètres physiquement impossibles
#'
#' @export
#'
#' @examples
#' # Exemple 1: Calcul de P2 (pression en aval)
#' bernoulli_standard(
#'   P1 = 101325,    # pression atmosphérique
#'   v1 = 2,         # vitesse initiale 2 m/s
#'   h1 = 5,         # hauteur initiale 5 m
#'   v2 = 5,         # vitesse finale 5 m/s
#'   h2 = 3,         # hauteur finale 3 m
#'   rho = 1000,     # eau
#'   solve_for = "P2"
#' )
#'
#' # Exemple 2: Calcul de v2 (vitesse en aval)
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
#' # Exemple 3: Calcul de h1 (hauteur initiale)
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
#' \code{\link{bernoulli_terms}} pour calculer les termes individuels de l'équation.
#'
bernoulli_standard <- function(P1 = NULL, v1 = NULL, h1 = NULL,
                               P2 = NULL, v2 = NULL, h2 = NULL,
                               rho, g = 9.81, solve_for = "P2") {

  # Validation des entrées
  params <- list(P1 = P1, v1 = v1, h1 = h1, P2 = P2, v2 = v2, h2 = h2)
  null_count <- sum(sapply(params, is.null))

  if (null_count != 1) {
    stop("Exactement un paramètre doit être NULL (variable à résoudre)")
  }

  if (is.null(rho) || rho <= 0) {
    stop("La masse volumique doit être positive")
  }

  # Équation de Bernoulli: P1 + ½rhov1² + rhogh1 = P2 + ½rhov2² + rhogh2
  bernoulli_eq <- function(P1, v1, h1, P2, v2, h2, rho, g) {
    P1 + 0.5 * rho * v1^2 + rho * g * h1 - (P2 + 0.5 * rho * v2^2 + rho * g * h2)
  }

  # Résolution selon la variable inconnue
  switch(solve_for,
         "P1" = {
           result <- P2 + 0.5 * rho * (v2^2 - v1^2) + rho * g * (h2 - h1)
         },
         "P2" = {
           result <- P1 + 0.5 * rho * (v1^2 - v2^2) + rho * g * (h1 - h2)
         },
         "v1" = {
           term <- (P2 - P1 + 0.5 * rho * v2^2 + rho * g * (h2 - h1)) / (0.5 * rho)
           if (term < 0) stop("Solution imaginaire - vérifiez les paramètres")
           result <- sqrt(term)
         },
         "v2" = {
           term <- (P1 - P2 + 0.5 * rho * v1^2 + rho * g * (h1 - h2)) / (0.5 * rho)
           if (term < 0) stop("Solution imaginaire - vérifiez les paramètres")
           result <- sqrt(term)
         },
         "h1" = {
           result <- h2 + (P2 - P1 + 0.5 * rho * (v2^2 - v1^2)) / (rho * g)
         },
         "h2" = {
           result <- h1 + (P1 - P2 + 0.5 * rho * (v1^2 - v2^2)) / (rho * g)
         },
         stop("Variable à résoudre non reconnue")
  )

  return(result)
}

#' Calcul des Termes de Bernoulli
#'
#' Calcule et retourne séparément chaque terme de l'équation de Bernoulli
#' ainsi que l'énergie totale par unité de volume.
#'
#' @param P Pression statique (Pascals, Pa).
#' @param v Vitesse d'écoulement (mètres par seconde, m/s).
#' @param h Hauteur (élévation) par rapport à une référence (mètres, m).
#' @param rho Masse volumique du fluide (kilogrammes par mètre cube, kg/m³).
#' @param g Accélération gravitationnelle (mètres par seconde carrée, m/s²).
#'   Valeur par défaut: 9.81.
#'
#' @return Une liste contenant:
#' \itemize{
#'   \item \code{pressure_term}: Terme de pression (Pa)
#'   \item \code{velocity_term}: Terme d'énergie cinétique (0.5 * ρ * v², en Pa)
#'   \item \code{elevation_term}: Terme d'énergie potentielle (ρ * g * h, en Pa)
#'   \item \code{total_energy}: Énergie totale par unité de volume (Pa)
#' }
#'
#' @details
#' L'énergie totale par unité de volume selon Bernoulli est la somme de trois termes:
#' \deqn{E_{tot} = P + \frac{1}{2} \rho v^2 + \rho g h}
#'
#' Cette fonction permet d'analyser la contribution relative de chaque terme.
#'
#' @export
#'
#' @examples
#' # Exemple 1: Écoulement d'eau
#' termes <- bernoulli_terms(
#'   P = 101325,   # pression atmosphérique
#'   v = 5,        # vitesse 5 m/s
#'   h = 10,       # hauteur 10 m
#'   rho = 1000    # eau
#' )
#'
#' # Afficher les résultats
#' str(termes)
#'
#' # Pourcentage de chaque terme
#' total <- termes$total_energy
#' pourcent_pression <- termes$pressure_term / total * 100
#' pourcent_cinetique <- termes$velocity_term / total * 100
#' pourcent_potentiel <- termes$elevation_term / total * 100
#'
#' cat(sprintf("Pression: %.1f%%\n", pourcent_pression))
#' cat(sprintf("Cinétique: %.1f%%\n", pourcent_cinetique))
#' cat(sprintf("Potentiel: %.1f%%\n", pourcent_potentiel))
#'
#' # Exemple 2: Analyse pour différents points d'un écoulement
#' point1 <- bernoulli_terms(P = 200000, v = 2, h = 5, rho = 1000)
#' point2 <- bernoulli_terms(P = 150000, v = 8, h = 2, rho = 1000)
#'
#' # Vérifier la conservation (sans pertes)
#' diff_energie <- point1$total_energy - point2$total_energy
#' cat(sprintf("Différence d'énergie: %.2f Pa\n", diff_energie))
#'
#' @seealso
#' \code{\link{bernoulli_standard}} pour résoudre l'équation complète.

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
  cat("=== TESTS COMPLETS DE L'ÉQUATION DE BERNOULLI DE BASE ===\n\n")

  # Test 1: Cas simple
  cat("Test 1 - Cas simple (P2):\n")
  res1 <- bernoulli_standard(P1 = 100000, v1 = 1, h1 = 0,
                             v2 = 2, h2 = 0, rho = 1000, solve_for = "P2")
  cat("  P2 =", round(res1, 2), "Pa\n")

  # Test 2: Avec hauteur
  cat("\nTest 2 - Avec différence de hauteur (v2):\n")
  res2 <- bernoulli_standard(P1 = 100000, v1 = 0, h1 = 10,
                             P2 = 100000, h2 = 0, rho = 1000, solve_for = "v2")
  cat("  v2 =", round(res2, 2), "m/s\n")

  # Test 3: Détail des termes
  cat("\nTest 3 - Détail des termes:\n")
  terms <- bernoulli_terms(P = 101325, v = 5, h = 10, rho = 1000)
  cat("  Pression:", round(terms$pressure_term/1000, 2), "kPa\n")
  cat("  Cinétique:", round(terms$velocity_term/1000, 2), "kPa\n")
  cat("  Hauteur:", round(terms$elevation_term/1000, 2), "kPa\n")
  cat("  Total:", round(terms$total_energy/1000, 2), "kPa\n")

  # Test 4: Conservation d'énergie
  cat("\nTest 4 - Conservation d'énergie:\n")
  P2_test <- bernoulli_standard(P1 = 150000, v1 = 2, h1 = 5,
                                v2 = 4, h2 = 3, rho = 1000, solve_for = "P2")
  E1 <- 150000 + 0.5*1000*2^2 + 1000*9.81*5
  E2 <- P2_test + 0.5*1000*4^2 + 1000*9.81*3
  cat("  Énergie point 1:", round(E1, 2), "J/m³\n")
  cat("  Énergie point 2:", round(E2, 2), "J/m³\n")
  cat("  Différence:", round(abs(E1-E2), 2), "J/m³\n")
  cat("  Test réussi:", abs(E1-E2) < 0.001, "\n")
}

# Exécuter tous les tests
test_bernoulli()
