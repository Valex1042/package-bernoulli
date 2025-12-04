#' Visualisation des Termes de Bernoulli
#'
#' Crée un graphique à barres montrant la contribution relative de chaque terme
#' de l'équation de Bernoulli et l'énergie totale.
#'
#' @param P Pression statique (Pascals, Pa).
#' @param v Vitesse d'écoulement (mètres par seconde, m/s).
#' @param h Hauteur (élévation) par rapport à une référence (mètres, m).
#' @param rho Masse volumique du fluide (kilogrammes par mètre cube, kg/m³).
#' @param g Accélération gravitationnelle (mètres par seconde carrée, m/s²).
#'   Valeur par défaut: 9.81.
#'
#' @return Un objet ggplot représentant les contributions des termes de Bernoulli.
#'   Le graphique montre sous forme de barres:
#'   \itemize{
#'     \item Le terme de pression (P)
#'     \item Le terme d'énergie cinétique (0.5 * ρ * v²)
#'     \item Le terme d'énergie potentielle (ρ * g * h)
#'     \item L'énergie totale (somme des trois termes)
#'   }
#'
#' @details
#' Cette fonction visualise l'équation de Bernoulli:
#' \deqn{E_{tot} = P + \frac{1}{2} \rho v^2 + \rho g h}
#'
#' Le graphique permet de:
#' \enumerate{
#'   \item Comparer l'importance relative de chaque terme
#'   \item Identifier le terme dominant dans un écoulement donné
#'   \item Visualiser l'énergie totale par unité de volume
#' }
#'
#' @section Avertissements:
#' - Nécessite le package `ggplot2` (installation automatique si manquant)
#' - Les valeurs sont converties en kPa pour une meilleure lisibilité
#' - Le terme "Total" est la somme algébrique des trois termes
#'
#' @export
#'
#' @examples
#' # Exemple 1: Écoulement typique d'eau
#' \dontrun{
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   plot_bernoulli_terms(
#'     P = 101325,   # pression atmosphérique
#'     v = 5,        # vitesse modérée
#'     h = 10,       # chute de 10 m
#'     rho = 1000    # eau
#'   )
#' }
#' }
#'
#' # Exemple 2: Analyse pour différents scénarios
#' \dontrun{
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   # Cas 1: Haute pression, basse vitesse (conduite forcée)
#'   p1 <- plot_bernoulli_terms(P = 500000, v = 1, h = 0, rho = 1000) +
#'     ggplot2::ggtitle("Conduite forcée - Haute pression")
#'
#'   # Cas 2: Basse pression, haute vitesse (jet libre)
#'   p2 <- plot_bernoulli_terms(P = 101325, v = 20, h = 0, rho = 1000) +
#'     ggplot2::ggtitle("Jet libre - Haute vitesse")
#'
#'   # Cas 3: Hauteur dominante (chute d'eau)
#'   p3 <- plot_bernoulli_terms(P = 101325, v = 1, h = 50, rho = 1000) +
#'     ggplot2::ggtitle("Chute d'eau - Hauteur dominante")
#'
#'   # Afficher les trois graphiques côte à côte
#'   if (requireNamespace("gridExtra", quietly = TRUE)) {
#'     gridExtra::grid.arrange(p1, p2, p3, ncol = 3)
#'   }
#' }
#' }
#'
#' # Exemple 3: Comparaison air/eau
#' \dontrun{
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   # Eau
#'   p_eau <- plot_bernoulli_terms(P = 101325, v = 10, h = 0, rho = 1000) +
#'     ggplot2::ggtitle("Eau (ρ=1000 kg/m³)")
#'
#'   # Air
#'   p_air <- plot_bernoulli_terms(P = 101325, v = 10, h = 0, rho = 1.225) +
#'     ggplot2::ggtitle("Air (ρ=1.225 kg/m³)")
#'
#'   if (requireNamespace("gridExtra", quietly = TRUE)) {
#'     gridExtra::grid.arrange(p_eau, p_air, ncol = 2)
#'   }
#' }
#' }
#'
#' @seealso
#' \code{\link{bernoulli_terms}} pour obtenir les valeurs numériques des termes.
#'
#' @importFrom ggplot2 ggplot aes geom_col labs theme_minimal
#' @importFrom graphics plot
plot_bernoulli_terms <- function(P, v, h, rho, g = 9.81) {
  # Vérification de ggplot2
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    warning("Le package ggplot2 n'est pas installé. Installation en cours...")
    install.packages("ggplot2", quiet = TRUE)
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
      stop("Échec de l'installation de ggplot2. Installez-le manuellement avec install.packages('ggplot2')")
    }
  }

  # Validation des entrées
  if (missing(P) || missing(v) || missing(h) || missing(rho)) {
    stop("Tous les paramètres P, v, h, rho doivent être spécifiés")
  }

  if (rho <= 0) {
    stop("La masse volumique rho doit être positive")
  }

  if (g <= 0) {
    stop("L'accélération gravitationnelle g doit être positive")
  }

  # Calcul des termes
  terms <- bernoulli_terms(P, v, h, rho, g)

  # Préparation des données
  data <- data.frame(
    Terme = factor(c("Pression", "Vitesse", "Élévation", "Total"),
                   levels = c("Pression", "Vitesse", "Élévation", "Total")),
    Valeur = c(terms$pressure_term, terms$velocity_term,
               terms$elevation_term, terms$total_energy) / 1000,  # Conversion en kPa
    Type = c("Composante", "Composante", "Composante", "Total")
  )

  # Création du graphique
  plot <- ggplot2::ggplot(data, ggplot2::aes(x = Terme, y = Valeur, fill = Type)) +
    ggplot2::geom_col(width = 0.7, color = "black", linewidth = 0.5) +
    ggplot2::labs(
      title = "Contribution des Termes de Bernoulli",
      subtitle = sprintf("P = %.0f Pa, v = %.1f m/s, h = %.1f m, ρ = %.1f kg/m³",
                         P, v, h, rho),
      x = "Terme",
      y = "Énergie spécifique (kPa)",
      fill = "Type"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = ggplot2::element_text(hjust = 0.5),
      axis.text.x = ggplot2::element_text(angle = 0, vjust = 0.5),
      legend.position = "top"
    ) +
    ggplot2::scale_fill_manual(values = c("Composante" = "steelblue", "Total" = "darkred"))

  # Ajout des valeurs sur les barres
  plot <- plot +
    ggplot2::geom_text(
      ggplot2::aes(label = sprintf("%.1f", Valeur)),
      vjust = -0.5,
      size = 3.5,
      fontface = "bold"
    )

  return(plot)
}

# Exemple de test avec vos paramètres (commenté pour éviter l'exécution automatique)
 plot_bernoulli_terms(P = 1000000, v = 20, h = 100, rho = 982.5)
