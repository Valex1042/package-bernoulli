#' @title Visualization of Bernoulli Terms
#'
#' @description
#' Creates a bar chart showing the relative contribution of each term
#' of the Bernoulli equation and the total energy.
#'
#' @param P Static pressure (Pascals, Pa).
#' @param v Flow velocity (meters per second, m/s).
#' @param h Elevation relative to a reference (meters, m).
#' @param rho Fluid density (kilograms per cubic meter, kg/m3).
#' @param g Gravitational acceleration (meters per second squared, m/s2).
#'   Default value: 9.81.
#'
#' @return A ggplot object representing the contributions of Bernoulli terms.
#'   The chart shows as bars:
#'   \itemize{
#'     \item Pressure term (P)
#'    \item Kinetic energy term (0.5 * rho * v2)
#'     \item Potential energy term (rho * g * h)
#'     \item Total energy (sum of the three terms)
#'   }
#'
#' @details
#' This function visualizes the Bernoulli equation:
#' \deqn{E_{tot} = P + \frac{1}{2} \rho v^2 + \rho g h}
#
#' The chart allows you to:
#' \enumerate{
#'   \item Compare the relative importance of each term
#'   \item Identify the dominant term in a given flow
#'   \item Visualize the total energy per unit volume
#' }
#'
#' @section Warnings:
#' - Requires the `ggplot2` package (automatic installation if missing)
#' - Values are converted to kPa for better readability
#' - The "Total" term is the algebraic sum of the three terms
#'
#' @export
#'
#' @examples
#' # Example 1: Typical water flow
#' \dontrun{
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   plot_bernoulli_terms(
#'     P = 101325,   # atmospheric pressure
#'     v = 5,        # moderate velocity
#'     h = 10,       # 10 m drop
#'     rho = 1000    # water
#'   )
#' }
#' }
#'
#' # Example 2: Analysis for different scenarios
#' \dontrun{
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   # Case 1: High pressure, low velocity (forced pipe)
#'   p1 <- plot_bernoulli_terms(P = 500000, v = 1, h = 0, rho = 1000) +
#'     ggplot2::ggtitle("Forced pipe - High pressure")
#'
#'   # Case 2: Low pressure, high velocity (free jet)
#'   p2 <- plot_bernoulli_terms(P = 101325, v = 20, h = 0, rho = 1000) +
#'     ggplot2::ggtitle("Free jet - High velocity")
#'
#'   #Case 3: Dominant height (waterfall)
#'   p3 <- plot_bernoulli_terms(P = 101325, v = 1, h = 50, rho = 1000) +
#'     ggplot2::ggtitle("Waterfall - Dominant height")
#'
#'    Display the three charts side by side
#'   if (requireNamespace("gridExtra", quietly = TRUE)) {
#'     gridExtra::grid.arrange(p1, p2, p3, ncol = 3)
#'   }
#' }
#' }
#'
#' Example 3: Comparison air/water
#' \dontrun{
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'  Water
#'  p_water <- plot_bernoulli_terms(P = 101325, v = 10, h = 0, rho = 1000) +
#'     ggplot2::ggtitle("Water (rho=1000 kg/m3)")
#'
#'  # Air
#'   p_air <- plot_bernoulli_terms(P = 101325, v = 10, h = 0, rho = 1.225) +
#'     ggplot2::ggtitle("Air (rho=1.225 kg/m3)")
#'
#'  if (requireNamespace("gridExtra", quietly = TRUE)) {
#'     gridExtra::grid.arrange(p_water, p_air, ncol = 2)
#'   }
#' }
#' }
#'
#' @seealso
#' \code{\link{bernoulli_terms}} to obtain the numerical values of the terms.
#'
#'
plot_bernoulli_terms <- function(P, v, h, rho, g = 9.81) {
  # Verification de ggplot2
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    warning("Le package ggplot2 nest pas installe. Installation en cours...")
    install.packages("ggplot2", quiet = TRUE)
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
      stop("Echec de linstallation de ggplot2. Installez-le manuellement avec install.packages(ggplot2)")
    }
  }

  # Validation des entrees
  if (missing(P) || missing(v) || missing(h) || missing(rho)) {
    stop("Tous les parametres P, v, h, rho doivent etre specifies")
  }

  if (rho <= 0) {
    stop("La masse volumique rho doit etre positive")
  }

  if (g <= 0) {
    stop("Lacceleration gravitationnelle g doit etre positive")
  }

  # Calcul des termes
  terms <- bernoulli_terms(P, v, h, rho, g)

  # Preparation des donnees
  data <- data.frame(
    Terme = factor(c("Pression", "Vitesse", "Elevation", "Total"),
                   levels = c("Pression", "Vitesse", "Elevation", "Total")),
    Valeur = c(terms$pressure_term, terms$velocity_term,
               terms$elevation_term, terms$total_energy) / 1000,  # Conversion en kPa
    Type = c("Composante", "Composante", "Composante", "Total")
  )

  # Creation du graphique
  plot <- ggplot2::ggplot(data, ggplot2::aes(x = data$Terme, y = data$Valeur, fill = Type)) +
    ggplot2::geom_col(width = 0.7, color = "black", linewidth = 0.5) +
    ggplot2::labs(
      title = "Contribution des Termes de Bernoulli",
      subtitle = sprintf("P = %.0f Pa, v = %.1f m/s, h = %.1f m, ? = %.1f kg/m3",
                         P, v, h, rho),
      x = "Terme",
      y = "Energie specifique (kPa)",
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


plot_bernoulli_terms(P =10090, v= 23, h=100, rho= 1000)
