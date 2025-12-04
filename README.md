---
title: "README"
author: "djiogap"
date: "02/12/2025"
output: html_document
---
# bernoulli : Package R pour l'Équation de Bernoulli

[![R-CMD-check](https://github.com/votre-username/bernoulli/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/votre-username/bernoulli/actions/workflows/R-CMD-check.yaml)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![CRAN Status](https://www.r-pkg.org/badges/version-ago/bernoulli)](https://cran.r-project.org/package=bernoulli)
[![Downloads](https://cranlogs.r-pkg.org/badges/bernoulli)](https://cran.r-project.org/package=bernoulli)

<img src="https://upload.wikimedia.org/wikipedia/commons/thumb/6/6c/BernoullisLawDerivationDiagram.svg/800px-BernoullisLawDerivationDiagram.svg.png" alt="Équation de Bernoulli" width="400" align="right">

Un package R pour implémenter l'équation de Bernoulli et ses applications en mécanique des fluides.

## 📋 Description

Le package `bernoulli` fournit des outils pour résoudre l'équation de Bernoulli et ses applications classiques en mécanique des fluides pour les fluides incompressibles. Il est conçu pour les étudiants, enseignants et ingénieurs en mécanique des fluides, hydraulique et aérodynamique.

**Équation de Bernoulli :**
\[
P_1 + \frac{1}{2}\rho v_1^2 + \rho g h_1 = P_2 + \frac{1}{2}\rho v_2^2 + \rho g h_2
\]

## ✨ Fonctionnalités

### 🔧 Fonctions principales
- **`bernoulli_standard()`** : Résout l'équation de Bernoulli complète
- **`velocity_torricelli()`** : Théorème de Torricelli (vitesse de sortie d'un réservoir)
- **`pressure_venturi()`** : Effet Venturi (différence de pression)
- **`velocity_pitot()`** : Tube de Pitot (mesure de vitesse)
- **`flow_rate_bernoulli()`** : Calcul de débit volumique
- **`bernoulli_terms()`** : Calcule les termes individuels de l'équation

### 📊 Visualisation
- **`plot_bernoulli_terms()`** : Graphique des contributions énergétiques
- **`plot_pressure_velocity()`** : Relation pression-vitesse

### 🛠️ Outils auxiliaires
- **`fluid_properties()`** : Propriétés des fluides communs (eau, air, huile...)
- **`validate_bernoulli()`** : Validation des hypothèses de Bernoulli

## 📥 Installation

### Depuis GitHub (développement)
```r
# Installer devtools si nécessaire
# install.packages("devtools")

# Installer depuis GitHub
devtools::install_github("votre-username/bernoulli")
