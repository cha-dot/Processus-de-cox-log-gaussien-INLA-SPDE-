# Processus-de-cox-log-gaussien-INLA-SPDE-
Les processus de cox log-gaussien permettent d'étudier la structure de points spatialisés. Ils sont ici utilisés à des fins de modélisation de la distribution spatiale d'une espèce.

## Table des matières

- [Installation des packages](#installation-des-packages)
- [Format des données](#format-des-données)
- [Code](#code)
- [Interprétations](#interprétations)


## Installation des packages

install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)

Si vous faites face à des problèmes de compatibilité à cause de la version du logiciel R, il est possible de télécharger des versions plus anciennes de INLA sur le site : 

[Site R-INLA Project](https://www.r-inla.org/download-install) dans la section **Version Compatibility**.

Télécharger les fonctions utilisées dans le code : `source("inlabookfunctions.R")`
