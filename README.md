# Processus de cox log-gaussien : méthode INLA-SPDE
Les processus de cox log-gaussien permettent d'étudier la structure de points spatialisés. Ils sont ici utilisés à des fins de modélisation de la distribution spatiale d'une espèce.

## Table des matières

- [Installation des packages](#installation-des-packages)
- [Données](#Données)
- [Code](#code)
- [Interprétations](#interprétations)


## Installation des packages

install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)

Si vous faites face à des problèmes de compatibilité à cause de la version du logiciel R, il est possible de télécharger des versions plus anciennes de INLA sur le site : 

[Site R-INLA Project](https://www.r-inla.org/download-install) dans la section **Version Compatibility**.

Télécharger les fonctions pour ce code : `source("inlabookfunctions.R")`

## Données

### Les observations

**Jeu de données :**

- **Lignes :**
  - Observations (un individu par ligne)

- **Colonnes :**
  - Coordonnées spatiales
  - Autres variables : année, numéro de passage (possibilité d'en mettre d'autres)


### Le domaine d'étude

Détenir les limites spatiales du domaine d'étude.

### Les zones de prospection (si elles ne recouvrent pas l'entiereté du domaine d'étude)

Détenir les limites spatiales des zones prospectées (ex : disques d'écoute, transects etc.)

### Les covariables environnementales (facultatives)

Pour des variables spatiales, détenir des fichiers géospatiaux.<br>
Pour des variables temporelles, détenir des tableaux avec les valeurs par année et par période de prospection.

Il est possible d'avoir les données sous un autre format mais cela nécessitera de retravailler les fonctions qui appellent les valeurs des covariables.
