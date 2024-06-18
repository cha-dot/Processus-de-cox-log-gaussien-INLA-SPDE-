# Processus de cox log-gaussien : méthode INLA-SPDE
Les processus de cox log-gaussien permettent d'étudier la structure de points spatialisés. Ils sont ici utilisés à des fins de modélisation de la distribution spatiale d'une espèce.

## Table des matières

- [Installation des packages](#installation-des-packages)
- [Données](#données)
  - [Les observations](#les-observations)
  - [Le domaine d'étude](#le-domaine-detude)
  - [Les zones de prospections](#les-zones-de-prospection)
  - [Les covariables environnementales](#les-covariables-environnementales)
- [Modèle](#modèle)
  - [Représentation des données](#représentation-des-données)
  - [Triangulation de l'espace](#triangulation-de-lespace)
  - [Paramétrage de la matrice Matérn](#paramétrage-de-la-matrice-matérn)
  - [Cellules de Voronoi](#cellules-de-voronoi)
  - [Vecteurs des observations et de pondération](#vecteurs-des-observations-et-de-pondération)
  - [Matrice d'observation A](#matrice-dobservation-A)
  - [Covariables environnementales](#covariables-environnementales)
  - [Stack INLA](#stack-inla)
  - [Formule du modèle](#formule-du-modèle)
  - [Représentation spatiale des prédictions](#représentation-spatiale-des-prédictions)
  - [Qualité d'ajustement du modèle](#qualité-dajustement-du-modèle)
  - [Courbes de densités a posteriori des effets fixes](#courbes-de-densités-a-posteriori-des-effets-fixes)
  - [Calcul des probabilités a posteriori](#calcul-des-probabilités-a-posteriori)
  - [Courbes de densités a posteriori des hyperparamètres](#courbes-de-densités-a-posteriori-des-hyperparamètres)


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
Pour des variables temporelles, détenir des tableaux avec les valeurs par année (lignes) et par période de prospection (colonnes).

Il est possible d'avoir les données sous un autre format mais cela nécessitera de retravailler les fonctions qui appellent les valeurs des covariables.

## Modèle

### Représentation des données

```r
ggplot() +
  geom_sf(data = st_as_sf(GB_PE_eau_fauche_LAMB), col = "blue", size = 0.1) +
  geom_sf(data = st_as_sf(contour_sp), alpha = 0, col = "black") +
  geom_sf(data = st_as_sf(Buffer_sp), alpha = 0.5, col = "orange", linewidth = 0.5)
```


