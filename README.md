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
  - [Paramétrage de la matrice de Matérn](#paramétrage-de-la-matrice-de-matérn)
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

`install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)`

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
load("dataprocess.rda")

ggplot() +
  geom_sf(data = st_as_sf(GB_PE_eau_fauche_LAMB), col = "blue", size = 0.1) +
  geom_sf(data = st_as_sf(contour_sp), alpha = 0, col = "black") +
  geom_sf(data = st_as_sf(Buffer_sp), alpha = 0.5, col = "orange", linewidth = 0.5)
```

### Triangulation de l'espace

Le domaine est triangulé afin de simuler et d'estimer les champs spatiaux aléatoires du modèle.

```r
tmp = spsample(contour_sp@polygons[[1]], n = 1000, type = "regular") # tirage régulier de 1000 points dans le domaine
coordLIM = rbind(contour_sp@polygons[[1]]@Polygons[[1]]@coords, tmp@coords) # coordonnées du domaine et des 1000 points
rm(tmp) # suppression de tmp pour libérer de la mémoire

bndint = inla.nonconvex.hull(coordLIM, convex = -.02) # limite autour des points de coordLIM (dite interne)
bndext = inla.nonconvex.hull(coordLIM, convex = -.1) # limite de l'extension du domaine (dite externe)

par(mar=rep(1,4), mfrow=c(1,1))
plot(rbind(bndext$loc, bndext$loc[1,]), type = "l", lwd=2)
lines(rbind(bndint$loc, bndint$loc[1,]), pch = 19, cex = .05, col = "orange",
      lwd = 1, lty = 2)
plot(contour_sp, add = T, border = 4)
```

Le domaine est étendu afin de limiter les effets de bord. La triangulation est réalisée à partir des points tirés régulièrement dans l'espace.

```r
mesh = inla.mesh.2d(loc = rbind(coordinates(GB_PE_eau_fauche_LAMB)), # contrainte = faire passer la mesh par les observations
                    boundary = list(int = bndint, out = bndext), # limites int et ext
                    max.edge = c(110, 1000), cutoff = 40, # taille des arrêtes des triangles, max. 110 dans limite int et max. 1000 dans limite ext ; min. 40
                    crs = "+proj=lcc +lat_0=46.5 +lon_0=3 +lat_1=49 +lat_2=44 
                    +x_0=700000 +y_0=6600000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")

mesh$n # nombre de triangles

par(mar=rep(1,4), mfrow=c(1,1))
plot(mesh, main = "", asp = 1)
plot(GB_PE_eau_fauche_LAMB, add = T, col = 2, pch = 16, cex = .3)
plot(contour_sp, add = T, border = 1, lwd = 2)
```
Nous faisons passer la mesh par les observations (les observations constituent donc des sommets de triangles). Cela permet de détenir des informations précises puisque les données à l'intérieur des triangles seront interpolées.

### Paramétrage de la matrice de Matérn

Les champs gaussiens suivent une loi normale multivariée d'espérance nulle et de fonction de variance-covariance, une matrice de type Matérn. Ce type de matrice est flexible et paramétrique permettant de modéliser efficacement la structure de dépendance spatiale entre les observations. Cela permet notamment de supposer que l'autocorrélation  spatiale entre les points diminue avec la distance.<br>
Le champ spatial est caractérise par 3 hyperparamètres :
- La rugosité, fixée à 2
- La variance dont la loi a priori repose sur un PC-prior
- La portée, correspondant à la distance à partir de laquelle les points ne sont plus significativement autocorrélés. La valeur de la portée dans les paramètres est choisie empiriquement : elle correspond approximativement à une distance qui est 3 fois inférieure à la distance entre les 2 points les plus éloignés du domaine.

```r
matern = inla.spde2.pcmatern(mesh,
                             alpha = 2, # rugosité du champ
                             prior.sigma = c(1, 0.5), # variance, P(sigma > 1) = 0.5 (en mètres ici)
                             prior.range = c(3000, 0.9)) # portée, P(range < 3000) = 0.9
```

### Cellules de Voronoi
