# Processus de cox log-gaussien : méthode INLA-SPDE
Les processus de cox log-gaussien permettent d'étudier la structure de points spatialisés. Ils sont ici utilisés à des fins de modélisation de la distribution spatiale d'une espèce.

## Table des matières

- [Installation des packages](#installation-des-packages)
- [Données](#données)
  - [Les observations](#les-observations)
  - [Le domaine d'étude](#le-domaine-détude)
  - [Les zones de prospections](#les-zones-de-prospection)
  - [Les covariables environnementales](#les-covariables-environnementales)
- [Modèle](#modèle)
  - [Représentation des données](#représentation-des-données)
  - [Triangulation de l'espace](#triangulation-de-lespace)
  - [Paramétrage de la matrice de Matérn](#paramétrage-de-la-matrice-de-matérn)
  - [Cellules de Voronoi](#cellules-de-voronoi)
  - [Vecteurs des observations et de pondération](#vecteurs-des-observations-et-de-pondération)
  - [Matrice d'observation A](#matrice-dobservation-a)
  - [Covariables environnementales](#covariables-environnementales)
  - [Stack INLA](#stack-inla)
  - [Formule du modèle](#formule-du-modèle)
  - [Représentation spatiale des prédictions](#représentation-spatiale-des-prédictions)
  - [Qualité d'ajustement du modèle](#qualité-dajustement-du-modèle)
  - [Courbes de densité a posteriori des effets fixes](#courbes-de-densité-a-posteriori-des-effets-fixes)
  - [Calcul des probabilités a posteriori](#calcul-des-probabilités-a-posteriori)
  - [Courbes de densité a posteriori des hyperparamètres](#courbes-de-densité-a-posteriori-des-hyperparamètres)
  - [Calcul des probabilités a posteriori des hyperparamètres](#calcul-des-probabilités-a-posteriori-des-hyperparamètres)


## Installation des packages

`install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)`

Si vous faites face à des problèmes de compatibilité à cause de la version du logiciel R, il est possible de télécharger des versions plus anciennes de INLA sur le site : 

[Site R-INLA Project](https://www.r-inla.org/download-install) dans la section **Version Compatibility**.

Charger les packages :

```r
library(INLA)
library(inlabru)
library(ggplot2)
library(raster)
library(sf)
library(sp)
library(tidyverse)
library(viridis)
library(hrbrthemes)
library(dplyr)
library(tidyr)
library(pROC)
```

Télécharger les fonctions utilisées dans ce code : `source("inlabookfunctions.R")`

## Données

### Les observations

**Jeu de données :**

- **Lignes :**
  - Observations (un individu par ligne)

- **Colonnes :**
  - Coordonnées spatiales
  - Autres variables : année, numéro de passage (possibilité d'en ajouter)


### Le domaine d'étude

Détenir les limites spatiales du domaine d'étude.

### Les zones de prospection

Si elles ne recouvrent pas l'entiereté du domaine, détenir les limites spatiales des zones prospectées (ex : disques d'écoute, transects etc.)

### Les covariables environnementales

Pour des variables spatiales, détenir des fichiers géospatiaux.<br>
Pour des variables temporelles, détenir des tableaux avec les valeurs par année (lignes) et par période de prospection (colonnes).

Il est possible d'avoir les données sous un autre format mais cela nécessitera de retravailler les fonctions qui appellent les valeurs des covariables.

L'ajout de covariables environnementales est facultatif.

## Modèle

### Représentation des données

```r
load("dataprocess.rda")

ggplot() +
  geom_sf(data = st_as_sf(GB_PE_eau_fauche_LAMB), col = "blue", size = 0.1) +
  geom_sf(data = st_as_sf(contour_sp), alpha = 0, col = "black") +
  geom_sf(data = st_as_sf(Buffer_sp), alpha = 0.5, col = "orange", linewidth = 0.5)
```

![Représentation des données](https://github.com/cha-dot/Processus-de-cox-log-gaussien-INLA-SPDE-/blob/images/carte_buffers.jpg?raw=true)

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

![Limites du domaine](https://github.com/cha-dot/Processus-de-cox-log-gaussien-INLA-SPDE-/blob/images/limites_complet.jpeg?raw=true)

Le domaine est étendu afin de limiter les effets de bord. La triangulation est réalisée à partir des points tirés régulièrement dans l'espace.

```r
mesh = inla.mesh.2d(loc = rbind(coordinates(GB_PE_eau_fauche_LAMB)), # contrainte (loc) = faire passer la mesh par les observations
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

![Triangulation](https://github.com/cha-dot/Processus-de-cox-log-gaussien-INLA-SPDE-/blob/images/mesh_complet.jpeg?raw=true)

Nous faisons passer la mesh par les observations (les observations constituent donc des sommets de triangles). Cela permet de détenir des informations précises puisque les données à l'intérieur des triangles seront interpolées.

### Paramétrage de la matrice de Matérn

Les champs gaussiens suivent une loi normale multivariée d'espérance nulle et de fonction de variance-covariance, une matrice de type Matérn. Ce type de matrice est flexible et paramétrique permettant de modéliser efficacement la structure de dépendance spatiale entre les observations (supposant notamment que l'autocorrélation  spatiale entre les points diminue avec la distance).<br>

Le champ spatial est caractérisé par 3 hyperparamètres :
- La rugosité, fixée à 2
- La variance dont la loi a priori repose sur un PC-prior
- La portée correspond à la distance à partir de laquelle les points ne sont plus significativement autocorrélés. La valeur de la portée dans les paramètres est choisie empiriquement : elle correspond approximativement à une distance qui est 3 fois inférieure à la distance entre les 2 points les plus éloignés du domaine.

```r
matern = inla.spde2.pcmatern(mesh,
                             alpha = 2, # rugosité du champ
                             prior.sigma = c(1, 0.5), # variance, P(sigma > 1) = 0.5 (en mètres ici)
                             prior.range = c(3000, 0.9)) # portée, P(range < 3000) = 0.9
```

### Cellules de Voronoi

L'intensité, c'est-à-dire le nombre d'individus observé par unité de surface, sera pondérée par la surface de prospection. Cela nécessite l'utilisation de cellules de Voronoi qui sont tracées autour de chacun des sommets des triangles de la première mesh.

```r
dmesh <- book.mesh.dual(mesh) # cellules de voronoi
plot(dmesh)
Buffer = st_as_sf(Buffer_sp) # conversion en format sf
dmesh_sf = st_as_sf(dmesh) # idem
st_crs(dmesh_sf) = st_crs(Buffer) # attribution des mêmes systèmes de référence
w <- sapply(1:length(dmesh), function(i) { # pour chaque cellule, calcul de l'aire d'intersection entre la cellule et les buffers
  if (nrow(st_intersection(dmesh_sf[i, ], Buffer)) > 0)
    return(st_area(st_intersection(dmesh_sf[i, ], Buffer)))
  else return(0)
})

wbis = unlist(lapply(w, max)) # attribution de l'aire d'intersection (w) max quand 2 buffers chevauchent une seule cellule
nv = mesh$n # nombre de cellules de voronoi
n = table(GB_PE_eau_fauche_LAMB$annee,GB_PE_eau_fauche_LAMB$num_passag) # nombre d'observations par année et par période

plot(dmesh,col=(wbis>0)*4) # cellules colorées = zones de prospection
plot(Buffer_sp, border="orange",lwd=2,add=T)
plot(GB_PE_eau_fauche_LAMB,pch=16,cex=0.6,col=2,add=T)
```

![Cellules de voronoi](https://github.com/cha-dot/Processus-de-cox-log-gaussien-INLA-SPDE-/blob/images/voronoi_complet.jpeg?raw=true)

La matrice d'observation A permet de faire le lien entre la mesh, les observations et les covariables. Nous allons construire manuellement cette matrice, ainsi que les vecteurs d'observations, de pondération et des covariables. Se référer à l'image ci-dessous pour une meilleure compréhension de l'élaboration des différentes parties.

<a id="matrice-a"></a>
![Matrice d'observation A](https://private-user-images.githubusercontent.com/173138382/340686454-5abee3ed-33f9-40f9-b746-c7e22400fb3e.PNG?jwt=eyJhbGciOiJIUzI1NiIsInR5cCI6IkpXVCJ9.eyJpc3MiOiJnaXRodWIuY29tIiwiYXVkIjoicmF3LmdpdGh1YnVzZXJjb250ZW50LmNvbSIsImtleSI6ImtleTUiLCJleHAiOjE3MTg3MTY0NjIsIm5iZiI6MTcxODcxNjE2MiwicGF0aCI6Ii8xNzMxMzgzODIvMzQwNjg2NDU0LTVhYmVlM2VkLTMzZjktNDBmOS1iNzQ2LWM3ZTIyNDAwZmIzZS5QTkc_WC1BbXotQWxnb3JpdGhtPUFXUzQtSE1BQy1TSEEyNTYmWC1BbXotQ3JlZGVudGlhbD1BS0lBVkNPRFlMU0E1M1BRSzRaQSUyRjIwMjQwNjE4JTJGdXMtZWFzdC0xJTJGczMlMkZhd3M0X3JlcXVlc3QmWC1BbXotRGF0ZT0yMDI0MDYxOFQxMzA5MjJaJlgtQW16LUV4cGlyZXM9MzAwJlgtQW16LVNpZ25hdHVyZT1kZjcyOGEyMjM4OTIxYzc3NzAzMDYzZTgwYmZkMzlhZjI2ZTkzYjE1ODZmOGJjM2I2OTE5OGQ3MjRlNTQ5MTQ0JlgtQW16LVNpZ25lZEhlYWRlcnM9aG9zdCZhY3Rvcl9pZD0wJmtleV9pZD0wJnJlcG9faWQ9MCJ9.ekKicQE8Z83QZ2M8rMqB-0TWtb1cFpvi4YmUFmLjgsw)

### Vecteurs des observations et de pondération

Le vecteur y.pp est composé de 0 et de 1, c'est le vecteur des observations. Le vecteur e.pp est composé des poids de chaque cellule de voronoi (c'est-à-dire la surface d'intersection entre les surfaces prospectées et les cellules de voronoi).

```r
y.pp = NULL
e.pp = NULL
for (i in 1:nrow(n)){ # pour chaque année
  for (j in 1:ncol(n)){ # pour chaque période de prospection
    y.pp = c(y.pp,rep(0:1, c(nv, n[i,j]))) # vecteur des observations (permettra de relier les observations à la mesh)
    e.pp = c(e.pp,c(wbis, rep(0, n[i,j]))) # vecteur des poids (permettra de pondérer l'intensité par la surface de prospection)
  }
}
```

### Matrice d'observation A

La [Matrice d'observation A](#matrice-a) est divisée en matrices I (représentant les cellules de voronoi), et en matrices L (représentant les observations d'une période et d'une année).

```r
imat <- Diagonal(nv, rep(1, nv)) # matrice I diagonale de dimensions nv (nombre de cellules de Voronoi), représente la mesh
lmat <- inla.spde.make.A(mesh, coordinates(GB_PE_eau_fauche_LAMB)) # matrice L qui représente les observations dans l'espace
A.pp = NULL
for (i in 2018:2023){
  for (j in 1:2){
    indice = (1:nrow(GB_PE_eau_fauche_LAMB))[(GB_PE_eau_fauche_LAMB$annee==i)&(GB_PE_eau_fauche_LAMB$num_passag==j)]
    A.pp = rbind(A.pp,rbind(imat,lmat[indice,])) # Matrice d'observation A différente pour chaque année (A11, A12...)
    # on ajoute les matrices I et L les unes à la suite des autres
  }
}
```

### Covariables environnementales

Chargement des données environnementales

```r
veg2 = raster("Rgpt_vg_3.tif") # raster végétation
vegSG2 = as(veg2, "SpatialGridDataFrame") # conversion du format

setwd("Fauche") # spécifier le chemin pour l'ouverture des fichiers plus tard
# Vecteur qui contient le nom des fichiers qui seront appelés dans la fonction
shapefiles = c("Fauche_2017_seule.shp", "Fauche_2018_seule.shp", "Fauche_2019_seule.shp", "Fauche_2020_seule.shp", "Fauche_2021_seule.shp", "Fauche_2022_seule.shp")

eau = read.csv("Hauteurs_eau_finales.csv", header = T, stringsAsFactors = T) # tableau des hauteurs d'eau
```
A présent, on détermine la fonction qui récupère la valeur de la végétation en chaque point du domaine d'étude. La fonction prend en argument les coordonnées spatiales, ainsi que les fichiers spatiaux de la fauche (pour séparer la végétation fauchée de celle qui ne l'est pas). 

```r
f.veg.moins = function(x, y, shape) {
  spp <- SpatialPoints(data.frame(x = x, y = y), proj4string = fm_sp_get_crs(vegSG2)) # création de points spatiaux à partir des coordonnées qui seront fournis
  proj4string(spp) <- fm_sp_get_crs(vegSG2) # attribution d'un système de référence
  v <- over(spp, vegSG2) # la valeur de végétation attribuée à chaque point
  v2 <- over(spp, shape) # la valeur de fauche attribuée à chaque point
  v$Rgpt_vg_3[is.na(v$Rgpt_vg_3)]=0 # si y a des NA, on transforme en 0
  v2$id_type_ge[is.na(v2$id_type_ge)]=0 # idem
  v$Rgpt_vg_3 = ifelse(v$Rgpt_vg_3 %in% c(13,1,11,8) & v2$id_type_ge == 1, 14, v$Rgpt_vg_3) # si la veg est herbacée haute (13, 1, 11, 8) et qu'elle est fauchée (1), alors on la met dans un groupe (14)
  v$Rgpt_vg_3 = ifelse(v$Rgpt_vg_3 %in% c(13,1,11,8) & v2$id_type_ge == 0, 15, v$Rgpt_vg_3) # si elle n'est pas fauchée (0), on la met dans un autre groupe (15)
  return(v$Rgpt_vg_3)
}
```
On applique la fonction sur les données pour chaque année et période de prospection. Dans la fonction, on renseigne les coordonnées spatiales de la mesh et des observations de l'année i et la période j, et le fichier spatial de la fauche à l'année i-1 (choix d'étude).

```r
vegTMP=NULL
for (i in 2018:2023){
  for (j in 1:2){
    shapefile = st_read(shapefiles[i - 2017]) # i-2017 pour ramener les indices de 1 à 6 : appelle les fichiers de fauche
    shapefile = as(shapefile, "Spatial") # conversion du format
    crs(shapefile) = "+proj=lcc +lat_0=46.5 +lon_0=3 +lat_1=49 +lat_2=44 +x_0=700000 +y_0=6600000 +ellps=GRS80
+towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
    indice = (1:nrow(GB_PE_eau_fauche_LAMB))[(GB_PE_eau_fauche_LAMB$annee==i)&(GB_PE_eau_fauche_LAMB$num_passag==j)] # indice qui détermine l'année et la période pour avoir les observations associées
    vegTMP = c(vegTMP,
               f.veg.moins(x=c(mesh$loc[,1], coordinates(GB_PE_eau_fauche_LAMB)[indice,1]),
                           y = c(mesh$loc[,2], coordinates(GB_PE_eau_fauche_LAMB)[indice,2]),
                           shape = shapefile)) # on applique la fonction
  }
}
```
Dans le cadre de notre étude, les groupements végétaux sont remaniés.

```r
veg = vegTMP
veg[vegTMP == 6 | vegTMP == 12 | vegTMP == 4 | vegTMP == 7| vegTMP == 9 | vegTMP == 10] = 1 # rase
veg[vegTMP == 14] = 2 # haute fauchee
veg[vegTMP == 15] = 3 # haute non fauchee
veg[vegTMP == 2] = 4 # arbustive
veg[vegTMP == 3] = 5 # roselières 
veg[vegTMP == 5] = 6 # friches
veg = as.factor(veg)
```

A présent, nous allons définir la fonction qui permettra d'associer les valeurs d'intensité de submersion à chaque observation. Il s'agit ici d'une variable temporelle, une même valeur sera donc associée à toutes les cellules de voronoi et toutes les observations d'une même période de prospection (voir la [Matrice d'observation A](#matrice-a)).

```r
f.eau.max = function(annee, prospection) {
  data = eau[eau$annee == annee,] # sélectionne la ligne correspondant à la bonne année dans le tableau des hauteurs d'eau
  if(prospection == 1) {return(data[,1])} else{return(data[,2])} # choix de la colonne où il y a la valeur de la prospection correspondante
}
```
On applique la fonction.

```r
eau_max=NULL
for (i in 2018:2023){ # pour chaque année
  for (j in 1:2){ # pour chaque période
    eau_max_1 = f.eau.max(annee = i, prospection = j)
    eau_max = c(eau_max, eau_max_1) # vecteur qui contient chaque valeur d'intensité de submersion
  }
}
eau_max = matrix(data = eau_max, ncol = 2, byrow = T)
```

On crée le vecteur de la covariable, de façon à pouvoir le relier à la [Matrice d'observation A](#matrice-a). C'est à cette étape que l'on associe une valeur à chaque cellule de Voronoi et à chaque observation.

```r
eau_maxTMP = NULL
for(i in 2018:2023) {
  for(j in 1:2) {
    eau_max_I = rep(eau_max[i-2017,j], nv) # on répète parce qu'il y a qu'une seule valeur par période pour toutes les positions
    eau_max_A = rep(eau_max[i-2017,j], n[i-2017,j]) # idem
    eau_maxTMP = c(eau_maxTMP, eau_max_I, eau_max_A) # on colle les vecteurs I et A au fur et à mesure
  }
}
```

De la même manière, on définit une fonction qui attribuera une valeur de durée de submersion (variable temporelle) pour toutes les cellules de voronoi et toutes les observations.

```r
f.eau.duree = function(annee, prospection) {
  data = eau[eau$annee == annee,]
  if(prospection == 1) {return(data[,3])} else{return(data[,5])}
}

eau_duree=NULL
for (i in 2018:2023){
  for (j in 1:2){
    eau_duree_1 = f.eau.duree(annee = i, prospection = j)
    eau_duree = c(eau_duree, eau_duree_1)
  }
}
eau_duree = matrix(data = eau_duree, ncol = 2, byrow = T)

eau_dureeTMP = NULL
for(i in 2018:2023) {
  for(j in 1:2) {
    eau_duree_I = rep(eau_duree[i-2017,j], nv)
    eau_duree_A = rep(eau_duree[i-2017,j], n[i-2017,j])
    eau_dureeTMP = c(eau_dureeTMP, eau_duree_I, eau_duree_A)
  }
}
```

### Stack INLA

Le stack INLA constitue un stock des matrices et vecteurs créés juste avant. Cela permet de combiner les différentes sources de données et les différents effets dans une seule entité.

```r
stk.pp <- inla.stack(                   
  data = list(y = y.pp, e = e.pp), # vecteurs d'observation et de pondération
  A = list(1, A.pp), # matrice d'observation
  effects = list(list(veg = veg, max_sub = eau_maxTMP, duree_sub = eau_dureeTMP), # covariables ; b0 = 1 quand y a juste le champ spatial
                 list(i = 1:nv)), # effet champ spatial
  tag = 'pp')
```

### Formule du modèle

Le modèle est composé d'effets fixes (les covariables environnementales) et d'un effet aléatoire (le champ spatial gaussien).<br>

Une première stratégie était de retirer l'intercept du modèle. Celui-ci serait donc compris dans chacune des modalités de la végétation (et pas uniquement confondu avec l'effet de la première modalité), permettant la comparaison de chaque effet.
Mais pour contourner des problèmes de corrélation entre les lois a posteriori des modalités de la variable "végétation", cette dernière est ajoutée comme une variable aléatoire au modèle en contraignant la somme des paramètres à 0 `constr = T`.

```r
ppVIV <- inla(y ~ 1 + f(i, model = matern) + max_sub + duree_sub  + f(veg, model = "iid", constr = T),  # modèle indépendant et identiquement distribué
              family = 'poisson', data = inla.stack.data(stk.pp),
              control.predictor = list(A = inla.stack.A(stk.pp)), # matrice d'observation extraite du stack
              E = inla.stack.data(stk.pp)$e, # vecteur de pondération extrait du stack
              control.compute=list(dic=TRUE, cpo=TRUE, config = TRUE, return.marginals.predictor=TRUE), # retourne les marginales des prédictions
              control.inla=list(strategy="simplified.laplace",int.strategy="eb")) # approximation de Laplace simplifiée ; eb = stratégie Empirical Bayes

summary(ppVIV)
```

### Représentation spatiale des prédictions

Afin de représenter spatialement les prédictions, il est nécessaire de créer une grille de pixels qui couvre la mesh (et donc le site d'étude). Un projecteur projettera donc les résultats du modèle sur la grille de pixels définie.

```r
toto = extent(contour_sp) # étendue totale des limites de la réserve (coordonnées minimales et maximales en x et y)
xrange = toto[2]-toto[1] # on fait ces manips pour avoir des pixels carrés et pas rectangulaires par la suite
yrange = toto[4]-toto[3]
resopred=300 # nombre de pixels (modulable)
st_crs(contour) = NA
inla.identical.CRS(mesh, contour)
contour = st_set_crs(contour, inla.CRS(mesh)) # mêmes systèmes de référence
pxl <- pixels(mesh, mask=contour, nx = resopred, # nx = nombre de pixels longitudinaux
              ny = round(resopred*yrange/xrange)) # ny = nombre de pixels latitudinaux
projgrid <- inla.mesh.projector(mesh, coordinates(pxl)) # projecteur
```

On stocke les prédictions de chaque observation dans un vecteur `predtot`. Cela nécessite de récupérer les médianes dans `ppVIV$summary.linear.predictor` pour chaque période et année, nécessitant l'utilisation d'indices `indN` et `indNV`.

```r
indN = 0 # indice
indNV = 1 # indice
predtot = NULL
for(i in 1:6){
  for(j in 1:2){
    predtot = rbind(predtot, ppVIV$summary.linear.predictor$'0.5quant'[((indNV-1)*nv+indN+1):(indNV*nv+indN)]) # toutes les prédictions
    indN = indN + n[i,j]
    indNV = indNV + 1
    
  }
}
```

Nous pouvons à présent représenter les prédictions spatiales (des prédicteurs linéaires ici). Toutes les prédictions sont sommées.

```
Wmean <- inla.mesh.project(projgrid, apply(predtot,2,sum)) # projette spatialement les valeurs des prédicteurs linéaires
ggplot()+gg(as(contour, "Spatial"))+gg(pxl, aes(fill=Wmean))+gg(GB_PE_eau_fauche_LAMB,size=0.5)
```

![Prédicteurs linéaires](https://github.com/cha-dot/Processus-de-cox-log-gaussien-INLA-SPDE-/blob/images/predicteurs_lineaires.jpg?raw=true)

Nous allons à présent représenter spatialement les densités prédites des individus de l'espèce. Nous récupérons toutes les prédictions de l'intensité (nombre d'individus moyen par unité de surface).

```r
indN = 0
indNV = 1
predtot = NULL
for(i in 1:6){
  for(j in 1:2){
    predtot = rbind(predtot, ppVIV$summary.fitted.values$'0.5quant'[((indNV-1)*nv+indN+1):(indNV*nv+indN)]) # donc fitted.values
    indN = indN + n[i,j]
    indNV = indNV + 1
    
  }
}

fitval = inla.mesh.project(projgrid,
                           apply(predtot,2,sum)) # projecteur des intensités prédites
```

Pour représenter la densité des individus, nous multiplions l'intensité par la surface du domaine.

```r
ggplot()+gg(as(contour, "Spatial"))+gg(pxl, aes(fill=fitval*as.numeric(st_area(rSF)))) + scale_fill_viridis(option = "H") + labs(fill = "Nombre d'individus/pixel") # espérance du nb de pts
```

![Carte densités individus](https://github.com/cha-dot/Processus-de-cox-log-gaussien-INLA-SPDE-/blob/images/derniere_carte_dalto.jpg?raw=true)

Afin d'évaluer la qualité d'ajustement du modèle, nous stockons les observations/pixel dans un vecteur.

```r
r = as(pxl, "SpatialPolygonsDataFrame") # conversion du format
rSF=st_as_sf(r) # idem
datvivSF = st_as_sf(GB_PE_eau_fauche_LAMB) # idem
st_crs(rSF) = st_crs(datvivSF) # homogénéisation des crs
test=st_intersects(rSF,datvivSF) # intersection entre les observations et chaque pixel
NobsVIV=unlist(lapply(test,length)) # nb d'observations par pixel
```

Il est alors possible d'explorer graphiquement l'ajustement des prédictions avec les observations. L'alignement des points sur la droite traduirait un ajustement parfait des prédicitions. Des points au-dessus de la droite correspondrait à un sur-ajustement, alors que des points en-dessous de la droite exprimeraient un sous-ajustement.

```r
plot(NobsVIV,fitval*st_area(rSF)) 
abline(a=0,b=1,col=4,lwd=2,lty=2) # graphe pour voir l'ajustement des pred
```

![Ajustement](https://github.com/cha-dot/Processus-de-cox-log-gaussien-INLA-SPDE-/blob/images/ajustement.png?raw=true)

### Qualité d'ajustement du modèle

L'AUC (Area Under the Curve), l'aire sous la courbe ROC (Receiver Operating Characteristic), est une métrique permettant d'évaluer la qualité d'ajustement d'un modèle. Cela repose sur la comparaison de paires de données (observations/prédictions). Il nous faut transformer les observations en données de présence/absence (données binaires), et les prédictions en probabilité de présence (données continues).

```r
NobsVIV_bis = as.numeric(NobsVIV > 0) # les pixels où il y a présence (observations)
predTEST = 1-exp(-as.numeric(fitval*st_area(rSF))) # estimation de la probabilité de présence
testroc=roc(NobsVIV_bis,predTEST)
testroc # AUC
plot.roc(NobsVIV_bis, predTEST, col="red", print.thres="best") # courbe roc
```



### Courbes de densité a posteriori des effets fixes

L'estimation bayésienne des paramètres permet d'analyser leur distribution a posteriori. Nous procédons alors à un échantillonnage du modèle.

```r
Nrep=2000 # nombre d'échantillons
pr.int.tot <- inla.posterior.sample(Nrep,ppVIV) # échantillonnage dans le modèle
ind=grep("veg",rownames(pr.int.tot[[1]]$latent)) # lignes de la composante latente qui ont les paramètres vg
m = grep("max_sub", rownames(pr.int.tot[[1]]$latent)) # lignes de la composante latente qui ont les paramètres max_sub
d = grep("duree_sub", rownames(pr.int.tot[[1]]$latent)) # lignes de la composante latente qui ont les paramètres duree_sub
post=matrix(unlist(lapply(pr.int.tot,function(x){x$latent[c(ind,m,d),]})),nrow=Nrep,byrow=T) # Nrep valeurs des paramètres pour chaque variable
post = as.data.frame(post)
colnames(post) = c("cultures", "rase", "haute fauchée", "haute non fauchée", "arbustive", "roselières/scirpaies", "friches", "intensite_sub", "duree_sub")
```

Les paramètres de chaque variable sont convertis à un format compatible avec la package `ggplot`. Afin d'afficher les intervalles de crédibilité sur les graphes, nous allons récupérer les valeurs des paramètres pour les quantiles à 2.5% et 97.5%.

```r
post_veg <- pivot_longer(post, cols = 1:7, names_to = "Variables", values_to = "Valeur") # format ggplot
post_duree = pivot_longer(post, cols = 9, names_to = "Variables", values_to = "Valeur")
post_max = pivot_longer(post, cols = 8, names_to = "Variables", values_to = "Valeur")

post.stat = apply(post,2,quantile,probs=c(0.025,0.5,0.975)) # valeurs des paramètres de chaque quantile pour chaque variable
post.stat.veg = post.stat[,1:7] # paramètres veg
post.stat.duree = post.stat[,9] # paramètre duree_sub
post.stat.max = post.stat[,8] # paramètre max_sub
```

Passons à la représentation graphique de ces densités a posteriori.

```r
# Végétation
ggplot(data = post_veg, aes(x = Valeur, color = Variables)) +
  geom_density(adjust = 1.5, fill = "transparent", size = 0.7) +
  labs(x = "estimation du paramètre", y = "densité")+
  theme_ipsum()
```

![Densités végétation](https://github.com/cha-dot/Processus-de-cox-log-gaussien-INLA-SPDE-/blob/images/densite_veg.jpg?raw=true)

```r
# Durée de submersion
ggplot(data = post_duree, aes(x = Valeur, fill = Variables)) +
  geom_density(adjust = 1.5, alpha = 0.4, fill = "gray") + # alpha = transparence
  geom_vline(aes(xintercept = post.stat[2, "duree_sub"], color = "Ligne 1"), linetype = "dotdash") + # médiane
  geom_vline(aes(xintercept = post.stat[1, "duree_sub"], color = "Ligne 2"), linetype = "solid") + # quantile à 2.5%
  geom_vline(aes(xintercept = post.stat[3, "duree_sub"], color = "Ligne 2"), linetype = "solid") + # quantile à 97.5%
  scale_color_manual(values = c("black", "brown"), 
                     labels = c("Médiane", "Intervalle de crédibilité à 95%")) +
  theme_ipsum() +
  labs(x = "estimation du paramètre", fill = "Variable", color = "Légendes")
```

[Densité durée](https://github.com/cha-dot/Processus-de-cox-log-gaussien-INLA-SPDE-/blob/images/densite_duree.jpg?raw=true)

```r
# Max
ggplot(data = post_max, aes(x = Valeur, fill = Variables)) +
  geom_density(adjust = 1.5, alpha = 0.4, fill = "gray") +
  geom_vline(aes(xintercept = post.stat[2, "intensite_sub"],color = "Ligne 1"), linetype = "dotdash")+ # médiane
  geom_vline(aes(xintercept = post.stat[1, "intensite_sub"],color = "Ligne 2"), linetype = "solid")+ # quantile à 2.5%
  geom_vline(aes(xintercept = post.stat[3, "intensite_sub"], color = "Ligne 2"), linetype = "solid")+ # quantile à 97.5%
  scale_color_manual(values = c("black", "brown"),
                     labels = c("Médiane", "Intervalle de crédibilité à 95%"))+
  theme_ipsum()+
  labs(x = "estimation du paramètre", fill = "Variables", color = "Légendes")
```

![Densité max](https://github.com/cha-dot/Processus-de-cox-log-gaussien-INLA-SPDE-/blob/images/max_veg.jpg?raw=true)

Une autre représentation peut s'avérer plus lisible pour l'affichage d'un grand nombre de paramètres ou de modalités sur un même graphe, comme c'est le cas ici avec la variable `végétation`.

```r
post.stat.veg.gg <- data.frame(
  Variables = 1:ncol(post.stat.veg),
  Mediane = post.stat.veg[2, ],
  Lower_CI = post.stat.veg[1, ], # quantile à 2.5%
  Upper_CI = post.stat.veg[3, ] # quantile à 97.5%
)

ggplot(post.stat.veg.gg, aes(x = Variables)) +
  geom_point(aes(y = Mediane), size = 3, shape = 16)+ # point représentant la médiane
  geom_errorbar(aes(ymin = Lower_CI, ymax = Upper_CI), width = 0.2)+ # barres d'erreur pour les intervalles de crédibilité
  geom_text(aes(y = Mediane, label = c("cultures", "rase", "haute fauchée", "haute non fauchée", "arbustive", "roselières/scirpaies", "friches")), vjust = -1, size = 4, hjust=-0.05)+
# vjust et hjust permettent d'ajuster la position des étiquettes
  scale_x_discrete(labels = NULL)+ # n'affiche pas les étiquettes sur l'axe des x
  labs(x = "Types de végétation", y = "Distributions a posteriori")
```

![A posteriori végétation](https://github.com/cha-dot/Processus-de-cox-log-gaussien-INLA-SPDE-/blob/images/a_posteriori_veg.jpg?raw=true)

### Calcul des probabilités a posteriori

```r
sum(post[,4]>post[,3])/Nrep # proba que le paramètre de la veg non fauchée > fauchée
sum(post[,9]<0)/Nrep # proba que la duree_sub < 0
sum(post[,8]>0)/Nrep # proba que max_sub > 0
```

### Courbes de densité a posteriori des hyperparamètres

Il peut être intéressant d'explorer les lois a posteriori des hyperparamètres, c'est-à-dire de la variance et de la portée. Pour cela, on échantillonne les hyperparamètres dans le modèle.

```r
pr.hyper.tot = inla.hyperpar.sample(Nrep, ppVIV) # échantillonnage
pr.hyper.tot.df = as.data.frame(pr.hyper.tot)
pr.hyper.tot.range = pivot_longer(pr.hyper.tot.df, cols = 1, names_to = "Portée",values_to = "Valeurs") # format ggplot
pr.hyper.tot.precision = pivot_longer(pr.hyper.tot.df, cols = 2, names_to = "Précision",values_to = "Valeurs")
hyper.stats = apply(pr.hyper.tot.df, 2, quantile, probs = c(0.025, 0.5, 0.975)) # quantiles des hyperparamètres

# Portée
ggplot(data = pr.hyper.tot.range, aes(x = Valeurs, fill = Portée)) +
  geom_density(adjust = 1.5, alpha = 0.4, fill = "gray") +
  geom_vline(xintercept = hyper.stats[2, "Range for i"], linetype = "dotdash", color = "black")+ # médiane
  geom_vline(xintercept = hyper.stats[1, "Range for i"], linetype = "solid", col = "brown")+ # quantile à 2.5%
  geom_vline(xintercept = hyper.stats[3, "Range for i"], linetype = "solid", col = "brown") + # quantile à 97.5%
  labs(x = "Portée (m)", y = "Densité")
```

![Densitée portée](https://github.com/cha-dot/Processus-de-cox-log-gaussien-INLA-SPDE-/blob/images/portee_3.jpg?raw=true)

```r
# Variance
ggplot(data = pr.hyper.tot.precision, aes(x = Valeurs, fill = Précision)) +
  geom_density(adjust = 1.5, alpha = 0.4, fill = "gray") +
  geom_vline(xintercept = hyper.stats[2, "Stdev for i"], linetype = "dashed", color = "black")+
  geom_vline(xintercept = hyper.stats[1, "Stdev for i"], linetype = "solid", col = "brown")+
  geom_vline(xintercept = hyper.stats[3, "Stdev for i"], linetype = "solid", col = "brown")+
  labs(x = "Variance", y = "Densité")
```

![Densité variance](https://github.com/cha-dot/Processus-de-cox-log-gaussien-INLA-SPDE-/blob/images/variance_3.jpg?raw=true)

### Calcul des probabilités a posteriori des hyperparamètres

Il est possible de vérifier les lois a priori des hyperparamètres :

```r
sum(pr.hyper.tot[,2]>1)/Nrep # variance
sum(pr.hyper.tot[,1]<3000)/Nrep # portée
```
![Lois des hyperparamètres](https://private-user-images.githubusercontent.com/173138382/341106209-41830df4-5b69-4ace-bbb0-7842ed9836ed.PNG?jwt=eyJhbGciOiJIUzI1NiIsInR5cCI6IkpXVCJ9.eyJpc3MiOiJnaXRodWIuY29tIiwiYXVkIjoicmF3LmdpdGh1YnVzZXJjb250ZW50LmNvbSIsImtleSI6ImtleTUiLCJleHAiOjE3MTg4MDY0NDAsIm5iZiI6MTcxODgwNjE0MCwicGF0aCI6Ii8xNzMxMzgzODIvMzQxMTA2MjA5LTQxODMwZGY0LTViNjktNGFjZS1iYmIwLTc4NDJlZDk4MzZlZC5QTkc_WC1BbXotQWxnb3JpdGhtPUFXUzQtSE1BQy1TSEEyNTYmWC1BbXotQ3JlZGVudGlhbD1BS0lBVkNPRFlMU0E1M1BRSzRaQSUyRjIwMjQwNjE5JTJGdXMtZWFzdC0xJTJGczMlMkZhd3M0X3JlcXVlc3QmWC1BbXotRGF0ZT0yMDI0MDYxOVQxNDA5MDBaJlgtQW16LUV4cGlyZXM9MzAwJlgtQW16LVNpZ25hdHVyZT1iOGM4NjdkMTE5YWU1NWFmZTc1ZmE2MzIyOTk5OGNiZTBiMWNmNjE2M2JkZmRiZmM4MTVkYjY5Y2IwMDQwNGU3JlgtQW16LVNpZ25lZEhlYWRlcnM9aG9zdCZhY3Rvcl9pZD0wJmtleV9pZD0wJnJlcG9faWQ9MCJ9.4fvxy4mRUI02cehSG-zWxRvQoiOLxKJ7rtCGp8R0_1A)

