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
La matrice d'observation A permet de faire le lien entre la mesh, les observations et les covariables. Nous allons construire manuellement cette matrice, ainsi que les vecteurs d'observations, de pondération et des covariables. Se référer à l'image ci-dessous pour une meilleure compréhension de l'élaboration des différentes parties.

<a id="matrice-a"></a>
![Matrice d'observation A](https://private-user-images.githubusercontent.com/173138382/340686454-5abee3ed-33f9-40f9-b746-c7e22400fb3e.PNG?jwt=eyJhbGciOiJIUzI1NiIsInR5cCI6IkpXVCJ9.eyJpc3MiOiJnaXRodWIuY29tIiwiYXVkIjoicmF3LmdpdGh1YnVzZXJjb250ZW50LmNvbSIsImtleSI6ImtleTUiLCJleHAiOjE3MTg3MTY0NjIsIm5iZiI6MTcxODcxNjE2MiwicGF0aCI6Ii8xNzMxMzgzODIvMzQwNjg2NDU0LTVhYmVlM2VkLTMzZjktNDBmOS1iNzQ2LWM3ZTIyNDAwZmIzZS5QTkc_WC1BbXotQWxnb3JpdGhtPUFXUzQtSE1BQy1TSEEyNTYmWC1BbXotQ3JlZGVudGlhbD1BS0lBVkNPRFlMU0E1M1BRSzRaQSUyRjIwMjQwNjE4JTJGdXMtZWFzdC0xJTJGczMlMkZhd3M0X3JlcXVlc3QmWC1BbXotRGF0ZT0yMDI0MDYxOFQxMzA5MjJaJlgtQW16LUV4cGlyZXM9MzAwJlgtQW16LVNpZ25hdHVyZT1kZjcyOGEyMjM4OTIxYzc3NzAzMDYzZTgwYmZkMzlhZjI2ZTkzYjE1ODZmOGJjM2I2OTE5OGQ3MjRlNTQ5MTQ0JlgtQW16LVNpZ25lZEhlYWRlcnM9aG9zdCZhY3Rvcl9pZD0wJmtleV9pZD0wJnJlcG9faWQ9MCJ9.ekKicQE8Z83QZ2M8rMqB-0TWtb1cFpvi4YmUFmLjgsw)

### Vecteur des observations et de pondération

Le vecteur y.pp est composé de 0 et de 1, c'est le vecteur des observations. Le vecteur e.pp est composé des poids pour chacune des cellules de voronoi (c'est-à-dire la surface d'intersection entre les surfaces prospectées et les cellules de voronoi).

```r
y.pp = NULL
e.pp = NULL
for (i in 1:nrow(n)){ # pour chaque année
  for (j in 1:ncol(n)){ # pour chaque période de prospection
    y.pp = c(y.pp,rep(0:1, c(nv, n[i,j]))) # vecteur qui permettra de relier les observations à la mesh
    e.pp = c(e.pp,c(wbis, rep(0, n[i,j]))) # vecteur qui permettra de pondérer l'intensité par la surface de prospection
  }
}
```

### Matrice d'observation A

La [Matrice d'observation A](#matrice-a) est divisée en matrices I (correspondant aux cellules de voronoi), et en matrices L (correspondant aux observations d'une période et d'une année).

```r
imat <- Diagonal(nv, rep(1, nv)) # matrice I diagonale de dimensions nv (nombre de cellules de Voronoi), représente la mesh
lmat <- inla.spde.make.A(mesh, coordinates(GB_PE_eau_fauche_LAMB)) # matrice L qui représente les observations dans l'espace
A.pp = NULL
for (i in 2018:2023){
  for (j in 1:2){
    indice = (1:nrow(GB_PE_eau_fauche_LAMB))[(GB_PE_eau_fauche_LAMB$annee==i)&(GB_PE_eau_fauche_LAMB$num_passag==j)]
    A.pp = rbind(A.pp,rbind(imat,lmat[indice,])) # Matrice d'observation A différente pour chaque année (A11, A12...) # on ajoute les matrices I et L les unes à la suite des autres
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
  v2$id_type_ge[is.na(v2$id_type_ge)]=0 #idem
  v$Rgpt_vg_3 = ifelse(v$Rgpt_vg_3 %in% c(13,1,11,8) & v2$id_type_ge == 1, 14, v$Rgpt_vg_3) # si la veg est herbacée haute (13, 1, 11, 8) et qu'elle est fauchée (1), alors on la met dans un groupe (14)
  v$Rgpt_vg_3 = ifelse(v$Rgpt_vg_3 %in% c(13,1,11,8) & v2$id_type_ge == 0, 15, v$Rgpt_vg_3) # si elle n'est pas fauchée (0), on la met dans un autre groupe (15)
  return(v$Rgpt_vg_3)
}
```
On applique la fonction sur les données pour chaque année et période de prospection. Dans la fonction, on renseigne les coordonnées spatiales de la mesh et des observations de l'année i et la période j, et le fichier spatiale de la fauche à l'année i-1 (choix d'étude).

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

