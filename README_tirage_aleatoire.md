# Suppression progressive de points d'écoute aléatoirement

Cela consiste à supprimer progressivement des points d'écoute aléatoirement sur le domaine d'étude afin de voir si les résultats obtenus sont les mêmes que ceux du modèle avec la totalité des points d'écoute.

## Table des matières

- [Packages](#packages)
- [Modification du nombre de points d'écoute](#modification-du-nombre-de-points-découte)
- [Nouveau jeu de données](#nouveau-jeu-de-données)
- [Modèle](#modèle)
  - [Sauvegarde du summary](#sauvegarde-du-summary)
  - [Sauvegarde de l'AUC, du MSE et du RMSE](#sauvegarde-de-lauc-du-mse-et-du-rmse)
  - [Sauvegarde des distributions a posteriori](#sauvegarde-des-distributions-a-posteriori)
- [Représentations graphiques des résultats](#représentations-graphiques-des-résultats)
  - [AUC](#auc)
  - [RMSE](#rmse)
  - [Paramètres du modèle](#paramètres-du-modèle)
    - [Intensité de submersion](#intensité-de-submersion)
    - [Durée de submersion](#durée-de-submersion)
    - [Végétation](#végétation)
      - [Cultures](#cultures)
      - [Végétation rase](#végétation-rase)
      - [Végétation haute fauchée](#végétation-haute-fauchée)
      - [Végétation haute non fauchée](#végétation-haute-non-fauchée)
      - [Végétation arbustive](#végétation-arbustive)
      - [Roselières/scirpaies](#roselières/scirpaies)
      - [Friches](#friches)

## Packages

Les packages à charger sont les suivants :

```r
library(INLA)
library(inlabru)
library(ggplot2)
library(raster)
library(sf)
library(sp)
library(tidyverse)
library(viridis)
library(rnaturalearth)
library(hrbrthemes)
library(dplyr)
library(tidyr)
library(pROC)
```

## Modification du nombre de points d'écoute

```r
buffers_a_supprimer = sample(unique(GB_PE_eau_fauche$numero_buf), points) # échantillonnage d'un nb (points) de points d'écoute
```

## Nouveau jeu de données

Nous supprimons les individus liés aux points d'écoute à supprimer (échantillonnés antérieurement) :

```r
GB_10_1 = GB_PE_eau_fauche[!(GB_PE_eau_fauche$numero_buf%in%buffers_a_supprimer),] # supprimer les individus des points d'écoute supprimés
GB_PE_eau_fauche_2 = SpatialPoints(coords = GB_10_1[,c("x_wgs84", "y_wgs84")]) # conversion de format
proj4string(GB_PE_eau_fauche_2) = CRS("+proj=longlat +datum=WGS84") # crs
points_lambert = spTransform(GB_PE_eau_fauche_2, CRSobj = CRS(contour_sp@proj4string@projargs)) # reprojection
GB_10_1_LAMB = SpatialPointsDataFrame(coords = coordinates(points_lambert), # conversion de format
                                               data = GB_10_1,
                                               proj4string = CRS(contour_sp@proj4string@projargs))
```

Nous supprimons les buffers tirés aléatoirement :

```r
coord_Buffer = GB_10_1[, c("X_PE", "Y_PE")] # coordonnées des buffers
coord_Buffer = coord_Buffer[!duplicated(coord_Buffer),] # liste des buffers et leurs coordonnées sans doublons
Buffer_reduit = SpatialPoints(coords = coord_Buffer, CRS("+proj=longlat +datum=WGS84")) # conversion de format
Buffer_reduit = st_as_sf(Buffer_reduit) # conversion de format
Buffer_reduit = st_transform(Buffer_reduit, crs = proj4string(GB_PE_eau_fauche_LAMB))
Buffer_reduit = st_buffer(Buffer_reduit, 150) # 150 = nb de traits qui servent à la formation du cercle
Buffer_reduit_sp = as(Buffer_reduit, "Spatial")
```
Nous pouvons représenter les buffers restants :

```r
ggplot() +
  geom_sf(data = st_as_sf(GB_10_1_LAMB), col = "blue", size = 0.1) +
  geom_sf(data = st_as_sf(contour_sp), alpha = 0, col = "black") +
  geom_sf(data = st_as_sf(Buffer_reduit_sp), alpha = 0.5, col = "orange", linewidth = 0.5)
```

![Aléatoire répartition des buffers](https://github.com/cha-dot/Processus-de-cox-log-gaussien-INLA-SPDE-/blob/images/buffers_aleatoires.jpg?raw=true)

## Modèle

Le reste du modèle ne change pas (pour plus d'explications du modèle, se référer à README_INLA.md)

```r
 tmp = spsample(contour_sp@polygons[[1]], n = 1000, type = "regular")
  coordLIM = rbind(contour_sp@polygons[[1]]@Polygons[[1]]@coords, tmp@coords)
  rm(tmp)
  
  bndint = inla.nonconvex.hull(coordLIM, convex = -.02)
  bndext = inla.nonconvex.hull(coordLIM, convex = -.1)

  mesh = inla.mesh.2d(loc = rbind(coordinates(GB_10_1_LAMB)),
                      boundary = list(int = bndint, out = bndext),
                      max.edge = c(125, 1000), cutoff = 50,
                      crs = "+proj=lcc +lat_0=46.5 +lon_0=3 +lat_1=49 +lat_2=44 
                    +x_0=700000 +y_0=6600000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")

  matern = inla.spde2.pcmatern(mesh,
                               alpha = 2,
                               prior.sigma = c(1, 0.5), # P(sigma > 1) = 0.5
                               prior.range = c(3000, 0.9)) # P(range < 100) = 0.9

  dmesh <- book.mesh.dual(mesh)
  Buffer_reduit = st_as_sf(Buffer_reduit_sp)
  dmesh_sf = st_as_sf(dmesh)
  st_crs(dmesh_sf) = st_crs(Buffer_reduit)
  w <- sapply(1:length(dmesh), function(i) { # calcul de l'aire d'intersection entre les cellules et les buffers pour chaque cellule
    if (nrow(st_intersection(dmesh_sf[i, ], Buffer_reduit)) > 0)
      return(st_area(st_intersection(dmesh_sf[i, ], Buffer_reduit)))
    else return(0)
  })
  
  wbis = unlist(lapply(w, max)) # attribution de l'aire d'intersection (w) max (quand 2 buffers)
  nv = mesh$n
  n = table(GB_10_1_LAMB$annee,GB_10_1_LAMB$num_passag) # nombre d'observations par année par période

 y.pp = NULL
  e.pp = NULL
  for (i in 1:nrow(n)){ # pour chaque année
    for (j in 1:ncol(n)){ # pour chaque période
      y.pp = c(y.pp,rep(0:1, c(nv, n[i,j])))
      e.pp = c(e.pp,c(wbis, rep(0, n[i,j])))
    }
  }
  
  imat <- Diagonal(nv, rep(1, nv))
  lmat <- inla.spde.make.A(mesh, coordinates(GB_10_1_LAMB))
  A.pp = NULL
  for (i in 2018:2023){
    for (j in 1:2){
      indice = (1:nrow(GB_10_1_LAMB))[(GB_10_1_LAMB$annee==i)&(GB_10_1_LAMB$num_passag==j)]
      A.pp = rbind(A.pp,rbind(imat,lmat[indice,]))
    }
  }

## Covariables
  veg2 = raster("Rgpt_vg_3.tif")
  vegSG2 = as(veg2, "SpatialGridDataFrame")
  
  setwd("Fauche") # à définir pour que la fonction fauche fonctionne
  shapefiles = c("Fauche_2017_seule.shp", "Fauche_2018_seule.shp", "Fauche_2019_seule.shp", "Fauche_2020_seule.shp", "Fauche_2021_seule.shp", "Fauche_2022_seule.shp")
  
  eau = read.csv("Hauteurs_eau_finales.csv", header = T, stringsAsFactors = T)
  
  f.veg.moins = function(x, y, shape) {
    spp <- SpatialPoints(data.frame(x = x, y = y), proj4string = fm_sp_get_crs(vegSG2))
    proj4string(spp) <- fm_sp_get_crs(vegSG2)
    v <- over(spp, vegSG2)
    v2 <- over(spp, shape)
    v$Rgpt_vg_3[is.na(v$Rgpt_vg_3)]=0
    v2$id_type_ge[is.na(v2$id_type_ge)]=0
    v$Rgpt_vg_3 = ifelse(v$Rgpt_vg_3 %in% c(13,1,11,8) & v2$id_type_ge == 1, 14, v$Rgpt_vg_3)
    v$Rgpt_vg_3 = ifelse(v$Rgpt_vg_3 %in% c(13,1,11,8) & v2$id_type_ge == 0, 15, v$Rgpt_vg_3) # haute fauchée sinon, non fauchée
    return(v$Rgpt_vg_3)
  }
  
  
  vegTMP=NULL
  for (i in 2018:2023){
    for (j in 1:2){
      shapefile = st_read(shapefiles[i - 2017]) # pour ramener les indices de 1 à 6
      shapefile = as(shapefile, "Spatial")
      crs(shapefile) = "+proj=lcc +lat_0=46.5 +lon_0=3 +lat_1=49 +lat_2=44 +x_0=700000 +y_0=6600000 +ellps=GRS80
+towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
      indice = (1:nrow(GB_10_1_LAMB))[(GB_10_1_LAMB$annee==i)&(GB_10_1_LAMB$num_passag==j)]
      vegTMP = c(vegTMP,
                 f.veg.moins(x=c(mesh$loc[,1], coordinates(GB_10_1_LAMB)[indice,1]),
                             y = c(mesh$loc[,2], coordinates(GB_10_1_LAMB)[indice,2]),
                             shape = shapefile))
    }
  }


# eau_max
  f.eau.max = function(annee, prospection) {
    data = eau[eau$annee == annee,]
    if(prospection == 1) {return(data[,1])} else{return(data[,2])}
  }
  
  eau_max=NULL
  for (i in 2018:2023){
    for (j in 1:2){
      eau_max_1 = f.eau.max(annee = i, prospection = j)
      eau_max = c(eau_max, eau_max_1)
    }
  }
  eau_max = matrix(data = eau_max, ncol = 2, byrow = T)
  
  eau_maxTMP = NULL
  for(i in 2018:2023) {
    for(j in 1:2) {
      eau_max_I = rep(eau_max[i-2017,j], nv)
      eau_max_A = rep(eau_max[i-2017,j], n[i-2017,j])
      eau_maxTMP = c(eau_maxTMP, eau_max_I, eau_max_A)
    }
  }
  
  veg = vegTMP
  veg[vegTMP == 6 | vegTMP == 12 | vegTMP == 4 | vegTMP == 7| vegTMP == 9 | vegTMP == 10] = 1 # rase
  veg[vegTMP == 14] = 2 # haute fauchee
  veg[vegTMP == 15] = 3 # haute non fauchee
  veg[vegTMP == 2] = 4 # arbustive
  veg[vegTMP == 3] = 5 # roselières 
  veg[vegTMP == 5] = 6 # friches
  veg = as.factor(veg)

 # eau durée
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

##stack inla
  stk.pp <- inla.stack(                   
    data = list(y = y.pp, e = e.pp),
    A = list(1, A.pp),
    effects = list(list(veg = veg, max_sub = eau_maxTMP, duree_sub = eau_dureeTMP), #b0 = 1 (quand y a juste le champ spatial)
                   list(i = 1:nv)),
    tag = 'pp')

ppVIV <- inla(y ~ 1 + f(i, model = matern) + max_sub + duree_sub  + f(veg, model = "iid", constr = T),
                family = 'poisson', data = inla.stack.data(stk.pp),
                control.predictor = list(A = inla.stack.A(stk.pp)),
                E = inla.stack.data(stk.pp)$e,
                control.compute=list(dic=TRUE, cpo=TRUE, config = TRUE, return.marginals.predictor=TRUE),
                control.inla=list(strategy="simplified.laplace",int.strategy="eb"))

modtest = summary(ppVIV)
```
### Sauvegarde du summary

Si le modèle est répété pour différents nombres de points d'écoute, il est possible de sauvegarder certains objets, tels que le summary, l'AUC ou encore les distributions a posteriori pour une comparaison de résultats.

Si le modèle est répété `repet` fois sur un cas de figure, `92-points` (par exemple, 7 points d'écoute supprimés donc 85 points d'écoute au total), alors il est possible de sauvegarder le `summary` de cette manière :

```r
chemin = "~/summaries"
fichier = paste0("nouv_summary_", repet, "_nPE_", 92-points, ".RData")
nom_complet = file.path(chemin, fichier)
save(modtest, file = nom_complet)
```
Si le script doit être répété de nombreuses fois, il peut être judicieux de réduire le nombre de pixels dans `resopred` pour réduire les temps de calcul (ici, `resopred = 100` semble être un bon compromis entre le temps de calculs et la précision des résultats).

```r
  toto = extent(contour_sp)
  xrange = toto[2]-toto[1]
  yrange = toto[4]-toto[3]
  resopred=100
  st_crs(contour) = NA
  inla.identical.CRS(mesh, contour)
  contour = st_set_crs(contour, inla.CRS(mesh))
  pxl <- pixels(mesh, mask=contour, nx = resopred,
                ny = round(resopred*yrange/xrange))
  projgrid <- inla.mesh.projector(mesh, coordinates(pxl))

  indN = 0
  indNV = 1
  predtot = NULL
  for(i in 1:6){
    for(j in 1:2){
      predtot = rbind(predtot, ppVIV$summary.fitted.values$'0.5quant'[((indNV-1)*nv+indN+1):(indNV*nv+indN)])
      indN = indN + n[i,j]
      indNV = indNV + 1
      
    }
  }
  
  r = as(pxl, "SpatialPolygonsDataFrame")
  rSF=st_as_sf(r)
  datvivSF = st_as_sf(GB_PE_eau_fauche_LAMB)
  st_crs(rSF) = st_crs(datvivSF)
  test=st_intersects(rSF,datvivSF)
  NobsVIV=unlist(lapply(test,length))
  fitval = inla.mesh.project(projgrid,
                             apply(predtot,2,sum))


# courbe ROC
  NobsVIV_bis = as.numeric(NobsVIV > 0) # les cellules où il y a présence
  predTEST = 1-exp(-as.numeric(fitval*st_area(rSF)))
  testroc=roc(NobsVIV_bis,predTEST)
  AUC = testroc$auc
  
  mse <- mean((predTEST - NobsVIV)^2)
  RMSE = sqrt(mse)
```

### Sauvegarde de l'AUC, du MSE et du RMSE

Il est possible de sauvegarder ces 3 métriques de la manière suivante :

```r
chemin = "~/metriques"
fichier = paste0("nouv_metriques_", repet, "_nPE_", 92-points,".RData")
nom_complet = file.path(chemin, fichier)
save(AUC, mse, RMSE, file = nom_complet)
```
Toujours dans l'objectif d'un gain de temps si le code doit être répliqué de nombreuses fois, il est possible de réduire le nombre d'échantillons `Nrep` du modèle (`Nrep=500` semble être un bon compromis entre les temps de calcul et la précision des résultats).

```r
  Nrep=500
  pr.int.tot <- inla.posterior.sample(Nrep,ppVIV) # Nrep échantillons du modèle
  ind=grep("veg",rownames(pr.int.tot[[1]]$latent)) # paramètres végétation
  m = grep("max_sub", rownames(pr.int.tot[[1]]$latent)) # paramètre intensité de submersion
  d = grep("duree_sub", rownames(pr.int.tot[[1]]$latent)) # paramètre durée de submersion
  post=matrix(unlist(lapply(pr.int.tot,function(x){x$latent[c(ind,m,d),]})),nrow=Nrep,byrow=T) # Nrep valeurs des paramètres pour chaque variable
  post = as.data.frame(post)
```

Lors de la suppression de points d'écoute, il est possible que certaines variables ou modalités spatiales, dans notre cas, certains types de végétation, ne soient plus représentés. Aucune information sur la distribution a posteriori du paramètre ne sera alors disponible.

`colnames(post) = c(rownames(pr.int.tot[[1]]$latent)[ind], "intensite_sub", "duree_sub")` à la place de `colnames(post)=c("cultures", "rase", "haute fauchée", "haute non fauchée", "arbustive", "roselières/scirpaies", "friches", "intensite_sub", "duree_sub")` permet de prendre en compte la possibilité d'une absence d'un type de végétation.

```r
  colnames(post) = c(rownames(pr.int.tot[[1]]$latent)[ind], "intensite_sub", "duree_sub") # nom des colonnes
  anciens = c("veg:1", "veg:2", "veg:3", "veg:4", "veg:5", "veg:6", "veg:7")
  nouv = c("cultures", "rase", "haute fauchée", "haute non fauchée", "arbustive", "roselières/scirpaies", "friches")
    
  for (p in 1:length(anciens)) { # remplacer "veg:x" par les noms des types de vg
      if (anciens[p] %in% colnames(post)) {
        colnames(post)[colnames(post) == anciens[p]] = nouv[p]
      }
    }
    
  post.stat = apply(post,2,quantile,probs=c(0.025,0.5,0.975)) # valeurs des paramètres pour chaque quantile de chaque variable
  post.stat.veg = post.stat[, which(colnames(post) %in% nouv)] # quantiles des paramètres veg
  post.stat.duree = post.stat[, "duree_sub"] # quantiles de la durée de submersion
  post.stat.max = post.stat[, "intensite_sub"] # quantiles de l'intensité de submersion
```
### Sauvegarde des distributions a posteriori

```r
chemin = "~/post_nouveaux"
fichier = paste0("post_", repet, "_nPE_", 92-points,".RData")
nom_complet = file.path(chemin, fichier)
save(post.stat.veg, post.stat.duree, post.stat.max, file = nom_complet)
```

## Représentations graphiques des résultats

Le modèle a été réalisé pour différents nombres de points d'écoute (de 0 à 85 par pas de 10). Chaque cas de figure a été répliqué 30 fois. Ainsi, nous ferons la moyenne des métriques pour chaque cas de figure.

### AUC

Dans un premier temps, nous chargeons les données d'AUC de chaque réplicat pour chaque nombre de points d'écoute.

```r
all_AUC = list() # stockera toutes les valeurs d'AUC

for(points in c(5, 15, 25, 35, 45, 55, 65, 75, 85)) {
  
  all_AUC[[as.character(points)]] <- list() # stockera les AUC d'un cas de figure
  
for (repet in 1:30) { # pour chaque réplicat
  
  chemin = "~/metriques"
  fichier = paste0("nouv_metriques_", repet, "_nPE_", points, ".RData") 
  nom_complet = file.path(chemin, fichier) # nom du fichier en fonction du nb de points d'écoute et du réplicat
  load(nom_complet) # charge le fichier
  
  all_AUC[[as.character(points)]][[paste0("repet_", repet)]] = AUC # valeur d'AUC (+indice pour les réplicats (x 30))
}}
```
A présent, nous allons faire la moyenne des AUC pour chaque cas de figure :

```r
mean_AUC <- list() # stockera les moyennes


for(points in c(5, 15, 25, 35, 45, 55, 65, 75, 85)) {
  auc_values <- unlist(all_AUC[[as.character(points)]]) # liste des AUC pour un cas de figure
  mean_AUC[[as.character(points)]] <- mean(auc_values, na.rm = TRUE) # moyenne de ces AUC
}
```

Nous établissons un dataframe à partir de ces données qui servira pour la représentation graphique :

```r
points <- c(5, 15, 25, 35, 45, 55, 65, 75, 85)
mean_values <- unlist(mean_AUC) # liste des moyennes convertie en vecteur

result_matrix <- cbind(points, mean_values) # tableau
result_matrix <- rbind(result_matrix, c(92, 0.8509)) # ajout de la valeur du modèle qui comprend tous les points d'écoute
AUC = as.data.frame(result_matrix)
```
Passons à la représentation graphique :

```r
 ggplot(data = AUC, aes(x = points, y = mean_values))+ # AUC en fonction du nombre de points d'écoute
  geom_line()+ # ligne qui relie les points moyens
  geom_point(color = "black", size = 1)+ # points moyens
  labs(x = "Nombre de points d'écoute", y = "AUC moyenne") + # légende
  theme_minimal()+
  scale_x_continuous(breaks = unique(AUC$points), labels = unique(AUC$points)) # affiche toutes les valeurs du vecteur "points" sur l'axe x

```

### RMSE

La racine de l'erreur quadratique moyenne (RMSE) est une métrique qui permet de comparer des modèles entre eux (plus le RMSE est faible, meilleur est le modèle).

La méthodologie pour charger et afficher les valeurs de RMSE est similaire à celle utilisée pour l'AUC.

Dans un premier temps, nous chargeons toutes les données de RMSE et les regroupons par nombre de points d'écoute (par cas de figure donc).

```r
all_RMSE = list() # stockera toutes les valeurs de RMSE

for(points in c(5, 15, 25, 35, 45, 55, 65, 75, 85)) {
  
  # Initialiser une sous-liste pour chaque point
  all_RMSE[[as.character(points)]] <- list() stockera les RMSE d'un cas de figure
  
  for (repet in 1:30) { # pour chaque réplicat
    
    chemin = "~/metriques"
    fichier = paste0("nouv_metriques_", repet, "_nPE_", points, ".RData") 
    nom_complet = file.path(chemin, fichier) # nom du fichier en fonction du nb de points d'écoute et du réplicat
    load(nom_complet) # charge le fichier
    
    # Ajouter les résultats à la liste avec une clé indiquant la répétition
    all_RMSE[[as.character(points)]][[paste0("repet_", repet)]] = RMSE # valeur de RMSE (+indice pour les réplicats (x30))
  }}
```

A présent, nous allons faire la moyenne des RMSE pour chaque cas de figure :

```r
mean_RMSE <- list() # stockera les moyennes

for(points in c(5, 15, 25, 35, 45, 55, 65, 75, 85)) {
  rmse_values <- unlist(all_RMSE[[as.character(points)]]) # liste des RMSE pour un cas de figure
  mean_RMSE[[as.character(points)]] <- mean(rmse_values, na.rm = TRUE) # moyenne de ces RMSE
}
```

Nous établissons un dataframe, à partir de ces données, qui servira pour la représentation graphique :

```r
points <- c(5, 15, 25, 35, 45, 55, 65, 75, 85)
mean_values2 <- unlist(mean_RMSE) # liste des moyennes convertie en vecteur

result_matrix2 <- cbind(points, mean_values2) # tableau
result_matrix2 <- rbind(result_matrix2, c(92, 2.3194)) # ajout de la valeur du modèle qui comprend tous les points d'écoute
RMSE = as.data.frame(result_matrix2)
```

Passons à la représentation graphique :

```r
ggplot(data = RMSE, aes(x = points, y = mean_values2))+
  geom_line()+ # ligne relier les points moyens
  geom_point(color = "black", size = 1)+ # points moyens
  labs(x = "Nombre de points d'écoute", y = "RMSE moyen") + # légende
  theme_minimal()+
  scale_x_continuous(breaks = unique(RMSE$points), labels = unique(RMSE$points)) # affiche toutes les valeurs du vecteur "points" sur l'axe x
```

### Paramètres du modèle

#### Intensité de submersion

Pour chaque nombre de points d'écoute, les 30 valeurs d'intensité de submersion sont chargées de la manière suivante :

```r
all_param = list() # stockera toutes les intensités de submersion
  
for(points in c(5, 15, 25, 35, 45, 55, 65, 75, 85)) { # pour chaque cas de figure
    
all_param[[as.character(points)]] <- list() # sous-liste correspondant à un cas de figure
    
for (repet in 1:30) { # 30 réplicats
      
chemin = "~/post_nouveaux"
fichier = paste0("post_", repet, "_nPE_", points, ".RData") 
nom_complet = file.path(chemin, fichier) # nom du fichier en fonction du nb de points d'écoute et du réplicat
load(nom_complet) # charge le fichier
      
all_param[[as.character(points)]][[paste0("repet_", repet)]] <- post.stat.max # intensité de submersion
}}
```

Nous allons à présent calculer les moyennes de ces intensités de submersion pour chaque nombre de points d'écoute, ainsi que les moyennes des intervalles de crédibilité associés :

```r
all_max = list() # stockera les moyennes d'intensité de submersion
all_max_IC1 = list() # stockera les moyennes des quantiles 2.5%
all_max_IC2 = list() # stockera les moyennes des quantiles 97.5%

for(points in c(5, 15, 25, 35, 45, 55, 65, 75, 85)) { # pour chaque cas de figure
  all_max[[as.character(points)]] = NULL
  all_max_IC1[[as.character(points)]] = NULL
  all_max_IC2[[as.character(points)]] = NULL

for(repet in 1:30) { # pour chaque réplicat
  modtest <- all_param[[as.character(points)]][[paste0("repet_", repet)]] # stocke l'intensité de submersion correspondante
      
  max = modtest[2] # médiane
  max_IC1 = modtest[1] # quantile 2.5%
  max_IC2 = modtest[3] # quantile 97.5%

  all_max[[as.character(points)]] = c(all_max[[as.character(points)]], max) # stocke toutes les médianes pour le nb de points d'écoute en question
  all_max_IC1[[as.character(points)]] = c(all_max_IC1[[as.character(points)]], max_IC1) # stocke quantiles 2.5% pour le nb de points d'écoute en question
  all_max_IC2[[as.character(points)]] = c(all_max_IC2[[as.character(points)]], max_IC2) # stocke quantiles 97.5% pour le nb de points d'écoute en question
    }
  }
  
  mean_max = sapply(all_max, function(x) mean(unlist(x), na.rm = TRUE)) # moyenne des médianes
  mean_max_IC1 = sapply(all_max_IC1, function(x) mean(unlist(x), na.rm = TRUE)) # moyenne des quantiles 2.5%
  mean_max_IC2 = sapply(all_max_IC2, function(x) mean(unlist(x), na.rm = TRUE)) # moyenne des quantiles 97.5%
```

Afin de représenter graphiquement les résultats, nous organisons les données dans un dataframe :

```r
  points <- c(5, 15, 25, 35, 45, 55, 65, 75, 85)
  
  mean_max_df <- data.frame( # on détermine les différentes colonnes
    points = points,
    mean_max = mean_max,
    max_IC1 = mean_max_IC1,
    max_IC2 = mean_max_IC2
  )
  
  ggplot(data=mean_max_df, aes(x=points, y=mean_max, ymin=max_IC1, ymax=max_IC2)) + 
    geom_line() + # ligne pour les médianes
    geom_ribbon(alpha=0.5)+ # ruban pour les intervalles de crédibilité
    labs(x = "Nombre de points d'écoute", y = "Estimation de l'intensité de submersion") # légende
```

#### Durée de submersion

Pour chaque nombre de points d'écoute, les 30 valeurs de durée de submersion sont chargées de la manière suivante :

```r
all_param = list() # stockera toutes les durées de submersion
  
for(points in c(5, 15, 25, 35, 45, 55, 65, 75, 85)) { # pour chaque cas de figure
    
  all_param[[as.character(points)]] <- list() # sous-liste correspondant à un cas de figure
    
  for (repet in 1:30) { # 30 réplicats
      
  chemin = "~/post_nouveaux"
  fichier = paste0("post_", repet, "_nPE_", points, ".RData")
  nom_complet = file.path(chemin, fichier) # nom du fichier en fonction du nb de points d'écoute et du réplicat
  load(nom_complet) # charge le fichier
      
  all_param[[as.character(points)]][[paste0("repet_", repet)]] <- post.stat.duree # durée de submersion
    }}
```
Nous allons à présent calculer les moyennes de ces durées de submersion pour chaque nombre de points d'écoute, ainsi que les moyennes des intervalles de crédibilité associés :

```r
all_duree = list() # stockera les moyennes des durées de submersion
all_duree_IC1 = list() # stockera les moyennes des quantiles 2.5%
all_duree_IC2 = list() # stockera les moyennes des quantiles 97.5%
  
for(points in c(5, 15, 25, 35, 45, 55, 65, 75, 85)) { # pour chaque cas de figure
  all_duree[[as.character(points)]] = NULL
  all_duree_IC1[[as.character(points)]] = NULL
  all_duree_IC2[[as.character(points)]] = NULL
    
for(repet in 1:30) { # pour chaque réplicat
  modtest <- all_param[[as.character(points)]][[paste0("repet_", repet)]] # stocke la durée de submersion correspondante
        
  duree = modtest[2] # médiane
  duree_IC1 = modtest[1] # quantile 2.5%
  duree_IC2 = modtest[3] # quantile 97.5%
        
  all_duree[[as.character(points)]] = c(all_duree[[as.character(points)]], duree) # stocke toutes les médianes pour le nb de points d'écoute en question
  all_duree_IC1[[as.character(points)]] = c(all_duree_IC1[[as.character(points)]], duree_IC1) # stocke quantiles 2.5% pour le nb de points d'écoute en question
  all_duree_IC2[[as.character(points)]] = c(all_duree_IC2[[as.character(points)]], duree_IC2) # stocke quantiles 97.5% pour le nb de points d'écoute en question
    }
  }
  
mean_duree = sapply(all_duree, function(x) mean(unlist(x), na.rm = TRUE)) # moyenne des médianes
mean_duree_IC1 = sapply(all_duree_IC1, function(x) mean(unlist(x), na.rm = TRUE)) # moyenne des quantiles 2.5%
mean_duree_IC2 = sapply(all_duree_IC2, function(x) mean(unlist(x), na.rm = TRUE)) # moyenne des quantiles 97.5%
```
Afin de représenter graphiquement les résultats, nous organisons les données dans un dataframe :

```r
points <- c(5, 15, 25, 35, 45, 55, 65, 75, 85)
  
mean_duree_df <- data.frame( # on détermine les différentes colonnes
  points = points,
  mean_duree = mean_duree,
  duree_IC1 = mean_duree_IC1,
  duree_IC2 = mean_duree_IC2
  )
  
    
ggplot(data=mean_duree_max, aes(x=points, y=mean_duree, ymin=duree_IC1, ymax=duree_IC2)) + 
  geom_line() + # ligne pour les médianes
  geom_ribbon(alpha=0.5)+ # ruban pour les intervalles de crédibilité
  labs(x = "Nombre de points d'écoute", y = "Estimation de la durée de submersion") # légende
  ```

#### Végétation

Nous allons procéder de la même manière pour tous les paramètres de la végétation.

Dans un premier temps, nous chargeons tous les fichiers :

```r
all_param = list() # stockera les paramètres de végétation

for(points in c(5, 15, 25, 35, 45, 55, 65, 75, 85)) { # pour chaque cas de figure
  
  all_param[[as.character(points)]] <- list() # sous-liste correspondant à un cas de figure
  
  for (repet in 1:30) { # 30 réplicats
    
    chemin = "~/post_nouveaux"
    fichier = paste0("post_", repet, "_nPE_", points, ".RData")
    nom_complet = file.path(chemin, fichier) # nom du fichier en fonction du nb de points d'écoute et du réplicat
    load(nom_complet) #charge le fichier
    
    all_param[[as.character(points)]][[paste0("repet_", repet)]] <- post.stat.veg # paramètres de végétation
  }}
```
A présent, nous calculons les moyennes des médianes et des quantiles (2.5% et 97.5%) de chacune des modalités :

```r
# Stockera les moyennes des paramètres
all_cultures = list()
all_rase = list()
all_haute_fauchee = list()
all_haute_non_fauchee = list()
all_arbustive = list()
all_roselieres = list()
all_friches = list()

# Stockera les moyennes des quantiles 2.5%
all_cultures_IC1 = list()
all_rase_IC1 = list()
all_haute_fauchee_IC1 = list()
all_haute_non_fauchee_IC1 = list()
all_arbustive_IC1 = list()
all_roselieres_IC1 = list()
all_friches_IC1 = list()

# Stockera les moyennes des quantiles 97.5%
all_cultures_IC2 = list()
all_rase_IC2 = list()
all_haute_fauchee_IC2 = list()
all_haute_non_fauchee_IC2 = list()
all_arbustive_IC2 = list()
all_roselieres_IC2 = list()
all_friches_IC2 = list()


for(points in c(5, 15, 25, 35, 45, 55, 65, 75, 85)) {
  all_cultures[[as.character(points)]] = NULL
  all_rase[[as.character(points)]] = NULL
  all_haute_fauchee[[as.character(points)]] = NULL
  all_haute_non_fauchee[[as.character(points)]] = NULL
  all_arbustive[[as.character(points)]] = NULL
  all_roselieres[[as.character(points)]] = NULL
  all_friches[[as.character(points)]] = NULL
  
  all_cultures_IC1[[as.character(points)]] = NULL
  all_rase_IC1[[as.character(points)]] = NULL
  all_haute_fauchee_IC1[[as.character(points)]] = NULL
  all_haute_non_fauchee_IC1[[as.character(points)]] = NULL
  all_arbustive_IC1[[as.character(points)]] = NULL
  all_roselieres_IC1[[as.character(points)]] = NULL
  all_friches_IC1[[as.character(points)]] = NULL  
  
  all_cultures_IC2[[as.character(points)]] = NULL
  all_rase_IC2[[as.character(points)]] = NULL
  all_haute_fauchee_IC2[[as.character(points)]] = NULL
  all_haute_non_fauchee_IC2[[as.character(points)]] = NULL
  all_arbustive_IC2[[as.character(points)]] = NULL
  all_roselieres_IC2[[as.character(points)]] = NULL
  all_friches_IC2[[as.character(points)]] = NULL
  
  for(repet in 1:30) {
    modtest <- all_param[[as.character(points)]][[paste0("repet_", repet)]]
    
    cultures = ifelse("cultures" %in% colnames(modtest), modtest[2, "cultures"], NA) # vérifie l'existence du paramètre, puis stocke la médiane
    rase = ifelse("rase" %in% colnames(modtest), modtest[2, "rase"], NA)
    haute_fauchee = ifelse("haute fauchée" %in% colnames(modtest), modtest[2, "haute fauchée"], NA)
    haute_non_fauchee = ifelse("haute non fauchée" %in% colnames(modtest), modtest[2, "haute non fauchée"], NA)
    arbustive = ifelse("arbustive" %in% colnames(modtest), modtest[2, "arbustive"], NA)
    roselieres = ifelse("roselières/scirpaies" %in% colnames(modtest), modtest[2, "roselières/scirpaies"], NA)
    friches = ifelse("friches" %in% colnames(modtest), modtest[2, "friches"], NA)
    
    cultures_IC1 = ifelse("cultures" %in% colnames(modtest), modtest[1, "cultures"], NA) # vérifie l'existence du paramètre, puis stocke le quantile 2.5%
    rase_IC1 = ifelse("rase" %in% colnames(modtest), modtest[1, "rase"], NA)
    haute_fauchee_IC1 = ifelse("haute fauchée" %in% colnames(modtest), modtest[1, "haute fauchée"], NA)
    haute_non_fauchee_IC1 = ifelse("haute non fauchée" %in% colnames(modtest), modtest[1, "haute non fauchée"], NA)
    arbustive_IC1 = ifelse("arbustive" %in% colnames(modtest), modtest[1, "arbustive"], NA)
    roselieres_IC1 = ifelse("roselières/scirpaies" %in% colnames(modtest), modtest[1, "roselières/scirpaies"], NA)
    friches_IC1 = ifelse("friches" %in% colnames(modtest), modtest[1, "friches"], NA)
    
    cultures_IC2 = ifelse("cultures" %in% colnames(modtest), modtest[3, "cultures"], NA) # vérifie l'existence du paramètre, puis stocke le quantile 97.5%
    rase_IC2 = ifelse("rase" %in% colnames(modtest), modtest[3, "rase"], NA)
    haute_fauchee_IC2 = ifelse("haute fauchée" %in% colnames(modtest), modtest[3, "haute fauchée"], NA)
    haute_non_fauchee_IC2 = ifelse("haute non fauchée" %in% colnames(modtest), modtest[3, "haute non fauchée"], NA)
    arbustive_IC2 = ifelse("arbustive" %in% colnames(modtest), modtest[3, "arbustive"], NA)
    roselieres_IC2 = ifelse("roselières/scirpaies" %in% colnames(modtest), modtest[3, "roselières/scirpaies"], NA)
    friches_IC2 = ifelse("friches" %in% colnames(modtest), modtest[3, "friches"], NA)
    
    # Ajoute les valeurs aux listes correspondantes
    all_cultures[[as.character(points)]] = c(all_cultures[[as.character(points)]], cultures)
    all_rase[[as.character(points)]] = c(all_rase[[as.character(points)]], rase)
    all_haute_fauchee[[as.character(points)]] = c(all_haute_fauchee[[as.character(points)]], haute_fauchee)
    all_haute_non_fauchee[[as.character(points)]] = c(all_haute_non_fauchee[[as.character(points)]], haute_non_fauchee)
    all_arbustive[[as.character(points)]] = c(all_arbustive[[as.character(points)]], arbustive)
    all_roselieres[[as.character(points)]] = c(all_roselieres[[as.character(points)]], roselieres)
    all_friches[[as.character(points)]] = c(all_friches[[as.character(points)]], friches)
    
    all_cultures_IC1[[as.character(points)]] = c(all_cultures_IC1[[as.character(points)]], cultures_IC1)
    all_rase_IC1[[as.character(points)]] = c(all_rase_IC1[[as.character(points)]], rase_IC1)
    all_haute_fauchee_IC1[[as.character(points)]] = c(all_haute_fauchee_IC1[[as.character(points)]], haute_fauchee_IC1)
    all_haute_non_fauchee_IC1[[as.character(points)]] = c(all_haute_non_fauchee_IC1[[as.character(points)]], haute_non_fauchee_IC1)
    all_arbustive_IC1[[as.character(points)]] = c(all_arbustive_IC1[[as.character(points)]], arbustive_IC1)
    all_roselieres_IC1[[as.character(points)]] = c(all_roselieres_IC1[[as.character(points)]], roselieres_IC1)
    all_friches_IC1[[as.character(points)]] = c(all_friches_IC1[[as.character(points)]], friches_IC1)
    
    all_cultures_IC2[[as.character(points)]] = c(all_cultures_IC2[[as.character(points)]], cultures_IC2)
    all_rase_IC2[[as.character(points)]] = c(all_rase_IC2[[as.character(points)]], rase_IC2)
    all_haute_fauchee_IC2[[as.character(points)]] = c(all_haute_fauchee_IC2[[as.character(points)]], haute_fauchee_IC2)
    all_haute_non_fauchee_IC2[[as.character(points)]] = c(all_haute_non_fauchee_IC2[[as.character(points)]], haute_non_fauchee_IC2)
    all_arbustive_IC2[[as.character(points)]] = c(all_arbustive_IC2[[as.character(points)]], arbustive_IC2)
    all_roselieres_IC2[[as.character(points)]] = c(all_roselieres_IC2[[as.character(points)]], roselieres_IC2)
    all_friches_IC2[[as.character(points)]] = c(all_friches_IC2[[as.character(points)]], friches_IC2)
  }
}


mean_cultures = sapply(all_cultures, function(x) mean(unlist(x), na.rm = TRUE)) # moyenne des médianes
mean_rase = sapply(all_rase, function(x) mean(unlist(x), na.rm = TRUE))
mean_haute_fauchee = sapply(all_haute_fauchee, function(x) mean(unlist(x), na.rm = TRUE))
mean_haute_non_fauchee = sapply(all_haute_non_fauchee, function(x) mean(unlist(x), na.rm = TRUE))
mean_arbustive = sapply(all_arbustive, function(x) mean(unlist(x), na.rm = TRUE))
mean_roselieres = sapply(all_roselieres, function(x) mean(unlist(x), na.rm = TRUE))
mean_friches = sapply(all_friches, function(x) mean(unlist(x), na.rm = TRUE))

mean_cultures_IC1 = sapply(all_cultures_IC1, function(x) mean(unlist(x), na.rm = TRUE)) # moyenne des quantiles 2.5%
mean_rase_IC1 = sapply(all_rase_IC1, function(x) mean(unlist(x), na.rm = TRUE))
mean_haute_fauchee_IC1 = sapply(all_haute_fauchee_IC1, function(x) mean(unlist(x), na.rm = TRUE))
mean_haute_non_fauchee_IC1 = sapply(all_haute_non_fauchee_IC1, function(x) mean(unlist(x), na.rm = TRUE))
mean_arbustive_IC1 = sapply(all_arbustive_IC1, function(x) mean(unlist(x), na.rm = TRUE))
mean_roselieres_IC1 = sapply(all_roselieres_IC1, function(x) mean(unlist(x), na.rm = TRUE))
mean_friches_IC1 = sapply(all_friches_IC1, function(x) mean(unlist(x), na.rm = TRUE))

mean_cultures_IC2 = sapply(all_cultures_IC2, function(x) mean(unlist(x), na.rm = TRUE)) # moyenne des quantiles 97.5%
mean_rase_IC2 = sapply(all_rase_IC2, function(x) mean(unlist(x), na.rm = TRUE))
mean_haute_fauchee_IC2 = sapply(all_haute_fauchee_IC2, function(x) mean(unlist(x), na.rm = TRUE))
mean_haute_non_fauchee_IC2 = sapply(all_haute_non_fauchee_IC2, function(x) mean(unlist(x), na.rm = TRUE))
mean_arbustive_IC2 = sapply(all_arbustive_IC2, function(x) mean(unlist(x), na.rm = TRUE))
mean_roselieres_IC2 = sapply(all_roselieres_IC2, function(x) mean(unlist(x), na.rm = TRUE))
mean_friches_IC2 = sapply(all_friches_IC2, function(x) mean(unlist(x), na.rm = TRUE))
```
Puis, nous créons les dataframes utiles aux représentations graphiques :

```r
points <- c(5, 15, 25, 35, 45, 55, 65, 75, 85)

# Cultures
mean_cultures_df <- data.frame(
  points = points,
  mean_cultures = mean_cultures,
  cultures_IC1 = mean_cultures_IC1,
  cultures_IC2 = mean_cultures_IC2
)
mean_cultures_df = rbind(mean_cultures_df, c(92, -0.7594380, -0.9763430, -0.09041827)) # ajout des valeurs correspondant au modèle avec tous les points d'écoute

# Végétation rase
mean_rase_df <- data.frame(
  points = points,
  mean_rase = mean_rase,
  rase_IC1 = mean_rase_IC1,
  rase_IC2 = mean_rase_IC2
)
mean_rase_df = rbind(mean_rase_df, c(92, -0.48107861, -0.93085301, -0.09041827))

# Végétation haute fauchée
mean_haute_fauchee_df = data.frame(
  points = points,
  mean_haute_fauchee = mean_haute_fauchee,
  haute_fauchee_IC1 = mean_haute_fauchee_IC1,
  haute_fauchee_IC2 = mean_haute_fauchee_IC2
)
mean_haute_fauchee_df = rbind(mean_haute_fauchee_df, c(92, -0.5879129, -0.8701863, -0.3092496))

# Végétation haute non fauchée
mean_haute_non_fauchee_df = data.frame(
  points = points,
  mean_haute_non_fauchee = mean_haute_non_fauchee,
  haute_non_fauchee_IC1 = mean_haute_non_fauchee_IC1,
  haute_non_fauchee_IC2 = mean_haute_non_fauchee_IC2
)
mean_haute_non_fauchee_df = rbind(mean_haute_non_fauchee_df, c(92, -0.04783005, -0.26903250, 0.13225568))

# Végétation arbustive basse
mean_arbustive_df = data.frame(
  points = points,
  mean_arbustive = mean_arbustive,
  arbustive_IC1 = mean_arbustive_IC1,
  arbustive_IC2 = mean_arbustive_IC2
)
mean_arbustive_df = rbind(mean_arbustive_df, c(92, 0.5366907, -0.2274477, 1.2822200))

# Roselières et scirpaies
mean_roselieres_df = data.frame(
  points = points,
  mean_roselieres = mean_roselieres,
  roselieres_IC1 = mean_roselieres_IC1,
  roselieres_IC2 = mean_roselieres_IC2
)
mean_roselieres_df = rbind(mean_roselieres_df, c(92, 0.28671340, -0.03498487, 0.60925905))

# Friches
mean_friches_df = data.frame(
  points = points,
  mean_friches = mean_friches,
  friches_IC1 = mean_friches_IC1,
  friches_IC2 = mean_friches_IC2
)
mean_friches_df = rbind(mean_friches_df, c(92, 1.0737676, 0.6223037, 1.5016847))
```

Et enfin les représentations graphiques :

```r
# Cultures
ggplot(data=mean_cultures_df, aes(x=points, y=mean_cultures, ymin=cultures_IC1, ymax=cultures_IC2)) + 
  geom_line() + 
  geom_ribbon(alpha=0.5)+
  labs(x = "Nombre de points d'écoute", y = "Estimation du paramètre -cultures-")

# Végétation rase
  ggplot(data=mean_rase_df, aes(x=points, y=mean_rase, ymin=rase_IC1, ymax=rase_IC2)) + 
  geom_line() + 
  geom_ribbon(alpha=0.5)+
  labs(x = "Nombre de points d'écoute", y = "Estimation de la végétation rase")
  
# Végétation haute fauchée
  ggplot(data=mean_haute_fauchee_df, aes(x=points, y=mean_haute_fauchee, ymin=haute_fauchee_IC1, ymax=haute_fauchee_IC2)) + 
  geom_line() + 
  geom_ribbon(alpha=0.5)+
  labs(x = "Nombre de points d'écoute", y = "Estimation de la végétation herbacée fauchée")

# Végétation haute non fauchée
  ggplot(data=mean_haute_non_fauchee_df, aes(x=points, y=mean_haute_non_fauchee, ymin=haute_non_fauchee_IC1, ymax=haute_non_fauchee_IC2)) + 
    geom_line() + 
    geom_ribbon(alpha=0.5)+
    labs(x = "Nombre de points d'écoute", y = "Estimation de la végétation herbacée non fauchée")

# Végétation arbustive basse
  ggplot(data=mean_arbustive_df, aes(x=points, y=mean_arbustive, ymin=arbustive_IC1, ymax=arbustive_IC2)) + 
    geom_line() + 
    geom_ribbon(alpha=0.5)+
    labs(x = "Nombre de points d'écoute", y = "Estimation de la végétation arbustive")

# Roselières et scirpaies
  ggplot(data=mean_roselieres_df, aes(x=points, y=mean_roselieres, ymin=roselieres_IC1, ymax=roselieres_IC2)) + 
    geom_line() + 
    geom_ribbon(alpha=0.5)+
    labs(x = "Nombre de points d'écoute", y = "Estimation du paramètre -roselières/scirpaies-")

# Friches
  ggplot(data=mean_friches_df, aes(x=points, y=mean_friches, ymin=friches_IC1, ymax=friches_IC2)) + 
    geom_line() + 
    geom_ribbon(alpha=0.5)+
    labs(x = "Nombre de points d'écoute", y = "Estimation du paramètre -friches-")
```
