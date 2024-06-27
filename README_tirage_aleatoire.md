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

### AUC

### RMSE

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
  
    
ggplot(data=tab_duree, aes(x=points, y=mean_duree, ymin=duree_IC1, ymax=duree_IC2)) + 
  geom_line() + # ligne pour les médianes
  geom_ribbon(alpha=0.5)+ # ruban pour les intervalles de crédibilité
  labs(x = "Nombre de points d'écoute", y = "Estimation de la durée de submersion") # légende
  ```

#### Végétation

Nous allons procéder de la même manière pour tous les paramètres de la végatation.
