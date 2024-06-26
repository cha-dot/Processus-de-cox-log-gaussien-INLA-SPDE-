# Suppression progressive de points d'écoute selon une grille

Cela consiste à définir une grille sur notre domaine d'étude et à ne garder qu'un seul point d'écoute par cellule, tiré aléatoirement.

## Table des matières

- [Packages](#packages)
- [Création de la grille](#création-de-la-grille)
- [Modification du nombre de points d'écoute](#modification-du-nombre-de-points-découte)
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
library(hrbrthemes)
library(dplyr)
library(tidyr)
library(pROC)
```

## Création de la grille

Nous créons une grille qui recouvre le domaine d'étude et dont chaque cellule a des dimensions spécifiées par l'objet `cellule`. Par exemple, `cellule = 200` permet de faire une grille divisée en cellules de dimensions 200x200.

```r
cellule = cel # dimensions des cellules
grille = st_make_grid(st_as_sf(contour_sp), cellsize = cellule) # crée une grille sur le domaine
grille_sf = st_sf(geometry = grille) # conversion du format
contour_sf = st_as_sf(contour_sp) # conversion du format
```

Représentons cette grille.

```r
 ggplot()+
   geom_sf(data = contour_sf, fill = NA, color = "blue") +
   geom_sf(data = grille_sf, fill = NA, color = "black") +
   geom_sf(data = Buffer, alpha = 0.5, color = "red")
```

![Grille complète](https://raw.githubusercontent.com/cha-dot/Processus-de-cox-log-gaussien-INLA-SPDE-/1c8361e0e5c9624d670f4f9609ccf18d63acce59/tirage_regulier_diapo.svg)

## Modification du nombre de points d'écoute

Pour limiter les chevauchements de buffers sur des cellules adjacentes, les buffers seront sélectionnés à partir de leur centroïde. Leur centroïde sont ainsi conservés dans l'objet `Buffer_centroids`.

```r
st_crs(Buffer) = st_crs(grille_sf)
Buffer = st_transform(Buffer, crs = st_crs(grille_sf)) # homogénéisation des crs
Buffer_centroids = st_centroid(Buffer) # centroïdes des buffers (disques d'écoute)
Buffer_centroids_sp = as(Buffer_centroids, "Spatial") # conversion de format
```

Tous les buffers sont attribués à une cellule de la grille. Un seul buffer est alors échantillonné aléatoirement.

```r
grille_buffers <- st_join(grille_sf, Buffer_centroids) # attribue les buffers à chaque cellule
buffer_unique <- grille_buffers %>%
  group_by(geometry) %>% # parmi les buffers dans une cellule
  sample_n(1) %>% # échantillonnage d'un buffer
  ungroup()
```
Nous dressons la liste des buffers sélectionnés. La commande `unique` permet de ne garder qu'un seul exemplaire des buffers sélectionnés (`buffer_unique` est un tableau dont les lignes représentent les cellules de la grille).

```r
liste = unique(buffer_unique$numero_buf[!is.na(buffer_unique$numero_buf)]) 
```

Nous allons à présent uniquement sélectionner les observations qui font partie de ces buffers.

```r
GB_10_1 = GB_PE_eau_fauche[GB_PE_eau_fauche$numero_buf%in%liste,] # garde les observations des buffers sélectionnés
GB_PE_eau_fauche_2 = SpatialPoints(coords = GB_10_1[,c("x_wgs84", "y_wgs84")]) # transformation des données en spatial points
proj4string(GB_PE_eau_fauche_2) = CRS("+proj=longlat +datum=WGS84")
points_lambert = spTransform(GB_PE_eau_fauche_2, CRSobj = CRS(contour_sp@proj4string@projargs)) # homogénéisation des crs
GB_10_1_LAMB = SpatialPointsDataFrame(coords = coordinates(points_lambert), # conversion du format
                                      data = GB_10_1,
                                      proj4string = CRS(contour_sp@proj4string@projargs))
```

On actualise la variable `Buffer` :

```r
coord_Buffer = GB_10_1[, c("X_PE", "Y_PE")] # coordonnées des buffers
coord_Buffer = coord_Buffer[!duplicated(coord_Buffer),] # liste des buffers et leurs coordonnées sans doublons
Buffer_reduit = SpatialPoints(coords = coord_Buffer, CRS("+proj=longlat +datum=WGS84")) # conversion du format
Buffer_reduit = st_as_sf(Buffer_reduit) # conversion du format
Buffer_reduit = st_transform(Buffer_reduit, crs = proj4string(GB_PE_eau_fauche_LAMB))
Buffer_reduit = st_buffer(Buffer_reduit, 150) # 150 = nb de traits qui servent à la formation du cercle
Buffer_reduit_sp = as(Buffer_reduit, "Spatial")
```

Nous pouvons représenter les buffers restant :

```r
 ggplot()+
     geom_sf(data = contour_sf, fill = NA, color = "blue") +
     geom_sf(data = grille_sf, fill = NA, color = "black") +
     geom_sf(data = Buffer_reduit, alpha = 0.5, color = "red")
```

![Grille incomplète](https://raw.githubusercontent.com/cha-dot/Processus-de-cox-log-gaussien-INLA-SPDE-/62253bf9fd400cf4a617092d9d79c5bda22c3131/tirage_regulier_diapo_1800.svg)

## Modèle

Le reste du modèle ne change pas (pour plus d'explications du modèle, se référer à README_INLA.md).

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

Si le modèle est répété pour différentes résolutions de la grille, il est possible de sauvegarder certains objets, tels que le summary, l'AUC ou encore les distributions a posteriori pour une comparaison de résultats.

Si le modèle est répété sur `cel` différentes résolutions de grille et `rep` fois pour détenir des réplicats de chaque cas de figure, alors il est possible de sauvegarder le summary de cette manière :

```r
  length(liste) # nb de points d'écoute
  chemin = "~/summaries_reg_nouveaux"
  fichier = paste0("nouv_summary_", cel, "_nPE_", length(liste), "_rep_", repetition, ".RData")
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

Il est possible de sauvegarder ces 3 dernières métriques :

```r
  chemin = "~/metriques_reg_nouveaux"
  fichier = paste0("nouv_metriques_", cel, "_nPE_", length(liste), "_rep_", repetition, ".RData")
  nom_complet = file.path(chemin, fichier)
  save(AUC, mse, RMSE, file = nom_complet)
```

Toujours dans l'optique de gain de temps si le code doit être répliqué de nombreuses fois, il est possible de réduire le nombre d'échantillons `Nrep` du modèle (`Nrep = 500` semble être un bon compromis entre le temps de calcul et la précision des résultats).

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

`colnames(post) = c(rownames(pr.int.tot[[1]]$latent)[ind], "intensite_sub", "duree_sub")` à la place de `colnames(post) = c("cultures", "rase", "haute fauchée", "haute non fauchée", "arbustive", "roselières/scirpaies", "friches", "intensite_sub", "duree_sub") permet de prendre en compte la possibilité d'une absence d'un type de végétation.

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
  chemin = "~/post_reg_nouveaux"
  fichier = paste0("reg_post_", cel, "_nPE_", length(liste), "_rep_", repetition, ".RData")
  nom_complet = file.path(chemin, fichier)
  save(post.stat.veg, post.stat.duree, post.stat.max, file = nom_complet)
```
## Représentations graphiques des résultats

Le modèle a été réalisé pour plusieurs résolutions de grille (cellules de dimensions 300x300 m² à 2500x2500 m², par pas de 100 m²). Chaque cas de figure a été répliqué 10 fois, puisque le buffer relevé par cellule est échantillonné aléatoirement. Ainsi, nous ferons la moyenne des métriques pour chaque cas de figure.

### AUC

Les commandes qui permettront de charger les bons fichiers :

```r
n_cell = seq(300, 2500, by = 100) # résolution de la grille
n_pe = c(92, 81, 69, 55, 48, 44, 37, 34, 31, 25, 27, 25, 18, 20, 22, 16, 15, 15, 14, 17, 15, 14, 12) # nb de points d'écoute (dans le nom des fichiers)
n_rep = 10 # nb de réplicats

generate_filename <- function(n_cell, n_pe, rep) { # génère le nom de fichier
  sprintf("~/nouv_metriques_%d_nPE_%d_rep_%d.RData", n_cell, n_pe, rep)
}
```

A présent, nous allons charger, regrouper les données d'AUC, et faire la moyenne pour chaque cas de figure :

```r
all_AUC_by_combination <- list() # vecteur qui stockera les AUC

# Boucle pour ouvrir chaque fichier, extraire les AUC et les stocker dans la liste
for (i in seq_along(n_cell)) { # pour chaque résolution de grille
  cell <- n_cell[i] # cas de figure i (dimensions des cellules)
  pe <- n_pe[i] # nb de points d'écoute associés au cas de figure i
  indice <- paste(cell, pe, sep = "_") # ex : "500_69"
  all_AUC_by_combination[[indice]] <- list() # permettra de regrouper par cas de figure
  for (rep in seq_len(n_rep)) { # pour chaque réplicat (x 10)
    filename <- generate_filename(cell, pe, rep) # charge le fichier
    if (file.exists(filename)) {
      load(filename)
      # Vérifier l'existence de la variable 'AUC'
      if (exists("AUC")) {
        all_AUC_by_combination[[indice]] <- c(all_AUC_by_combination[[indice]], AUC) # stocke l'AUC
      } else {
        warning(sprintf("La variable 'AUC' n'existe pas dans le fichier %s.", filename))
      }
    } else {
      warning(sprintf("Le fichier %s n'existe pas.", filename))
    }
  }
}

# Initialiser une liste pour stocker les moyennes des AUCs
moy_AUC <- list()

for (indice in names(all_AUC_by_combination)) { # pour chaque cas de figure
  if (length(all_AUC_by_combination[[indice]]) > 0) { # vérifie l'existence de la valeur
    moy_AUC[[indice]] <- mean(unlist(all_AUC_by_combination[[indice]])) # moyenne AUC
  } else {
    moy_AUC[[indice]] <- NA
  }
}
```

Passons à la représentation graphique. L'axe des abscisses est inversé. En effet, plus les cellules sont grandes, moins il y a de points d'écoute. Je trouvais plus cohérent de représenter les valeurs des métriques en fonction du nombre croissant de points d'écoute (permettant une comparaison plus facile avec le tirage aléatoire).

```r
AUC_reg_final = unlist(moy_AUC) # format vecteur
AUC_reg_final = as.data.frame(cbind(n_cell, AUC_reg_final)) # tableau des AUC moyennes pour chaque résolution de grille
colnames(AUC_reg_final) = c("cellules", "AUC") # nom des colonnes

ggplot(AUC_reg_final, aes(x = cellules, y = AUC))+ # AUC en fonction de la résolution de la grille
  geom_line()+
  scale_y_continuous(limits = c(0,1))+ # limites de l'axe y
  scale_x_reverse()+ # inverse le sens de l'axe x (facultatif)
  labs(x = "Taille des cellules (en m²)", y = "AUC")+ # légende
  theme_minimal()
```

![Régulier AUC](https://github.com/cha-dot/Processus-de-cox-log-gaussien-INLA-SPDE-/blob/images/AUC_reg_finale.jpg?raw=true)

### RMSE

Les commandes qui permettront de charger les bons fichiers :

```r
n_cell = seq(300, 2500, by = 100) # résolution de la grille
n_pe = c(92, 81, 69, 55, 48, 44, 37, 34, 31, 25, 27, 25, 18, 20, 22, 16, 15, 15, 14, 17, 15, 14, 12) # nb de points d'écoute (dans le nom des fichiers)
n_rep = 10 # nb de réplicats

generate_filename <- function(n_cell, n_pe, rep) { # génère le nom de fichier
  sprintf("~/nouv_metriques_%d_nPE_%d_rep_%d.RData", n_cell, n_pe, rep)
}
```

Le RMSE est une métrique utile pour comparer des modèles entre eux (plus il est faible, meilleur est le modèle). Nous procédons exactement de la même manière que l'AUC pour sa représentation graphique.

```r
all_RMSE_by_combination <- list() # vecteur qui stockera les RMSE

# Boucle pour ouvrir chaque fichier, extraire les AUC et les stocker dans la liste
for (i in seq_along(n_cell)) { # pour chaque résolution de grille
  cell <- n_cell[i] # cas de figure i (dimensions des cellules)
  pe <- n_pe[i] # nb de points d'écoute associé au cas de figure i
  indice <- paste(cell, pe, sep = "_") # ex : "500_69"
  all_RMSE_by_combination[[indice]] <- list() # permettra de regrouper par cas de figure
  for (rep in seq_len(n_rep)) { # pour chaque réplicat (x 10)
    filename <- generate_filename(cell, pe, rep) # charge le fichier
    if (file.exists(filename)) {
      load(filename)
      # Vérifier l'existence de la variable 'RMSE'
      if (exists("RMSE")) {
        all_RMSE_by_combination[[indice]] <- c(all_RMSE_by_combination[[indice]], RMSE) # stocke le RMSE
      } else {
        warning(sprintf("La variable 'RMSE' n'existe pas dans le fichier %s.", filename))
      }
    } else {
      warning(sprintf("Le fichier %s n'existe pas.", filename))
    }
  }
}

# Initialiser une liste pour stocker les moyennes des RMSE
moy_RMSE <- list()

for (indice in names(all_RMSE_by_combination)) { # pour chaque cas de figure
  if (length(all_RMSE_by_combination[[indice]]) > 0) { # vérifier l'existence de la valeur
    moy_RMSE[[indice]] <- mean(unlist(all_RMSE_by_combination[[indice]])) # moyenne RMSE
  } else {
    moy_RMSE[[indice]] <- NA
  }
}
```

Passons à la représentation graphique :

```r
RMSE_reg_final = unlist(moy_RMSE) # format vecteur
RMSE_reg_final = as.data.frame(cbind(n_cell, RMSE_reg_final)) # tableau des RMSE moyens pour chaque résolution de grille
colnames(RMSE_reg_final) = c("cellules", "RMSE") # nom des colonnes

ggplot(RMSE_reg_final, aes(x = cellules, y = RMSE))+ # RMSE en fonction de la résolution de la grille
  geom_line()+
  scale_y_continuous(limits = c(0,1.5))+ # limites de l'axe y
  scale_x_reverse()+ # inverse l'axe des x (facultatif)
  labs(x = "Taille des cellules (en m²)", y = "RMSE")+ # légende
  theme_minimal()
```

![Régulier RMSE](https://github.com/cha-dot/Processus-de-cox-log-gaussien-INLA-SPDE-/blob/images/RMSE_reg_finale.jpg?raw=true)

### Paramètres du modèle

Concernant les paramètres du modèle, la médiane et les intervalles de crédibilité moyens de chaque cas de figure seront représentés en fonction de la résolution de la grille.

Les fichiers contenant les données sur les paramètres du modèle peuvent être chargés de cette manière :

```r
n_cell = seq(300, 2500, by = 100)
n_pe = c(92, 81, 69, 55, 48, 44, 37, 34, 31, 25, 27, 25, 18, 20, 22, 16, 15, 15, 14, 17, 15, 14, 12)
n_rep = 10

generate_filename <- function(n_cell, n_pe, rep) { # génère le nom de fichier
  sprintf("C:/Users/Charlotte.Marques/OneDrive - LPO/Documents/Données/Resultats_Florian/post_reg_nouveaux/reg_post_%d_nPE_%d_rep_%d.RData", n_cell, n_pe, rep)
}
```

La méthode utilisée pour représenter les données est similaire à celle utilisée pour l'AUC et le RMSE.

#### Intensité de submersion

Calcul de la valeur moyenne du paramètre et de ses intervalles de crédibilité, pour chaque cas de figure (résolution de la grille) :

```r
all_max_by_combination <- list() vecteur qui stockera les intensités de submersion

# Boucle pour ouvrir chaque fichier, extraire les quantiles et calculer les moyennes
for (i in seq_along(n_cell)) { # pour chaque résolution de grille
  cell <- n_cell[i] # cas de figure i (dimensions des cellules)
  pe <- n_pe[i] # nb de points d'écoute associé au cas de figure i
  indice <- paste(cell, pe, sep = "_") # ex : "500_69"
  temp_quantiles_Q2.5 <- c() # stockera les quantiles 2.5%
  temp_quantiles_Q50 <- c() # stockera les médianes 
  temp_quantiles_Q97.5 <- c() # stockera les quantiles 97.5%
  
  for (rep in seq_len(n_rep)) { # pour chaque réplicat (x 10)
    filename <- generate_filename(cell, pe, rep) # charge le fichier
    if (file.exists(filename)) {
      load(filename)
      # Vérifier l'existence des variables de quantiles
      if (exists("post.stat.max")) {
        # Stocker les résultats dans les vecteurs temporaires
        temp_quantiles_Q2.5 <- c(temp_quantiles_Q2.5, post.stat.max["2.5%"])
        temp_quantiles_Q50 <- c(temp_quantiles_Q50, post.stat.max["50%"])
        temp_quantiles_Q97.5 <- c(temp_quantiles_Q97.5, post.stat.max["97.5%"])
      } else {
        warning(sprintf("Les quantiles n'existent pas dans le fichier %s.", filename))
      }
    } else {
      warning(sprintf("Le fichier %s n'existe pas.", filename))
    }
  }
  
  # Calculer les moyennes des quantiles
  if (length(temp_quantiles_Q2.5) > 0 & length(temp_quantiles_Q50) > 0 & length(temp_quantiles_Q97.5) > 0) { # vérifie l'existence des quantiles
    all_max_by_combination[[indice]] <- c( # pour chaque cas de figure
      moy_Q2.5 = mean(temp_quantiles_Q2.5), # moyenne du quantile 2.5%
      moy_Q50 = mean(temp_quantiles_Q50), # moyenne de la médiane
      moy_Q97.5 = mean(temp_quantiles_Q97.5) # moyenne du quantile 97.5%
    )
  } else {
    all_max_by_combination[[indice]] <- NA
  }
}
````

Représentation graphique :

```r
# Créer un dataframe à partir de la liste de moyennes des quantiles
all_max_by_combination_df <- do.call(rbind, lapply(all_max_by_combination, as.data.frame))
all_max_by_combination_df = data.frame(t(sapply(all_max_by_combination, unlist)))
all_max_by_combination_df = cbind(all_max_by_combination_df, n_cell) # ajout de la colonne des dimensions des cellules de la grille
# "sapply" convertit chaque élément de la liste en un vecteur et "t" transpose pour que chaque ligne corresponde à une combinaison
rownames(all_max_by_combination_df) = NULL # suppression des noms des lignes

ggplot(all_max_by_combination_df, aes(x = n_cell, y = moy_Q50))+
  geom_line()+ # ligne pour les médianes
  geom_ribbon(aes(ymin = moy_Q2.5, ymax = moy_Q97.5), alpha=0.5)+ # ruban pour les intervalles de crédibilité
  labs(x = "Taille des cellules (en m²)", y = "Estimation de l'intensité de submersion")+ # légende
  scale_y_continuous(limits = c(-1,5))+ # limites de l'axe y
  scale_x_reverse()+ # inverse le sens de l'axe x (facultatif)
  theme_minimal()
```

![Régulier intensité de submersion](https://github.com/cha-dot/Processus-de-cox-log-gaussien-INLA-SPDE-/blob/images/max_reg_pte_final.jpg?raw=true)

#### Durée de submersion

Calcul de la valeur moyenne du paramètre et de ses intervalles de crédibilité, pour chaque cas de figure (résolution de la grille) :

```r
all_duree_by_combination <- list() vecteur qui stockera les durées de submersion

# Boucle pour ouvrir chaque fichier, extraire les quantiles et calculer les moyennes
for (i in seq_along(n_cell)) { # pour chaque résolution de grille
  cell <- n_cell[i] # cas de figure i (dimensions des cellules)
  pe <- n_pe[i] # nb de points d'écoute associé au cas de figure i
  indice <- paste(cell, pe, sep = "_") # ex : "500_69"
  temp_quantiles_Q2.5 <- c() # stockera les quantiles 2.5%
  temp_quantiles_Q50 <- c() # stockera les médianes
  temp_quantiles_Q97.5 <- c() # stockera les quantiles 97.5%
  
  for (rep in seq_len(n_rep)) { # pour chaque réplicat (x 10)
    filename <- generate_filename(cell, pe, rep) # charge le fichier
    if (file.exists(filename)) {
      load(filename)
      # Vérifier l'existence des variables de quantiles
      if (exists("post.stat.duree")) {
        # Stocker les résultats dans les vecteurs temporaires
        temp_quantiles_Q2.5 <- c(temp_quantiles_Q2.5, post.stat.duree["2.5%"])
        temp_quantiles_Q50 <- c(temp_quantiles_Q50, post.stat.duree["50%"])
        temp_quantiles_Q97.5 <- c(temp_quantiles_Q97.5, post.stat.duree["97.5%"])
      } else {
        warning(sprintf("Les quantiles n'existent pas dans le fichier %s.", filename))
      }
    } else {
      warning(sprintf("Le fichier %s n'existe pas.", filename))
    }
  }
  
  # Calculer les moyennes des quantiles
  if (length(temp_quantiles_Q2.5) > 0 & length(temp_quantiles_Q50) > 0 & length(temp_quantiles_Q97.5) > 0) { # vérifie l'existence des quantiles
    all_duree_by_combination[[indice]] <- c( # pour chaque cas de figure
      moy_Q2.5 = mean(temp_quantiles_Q2.5), # moyenne du quantile 2.5%
      moy_Q50 = mean(temp_quantiles_Q50), # moyenne de la médiane
      moy_Q97.5 = mean(temp_quantiles_Q97.5) # moyenne du quantile 97.5%
    )
  } else {
    all_duree_by_combination[[indice]] <- NA
  }
}

```
Représentation graphique :

```r
# Créer un dataframe à partir de la liste de moyennes des quantiles
all_duree_by_combination_df <- do.call(rbind, lapply(all_duree_by_combination, as.data.frame))
all_duree_by_combination_df = data.frame(t(sapply(all_duree_by_combination, unlist)))
all_duree_by_combination_df = cbind(all_duree_by_combination_df, n_cell) # ajout de la colonne des dimensions des cellules de la grille
# "sapply" convertit chaque élément de la liste en un vecteur et "t" transpose pour que chaque ligne corresponde à une combinaison
rownames(all_duree_by_combination_df) = NULL # suppression des noms des lignes

ggplot(all_duree_by_combination_df[-c(11:23),], aes(x = n_cell, y = moy_Q50))+
  geom_line()+ # ligne pour les médianes
  geom_ribbon(aes(ymin = moy_Q2.5, ymax = moy_Q97.5), alpha=0.5)+ # ruban pour les intervalles de crédibilité
  labs(x = "Taille des cellules (en m²)", y = "Estimation de la durée de submersion")+ # légende
  scale_y_continuous(limits = c(-0.0075,0.0075))+ # limites de l'axe y
  scale_x_reverse(breaks = all_duree_by_combination_df$n_cell)+ # inverse le sens de l'axe x (facultatif)
  theme_minimal()
```

![Régulier durée de submersion](https://github.com/cha-dot/Processus-de-cox-log-gaussien-INLA-SPDE-/blob/images/duree_reg_pte_final.jpg?raw=true)

#### Végétation

##### Cultures

Calcul de la valeur moyenne du paramètre et de ses intervalles de crédibilité, pour chaque cas de figure (résolution de la grille) :

```r
all_cultures_by_combination <- list() # vecteur qui stockera les paramètres de la modalité "cultures"

# Boucle pour ouvrir chaque fichier, extraire les quantiles et calculer les moyennes
for (i in seq_along(n_cell)) { # pour chaque résolution de grille
  cell <- n_cell[i] # cas de figure i (dimensions des cellules)
  pe <- n_pe[i] # nb de points d'écoute associé au cas de figure i
  indice <- paste(cell, pe, sep = "_") # ex : "500_69"
  temp_quantiles_Q2.5 <- c() # stockera les quantiles 2.5%
  temp_quantiles_Q50 <- c() # stockera les médianes
  temp_quantiles_Q97.5 <- c() # stockera les quantiles 97.5%
  
  for (rep in seq_len(n_rep)) { pour chaque réplicat (x 10)
    filename <- generate_filename(cell, pe, rep) # charge le fichier
    if (file.exists(filename)) {
      load(filename)
      # Vérifier l'existence des variables de quantiles
      if (exists("post.stat.veg")) {
        # Stocker les résultats dans les vecteurs temporaires
        temp_quantiles_Q2.5 <- c(temp_quantiles_Q2.5, post.stat.veg["2.5%", "cultures"])
        temp_quantiles_Q50 <- c(temp_quantiles_Q50, post.stat.veg["50%", "cultures"])
        temp_quantiles_Q97.5 <- c(temp_quantiles_Q97.5, post.stat.veg["97.5%", "cultures"])
      } else {
        warning(sprintf("Les quantiles n'existent pas dans le fichier %s.", filename))
      }
    } else {
      warning(sprintf("Le fichier %s n'existe pas.", filename))
    }
  }
  
  # Calculer les moyennes des quantiles
  if (length(temp_quantiles_Q2.5) > 0 & length(temp_quantiles_Q50) > 0 & length(temp_quantiles_Q97.5) > 0) { # vérifie l'existence des quantiles
    all_cultures_by_combination[[indice]] <- c( # pour chaque cas de figure
      moy_Q2.5 = mean(temp_quantiles_Q2.5), # moyenne du quantile 2.5%
      moy_Q50 = mean(temp_quantiles_Q50), # moyenne de la médiane
      moy_Q97.5 = mean(temp_quantiles_Q97.5) # moyenne du quantile 97.5%
    )
  } else {
    all_cultures_by_combination[[indice]] <- NA
  }
}
```

Représentation graphique :

```r
# Créer un dataframe à partir de la liste de moyennes des quantiles
all_cultures_by_combination_df <- do.call(rbind, lapply(all_cultures_by_combination, as.data.frame))
all_cultures_by_combination_df = data.frame(t(sapply(all_cultures_by_combination, unlist)))
all_cultures_by_combination_df = cbind(all_cultures_by_combination_df, n_cell) # ajout de la colonne des dimensions des cellules de la grille
rownames(all_cultures_by_combination_df) = NULL # suppression des noms des lignes

ggplot(all_cultures_by_combination_df[-c(11:23),], aes(x = n_cell, y = moy_Q50))+
  geom_line()+ # ligne pour les médianes
  geom_ribbon(aes(ymin = moy_Q2.5, ymax = moy_Q97.5), alpha=0.5)+ # ruban pour les intervalles de crédibilité
  labs(x = "Taille des cellules (en m²)", y = "Estimation du paramètre -cultures-")+ # légende
  #scale_y_continuous(limits = c(-2,1))+ # limites de l'axe y
  scale_x_reverse(breaks = all_cultures_by_combination_df$n_cell)+ # inverse le sens de l'axe x (facultatif)
  theme_minimal()
```

![Régulier cultures](https://github.com/cha-dot/Processus-de-cox-log-gaussien-INLA-SPDE-/blob/images/cultures_reg_pte_final.jpg?raw=true)

##### Végétation rase

Calcul de la valeur moyenne du paramètre et de ses intervalles de crédibilité, pour chaque cas de figure (résolution de la grille) :

```r
all_rase_by_combination <- list()

# Boucle pour ouvrir chaque fichier, extraire les quantiles et calculer les moyennes
for (i in seq_along(n_cell)) {
  cell <- n_cell[i]
  pe <- n_pe[i]
  indice <- paste(cell, pe, sep = "_")
  temp_quantiles_Q2.5 <- c()
  temp_quantiles_Q50 <- c()
  temp_quantiles_Q97.5 <- c()
  
  for (rep in seq_len(n_rep)) {
    filename <- generate_filename(cell, pe, rep)
    if (file.exists(filename)) {
      load(filename)
      # Vérifier l'existence des variables de quantiles
      if (exists("post.stat.veg")) {
        # Ajouter les résultats aux vecteurs temporaires
        temp_quantiles_Q2.5 <- c(temp_quantiles_Q2.5, post.stat.veg["2.5%", "rase"])
        temp_quantiles_Q50 <- c(temp_quantiles_Q50, post.stat.veg["50%", "rase"])
        temp_quantiles_Q97.5 <- c(temp_quantiles_Q97.5, post.stat.veg["97.5%", "rase"])
      } else {
        warning(sprintf("Les quantiles n'existent pas dans le fichier %s.", filename))
      }
    } else {
      warning(sprintf("Le fichier %s n'existe pas.", filename))
    }
  }
  
  # Calculer les moyennes des quantiles pour cette combinaison
  if (length(temp_quantiles_Q2.5) > 0 & length(temp_quantiles_Q50) > 0 & length(temp_quantiles_Q97.5) > 0) {
    all_rase_by_combination[[indice]] <- c(
      moy_Q2.5 = mean(temp_quantiles_Q2.5),
      moy_Q50 = mean(temp_quantiles_Q50),
      moy_Q97.5 = mean(temp_quantiles_Q97.5)
    )
  } else {
    all_rase_by_combination[[indice]] <- NA
  }
}
```

Représentation graphique :

```r
# Créer un dataframe à partir de la liste de moyennes des quantiles
all_rase_by_combination_df <- do.call(rbind, lapply(all_rase_by_combination, as.data.frame))
all_rase_by_combination_df = data.frame(t(sapply(all_rase_by_combination, unlist)))
all_rase_by_combination_df = cbind(all_rase_by_combination_df, n_cell)
rownames(all_rase_by_combination_df) = NULL

ggplot(all_rase_by_combination_df[-c(11:23),], aes(x = n_cell, y = moy_Q50))+
  geom_line()+
  geom_ribbon(aes(ymin = moy_Q2.5, ymax = moy_Q97.5), alpha=0.5)+
  labs(x = "Taille des cellules (en m²)", y = "Estimation de la végétation rase")+
  #scale_y_continuous(limits = c(-0.01,0.01))+
  scale_x_reverse(breaks = all_rase_by_combination_df$n_cell)+
  theme_minimal()
```

![Régulier végétation rase](https://github.com/cha-dot/Processus-de-cox-log-gaussien-INLA-SPDE-/blob/images/rase_reg_pte_final.jpg?raw=true)

##### Végétation haute fauchée

Calcul de la valeur moyenne du paramètre et de ses intervalles de crédibilité, pour chaque cas de figure (résolution de la grille) :

```r
all_fauchee_by_combination <- list()

# Boucle pour ouvrir chaque fichier, extraire les quantiles et calculer les moyennes
for (i in seq_along(n_cell)) {
  cell <- n_cell[i]
  pe <- n_pe[i]
  indice <- paste(cell, pe, sep = "_")
  temp_quantiles_Q2.5 <- c()
  temp_quantiles_Q50 <- c()
  temp_quantiles_Q97.5 <- c()
  
  for (rep in seq_len(n_rep)) {
    filename <- generate_filename(cell, pe, rep)
    if (file.exists(filename)) {
      load(filename)
      # Vérifier l'existence des variables de quantiles
      if (exists("post.stat.veg")) {
        # Ajouter les résultats aux vecteurs temporaires
        temp_quantiles_Q2.5 <- c(temp_quantiles_Q2.5, post.stat.veg["2.5%", "haute fauchée"])
        temp_quantiles_Q50 <- c(temp_quantiles_Q50, post.stat.veg["50%", "haute fauchée"])
        temp_quantiles_Q97.5 <- c(temp_quantiles_Q97.5, post.stat.veg["97.5%", "haute fauchée"])
      } else {
        warning(sprintf("Les quantiles n'existent pas dans le fichier %s.", filename))
      }
    } else {
      warning(sprintf("Le fichier %s n'existe pas.", filename))
    }
  }
  
  # Calculer les moyennes des quantiles pour cette combinaison
  if (length(temp_quantiles_Q2.5) > 0 & length(temp_quantiles_Q50) > 0 & length(temp_quantiles_Q97.5) > 0) {
    all_fauchee_by_combination[[indice]] <- c(
      moy_Q2.5 = mean(temp_quantiles_Q2.5),
      moy_Q50 = mean(temp_quantiles_Q50),
      moy_Q97.5 = mean(temp_quantiles_Q97.5)
    )
  } else {
    all_fauchee_by_combination[[indice]] <- NA
  }
}
```

Représentation graphique :

```r
# Créer un dataframe à partir de la liste de moyennes des quantiles
all_fauchee_by_combination_df <- do.call(rbind, lapply(all_fauchee_by_combination, as.data.frame))
all_fauchee_by_combination_df = data.frame(t(sapply(all_fauchee_by_combination, unlist)))
all_fauchee_by_combination_df = cbind(all_fauchee_by_combination_df, n_cell)
rownames(all_fauchee_by_combination_df) = NULL

ggplot(all_fauchee_by_combination_df[-c(11:23),], aes(x = n_cell, y = moy_Q50))+
  geom_line()+
  geom_ribbon(aes(ymin = moy_Q2.5, ymax = moy_Q97.5), alpha=0.5)+
  labs(x = "Taille des cellules (en m²)", y = "Estimation de la végétation haute fauchée")+
  #scale_y_continuous(limits = c(-0.01,0.01))+
  scale_x_reverse(breaks = all_fauchee_by_combination_df$n_cell)+
  theme_minimal()
```

![Régulier végétation haute fauchée](https://github.com/cha-dot/Processus-de-cox-log-gaussien-INLA-SPDE-/blob/images/fauchee_reg_pte_final.jpg?raw=true)

##### Végétation haute non fauchée

Calcul de la valeur moyenne du paramètre et de ses intervalles de crédibilité, pour chaque cas de figure (résolution de la grille) :


```r
all_non_fauchee_by_combination <- list()

# Boucle pour ouvrir chaque fichier, extraire les quantiles et calculer les moyennes
for (i in seq_along(n_cell)) {
  cell <- n_cell[i]
  pe <- n_pe[i]
  indice <- paste(cell, pe, sep = "_")
  temp_quantiles_Q2.5 <- c()
  temp_quantiles_Q50 <- c()
  temp_quantiles_Q97.5 <- c()
  
  for (rep in seq_len(n_rep)) {
    filename <- generate_filename(cell, pe, rep)
    if (file.exists(filename)) {
      load(filename)
      # Vérifier l'existence des variables de quantiles
      if (exists("post.stat.veg")) {
        # Ajouter les résultats aux vecteurs temporaires
        temp_quantiles_Q2.5 <- c(temp_quantiles_Q2.5, post.stat.veg["2.5%", "haute non fauchée"])
        temp_quantiles_Q50 <- c(temp_quantiles_Q50, post.stat.veg["50%", "haute non fauchée"])
        temp_quantiles_Q97.5 <- c(temp_quantiles_Q97.5, post.stat.veg["97.5%", "haute non fauchée"])
      } else {
        warning(sprintf("Les quantiles n'existent pas dans le fichier %s.", filename))
      }
    } else {
      warning(sprintf("Le fichier %s n'existe pas.", filename))
    }
  }
  
  # Calculer les moyennes des quantiles pour cette combinaison
  if (length(temp_quantiles_Q2.5) > 0 & length(temp_quantiles_Q50) > 0 & length(temp_quantiles_Q97.5) > 0) {
    all_non_fauchee_by_combination[[indice]] <- c(
      moy_Q2.5 = mean(temp_quantiles_Q2.5),
      moy_Q50 = mean(temp_quantiles_Q50),
      moy_Q97.5 = mean(temp_quantiles_Q97.5)
    )
  } else {
    all_non_fauchee_by_combination[[indice]] <- NA
  }
}
```

Représentation graphique :

```r
# Créer un dataframe à partir de la liste de moyennes des quantiles
all_non_fauchee_by_combination_df <- do.call(rbind, lapply(all_non_fauchee_by_combination, as.data.frame))
all_non_fauchee_by_combination_df = data.frame(t(sapply(all_non_fauchee_by_combination, unlist)))
all_non_fauchee_by_combination_df = cbind(all_non_fauchee_by_combination_df, n_cell)
rownames(all_non_fauchee_by_combination_df) = NULL

ggplot(all_non_fauchee_by_combination_df[-c(11:23),], aes(x = n_cell, y = moy_Q50))+
  geom_line()+
  geom_ribbon(aes(ymin = moy_Q2.5, ymax = moy_Q97.5), alpha=0.5)+
  labs(x = "Taille des cellules (en m²)", y = "Estimation de la végétation haute non fauchée")+
  #scale_y_continuous(limits = c(-0.01,0.01))+
  scale_x_reverse(breaks = all_non_fauchee_by_combination_df$n_cell)+
  theme_minimal()
```

![Régulier végétation haute non fauchée](https://github.com/cha-dot/Processus-de-cox-log-gaussien-INLA-SPDE-/blob/images/non_fauchee_reg_pte_final.jpg?raw=true)

##### Végétation arbustive

Calcul de la valeur moyenne du paramètre et de ses intervalles de crédibilité, pour chaque cas de figure (résolution de la grille) :

```r
all_arbustive_by_combination <- list()

# Boucle pour ouvrir chaque fichier, extraire les quantiles et calculer les moyennes
for (i in seq_along(n_cell)) {
  cell <- n_cell[i]
  pe <- n_pe[i]
  indice <- paste(cell, pe, sep = "_")
  temp_quantiles_Q2.5 <- c()
  temp_quantiles_Q50 <- c()
  temp_quantiles_Q97.5 <- c()
  
  for (rep in seq_len(n_rep)) {
    filename <- generate_filename(cell, pe, rep)
    if (file.exists(filename)) {
      load(filename)
      # Vérifier l'existence des variables de quantiles
      if (exists("post.stat.veg")) {
        # Ajouter les résultats aux vecteurs temporaires
        temp_quantiles_Q2.5 <- c(temp_quantiles_Q2.5, post.stat.veg["2.5%", "arbustive"])
        temp_quantiles_Q50 <- c(temp_quantiles_Q50, post.stat.veg["50%", "arbustive"])
        temp_quantiles_Q97.5 <- c(temp_quantiles_Q97.5, post.stat.veg["97.5%", "arbustive"])
      } else {
        warning(sprintf("Les quantiles n'existent pas dans le fichier %s.", filename))
      }
    } else {
      warning(sprintf("Le fichier %s n'existe pas.", filename))
    }
  }
  
  # Calculer les moyennes des quantiles pour cette combinaison
  if (length(temp_quantiles_Q2.5) > 0 & length(temp_quantiles_Q50) > 0 & length(temp_quantiles_Q97.5) > 0) {
    all_arbustive_by_combination[[indice]] <- c(
      moy_Q2.5 = mean(temp_quantiles_Q2.5),
      moy_Q50 = mean(temp_quantiles_Q50),
      moy_Q97.5 = mean(temp_quantiles_Q97.5)
    )
  } else {
    all_arbustive_by_combination[[indice]] <- NA
  }
}
```

Représentation graphique :

```r
# Créer un dataframe à partir de la liste de moyennes des quantiles
all_arbustive_by_combination_df <- do.call(rbind, lapply(all_arbustive_by_combination, as.data.frame))
all_arbustive_by_combination_df = data.frame(t(sapply(all_arbustive_by_combination, unlist)))
all_arbustive_by_combination_df = cbind(all_arbustive_by_combination_df, n_cell)
rownames(all_arbustive_by_combination_df) = NULL

ggplot(all_arbustive_by_combination_df[-c(11:23),], aes(x = n_cell, y = moy_Q50))+
  geom_line()+
  geom_ribbon(aes(ymin = moy_Q2.5, ymax = moy_Q97.5), alpha=0.5)+
  labs(x = "Taille des cellules (en m²)", y = "Estimation de la végétation arbustive")+
  #scale_y_continuous(limits = c(-0.01,0.01))+
  scale_x_reverse(breaks = all_arbustive_by_combination_df$n_cell)+
  theme_minimal()
```

![Régulier végétation arbustive basse](https://github.com/cha-dot/Processus-de-cox-log-gaussien-INLA-SPDE-/blob/images/arbustive_reg_pte_final.jpg?raw=true)

##### Roselières/scirpaies

Calcul de la valeur moyenne du paramètre et de ses intervalles de crédibilité, pour chaque cas de figure (résolution de la grille) :

```r
all_roselieres_by_combination <- list()

# Boucle pour ouvrir chaque fichier, extraire les quantiles et calculer les moyennes
for (i in seq_along(n_cell)) {
  cell <- n_cell[i]
  pe <- n_pe[i]
  indice <- paste(cell, pe, sep = "_")
  temp_quantiles_Q2.5 <- c()
  temp_quantiles_Q50 <- c()
  temp_quantiles_Q97.5 <- c()
  
  for (rep in seq_len(n_rep)) {
    filename <- generate_filename(cell, pe, rep)
    if (file.exists(filename)) {
      load(filename)
      # Vérifier l'existence des variables de quantiles
      if (exists("post.stat.veg")) {
        # Ajouter les résultats aux vecteurs temporaires
        temp_quantiles_Q2.5 <- c(temp_quantiles_Q2.5, post.stat.veg["2.5%", "roselières/scirpaies"])
        temp_quantiles_Q50 <- c(temp_quantiles_Q50, post.stat.veg["50%", "roselières/scirpaies"])
        temp_quantiles_Q97.5 <- c(temp_quantiles_Q97.5, post.stat.veg["97.5%", "roselières/scirpaies"])
      } else {
        warning(sprintf("Les quantiles n'existent pas dans le fichier %s.", filename))
      }
    } else {
      warning(sprintf("Le fichier %s n'existe pas.", filename))
    }
  }
  
  # Calculer les moyennes des quantiles pour cette combinaison
  if (length(temp_quantiles_Q2.5) > 0 & length(temp_quantiles_Q50) > 0 & length(temp_quantiles_Q97.5) > 0) {
    all_roselieres_by_combination[[indice]] <- c(
      moy_Q2.5 = mean(temp_quantiles_Q2.5),
      moy_Q50 = mean(temp_quantiles_Q50),
      moy_Q97.5 = mean(temp_quantiles_Q97.5)
    )
  } else {
    all_roselieres_by_combination[[indice]] <- NA
  }
}
```

Représentation graphique :

```r
# Créer un dataframe à partir de la liste de moyennes des quantiles
all_roselieres_by_combination_df <- do.call(rbind, lapply(all_roselieres_by_combination, as.data.frame))
all_roselieres_by_combination_df = data.frame(t(sapply(all_roselieres_by_combination, unlist)))
all_roselieres_by_combination_df = cbind(all_roselieres_by_combination_df, n_cell)
rownames(all_roselieres_by_combination_df) = NULL

ggplot(all_roselieres_by_combination_df[-c(11,23),], aes(x = n_cell, y = moy_Q50))+
  geom_line()+
  geom_ribbon(aes(ymin = moy_Q2.5, ymax = moy_Q97.5), alpha=0.5)+
  labs(x = "Taille des cellules (en m²)", y = "Estimation des roselières/scirpaies")+
  #scale_y_continuous(limits = c(-0.01,0.01))+
  scale_x_reverse(breaks = all_roselieres_by_combination_df$n_cell)+
  theme_minimal()
```

![Régulier roselières et scirpaies](https://github.com/cha-dot/Processus-de-cox-log-gaussien-INLA-SPDE-/blob/images/roselieres_reg_pte_final.jpg?raw=true)

##### Friches

Calcul de la valeur moyenne du paramètre et de ses intervalles de crédibilité, pour chaque cas de figure (résolution de la grille) :

```r
all_friches_by_combination <- list()

# Boucle pour ouvrir chaque fichier, extraire les quantiles et calculer les moyennes
for (i in seq_along(n_cell)) {
  cell <- n_cell[i]
  pe <- n_pe[i]
  indice <- paste(cell, pe, sep = "_")
  temp_quantiles_Q2.5 <- c()
  temp_quantiles_Q50 <- c()
  temp_quantiles_Q97.5 <- c()
  
  for (rep in seq_len(n_rep)) {
    filename <- generate_filename(cell, pe, rep)
    if (file.exists(filename)) {
      load(filename)
      # Vérifier l'existence des variables de quantiles
      if (exists("post.stat.veg")) {
        # Ajouter les résultats aux vecteurs temporaires
        temp_quantiles_Q2.5 <- c(temp_quantiles_Q2.5, post.stat.veg["2.5%", "friches"])
        temp_quantiles_Q50 <- c(temp_quantiles_Q50, post.stat.veg["50%", "friches"])
        temp_quantiles_Q97.5 <- c(temp_quantiles_Q97.5, post.stat.veg["97.5%", "friches"])
      } else {
        warning(sprintf("Les quantiles n'existent pas dans le fichier %s.", filename))
      }
    } else {
      warning(sprintf("Le fichier %s n'existe pas.", filename))
    }
  }
  
  # Calculer les moyennes des quantiles pour cette combinaison
  if (length(temp_quantiles_Q2.5) > 0 & length(temp_quantiles_Q50) > 0 & length(temp_quantiles_Q97.5) > 0) {
    all_friches_by_combination[[indice]] <- c(
      moy_Q2.5 = mean(temp_quantiles_Q2.5),
      moy_Q50 = mean(temp_quantiles_Q50),
      moy_Q97.5 = mean(temp_quantiles_Q97.5)
    )
  } else {
    all_friches_by_combination[[indice]] <- NA
  }
}
```

Représentation graphique :

```
# Créer un dataframe à partir de la liste de moyennes des quantiles
all_friches_by_combination_df <- do.call(rbind, lapply(all_friches_by_combination, as.data.frame))
all_friches_by_combination_df = data.frame(t(sapply(all_friches_by_combination, unlist)))
all_friches_by_combination_df = cbind(all_friches_by_combination_df, n_cell)
rownames(all_friches_by_combination_df) = NULL

ggplot(all_friches_by_combination_df[-c(11:23),], aes(x = n_cell, y = moy_Q50))+
  geom_line()+
  geom_ribbon(aes(ymin = moy_Q2.5, ymax = moy_Q97.5), alpha=0.5)+
  labs(x = "Taille des cellules (en m²)", y = "Estimation des friches")+
  #scale_y_continuous(limits = c(-0.01,0.01))+
  scale_x_reverse(breaks = all_friches_by_combination_df$n_cell)+
  theme_minimal()
```

![Régulier friches](https://github.com/cha-dot/Processus-de-cox-log-gaussien-INLA-SPDE-/blob/images/friches_reg_pte_final.jpg?raw=true)
