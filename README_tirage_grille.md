# Suppression progressive de points d'écoute selon une grille

Cela consiste à définir une grille sur notre domaine d'étude et à ne garder qu'un seul point d'écoute, tiré aléatoirement, par cellule.

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

Le reste du modèle en change pas (pour plus d'explications du modèle, se référer à README_INLA.md).

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

Sauvegardons le `summary` :

```r
  length(liste)
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

Sauvegardons ces 3 dernières métriques :

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
  ind=grep("veg",rownames(pr.int.tot[[1]]$latent))
  m = grep("max_sub", rownames(pr.int.tot[[1]]$latent))
  d = grep("duree_sub", rownames(pr.int.tot[[1]]$latent))
  post=matrix(unlist(lapply(pr.int.tot,function(x){x$latent[c(ind,m,d),]})),nrow=Nrep,byrow=T)
  post = as.data.frame(post)
  colnames(post) = c(rownames(pr.int.tot[[1]]$latent)[ind], "intensite_sub", "duree_sub")

  anciens = c("veg:1", "veg:2", "veg:3", "veg:4", "veg:5", "veg:6", "veg:7")
  nouv = c("cultures", "rase", "haute fauchée", "haute non fauchée", "arbustive", "roselières/scirpaies", "friches")
    
  for (p in 1:length(anciens)) {
      if (anciens[p] %in% colnames(post)) {
        colnames(post)[colnames(post) == anciens[p]] <- nouv[p]
      }
    }
    
  post.stat = apply(post,2,quantile,probs=c(0.025,0.5,0.975)) # valeurs des paramètres pour chaque quantile de chaque variable
  post.stat.veg = post.stat[, which(colnames(post) %in% nouv)]
  post.stat.duree = post.stat[, "duree_sub"]
  post.stat.max = post.stat[, "intensite_sub"]
  
  chemin = "~/post_reg_nouveaux"
  fichier = paste0("reg_post_", cel, "_nPE_", length(liste), "_rep_", repetition, ".RData")
  nom_complet = file.path(chemin, fichier)
  save(post.stat.veg, post.stat.duree, post.stat.max, file = nom_complet)
```
