
# remotes::install_github("ropensci/rnaturalearthhires")
# remotes::install_github("kaskr/adcomp/TMB")
library(raster)
library(sf)
library(rnaturalearth)
library(splines2)
library(TMB)

setwd("C:/Users/James.Thorson/Desktop/Work files/AFSC/2022-02 -- Nielsen replace imager with CTMC")
#compile( "hmm_TMB.cpp", framework='TMBad' )

load("hmm_data_178690.RData")
load("ai_bathy_3km.RData") #smaller spatial extent and 3 km grid size
load("Pcod178690_likelihood_3km.RData") #smaller spatial extent and 3 km grid size

# Convert bathymetry
bathy_sp = as( bathy.maps.3k[[1]], 'SpatialPolygonsDataFrame')
bathy_sf = st_as_sf(bathy_sp)   # Goes columns -> rows
bathy_sf$northings = st_coordinates(st_centroid(bathy_sf))[,'Y']

# Get land layer
sf_states = ne_states( "United States of America", return="sf")
sf_states = sf_states[pmatch("Alas", sf_states$name_en),]
sf_states = st_transform( sf_states, crs=st_crs(bathy_sf) )
sf_states = st_geometry(sf_states)

# Format likelihood
L_gt = apply( L, MARGIN=3, FUN=function(mat){as.vector(t(mat))} )  
sf_L_gt = st_sf( st_geometry(bathy_sf), L_gt )

# Plot likelihood as sanity check
cols = round( seq(1,ncol(st_drop_geometry(sf_L_gt)),length=9) )
plotgrid = st_sf( st_geometry(bathy_sf), "density"=st_drop_geometry(sf_L_gt)[,cols] )
#plot( plotgrid, border=NA)

# Get adjacency matrix using sf
# Doesn't work because sf results in non-connected graph
# grid_A is same order as bathy_sf, i.e., column -> row
#st_rook = function(a, b = a, ...) st_relate(a, b, pattern = "F***1****", ... )
#st_queen = function(a, b = a, ...) st_relate(a, b, pattern = "F***T****", ... )
#grid_A = st_rook( bathy_sf, sparse=TRUE )
#A_gg = as(grid_A, "sparseMatrix")
#  A_gg = as(A_gg, "TsparseMatrix")

# Get adjacency matrix using Raster
# A_gg is column -> row  ("left to right and then top to bottom")
A_gg = adjacent( bathy.maps.3k[[1]], cells=1:prod(dim(bathy.maps.3k[[1]])) )
A_gg = Matrix::sparseMatrix( i=A_gg[,1], j=A_gg[,2], x=rep(1,nrow(A_gg)) )
A_gg = as(A_gg, "TsparseMatrix")

# Drop geometry and likelihood on land
# as.vector goes rows -> column so must transpose to conform to st_as_sf
which_intersect = which( as.integer(st_intersects(sf_L_gt, sf_states))==TRUE )
#sf_L_gt = st_difference( sf_L_gt, sf_states )
#bathy_sf = st_difference( bathy_sf, sf_states )
sf_L_gt = sf_L_gt[-which_intersect,]
bathy_sf = bathy_sf[-which_intersect,]
A_gg = A_gg[-which_intersect,-which_intersect]
#plot( bathy_sf[,c('northings','ai_bathy_fill')], border=NA )

# Assemble inputs
colsumA_g = colSums(A_gg)
I_gg = Matrix::sparseMatrix( i=1:nrow(bathy_sf), j=1:nrow(bathy_sf), x=rep(1,nrow(bathy_sf)) )
  I_gg = as(I_gg, "TsparseMatrix")
D_gg = Matrix::sparseMatrix( i=1:nrow(bathy_sf), j=1:nrow(bathy_sf), x=colsumA_g )
  D_gg = as(D_gg, "TsparseMatrix")
At_zz = cbind( attr(A_gg,"i"), attr(A_gg,'j') )

# Check whether fully connected -- GOOD
if( FALSE ){
  graph = igraph::graph_from_adjacency_matrix( A_gg )
  igraph::is_connected( graph, mode = "weak" )
  igraph::is_connected( graph, mode = "strong" )
}

# Make covariates
#X_gz = model.matrix( ~ 0, data=bathy_sf )
#X_gz = model.matrix( ~ scale(ai_bathy_fill), data=bathy_sf )
X_gz = bSpline( x=scale(bathy_sf$ai_bathy_fill), df=1, degree=1 )

#####################
# Diffusion-only model
#####################

# Compile
compile( "hmm_TMB.cpp", framework='TMBad' )
dyn.load( dynlib("hmm_TMB") )

# Build inputs
Params = list( "ln_D" = 1,
               "gamma_z" = 0.01*rnorm(ncol(X_gz)) )
Data = list( "CTMC_version" = 1,
             "Nmax" = 300,
             "colsumA_g" = colsumA_g,
             "L_gt" = as.matrix(st_drop_geometry(sf_L_gt)),
             "X_gz" = X_gz,
             "At_zz" = At_zz,
             "gproj" = which(st_drop_geometry(sf_L_gt)[,1]>0)-1 ) #,
             #"A_gg" = A_gg,
             #"D_gg" = D_gg )
Map = list()
  #Map$gamma_z = factor(rep(NA,ncol(X_gz)))

# Build object
obj = MakeADFun( data=Data, parameters=Params, map=Map, DLL="hmm_TMB" )
#report = obj$report()

# Sanity checks
if( FALSE ){
  # should be nonzero
  report$nll_t
  # Should equal 0 -- GOOD
  range(rowSums( report$Mrate_gg ))
  # Check Metzler
  all( attr(report$Mrate_gg,"x")[attr(report$Mrate_gg,"i")==attr(report$Mrate_gg,"j")] >= 0 )
  # Should equal 1 -- GOOD
  colSums( report$prob_gt )
  # Should all be positive -- GOOD
  all( report$pred_gt>=0 )
  all( report$proj_gt>=0 )
  all( report$prob_gt>=0 )
  # Should equal 1 -- NOT WORKING FOR SOME REASON
  range(colSums( report$pred_gt )[-1])
  # Should equal 1 -- NOT WORKING FOR SOME REASON
  range(colSums( report$proj_gt )[-1])
}

# Optimize
opt = TMBhelper::fit_tmb( obj=obj )
#report = obj$report()
cols = round( seq(2,ncol(report$prob_gt)-1,length=8) )

# Plot covariates
if( FALSE ){
  Xval = seq( min(scale(bathy_sf$ai_bathy_fill)),max(scale(bathy_sf$ai_bathy_fill)),length=1000)
  Xpred_xz = predict(X_gz, Xval )
  matplot( y=Xpred_xz, x=Xval, type="l", lwd=2 )
  Yval = Xpred_xz %*% obj$env$parList()$gamma_z
  plot( x=Xval, y=Yval, type="l", lwd=2 )
}

# filtered probability
plotgrid = st_sf( st_geometry(bathy_sf), "density"=report$pred_gt[,cols] )
plot( plotgrid, border=NA)

# one-step-prediction
plotgrid = st_sf( st_geometry(bathy_sf), "density"=report$prob_gt[,cols] )
plot( plotgrid, border=NA)

# Unconditional projection
plotgrid = st_sf( st_geometry(bathy_sf), 
                  "bathy" = bathy_sf$ai_bathy_fill, 
                  "pref" = exp(report$h_g), 
                  "density"=report$proj_gt[,cols] )
plot( plotgrid, border=NA, logz=TRUE )
