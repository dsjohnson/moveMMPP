#' @title Produce design data for use in fitting MMPP movement models
#' @param sighting_data ---. 
#' @param cell_data ---.
#' @param rast_mask Raster mask for inaccessible cells when \code{cell_data} is of type \code{SpatRaster} from the \code{terra} package. This is ignored
#' if \code{cell_data} is an \code{POLYGON} data frame from the \code{sf} package.
#' @param debug Debugging level: 1-3 mainly for package developers.
#' @param ... Ignored arguments.
#' @import dplyr 
#' @importFrom units set_units
#' @export
make_design_data <- function(sighting_data=NULL, cell_data, rast_mask=NULL, debug=0,...){
  
  cells <- obs <- timestamp <- quad <- NULL
  
  if(debug==1) browser()
  
  have_sight_data <- !is.null(sighting_data)
  if(have_sight_data){
    if(!attr(sighting_data, "proc_data")) stop("Sighting data must first be processed through the 'process_data(...) function!")
    cell_name <- attr(sighting_data,"cell_name")
    obs_cells <- unique(sighting_data$cell)
    obs_cells <- obs_cells[!is.na(obs_cells)]
    quad_data <- select(sighting_data, timestamp, quad, period) %>% 
      filter(quad==1) %>% distinct() %>% arrange(timestamp) %>% select(-quad)
  } else{
    cell_name <- NULL
    quad_data <- NULL
  }
  
  if(debug==2) browser()
  
  if(inherits(cell_data, "SpatRaster")){
    q_list <- make_q_data_rast(cell_data, rast_mask=rast_mask)
  } else if(inherits(cell_data, "sf")){
    q_list <- make_q_data_sf(cell_data, cell_name)
  }
  
  lambda_data <- q_list$q_r
  
  if(have_sight_data){
    period_data <- bind_cols(
      quad_data[rep(seq_len(nrow(quad_data)-1), each=nrow(lambda_data)),],
      cell=lambda_data$cell[rep(seq_len(nrow(lambda_data)), nrow(quad_data)-1)]
    )
    lambda_data <- left_join(lambda_data, period_data, by="cell")
    lambda_data$obs <- ifelse(lambda_data$cell%in%obs_cells, 1, 0)
    avail_cells <- unique(lambda_data$cell)
    if(!all(obs_cells %in% avail_cells)) stop("There are observed cells in the sighting data not in the cell data!")
  }
  
  out <- list(
    lambda = lambda_data,
    q_r = q_list$q_r,
    q_m = q_list$q_m,
    quad_pts = quad_data
  )
  
  if(debug==3) browser()
  
  return(out)
  
}
  
  
  
  
   