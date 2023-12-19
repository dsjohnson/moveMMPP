
###-----------------------------------------------------------------------------
### Q data for {terra} rasters
###-----------------------------------------------------------------------------
#' @importFrom terra compareGeom xFromCell yFromCell adjacent ncell values
#' @import dplyr
make_q_data_rast <- function(cell_data, rast_mask=NULL,...){
  
  mask <- cell <- from_cell <- NULL
  
  ### data checks
  if(!is.null(rast_mask)){
    if(!terra::compareGeom(cell_data, rast_mask)) stop("'cell_data' and 'rast_mask' must have the same geometry!")
  }
  
  N <- terra::ncell(cell_data)
  
  ### Make Q_r data
  q_r_data <- dplyr::tibble(
    cell = 1:N,
    x = terra::xFromCell(cell_data, 1:N),
    y = terra::yFromCell(cell_data, 1:N),
    mask = terra::values(rast_mask, dataframe=TRUE)[,1]
  ) 
  
  q_r_data <- q_r_data %>% bind_cols(terra::values(cell_data, dataframe=TRUE)) %>% 
    filter(!is.na(mask)) %>% select(-mask)
  
  q_r_data <- q_r_data %>% mutate(
    cellx =  as.numeric(factor(cell))
  )
  
  ### Make Q_m data
  q_m_data <- tibble(
    from_cell = rep(1:N, each=4), 
    cell=as.vector(t(terra::adjacent(cell_data, 1:N)))
  ) %>% filter(!is.na(cell)) %>% filter(from_cell %in% q_r_data$cell) %>% 
    filter(cell %in% q_r_data$cell)
  
  covs <- q_r_data 
  q_m_data <- q_m_data %>% left_join(covs, by="cell")
  colnames(covs) <- paste0("from_", colnames(covs))
  q_m_data <- q_m_data %>% left_join(covs, by="from_cell")
  rm(covs)
  
  return(list(q_r = q_r_data, q_m = q_m_data))

}


###-----------------------------------------------------------------------------
### Q data for {sf} polygons
###-----------------------------------------------------------------------------

#' @importFrom sf st_drop_geometry st_geometry_type
#' @importFrom spdep poly2nb
#' @import dplyr
make_q_data_sf <- function(cell_data, cell_name,...){
  
  cellx <- from_cellx <- neighborhood <- NULL
  
  ### checks
  if(!all(st_geometry_type(cell_data) %in% c("POLYGON","MULTIPOLYGON"))){stop("'cell_data' all geometry types must be POLYGON or MULTIPOLYGON!")}
  if(is.null(cell_name)){
    cell_name <- "cell"
    cell_data$cell <- 1:nrow(cell_data)
  }
  cell_data$cellx <- as.integer(factor(cell_data[[cell_name]]))
  cell_data <- cell_data %>% rename(cell=.data[[cell_name]])
  ### Make Q_r data
  q_r_data <- cell_data %>% st_drop_geometry()
  
  ### Make Q_m data
  nb <- spdep::poly2nb(cell_data, queen=FALSE)
  q_m_data <- cell_data %>% select(cellx) %>% st_drop_geometry() %>% rename(from_cellx = cellx) %>%
    rowwise() %>%
    mutate(
      neighborhood = list(
        tibble(cellx = nb[[from_cellx]], num_neigh=length(nb[[from_cellx]]))
      )
    ) %>% ungroup() %>% arrange(from_cellx) %>% unnest(cols=neighborhood) 
  q_m_data <- q_m_data %>% arrange(from_cellx, cellx)
  covs <- q_r_data 
  q_m_data <- q_m_data %>% left_join(covs, by="cellx")
  colnames(covs) <- paste0("from_", colnames(covs))
  q_m_data <- q_m_data %>% left_join(covs, by="from_cellx")
  rm(covs)

  return(list(q_r = q_r_data, q_m = q_m_data))
}