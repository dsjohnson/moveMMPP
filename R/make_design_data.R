#' @title Produce design data for use in fitting MMPP movement models
#' @param sighting_data ---. 
#' @param cell_data ---.
#' @param add_fix ---.
#' @param add_boundary_length Logical. ---. Defaults to \code{FALSE}.
#' @param debug Debugging level: 1-3 mainly for package developers.
#' @param ... Ignored arguments.
#' @import dplyr sf
#' @importFrom spdep poly2nb
#' @importFrom rmapshaper ms_innerlines
#' @importFrom units set_units
#' @export
make_design_data <- function(sighting_data, cell_data, add_fix=TRUE, 
                             add_boundary_length=FALSE, debug=0,...){
  
  cell <- cellx <- timestamp <- quad <- period <- fix <- from_cellx <- NULL
  neighborhood <- to_cellx <- boundary <- NULL
  
  if(debug==1) browser()
  
  obs_cells <- unique(sighting_data$cell)
  obs_cells <- obs_cells[!is.na(obs_cells)]
  cell_data$cellx <- as.integer(factor(cell_data$cell))
  cells <- select(cell_data, cell, cellx) %>% st_drop_geometry()
  if(!all(obs_cells %in% cells$cell)) stop("There are observed cells not in the cell data set!")
  
  quad_data <- select(sighting_data, timestamp, quad, period) %>% 
    filter(quad==1) %>% distinct() %>% arrange(timestamp) %>% select(-quad)
  
  # time*cell data frame for lambda model
  lambda_data <- bind_cols(
    quad_data[rep(seq_len(nrow(quad_data)-1), each=nrow(cells)),],
    cells[rep(seq_len(nrow(cells)), nrow(quad_data)-1),]
  )
  if(add_fix){
    fix_data <- tibble(cell=obs_cells, fix=1)
    lambda_data <- left_join(lambda_data, fix_data, by='cell')
    lambda_data <- lambda_data %>% mutate(fix = ifelse(is.na(fix), 0, NA))
  }
  
  q_r_data <- cell_data %>% st_drop_geometry()
  lambda_data <- left_join(lambda_data, q_r_data, by = c("cell", "cellx"))
  

  ### Create Q matrix data frame
  if(debug==2) browser()
  
  nb <- spdep::poly2nb(cell_data, queen=FALSE)
  q_m_data <- cell_data %>% select(cellx) %>% st_drop_geometry() %>% rename(from_cellx = cellx) %>%
    rowwise() %>%
    mutate(
      neighborhood = list(
        tibble(cellx = nb[[from_cellx]], num_neigh=length(nb[[from_cellx]]))
      )
    ) %>% ungroup() %>% arrange(from_cellx) %>% unnest(cols=neighborhood) 
  
  q_m_data <- q_m_data %>% mutate(
    from_area = st_area(cell_data[from_cellx,]) %>% set_units("km^2"),
    area = st_area(cell_data[cellx,]) %>% set_units("km^2")
  ) %>% arrange(from_cellx, cellx)
  
  if(add_boundary_length) warning("'add_boudary_length' not currently implemented.")
  # if(add_boundary_length){
  #   move_data <- move_data %>% rowwise() %>%
  #     mutate(
  #       boundary = list(ms_innerlines(cell_data[c(from_cellx,cellx),])),
  #       boundary_length = st_length(boundary) %>% set_units("km")
  #     ) %>% ungroup() %>% select(-boundary)
  # }
  
  cov_nms <- colnames(q_r_data)
  colnames(q_r_data) <- paste0("from_",cov_nms)
  q_m_data <- left_join(q_m_data, q_r_data, by="from_cellx")
  colnames(q_r_data) <- cov_nms
  q_m_data <- left_join(q_m_data, q_r_data, by="cellx")
  
  q_m_data$fix <- NA
  
  # if(dynamic_movement){
  #   # add code here to expand move_data by quad_data if time indexed movement is 
  #   # desired.
  # }
  
  obs_cell <- cells %>% filter(cell %in% obs_cells)
  
  out <- list(
    lambda = lambda_data,
    q_r = q_r_data,
    q_m = q_m_data,
    quad_pts = quad_data,
    obs_cell = obs_cell
  )
  
  if(debug==3) browser()
  
  return(out)
  
}