#' @title Process MMPP movement detection data for model fitting
#' @param data Detection data. Must contain columns \code{timestamp} and \code{cell} to 
#' indicate time and place of detection. Also, must contain column named
#' \code{id} to indicate different individuals. 
#' @param cell_name Character. The column name in the data representing the 
#' discrete spatial cells on which the animal moves.
#' @param aug_timestamp Augmented times for modeling temporally dynamic movement or
#' detection. See 'Details'
#' @param time_unit ---.
#' @details (\code{aug_timestamp} can be either passed as a separate vector of 
#' POSIXct or \code{process_data} can derive a sequence of regularly spaced 
#' augmentation times from the original data. This is specified by providing a 
#' character string that corresponds to the by argument of the seq.POSIXt 
#' function (e.g. '1 year', '1 month'). 
#' @import dplyr tidyr lubridate
#' @export
process_data <- function(data, cell_name=NULL, aug_timestamp=NULL, time_unit="days"){
  
  timestamp <- period <- time <- elapse <- cell <- NULL
  
  # browser()
  
  if(is.null(cell_name)) stop("The argument 'cell_name' must be specified as a column in the data!")
  names(data)[names(data) == cell_name] <- "cell"
  # data$cellx <- as.numeric(factor(data[[cell_name]]))
  
  min_dt <- min(data$timestamp,na.rm=TRUE)
  max_dt <- max(data$timestamp,na.rm=TRUE)
  if(!is.null(aug_timestamp)){
    if(is.character(aug_timestamp)){
      t_int <- unlist(strsplit(aug_timestamp, " "))
      unit_idx <- ifelse(length(t_int)==2, 2, 1)
      min_seq <- lubridate::floor_date(min_dt,t_int[unit_idx])
      max_seq <- lubridate::ceiling_date(max_dt,t_int[unit_idx])
      aug_times <- seq(min_seq, max_seq, by = aug_timestamp)
    } else if(length(aug_timestamp)>1 & any(class(aug_timestamp)%in%c("POSIXct","POSIXt"))){
      t_int <- attr(diff(aug_timestamp), "units")
      min_seq <- lubridate::floor_date(min_dt,t_int)
      max_seq <- lubridate::ceiling_date(max_dt,t_int)
      aug_times <- c(min_seq, aug_timestamp[aug_timestamp>min_dt & aug_timestamp<max_dt], max_seq)
    }
    aug_df <- data.frame(timestamp=aug_times, period=1:length(aug_times), quad=1)
  } else{
    aug_df <-  data.frame(timestamp=floor_date(min_dt, time_unit), period=1, quad=1)
  }
  # %>% mutate(start = timestamp, end = c(tail(timestamp,-1),NA))
  data <- data %>% group_by(id) %>% arrange(timestamp) %>% nest() %>% 
    mutate(ns = sapply(data, nrow) %>% as.integer)
  if(any(data$ns<=1)){
    ns <- sum(data$ns<=1)
    warning(paste0(ns, " Individuals with <= 1 sighting have been removed!"))
    data <- filter(data, ns>1) 
  }
  out <- NULL
  for(i in 1:nrow(data)){
    min_t <- min(data$data[[i]]$timestamp)
    max_t <- max(data$data[[i]]$timestamp)
    tmp <- bind_rows(data$data[[i]], aug_df) %>% arrange(timestamp) %>% 
      fill(period) %>% filter(timestamp>= min_t, timestamp<= max_t) %>% 
      mutate(id = data$id[[i]])
    tmp <- tmp[,c('id','timestamp','period','cell','quad')]
    out <- bind_rows(out,tmp)
  }
  out <- arrange(out, id, timestamp)
  out$quad <- ifelse(is.na(out$quad), 0, out$quad)
  out$quad[out$timestamp==min(out$timestamp)] <- 1
  out$quad[out$timestamp==max(out$timestamp)] <- 1
  
  out$idx <- as.integer(factor(out$id))
  out$time <- out$timestamp-min(out$timestamp)
  units(out$time) <- time_unit
  out$time <- as.numeric(out$time)
  out <- group_by(out, id) %>% 
    mutate(
      elapse = time-min(time),
      delta = c(NA, diff(elapse))
    ) %>% ungroup()
  
  # out$period <- ifelse(is.na(out$cell), out$period-1, out$period)
  attr(out, "cell_name") <- cell_name
  attr(out, 'proc_data') <- TRUE
  
  return(out)
}
