#' @title Define specifics for model structure and fitting control
#' @name arg_funs
#' @aliases mmpp_model
#' @aliases mmpp_control
#' @param form Model formula for MMPP model components, residency, movement, and detection
#' @param link A link function for the model component. One of `"soft_plus"`, `"log"`, or `"logit"`.
#' @param a Scale parameter for the `"soft_plus"` link function. Ignored for `link="log"` or `"logit"`. The `a` parameter determines the approximation to 
#' a hard plus function, i.e., as `a` becomes large the soft plus function converges to `g^{-1}(x) = max(0,x)`. For this specification, `a` must be greater than or equal to 1.
#' @param L Lover bound for general logit link.
#' @param U Upper bound for general logit link. 
#' @param q_r A named list for the residency model. Must contain elements: `form`, `link`, and `a`. 
#' The easiest way to construct this is the `\link{mmpp_model}` function.
#' @param q_m A named list for the movement portion of the model. See `q_r` and `\link{mmpp_model}`.
#' @param lambda A named list for the detection portion of the model. See `\link{mmpp_model}`
#' @param struc The form of the rate matrix entries. Can be one of: `"mult"` for residency times movement, i.e., q_{ij} = q_r(i) * pi_m(i,j) or `"add"` for a model of 
#' the form q_{ij} = q_r(i) + pi_m(i,j). 
#' @param norm Should the movement portion be adjusted to sum to 1. This is the parameterization suggested by Hewitt et al. (2023). 
#' @references Hewitt, J., Gelfand, A. E., & Schick, R. S. (2023). Time-discretization approximation enriches continuous-time discrete-space models for animal movement. The Annals of Applied Statistics, 17:740-760.

#' @export
mmpp_model <- function(form=~1, link="soft_plus", a=1, L=0, U=0){
  if(link=="logit" & (U==0 | U<=L)) stop("U>0 and U>L must be specified when using the logit link.") 
  out <- list(
    form = form, link = link, a = a, L=L, U=U
  )
  class(out) <- c("list", "mmpp_mod")
  return(out)
}

#' @name arg_funs
#' @export
mmpp_control <- function(lambda=mmpp_model(), q_r = mmpp_model(), q_m=mmpp_model(), struc="mult", norm=TRUE){
  out <- list(
    lambda=lambda, q_r = q_r, q_m=q_m, struc=struc, norm=norm
  )
  class(out) <- c("list","mmpp_cont")
  return(out)
}
