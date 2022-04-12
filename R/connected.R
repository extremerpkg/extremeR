#' Test graph connected
#'
#' Test if nb object is fully connected. Created by Mitzi Morris.
#'  @seealso \url{https://github.com/stan-dev/example-models/blob/master/knitr/car-iar-poisson/nb_data_funs.R}
#'
#' @param nb_object nb object in prepared data list
#'
#' @importFrom spdep n.comp.nb
#' @return
#' @export
#'
#' @examples
testconnected <- function(nb_object) {
  if (!isdisconnected(nb_object)) {
    print("Success! The Graph is Fully Connected")
  } else{
    warning("Failure... Some parts of the graph are still disconnected...")
  }
}

#' Check graph disconnected
#'
#' Test if nb object is disconnected
#'
#' @param x nb object
#'
#' @importFrom spdep n.comp.nb
#' @return
#' @export
#'
#' @examples
isdisconnected = function(x) {
  return(n.comp.nb(x)[[1]] > 1);
}
