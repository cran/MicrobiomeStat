#' Microbiome data from the human upper respiratory tract (left and right throat)
#'
#' A dataset containing "otu", "tax", meta", "genus", family"
#'
#' @format A  list with elements
#' \describe{
#'   \item{otu}{otu table, 2156 taxa by 290 samples}
#'   \item{tax}{taxonomy table, 2156 taxa by 7 taxonomic ranks}
#'   \item{meta}{meta table, 290 samples by 53 sample variables}
#'   \item{genus}{304 by 290}
#'   \item{family}{113 by 290}
#' }
#'
#' @source
#' %%  ~~ reference to a publication or URL from which the data were obtained ~~
#' \url{https://qiita.ucsd.edu/} study ID:524
#' Reference: Charlson ES, Chen J, Custers-Allen R, Bittinger K, Li H, et al. (2010)
#' Disordered Microbial Communities in the Upper Respiratory Tract of
#' Cigarette Smokers. PLoS ONE 5(12): e15216.
#'
#' @usage data(smokers)
#' @docType data
#' @keywords datasets
"smokers"
