#' Apply PeacoQC quality control to a GatingSet
#' 
#' @description
#' Applies PeacoQC quality control to all samples in a GatingSet object. The function
#' processes each sample individually and adds a new gate containing the QC results
#' to the gating hierarchy.
#'
#' @usage
#' PeacoQC_GatingSet(gs, channels, ...)
#'
#' @param gs A \code{GatingSet} object containing the flow cytometry data to be
#'   quality controlled
#' @param channels Character vector specifying which channels should be used for
#'   quality control
#' @param ... Additional arguments passed to \code{PeacoQC()}
#'
#' @return A \code{GatingSet} object with an additional population node named
#'   "PeacoQC_good" containing the events that passed quality control
#'
#' @details
#' This function is a wrapper around the original PeacoQC algorithm that makes it
#' compatible with GatingSet objects. For each sample in the GatingSet:
#' 1. Extracts the underlying flow data
#' 2. Converts it to a flowFrame
#' 3. Runs PeacoQC
#' 4. Adds the results as a new gate in the hierarchy
#'
#' @importFrom flowCore read.FCS
#' @importFrom flowWorkspace load_cytoset_from_fcs GatingSet
#' @importFrom flowWorkspace sampleNames cytoframe_to_flowFrame 
#' @importFrom flowWorkspace gh_pop_get_data gs_get_pop_paths gs_pop_add
#' @importFrom utils tail
#'
#' @examples
#' 
#' # prepare GatingSet
#' fileName <- system.file("extdata", "111_Comp_Trans.fcs", package="PeacoQC")
#' cs <- flowWorkspace::load_cytoset_from_fcs(fileName)
#' gs <- flowWorkspace::GatingSet(cs)
#' 
#' channels <- c("FSC-A", "SSC-A", "B515-A", "G780-A", "G710-A", 
#'               "G660-A", "G610-A", "G560-A", "R780-A", "R710-A", 
#'               "R660-A", "V800-A", "V585-A", "V450-A")
#' 
#' 
#' # Run QC
#' PeacoQC_GatingSet(gs, 
#'                   channels = channels,
#'                   save_fcs = FALSE,
#'                   report = FALSE,
#'                   output_directory=NULL)
#' 
#' # check that the filter "gates" were added to the GatingSet
#' flowWorkspace::gs_get_pop_paths(gs)
#'
#' @seealso
#' \code{\link[PeacoQC]{PeacoQC}} for the underlying QC algorithm
#' \code{\link[flowWorkspace]{GatingSet}} for GatingSet operations
#'
#' @export
PeacoQC_GatingSet <- function(gs, channels, ...) {
  
  sample_ids <- sampleNames(gs) 
  
  # generate a list of filter "gates"
  lg <- lapply(sample_ids, function(sample_id) {
    ff <- gh_pop_get_data(gs[[sample_id]])
    qc_results <- PeacoQC(cytoframe_to_flowFrame(ff), 
                          channels=channels,
                          ...)
    good_events <- qc_results$GoodCells
    return(good_events)
  })
  names(lg) <- sample_ids
  
  ## add directly a list of filter "gates" to gs
  return(
    gs_pop_add(gs,
               gate = lg,
               parent =  tail(gs_get_pop_paths(gs), 1),
               name = "PeacoQC_good")
  )
}