#' Check Signal peptide fidelity
#'
#' A convenience function to compute atomic contacts between a predicted signal peptide and the rest of the protein in AlphaFold2 structures. This implements the FORTRAN code from Elcock-Lab/AlphaFold2-Signal, with added features to calculate the median pLDDT and count residues surpassing the pLDDT threshold for the signal peptide. Unlike the FORTRAN version, this function highlights cases where low-confidence signal peptide structures filter all residues, suggesting they may not represent true outward-facing signal peptides.
#'
#' @importFrom bio3d read.pdb basename.pdb
#' @import dplyr purrr
#' @export
#'
#' @param pdb A link to PDB file or path to PDB file if the file is locally present.
#' @param cut_dist Distance cutoff for defining atomic contacts (Default: 4 Angstroms).
#' @param nsignal Number of residues that comprise of N-terminal signal peptide. If SignalP prediction is available for your protein, the predicted length should be used. (Default: 25).
#' @param bfac_thresh pLDDT threshold for including residues in the atomic contact calculation.
#' @param nskip No. of residues immediately next to the cleavage site to exclude from atomic contact calculations.
#'
#' @return A data frame, containing statistics about signal peptide and rest of the protein. Zero atomic or residue-residue contact is indicative of True positives while non-zero values are putative False Positives.
#' @examples
#' \dontrun{
#' df <- easyAF2Signal("https://alphafold.ebi.ac.uk/files/AF-Q9TY95-F1-model_v4.pdb")
#' }
#'
easyAF2Signal <- function(pdb, cut_dist = 4, nsignal = 25, bfac_thresh = 90, nskip = 1) {
  # Read PDB file and extract atom data
  atom_data <- bio3d::read.pdb(pdb)$atom

  # Calculate number of residues and atoms
  ires <- max(atom_data$resno)

  # Calculate average bfact for each residue
  bfac <- dplyr::tibble(resno = atom_data$resno, b = atom_data$b) %>%
    dplyr::group_by(resno) %>%
    dplyr::summarize(bfac = mean(b)) %>%
    dplyr::pull(bfac)

  # Signal peptide and protein region calculations
  signalp_res_with_significant_bfac <- sum(bfac[1:nsignal] > bfac_thresh)
  median_signalp_bfac <- median(bfac[1:nsignal])
  cleavage_site_bfac <- bfac[nsignal + 1]
  rest_res_with_significant_bfac <- sum(bfac[(nsignal + 2):length(bfac)] > bfac_thresh)
  median_rest_bfac <- median(bfac[(nsignal + 2):length(bfac)])

  # Initialize contact counters
  ncont_atm_sig_rest <- 0
  ncont_res_sig_rest <- 0

  # Contact calculations
  purrr::walk(1:(nsignal - nskip), function(i1) {
    if (bfac[i1] < bfac_thresh) {
      return()
    }
    purrr::walk((nsignal + nskip + 1):ires, function(j1) {
      if (bfac[j1] < bfac_thresh) {
        return()
      }

      atoms_i <- atom_data[atom_data$resno == i1, c("x", "y", "z")]
      atoms_j <- atom_data[atom_data$resno == j1, c("x", "y", "z")]
      distances <- as.matrix(dist(rbind(atoms_i, atoms_j)))[1:nrow(atoms_i), (nrow(atoms_i) + 1):ncol(distances)]
      contacts <- sum(distances <= cut_dist)
      ncont_atm_sig_rest <<- ncont_atm_sig_rest + contacts

      if (contacts > 0) ncont_res_sig_rest <<- ncont_res_sig_rest + 1
    })
  })

  # Prepare output
  data.frame(
    Name = basename.pdb(pdb)[[1]],
    length_signalpeptide = nsignal,
    length_protein = ires,
    postbfacFiltered_signalP_res = signalp_res_with_significant_bfac,
    medianSignalP_bfac = median_signalp_bfac,
    clevagesite_bfac = cleavage_site_bfac,
    postbfacFiltered_rest_res = rest_res_with_significant_bfac,
    medianRest_bfac = median_rest_bfac,
    res_res_conts = ncont_res_sig_rest,
    atm_atm_conts = ncont_atm_sig_rest
  )
}
