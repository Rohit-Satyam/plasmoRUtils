
#' Check Signal Peptide Fidelity
#'
#' A convenience function to compute atomic contacts between a predicted
#' signal peptide and the rest of the protein in AlphaFold2 structures.
#' This implements the FORTRAN code from Elcock-Lab/AlphaFold2-Signal,
#' with added features to calculate the median pLDDT and count residues
#' surpassing the pLDDT threshold for the signal peptide.
#'
#' Unlike the FORTRAN version, this function highlights cases where
#' low-confidence signal peptide structures filter all residues,
#' suggesting they may not represent true outward-facing signal peptides.
#'
#' @importFrom bio3d read.pdb
#' @import dplyr purrr
#' @export
#'
#' @param pdb URL or local path to a PDB file.
#' @param cut_dist Distance cutoff for defining atomic contacts
#'   (Default: 4 Angstroms).
#' @param nsignal Number of residues comprising the N-terminal signal
#'   peptide (Default: 25).
#' @param bfac_thresh pLDDT threshold for including residues in atomic
#'   contact calculations (Default: 90).
#' @param nskip Number of residues adjacent to the cleavage site to
#'   exclude from atomic contact calculations (Default: 1).
#'
#' @return A data frame containing signal peptide confidence statistics,
#'   atomic contacts and residue-residue contacts.
#'
#' @examples
#' \dontrun{
#' df <- easyAF2Signal(
#'   "https://alphafold.ebi.ac.uk/files/AF-A0A140KXW6-F1-model_v6.pdb"
#' )
#' }

easyAF2Signal <- function(
    pdb,
    cut_dist = 4,
    nsignal = 25,
    bfac_thresh = 90,
    nskip = 1
) {
  
  # 1. Validate input
  if (!is.character(pdb) ||
      length(pdb) != 1L ||
      is.na(pdb) ||
      !nzchar(trimws(pdb))) {
    stop("pdb must be a valid URL or local file path.", call. = FALSE)
  }
  
  # Preserve the original filename before parsing the PDB
  pdb_name <- sub(
    "\\.pdb$",
    "",
    basename(sub("[?#].*$", "", pdb)),
    ignore.case = TRUE
  )
  
  pdb_path <- pdb
  
  # 2. Download remote PDB files before parsing
  if (grepl("^https?://", pdb, ignore.case = TRUE)) {
    
    pdb_tmp <- tempfile(fileext = ".pdb")
    
    # Remove temporary file when function exits
    on.exit(unlink(pdb_tmp), add = TRUE)
    
    httr2::request(pdb) |>
      httr2::req_user_agent("plasmoRUtils easyAF2Signal") |>
      httr2::req_timeout(60) |>
      httr2::req_retry(
        max_tries = 3,
        retry_on_failure = TRUE
      ) |>
      httr2::req_perform(path = pdb_tmp)
    
    # Check downloaded file
    if (!file.exists(pdb_tmp) ||
        is.na(file.info(pdb_tmp)$size) ||
        file.info(pdb_tmp)$size == 0) {
      stop(
        "PDB download failed or returned an empty file.",
        call. = FALSE
      )
    }
    
    pdb_path <- pdb_tmp
    
  } else {
    
    if (!file.exists(pdb_path)) {
      stop(
        "PDB file does not exist: ", pdb_path,
        call. = FALSE
      )
    }
  }
  
  # 3. Parse the local PDB structure
  pdb_data <- bio3d::read.pdb(pdb_path)
  
  atom_data <- pdb_data$atom
  
  # 4. Calculate number of residues and atoms
  ires <- max(atom_data$resno)
  
  # 5. Calculate average pLDDT for each residue
  bfac <- dplyr::tibble(
    resno = atom_data$resno,
    b = atom_data$b
  ) %>%
    dplyr::group_by(resno) %>%
    dplyr::summarise(
      bfac = mean(b),
      .groups = "drop"
    ) %>%
    dplyr::pull(bfac)
  
  # 6. Signal peptide and protein region calculations
  signalp_res_with_significant_bfac <- sum(
    bfac[1:nsignal] > bfac_thresh
  )
  
  median_signalp_bfac <- stats::median(
    bfac[1:nsignal]
  )
  
  cleavage_site_bfac <- bfac[nsignal + 1]
  
  rest_res_with_significant_bfac <- sum(
    bfac[(nsignal + 2):length(bfac)] > bfac_thresh
  )
  
  median_rest_bfac <- stats::median(
    bfac[(nsignal + 2):length(bfac)]
  )
  
  # 7. Initialize contact counters
  ncont_atm_sig_rest <- 0
  ncont_res_sig_rest <- 0
  
  # 8. Calculate atomic contacts
  purrr::walk(1:(nsignal - nskip), function(i1) {
    
    # Skip signal peptide residues below pLDDT threshold
    if (bfac[i1] < bfac_thresh) {
      return()
    }
    
    purrr::walk((nsignal + nskip + 1):ires, function(j1) {
      
      # Skip protein residues below pLDDT threshold
      if (bfac[j1] < bfac_thresh) {
        return()
      }
      
      # Extract atomic coordinates for both residues
      atoms_i <- atom_data[
        atom_data$resno == i1,
        c("x", "y", "z")
      ]
      
      atoms_j <- atom_data[
        atom_data$resno == j1,
        c("x", "y", "z")
      ]
      
      # Skip residues without atomic coordinates
      if (nrow(atoms_i) == 0 || nrow(atoms_j) == 0) {
        return()
      }
      
      # Calculate pairwise atomic distances
      dist_mat <- as.matrix(
        stats::dist(rbind(atoms_i, atoms_j))
      )
      
      distances <- dist_mat[
        1:nrow(atoms_i),
        (nrow(atoms_i) + 1):ncol(dist_mat)
      ]
      
      # Count atomic contacts within the distance cutoff
      contacts <- sum(distances <= cut_dist)
      
      ncont_atm_sig_rest <<-
        ncont_atm_sig_rest + contacts
      
      # Count residue-residue contacts
      if (contacts > 0) {
        ncont_res_sig_rest <<-
          ncont_res_sig_rest + 1
      }
      
    })
  })
  
  # 9. Prepare output
  df <- data.frame(
    Name = pdb_name,
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
  
  return(df)
}
