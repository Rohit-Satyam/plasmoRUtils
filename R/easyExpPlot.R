#' Plot Normalised expression quickly
#'
#' A convenience function to make line plots and bubble plots to showcase gene expression trends across time points or sample groups
#'
#' @import ggplot2
#'
#' @param df A 3 column data frame that has bee transformed using reshape2::melt function.
#' @param x,y x and y variables for drawing.
#' @param value Remaining column after choosing x and y. This column should be numeric if type="line" and character if type="bubble".
#' @param type Type of plot. Default ("line").
#' @return A line plot or a bubble plot that represents trends across different time points and sample types.
#' @export
#' @examples
#' \dontrun{
#'   # Load sample data that contain Z-score transformed TPM values (randomly generated) for some genes
#'   data(pf3d7TPMs)
#'
#'
#'   df <- reshape2::melt(pf3d7TPMs,"Probe",na.rm = T)
#'
#'   easyExpPlot(df,x="variable",y="Probe",value="value", type = "bubble")
#'   easyExpPlot(df,x="variable",y="value",value="Probe")
#' }
#'
#'



easyExpPlot <- function(df, x, y, value, type = "line", scaleBubbles = c(2, 6)) {
  if (ncol(df) > 3) {
    message(
      "Provided data frame contains more than 3 columns. Consider reshaping it first ",
      "(e.g., with reshape2::melt)."
    )
  }
  
  base_theme <- ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(face = "bold"),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank()
    )
  
  if (type == "line") {
    ggplot2::ggplot(
      df,
      ggplot2::aes(
        x = .data[[x]],
        y = .data[[y]],
        group = .data[[value]]
      )
    ) +
      ggplot2::geom_line(color = "grey") +
      ggplot2::geom_hline(
        yintercept = 0,
        color = "black",
        linewidth = 0.5
      ) +
      ggplot2::stat_summary(
        ggplot2::aes(group = 1),
        fun = mean,
        geom = "line",
        color = "red",
        linewidth = 1
      ) +
      base_theme +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(
          angle = 90,
          vjust = 0.5,
          hjust = 1,
          face = "bold"
        )
      )
    
  } else if (type == "bubble") {
    ggplot2::ggplot(
      df,
      ggplot2::aes(
        x = .data[[x]],
        y = .data[[y]],
        size = .data[[value]]
      )
    ) +
      ggplot2::geom_point(
        shape = 21,
        fill = "#9FB9E3",
        color = "white",
        stroke = 1.5
      ) +
      base_theme +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(
          angle = 90,
          vjust = 0.5,
          hjust = 1,
          face = "bold"
        ),
        panel.spacing = grid::unit(0.1, "lines"),
        axis.ticks = ggplot2::element_blank(),
        axis.title.x = ggplot2::element_text(face = "bold"),
        axis.title.y = ggplot2::element_text(face = "bold")
      ) +
      ggplot2::scale_size(range = scaleBubbles)
    
  } else {
    stop("Invalid plot type. Use 'line' or 'bubble'.")
  }
}
