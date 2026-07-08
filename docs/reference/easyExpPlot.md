# Plot Normalised expression quickly

A convenience function to make line plots and bubble plots to showcase
gene expression trends across time points or sample groups

## Usage

``` r
easyExpPlot(df, x, y, value, type = "line", scaleBubbles = c(2, 6))
```

## Arguments

- df:

  A 3 column data frame that has bee transformed using reshape2::melt
  function.

- x, y:

  x and y variables for drawing.

- value:

  Remaining column after choosing x and y. This column should be numeric
  if type="line" and character if type="bubble".

- type:

  Type of plot. Default ("line").

## Value

A line plot or a bubble plot that represents trends across different
time points and sample types.

## Examples

``` r
if (FALSE) { # \dontrun{
  # Load sample data that contain Z-score transformed TPM values (randomly generated) for some genes
  data(pf3d7TPMs)


  df <- reshape2::melt(pf3d7TPMs,"Probe",na.rm = T)

  easyExpPlot(df,x="variable",y="Probe",value="value", type = "bubble")
  easyExpPlot(df,x="variable",y="value",value="Probe")
} # }

```
