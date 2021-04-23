# Load example data
data("spikeindata", package = "MSstatsLOBD")
# linear_quantlim() produces expected output
linear_result = MSstatsLOBD::linear_quantlim(spikeindata)
expect_inherits(linear_result, "data.frame")
expect_equal(colnames(linear_result), 
                       c("CONCENTRATION", "MEAN", "LOW", "UP", "LOB", "LOD",
                         "SLOPE", "INTERCEPT", "NAME", "METHOD"))
expect_equal(as.character(linear_result$METHOD[1]), "LINEAR")
expect_equal(as.character(linear_result$NAME[1]), 
                       as.character(spikeindata$NAME[1]))
expect_inherits(linear_result$MEAN, "numeric")
expect_inherits(linear_result$LOW, "numeric")
expect_inherits(linear_result$UP, "numeric")
expect_inherits(linear_result$LOB, "numeric")
expect_inherits(linear_result$LOD, "numeric")
expect_inherits(linear_result$SLOPE, "numeric")
expect_inherits(linear_result$INTERCEPT, "numeric")
# nonlinear_quantlim() produces expected output
nonlinear_result = MSstatsLOBD::nonlinear_quantlim(spikeindata)
expect_inherits(nonlinear_result, "data.frame")
expect_equal(colnames(nonlinear_result), 
                       c("CONCENTRATION", "MEAN", "LOW", "UP", "LOB", "LOD",
                         "SLOPE", "INTERCEPT", "NAME", "METHOD"))
expect_equal(as.character(nonlinear_result$METHOD[1]), "NONLINEAR")
expect_equal(as.character(nonlinear_result$NAME[1]), 
                       as.character(spikeindata$NAME[1]))
expect_inherits(nonlinear_result$MEAN, "numeric")
expect_inherits(nonlinear_result$LOW, "numeric")
expect_inherits(nonlinear_result$UP, "numeric")
expect_inherits(nonlinear_result$LOB, "numeric")
expect_inherits(nonlinear_result$LOD, "numeric")
expect_inherits(nonlinear_result$SLOPE, "numeric")
expect_inherits(nonlinear_result$INTERCEPT, "numeric")
# plot_quantlim produces expected output - two plots
linear_plot = MSstatsLOBD::plot_quantlim(spikeindata, linear_result, address = FALSE)
expect_inherits(linear_plot, "list")
expect_equal(length(linear_plot), 2)
expect_inherits(linear_plot[[1]], "gg")
expect_inherits(linear_plot[[2]], "gg")
# Theme is OK
expect_equal(linear_plot[[1]]$theme$panel.background$fill, "white")
nonlinear_plot = MSstatsLOBD::plot_quantlim(spikeindata, linear_result, address = FALSE)
expect_inherits(nonlinear_plot, "list")
expect_equal(length(nonlinear_plot), 2)
expect_inherits(nonlinear_plot[[1]], "gg")
expect_inherits(nonlinear_plot[[2]], "gg")
