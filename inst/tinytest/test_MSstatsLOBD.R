# Load example data
data("spikeindata", package = "MSstatsLOBD")
# linear_quantlim() produces expected output
linear_result = MSstatsLOBD::linear_quantlim(spikeindata)
tinytest::expect_inherits(linear_result, "data.frame")
tinytest::expect_equal(colnames(linear_result), 
                       c("CONCENTRATION", "MEAN", "LOW", "UP", "LOB", "LOD",
                         "SLOPE", "INTERCEPT", "NAME", "METHOD"))
tinytest::expect_equal(as.character(linear_result$METHOD[1]), "LINEAR")
tinytest::expect_equal(as.character(linear_result$NAME[1]), 
                       as.character(spikeindata$NAME[1]))
tinytest::expect_inherits(linear_result$MEAN, "numeric")
tinytest::expect_inherits(linear_result$LOW, "numeric")
tinytest::expect_inherits(linear_result$UP, "numeric")
tinytest::expect_inherits(linear_result$LOB, "numeric")
tinytest::expect_inherits(linear_result$LOD, "numeric")
tinytest::expect_inherits(linear_result$SLOPE, "numeric")
tinytest::expect_inherits(linear_result$INTERCEPT, "numeric")
# nonlinear_quantlim() produces expected output
nonlinear_result = MSstatsLOBD::nonlinear_quantlim(spikeindata)
tinytest::expect_inherits(nonlinear_result, "data.frame")
tinytest::expect_equal(colnames(nonlinear_result), 
                       c("CONCENTRATION", "MEAN", "LOW", "UP", "LOB", "LOD",
                         "SLOPE", "INTERCEPT", "NAME", "METHOD"))
tinytest::expect_equal(as.character(nonlinear_result$METHOD[1]), "NONLINEAR")
tinytest::expect_equal(as.character(nonlinear_result$NAME[1]), 
                       as.character(spikeindata$NAME[1]))
tinytest::expect_inherits(nonlinear_result$MEAN, "numeric")
tinytest::expect_inherits(nonlinear_result$LOW, "numeric")
tinytest::expect_inherits(nonlinear_result$UP, "numeric")
tinytest::expect_inherits(nonlinear_result$LOB, "numeric")
tinytest::expect_inherits(nonlinear_result$LOD, "numeric")
tinytest::expect_inherits(nonlinear_result$SLOPE, "numeric")
tinytest::expect_inherits(nonlinear_result$INTERCEPT, "numeric")
# plot_quantlim produces expected output - two plots
linear_plot = MSstatsLOBD::plot_quantlim(spikeindata, linear_result, address = FALSE)
tinytest::expect_inherits(linear_plot, "list")
tinytest::expect_equal(length(linear_plot), 2)
tinytest::expect_inherits(linear_plot[[1]], "gg")
tinytest::expect_inherits(linear_plot[[2]], "gg")
# Theme is OK
tinytest::expect_equal(linear_plot[[1]]$theme$panel.background$fill, "white")
nonlinear_plot = MSstatsLOBD::plot_quantlim(spikeindata, linear_result, address = FALSE)
tinytest::expect_inherits(nonlinear_plot, "list")
tinytest::expect_equal(length(nonlinear_plot), 2)
tinytest::expect_inherits(nonlinear_plot[[1]], "gg")
tinytest::expect_inherits(nonlinear_plot[[2]], "gg")
