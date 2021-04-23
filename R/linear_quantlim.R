#' Calculation of the LOB and LOD with a  linear fit
#' 
#' This function calculates the value of the LOB (limit of blank) and LOD (limit of detection) 
#' from the (Concentration, Intensity) spiked in data.
#' The function also returns the values of the linear curve fit that allows it to be plotted. 
#' At least 2 blank samples (characterized by Intensity = 0) are required by this function
#' which are used to calculate the background noise. T
#' he LOB is defined as the concentration at which the value of the linear fit is equal to the 95\% upper bound of the noise. 
#' The LOD is the concentration at which the latter is equal to the 90\% lower bound of the prediction interval (5\% quantile) of the linear fit. 
#' A weighted linear fit is used with weights for every unique concentration proportional to the inverse of variance between replicates.
#' 
#' @export
#' @import minpack.lm 
#' @param datain Data frame that contains the input data. 
#' The input data frame has to contain the following columns : CONCENTRATION, 
#' INTENSITY (both of which are measurements from the spiked in experiment) 
#' and NAME which designates the name of the assay (e.g. the name of the peptide or protein)
#' @param alpha Probability level to estimate the LOB/LOD
#' @param Npoints Number of points to use to discretize the concentration line between 0 and the maximum spiked concentration
#' @param Nbootstrap Number of bootstrap samples to use to calculate the prediction interval of the fit. 
#' This number has to be increased for very low alpha values or whenever very accurate assay characterization is required.
#' @param num_changepoint_samples Number of bootstrap samples for the prediction 
#' inteval for the changepoint. Large values can make calculations very expensive
#' @param num_prediction_samples Number of prediction samples to generate
#' @param max_iter Number of trials for convergence of every curvefit algorithm
#' @return data.frame, It contains the following columns: 
#' i) CONCENTRATION: Concentration values at which the value of the fit is calculated 
#' ii) MEAN: The value of the curve fit 
#' iii) LOW: The value of the lower bound of the 95\% prediction interval 
#' iv) UP: The value of the upper bound of the 95\% prediction interval 
#' v) LOB: The value of the LOB (one column with identical values) 
#' vi) LOD: The value of the LOD (one column with identical values) 
#' vii) SLOPE: Value of the slope of the linear curve fit where only the spikes above LOD are considered 
#' viii) INTERCEPT: Value of the intercept of the linear curve fit where only the spikes above LOD are considered 
#' ix) NAME: The name of the assay (identical to that provided in the input) 
#' x) METHOD which is always set to LINEAR when this function is used. 
#' Each line of the data frame corresponds to a unique concentration value 
#' at which the value of the fit and prediction interval are evaluated. 
#' More unique concentrations values than in the input data frame are used to increase the accuracy of the LOB/D calculations.
#' @details The LOB and LOD can only be calculated when more than 2 blank samples are included.
#' The data should ideally be plotted using the companion function plot_quantlim to ensure that a linear fit is suited to the data.	
#' @examples
#' # Consider data from a spiked-in contained in an example dataset. This dataset contains 
#' # a significant threshold at low concentrations that is not well captured by a linear fit
#' 
#' head(spikeindata)
#' 
#' linear_quantlim_out <- linear_quantlim(spikeindata)


linear_quantlim = function(datain, alpha = 0.05, Npoints = 100, 
                           Nbootstrap = 500, num_changepoint_samples = 30,
                           num_prediction_samples = 200, max_iter = 30) {
  if(alpha >= 1 | alpha  <= 0) {
    stop("incorrect specified value for alpha,  0 < alpha < 1")
  }
  switch(Sys.info()[['sysname']],
         Windows = {null_output = "NUL"},
         Linux  = {null_output = "/dev/null"},
         Darwin = {null_output = "/dev/null"})
  
  names(datain)[names(datain) == 'CONCENTRATION'] = 'C'
  names(datain)[names(datain) == 'INTENSITY'] = 'I' 
  datain = datain[!is.na(datain$I) & !is.na(datain$C), ]  
  datain = datain[!is.infinite(datain$I) & !is.infinite(datain$C), ]  
  datain = datain[order(datain$C), ]
  tmp_nob = datain[datain$C > 0, ] 
  tmp_all = datain
  tmp_blank = datain[datain$C == 0, ]
  
  # Calculate the value of the noise:
  # Use the zero concentration to calculate the LOD:
  pb = txtProgressBar(min = 0, max = 1, initial = 0, char = "%",
                      width = 40, title, label, style = 1, file = "")
  # Prediction interval: gaussian, unknown parameters.
  n_blank = length(unique(tmp_blank$I))
  noise = mean(tmp_blank$I)
  var_noise = var(tmp_blank$I)
  if(nrow(tmp_blank) <= 1 || var_noise  <= 0) {
    stop("Not enough blank samples!!!")
  }
  fac = qt(1-alpha,n_blank - 1)*sqrt(1+1/n_blank)
  # upper bound of noise prediction interval
  up_noise = noise   +  fac* sqrt(var_noise)
  unique_c = sort(unique(tmp_all$C))   
  var_v = rep(0.0, length(unique_c))
  weights  = rep(0.0, length(tmp_all$C))
  weights_nob  = rep(0.0, length(tmp_nob$C))
  
  # Calculate variance for all concentrations:
  ii = 1
  for (j in unique_c){
    data_f = subset(tmp_all, C == j)
    var_v[ii] = var(data_f$I)
    ii = ii +1
  }
  
  #Log scale discretization:
  xaxis_orig_2 = exp(c(seq(from = log(10 + 0), to = log(1 + max(unique_c)), 
                           by = log(1 + max(unique_c)) / Npoints))) - 10 # 0.250 to go fast here
  xaxis_orig_2 = unique(sort(c(xaxis_orig_2, unique_c))) 
  ###Loop here to make sure that the discretization is fine enough
  # Create a  piecewise linear approximation:
  var_v_lin = approx(unique_c[!is.na(var_v)], 
                     var_v[!is.na(var_v)], xout = xaxis_orig_2)$y
  var_v_lin_unique = approx(unique_c[!is.na(var_v)], 
                            var_v[!is.na(var_v)], xout = unique_c)$y
  # In the following, consider that the smoothed variance is that from the log:
  var_v_s = var_v_lin
  var_v_s_unique = var_v_lin_unique
  var_v_s_log = var_v_s
  var_v_s_log_unique = var_v_s_unique
  
  outB = matrix(NA_real_, nrow=Nbootstrap, ncol=length(xaxis_orig_2))
  outBB_pred = matrix(NA_real_, nrow=Nbootstrap*num_prediction_samples, ncol=length(xaxis_orig_2))
  change_B = rep(NA, Nbootstrap)
  for (j in 1:Nbootstrap) {
    setTxtProgressBar(pb, j / Nbootstrap, title = NULL, label = NULL)
    lin.blank_B = NULL
    tmpB = tmp_all[sample(1:length(tmp_all$C), replace=TRUE), ] 
    # Pick **  observations with replacement among all that are available.
    # Blank samples are included
    weights = rep(0,length(tmpB$C))
    for (kk in 1:length(tmpB$C)) {
      weights[kk] = 1 / var_v_s_unique[which( unique_c == tmpB$C[kk])]
    }   
    noise_B = mean(sample(tmp_blank$I, length(tmp_blank$I), replace=TRUE))
    
    ll = 0
    while(ll < max_iter) { 
      ll = ll + 1
      slope = median(tmpB$I)/median(tmpB$C)*runif(1) 
      intercept = noise * runif(1)
      sink(null_output)
      lin.blank_B <- tryCatch({
        nlsLM(I ~ .linear(C , intercept, slope),
              data = tmpB, trace = TRUE,
              start = c(intercept = intercept, slope = slope), 
              weights = weights,
              control = nls.lm.control(nprint = 1,
                                       ftol = sqrt(.Machine$double.eps) / 2, 
                                       maxiter = 50))
      }, error = function(e) NULL)
      sink()
      if(!is.null(lin.blank_B)) break
    }
    
    # Store the curve fits obtained via bootstrap with bilinear and linear:
    if(!is.null(lin.blank_B)){ #If linear fit, change = 0 anyway
      outB[j, ]  = .linear(xaxis_orig_2, summary(lin.blank_B)$coefficient[1], 
                           summary(lin.blank_B)$coefficient[2])
      change_B[j] = 0
    } else {
      outB[j, ]  = rep(NA, length(xaxis_orig_2))
      change_B[j] = NA
    }
    
    if(!is.null(lin.blank_B)) {
      for (jj in 1:num_prediction_samples){ # Calculate predictions
        outBB_pred[(j - 1) * num_prediction_samples + jj, ] = outB[j, ] + 
          rnorm(length(xaxis_orig_2), 0, sqrt(var_v_s_log))
      }
    }  else {
      for (jj in 1:num_prediction_samples){ # Calculate predictions
        outBB_pred[(j - 1) * num_prediction_samples + jj, ] = NA
      }
    }
  }
  
  #Calculate the variance of the fits: 
  var_bilinear = apply(outB, 2, var,na.rm = TRUE)
  mean_bilinear = apply(outB, 2, mean,na.rm = TRUE)
  #Calculate confidence interval based on quantiles:
  mean_pred = apply(outBB_pred, 2, mean, na.rm = TRUE) 
  var_pred = apply(outBB_pred, 2, var, na.rm = TRUE) 
  lower_Q_pred = apply(outBB_pred, 2, quantile, probs = alpha, na.rm = TRUE)
  upper_Q_pred = apply(outBB_pred, 2, quantile, probs = 1 - alpha, na.rm = TRUE)      
  
  #Calculate the LOD/LOQ from the prediction interval:
  i_before = which(diff(sign(up_noise - mean_bilinear)) != 0) #before sign change
  if (length(i_before) > 0) {
    i_after = i_before + 1
    x1 = xaxis_orig_2[i_before] 
    f1 = (up_noise - mean_bilinear)[i_before]
    x2 = xaxis_orig_2[i_after]
    f2 = (up_noise - mean_bilinear)[i_after]
    #Linear interpolation to find where the function changes sign:
    LOD_pred_mean =  x1 - f1 * (x2 - x1) / (f2 - f1)
    LOD_pred = LOD_pred_mean
    y_LOD_pred = up_noise
  } else { 
    LOD_pred = 0
    y_LOD_pred = up_noise
  }
  
  i_before = which(diff(sign(up_noise - lower_Q_pred)) != 0)
  if (length(i_before) > 0) {
    i_after = i_before+1
    x1 = xaxis_orig_2[i_before] 
    f1 = (up_noise - lower_Q_pred)[i_before]
    x2 = xaxis_orig_2[i_after] 
    f2 = (up_noise - lower_Q_pred)[i_after]
    x_inter = x1 - f1 * (x2 - x1) / (f2 - f1)
    LOQ_pred = x_inter
    y_LOQ_pred = up_noise
  } else {
    LOQ_pred = 0
    y_LOQ_pred = up_noise
  }
  
  data_linear = datain[datain$C > LOQ_pred,]
  ii = 0
  while(ii < max_iter){
    ii = ii + 1
    { 
      intercept = noise*runif(1)
      slope = median(data_linear$I) / median(data_linear$C) * runif(1) 
      weights = rep(0, length(data_linear$C))
      for (kk in 1:length(data_linear$C)) {
        weights[kk] = 1 / var_v_s[which(unique_c == data_linear$C[kk])] 
      } 
      sink(null_output)
      fit.blank_lin = NULL 
      fit.blank_lin = tryCatch({
        nlsLM(I ~ .linear(C, intercept, slope),
              data = data_linear, trace = TRUE,
              start = c(slope = slope, intercept = intercept), 
              weights = weights,
              control = nls.lm.control(nprint = 1, ftol = sqrt(.Machine$double.eps) / 2, 
                                       maxiter = 50))
      }, error = function(e) NULL)
      sink()
    }
    if(!is.null(fit.blank_lin)) break
  }
  
  if(!is.null(fit.blank_lin)){
    slope_lin = summary(fit.blank_lin)$coefficients[[1]]
    intercept_lin = summary(fit.blank_lin)$coefficients[[2]]
  } else {
    slope_lin = NA
    intercept_lin = NA
    message("Linear assay characterization above LOD did not converge")  
  }

  data.frame(CONCENTRATION = xaxis_orig_2, 
             MEAN = mean_bilinear,
             LOW = lower_Q_pred, 
             UP = upper_Q_pred, 
             LOB = rep(LOD_pred, length(upper_Q_pred)),  
             LOD = rep(LOQ_pred, length(upper_Q_pred)),
             SLOPE = slope_lin, 
             INTERCEPT = intercept_lin, 
             NAME = rep(datain$NAME[1], length(upper_Q_pred)),
             METHOD = rep("LINEAR", length(upper_Q_pred)))
}
                               
