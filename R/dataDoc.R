
#' Example of dataset that contains spike in data for 43 distinct peptides.
#'
#' This example dataset is from CPTAC(Clinical Proteomic Tumor Analysis Consortium) assay portan (Thomas and others 2015).
#' The dataset contains spike in data for 43 distinct peptides. 
#' For each peptide, 8 distinct concentration spikes for 3 different replicates are measured.
#' The Skyline files for the assay along with details about the experiment can be obtained from this webpage:
#' https://assays.cancer.gov. The particular dataset examined here (called JHU_DChan_HZhang_ZZhang) can
#' be found at https://panoramaweb.org/labkey/project/CPTAC%20Assay%20Portal/JHU_DChan_HZhang_
#' ZZhang/Serum_QExactive_GlycopeptideEnrichedPRM/begin.view?. It should be downloaded from the
#' MSStats website http://msstats.org/?smd_process_download=1&download_id=548. The data is then
#' exported in a csv file (calibration_data_raw.csv) from Skyline.
#' The csv file contains the measured peak area for each fragment of each light 
#' and heavy version of each peptide. Depending on the format of the Skyline
#' file and depending on whether standards were used, the particular outputs obtained in the csv file may vary.
#' In this particular case the following variables are obtained in the output file calibration_data_raw.csv:
#' File Name, Sample Name, Replicate Name, Protein Name, Peptide Sequence, Peptide Modified Sequence,
#' Precursor Charge, Product Charge, Fragment Ion, Average Measured Retention Time, SampleGroup, IS Spike,
#' Concentration, Replicate, light Area, heavy Area. A number of variables are byproducts of the acquisition
#' process and will not be considered for the following, i.e. File Name, Sample Name, Replicate Name, SampleGroup, IS Spike.
#' Variables that are important for the assay characterization are detailed below 
#' (others are assumed to be self explanatory):
#'
#' \itemize{
#'   \item Pepdidesequence Name of the peptide sequence
#'   \item Concentration Value of the known spiked concentration in pmol.
#'   \item Replicate Number of the technical replicate
#'   \item light Area Peak area of the light (measured)
#'   \item heavy Area Peak area of the heavy (reference) peptide
#' }
#'
#' @format A data frame with 3870 rows and 16 variables.
#' @examples
#' head(raw_data)
#'
"raw_data"


#' Example of normalized datasets from raw_data,
#' 
#' We normalize the intensity of the light peptides using that of the heavy peptides. 
#' This corrects any systematic errors that can occur during a run or across replicates. 
#' The calculation is greatly simplified by the use of the tidyr and dplyr packages. 
#' The area from all the different peptide fragments is first summed then log transformed. 
#' The median intensity of the reference heavy peptides medianlog2heavy is calculated.
#' Their intensities should ideally remain constant across runs 
#' since the spiked concentration of the heavy peptide is constant. 
#' The difference between the median for all the heavy peptide spikes is calculated. 
#' It is then used to correct (i.e. to normalize) the intensity of the light peptides log2light 
#' to obtain the adjusted intensity log2light_norm. 
#' The intensity is finally converted back to original space.
#' Details are available in vignette.
#' The variables are as follows:
#'
#' \itemize{
#'   \item CONCENTRATION: Concentration values at which the value of the fit is calculated 
#'   \item MEAN: The value of the curve fit 
#'   \item LOW: The value of the lower bound of the 95\% prediction interval 
#'   \item UP: The value of the upper bound of the 95\% prediction interval 
#'   \item LOB: The value of the LOB (one column with identical values) 
#'   \item LOD: The value of the LOD (one column with identical values) 
#'   \item SLOPE: Value of the slope of the linear curve fit where only the spikes above LOD are considered 
#'   \item INTERCEPT: Value of the intercept of the linear curve fit where only the spikes above LOD are considered 
#'   \item NAME: The name of the assay (identical to that provided in the input) 
#'   \item METHOD which is always set to NONLINEAR when this function is used. 
#'   \item Each line of the data frame corresponds to a unique concentration value 
#' at which the value of the fit and prediction interval are evaluated. 
#'   \item More unique concentrations values than in the input data frame are used to increase the accuracy of the LOB/D calculations.
#' }
#'
#' @format A data frame with 30 rows and 4 variables.
#' @examples
#' head(spikeindata)
#'
"spikeindata"
