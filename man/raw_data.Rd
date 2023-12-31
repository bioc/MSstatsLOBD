% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dataDoc.R
\docType{data}
\name{raw_data}
\alias{raw_data}
\title{Example of dataset that contains spike in data for 43 distinct peptides.}
\format{A data frame with 3870 rows and 16 variables.}
\usage{
raw_data
}
\description{
This example dataset is from CPTAC(Clinical Proteomic Tumor Analysis Consortium) assay portan (Thomas and others 2015).
The dataset contains spike in data for 43 distinct peptides.
For each peptide, 8 distinct concentration spikes for 3 different replicates are measured.
The Skyline files for the assay along with details about the experiment can be obtained from this webpage:
https://assays.cancer.gov. The particular dataset examined here (called JHU_DChan_HZhang_ZZhang) can
be found at https://panoramaweb.org/labkey/project/CPTAC\%20Assay\%20Portal/JHU_DChan_HZhang_
ZZhang/Serum_QExactive_GlycopeptideEnrichedPRM/begin.view?. It should be downloaded from the
MSStats website http://msstats.org/?smd_process_download=1&download_id=548. The data is then
exported in a csv file (calibration_data_raw.csv) from Skyline.
The csv file contains the measured peak area for each fragment of each light
and heavy version of each peptide. Depending on the format of the Skyline
file and depending on whether standards were used, the particular outputs obtained in the csv file may vary.
In this particular case the following variables are obtained in the output file calibration_data_raw.csv:
File Name, Sample Name, Replicate Name, Protein Name, Peptide Sequence, Peptide Modified Sequence,
Precursor Charge, Product Charge, Fragment Ion, Average Measured Retention Time, SampleGroup, IS Spike,
Concentration, Replicate, light Area, heavy Area. A number of variables are byproducts of the acquisition
process and will not be considered for the following, i.e. File Name, Sample Name, Replicate Name, SampleGroup, IS Spike.
Variables that are important for the assay characterization are detailed below
(others are assumed to be self explanatory):
}
\details{
\itemize{
\item Pepdidesequence Name of the peptide sequence
\item Concentration Value of the known spiked concentration in pmol.
\item Replicate Number of the technical replicate
\item light Area Peak area of the light (measured)
\item heavy Area Peak area of the heavy (reference) peptide
}
}
\examples{
head(raw_data)

}
\keyword{datasets}
