# DynamicModelUpdate

R Code to accompany the manuscript "Dynamic Updating of Clinical Survival Prediction Models in a Changing Environment" by Kamaryn T. Tanner, Ruth H. Keogh, Carol AC Coupland, Julia Hippisley-Cox and Karla Diaz-Ordaz.

For questions and comments regarding the code, please contact Kamaryn Tanner (kamaryn.tanner1@lshtm.ac.uk). As code is updated, it will be posted to https://github.com/KamTan/DynamicModelUpdate.

This code was written in R v4.0.2 running on a Windows 10 PC.

We ilustrate the use of this code with the Mayo Clinic Primary Biliary Cholangitis dataset publicly available via the R package survival. Although we are not permitted to make the QResearch data used in the manuscript public, interested researchers may get more information about the data at: https://www.qresearch.org/information/information-for-researchers/

The code is organised in two files:

expSurvivalModel.stan contains the stan model file for the Bayesian dynamic update.

In sampleUpdateCode.R, we first format the pbc data.  We assume that the first 200 records comprise the original development dataset and the remaining records represent new data for updating the original model.  The original model is fit with a Cox proportional hazards model.  Using the "new" data, we illustrate intercept recalibration, refitting and Bayesian dynamic updating.

