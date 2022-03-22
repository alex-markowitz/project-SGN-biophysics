### PREAMBLE ######################################################################################
library(stats);
setwd('/Users/alex/Library/Mobile Documents/com~apple~CloudDocs/Documents/project-SGN-biophysics');
### LINEAR REGRESSION MODELS ######################################################################
run.lm.models <- function(model.num, impute.missing = FALSE) {

    ### LOAD DATA #################################################################################
    SGN.data        <- read.delim('Derived/2022-03-21_SGN_ephys_morphology_spreadsheet.txt');

    ### FORMAT DATA ###############################################################################
    if (model.num == 1) {
        # Subset data for Model 1
        model.features <- c('CDR', 'Current.Threshold', 'gmax', 'Outward.inactivation.tau');
        predictor.features <- c('Current.Threshold', 'gmax', 'Outward.inactivation.tau');
        model.data <- subset(SGN.data, SGN.data$Age >= 3 & SGN.data$Age <=10 & SGN.data$Spike.Type != 'Graded', select = colnames(SGN.data) %in% model.features);
    }

    else if (model.num == 2) {
        # Subset data for Model 2
        model.features <- c('CDR', 'AP.Latency.threshold', 'Outward.inactivation.tau');
        predictor.features <- c('AP.Latency.threshold', 'Outward.inactivation.tau');
        model.data <- subset(SGN.data, SGN.data$Age >= 3 & SGN.data$Age <=10, select = colnames(SGN.data) %in% model.features);
    }

    else if (model.num == 3) {
        # Subset data for Model 3
        model.features <- c('CDR', 'AP.Latency.threshold');
        predictor.features <- c('AP.Latency.threshold');
        model.data <- subset(SGN.data, SGN.data$Age >= 3 & SGN.data$Age <=5, select = colnames(SGN.data) %in% model.features);
    }

    else if (model.num == 4) {
        # Subset data for Model 4
        model.features <- c('CDR', 'AP.Latency.threshold');
        predictor.features <- c('AP.Latency.threshold');
        model.data <- subset(SGN.data, SGN.data$Age >= 6 & SGN.data$Age <= 8, select = colnames(SGN.data) %in% model.features);
    }

    # Remove rows with missing CDR values
    model.data <- model.data[!(is.na(model.data$CDR)),];

    # Imputation of missing values
    if (impute.missing == TRUE) {
        model.data <- as.data.frame(apply(model.data, 2, function(x) ifelse(is.na(x), median(x), x)));
        }
    else if (impute.missing == FALSE) {
        model.data <- na.omit(model.data);
        }

    ### RUN LINEAR MODELS #########################################################################
    if (model.num == 1) {
        # Fit model
        linear.reg.model <- lm(CDR ~ Current.Threshold + gmax + Outward.inactivation.tau, data = model.data);
        # Extract Beta Coefficients
        model.summary <- as.data.frame(summary(linear.reg.model)$coefficients);
        }

    else if (model.num == 2) {
        # Fit model
        linear.reg.model <- lm(CDR ~ AP.Latency.threshold + Outward.inactivation.tau, data = model.data);
        # Extract Beta Coefficients
        model.summary <- as.data.frame(summary(linear.reg.model)$coefficients);
        }

    else if (model.num == 3) {
        # Fit model
        linear.reg.model <- lm(CDR ~ AP.Latency.threshold, data = model.data);
        # Extract Beta Coefficients
        model.summary <- as.data.frame(summary(linear.reg.model)$coefficients);
        }

    else if (model.num == 4) {
        # Fit model
        linear.reg.model <- lm(CDR ~ AP.Latency.threshold, data = model.data);
        # Extract Beta Coefficients
        model.summary <- as.data.frame(summary(linear.reg.model)$coefficients);
        }

    ### SAVE FILE #################################################################################
    write.csv(model.summary, paste0(Sys.Date(), '_beta.coefficients.model.num.', model.num, '.csv'));
    }
### END ###########################################################################################