### DESCRIPTION ###################################################################################
# This function predicts contact distance ratio (CDR) of Type I SGN based on a list of features
# explored in Markowitz & Kalluri (2020).
# This functions uses the R package mlr3 and mlr3viz for machine learning protocols.
# Documentation for mlr3 can be found at https://mlr3book.mlr-org.com/index.html

# INSTRUCTIONS ON INPUTS:
# Enter model number ID (1-4) based on Markowitz & Kalluri (2020)
# Enter directory path to input csv
# Please format input csv of unlabelled data into four separate csv files based on the model number (model.num)

# Model 1 (model.num = 1) :
# Column names of input csv: 'cell.ID', Current.Threshold', 'gmax', 'Outward.inactivation.tau'

# Model 2 (model.num = 2) :
# Column names of input csv: 'cell.ID', 'AP.Latency.threshold', 'Outward.inactivation.tau'

# Model 3 (3 <= Age <= 5) (model.num = 3) :
# Column names of input csv: 'cell.ID', 'AP.Latency.threshold'

# Model 4 (6 <= Age <= 8)(model.num = 4) :
# Column names of input csv: 'cell.ID', 'AP.Latency.threshold'

# EXPECTED OUTPUT
# The output of this function will present the RSME of the linear regression in the terminal
# An autoplot of the pair-wise correlations of predictor variables per model will be generated
# An output csv file of predicted CDR will be generated with CellIDs as row names and colname Pred.CDR

# EXAMPLE USAGE
# predict.SGN.CDR.lm.reg(model1.num = 1, 'Derived/model1.unlabelled.data.sample.csv')

### PREAMBLE ######################################################################################
library(mlr3verse);
library(mlr3);
setwd('/Users/alex/Library/Mobile Documents/com~apple~CloudDocs/Documents/project-SGN-biophysics');
### PLOT SGN LANDSCAPE ############################################################################
predict.SGN.CDR.lm.reg <- function(model.num, path.to.unlabelled.csv, impute.missing = FALSE) {

    ### LOAD DATA #################################################################################
    SGN.data        <- read.delim('Derived/2022-03-21_SGN_ephys_morphology_spreadsheet.txt');
    unlabeled.data <- read.csv(unlabelled.data.path);

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

    ## MACHINE LEARNING ###########################################################################
    # Using Regression to predict CDR
    task        <- as_task_regr(model.data, target = 'CDR');
    learner     <- lrn('regr.lm');
    resampling  <- rsmp('cv', folds = 3); # 3-folds is standard
    rr          <- resample(task, learner, resampling, store_models = TRUE);
    pred <- learner$train(task)$predict(task);

    # Extracted model results
    model.results <- data.frame(
        RMSE = rr$aggregate(msr('regr.rmse'))
        );

    # RMSE is your best measurement of model accuracy
    cat(paste0('Root mean squared error (RMSE): ', model.results$RMSE));

    # Visualize Pair correlations
    autoplot(task, type = 'pairs');

    ### PREDICT WELLS #############################################################################
    contact.predictions <- predict(learner, unlabeled.data[, predictor.features]);

    predictions.dataframe <- data.frame(
        cell.ID = row.names(unlabeled.data),
        Pred.CDR = contact.predictions
        );

    ### SAVE FILE #################################################################################
    write.csv(predictions.dataframe, paste0(Sys.Date(), '_Predicted.CDR_model.num.', model.num, '.csv'));
    write.csv(beta.coefficients, paste0(Sys.Date(), '_beta.coefficients.model.num.', model.num, '.csv'))
    }
### END ###########################################################################################
