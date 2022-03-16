### PREAMBLE ######################################################################################

setwd('/Users/alex/Library/Mobile Documents/com~apple~CloudDocs/Documents/USC/PythonCoding/');

### FORMAT SGN SPREADSHEET FOR PLOTTING ###########################################################
format.SGN.data <- function() {

    ### LOAD DATA #################################################################################
    SGN.data <- read.csv('SGN_Spreadsheet.csv');
    ### FORMAT DATA ###############################################################################
    
    # Extract features of interest for landscape plot and analysis
    features.of.interest <- c(
        'Classification',
        'SGN.Type',
        'Age',
        'AP.Latency..ms.',
        'Code..Graded.Index',
        'Code..Cthres',
        'Spike.Type',
        'Branch.Number',
        'Peak.Current',
        'SS.Current', 
        'CDR..Type.I.only.', 
        'AHP.Time.Constant',
        'Code..AP.Height',
        'Code..Rin', 
        'Code..Vthres',
        'Code..Vm',
        'ss.30NC',
        'ssmax_NC',
        'g.30NC',
        'gmax_NC',
        'g30.gmax_NC',
        'Vmid...NC',
        'Code..Lat.at.Thresh'
        );
    
    # Subset features of interest
    SGN.data <- subset(SGN.data, select = colnames(SGN.data) %in% features.of.interest);
    
    # Rename colnames
    colnames(SGN.data) <- c(
        'Graded.Index',
        'Peak.Current',
        'SS.Current',
        'Branch.Number',
        'SGN.Type',
        'Spike.Type',
        'MP.Classification',
        'Age',
        'CDR',
        'AP.Latency.threshold',
        'AHP.tau',
        'Current.Threshold',
        'Vm',
        'AP.Height',
        'Rin',
        'ss.30',
        'ss.max',
        'g.30',
        'gmax',
        'Vmid',
        'g30.gmax',
        'AP.Latency'
        );
    
    ### SAVE FILE #################################################################################
    write.table(SGN.data, '2022-03-16_SGN_ephys_morphology_spreadsheet.txt', sep = '\t');
    }
### END ###########################################################################################
