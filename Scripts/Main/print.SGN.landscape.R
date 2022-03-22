### PREAMBLE ######################################################################################
library(BoutrosLab.plotting.general);
library(BoutrosLab.statistics.general);
library(BoutrosLab.utilities);
setwd('/Users/alex/Library/Mobile Documents/com~apple~CloudDocs/Documents/project-SGN-biophysics');
### PLOT SGN LANDSCAPE ############################################################################
print.SGN.landscape <- function() {
    
    ### LOAD DATA #################################################################################
    SGN.data <- read.delim('Derived/2022-03-16_SGN_ephys_morphology_spreadsheet.txt');

    ### FORMAT DATA ###############################################################################
    SGN.data$CDR[SGN.data$MP.Classification == 'T'] <- -0.6;
    SGN.data <- subset(SGN.data, !(is.na(SGN.data$CDR)));

    SGN.data <- SGN.data[order(desc(SGN.data$CDR)),];

    VClamp.features <- c(
        'Peak.Current',
        'SS.Current',
        'Vmid',
        'g.30',
        'gmax',
        'g30.gmax'
        );

    CClamp.features <- c(
        'Graded.Index',
        'AP.Latency.threshold',
        'AHP.tau',
        'Current.Threshold',
        'Voltage.Threshold',
        'Vm',
        'AP.Height'
        );

    other.features <- c(
        'SGN.Type',
        'Spike.Type',
        'MP.Classification',
        'Age',
        'CDR'
        );
    # Subset data frames by patch clamp configuration
    VClamp.data <- SGN.data[, VClamp.features];
    CClamp.data <- SGN.data[, CClamp.features];

    # Transform data into z-scores
    VClamp.data <- apply(VClamp.data, 2, zscore);
    CClamp.data <- apply(CClamp.data, 2, zscore);

    ### FORMAT COVARIATE INFO ###########################################################################
    covariate.data <- SGN.data[, other.features];
    covariate.data$z.age <- zscore(covariate.data$Age);

    covariate.data$SGN.type.col <- recode.vector(covariate.data$SGN.Type, list(
        'orange' = 'Type I',
        'darkorchid4' = 'Type II'
        ));

    covariate.data$Spike.Type.col <- recode.vector(covariate.data$Spike.Type, list(
        'darkorange1' = 'Slow',
        'dodgerblue' = 'Transient',
        'goldenrod1' = 'Graded',
        'darkgreen' = 'Intermediate'
    ));
    covariate.data$Spike.Type.col[is.na(covariate.data$Spike.Type.col)] <- 'white';

    covariate.data$MP.Classification <- recode.vector(covariate.data$MP.Classification, list(
        'firebrick3' = 'M',
        'lightskyblue' = 'P',
        'lightgoldenrod' = 'T'
    ));

    ### CREATE CDR WATERFALL ###########################################################################
    CDR.barplot <- create.barplot(
        SGN.data$CDR ~ c(1:sum(!is.na(SGN.data$CDR))),
        SGN.data,
        xlab.cex = 0,
        ylab.label = 'CDR',
        ylimits = c(-0.8, 0.8),
        xaxis.tck = 0,
        xaxis.cex = 0,
        use.legacy.settings = TRUE
        );

    ### CREATE COVARIATE BARS ###########################################################################
    covariate.heatmap <- create.heatmap(
        covariate.data[, c('SGN.type.col', 'Spike.Type.col', 'MP.Classification')],
        clustering.method = 'none',
        yaxis.lab = c('SGN Type', 'Spike Type', 'Classification'),
        use.legacy.settings = TRUE,
        input.colours = TRUE,
        print.colour.key = FALSE,
        xaxis.tck = 0,
        yaxis.tck = 0
        );

    age.cdr.heatmap <- create.heatmap(
        covariate.data[, c('z.age', 'z.age')],
        clustering.method = 'none',
        print.colour.key = FALSE,
        yaxis.lab = c('Age', ''),
        use.legacy.settings = TRUE,
        xaxis.tck = 0,
        yaxis.tck = 0
        );

    ### CREATE EPHYS HEATMAPS ###########################################################################

    CC.heatmap <- create.heatmap(
        CClamp.data,
        clustering.method = 'none',
        #cluster.dimensions = 'both',
        yaxis.lab = colnames(CClamp.data),
        xaxis.tck = 0,
        use.legacy.settings = TRUE,
        plot.dendrograms = FALSE,
        print.colour.key = FALSE
        );

    VC.heatmap <- create.heatmap(
        VClamp.data[,rev(colnames(VClamp.data))],
        clustering.method = 'none',
        #cluster.dimensions = 'both',
        yaxis.lab = rev(colnames(VClamp.data)),
        xaxis.tck = 0,
        use.legacy.settings = TRUE,
        plot.dendrograms = TRUE,
        print.colour.key = FALSE
        );

    ### CREATE MULTIPANEL PLOT ###########################################################################
    create.multipanelplot(
        plot.objects = list(
            CDR.barplot,
            covariate.heatmap,
            age.cdr.heatmap,
            VC.heatmap,
            CC.heatmap
        ),
        layout.width = 1,
        layout.height = 5,
        plot.objects.heights = c(
            2,
            2,
            1,
            4,
            4
        ),
        plot.objects.widths = 10
        );
    }
### END ###########################################################################################