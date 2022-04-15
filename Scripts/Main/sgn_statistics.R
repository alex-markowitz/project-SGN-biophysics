### DESCRIPTION ###############################################################
# Statistical assessment of biophysical features of SGN spike subtypes
### PREAMBLE ##################################################################
library(stats);
library(GGally);
library(lme4);
library(plyr);
setwd('/Users/alex/Library/Mobile\ Documents/com~apple~CloudDocs/Documents/project-SGN-biophysics');

### SGN ANALYSIS ##############################################################
### LOAD DATA #################################################################
sgn_data <- read.delim(
    "Derived/2022-03-21_SGN_ephys_morphology_spreadsheet.txt"
    );
### FORMAT DATA ###############################################################
colnames(sgn_data) # Show output
sgn_cc_features <- c("Graded.Index", "Current.Threshold",
                          "Voltage.Threshold", "Vm", "AP.Height",
                          "AP.Latency.threshold");
sgn_vc_features <- c("Peak.Current", "SS.Current", "ss.30", "ss.max", "g.30",
                     "gmax", "Vmid", "g30.gmax", "Outward.inactivation.tau");

sgn_data$SGN.Type[sgn_data$SGN.Type == ""]     <- NA;
sgn_data$Spike.Type[sgn_data$Spike.Type == ""] <- NA;

sgn_data$SGN.Type   <- factor(sgn_data$SGN.Type);
sgn_data$Spike.Type <- factor(sgn_data$Spike.Type);

### RUN STATISTICS ############################################################
## Hypothesis 1: Type I and Type II have different biophysical properties
ggpairs(sgn_data,
        columns = sgn_cc_features,
        aes(color = sgn_data$SGN.Type,
            alpha = 0.5)); ## Show output

# Normality Test
shapiro.test(sgn_data$Vm);
shapiro.test(sgn_data$Vm[sgn_data$SGN.Type == 'Type I']);
shapiro.test(sgn_data$Vm[sgn_data$SGN.Type == 'Type II']);

# Pairwise statistics
wilcox.test(formula = Vm ~ SGN.Type,
            data = sgn_data,
            na.action = "na.omit"); # Show output

# Write function
sgn_type_wilcox_stats <- function(feature) {
    wilcox.stat <- wilcox.test(formula = sgn_data[, feature] ~ SGN.Type,
                               data = sgn_data,
                               na.action = "na.omit");
    stat_summary <- data.frame(
        feature = feature,
        wilcox.p = wilcox.stat$p.value,
        wilcox.statistic = wilcox.stat$statistic
        );
    return(stat_summary);
    }

## Hypothesis 2: Type I Spike Types have differential properties
type_I_sgn_data <- subset(sgn_data, sgn_data$SGN.Type == 'Type I');
type_I_sgn_data <- subset(sgn_data, sgn_data$Spike.Type != 'Graded');

aov_object <- aov(Current.Threshold ~ Spike.Type, data = type_I_sgn_data);

TukeyHSD(aov_object); # Shot output

# Write function
spike_type_tukeyhsd_stats <- function(feature) {
    
    type_I_sgn_data <- subset(sgn_data, sgn_data$SGN.Type == 'Type I');
    type_I_sgn_data <- subset(sgn_data, sgn_data$Spike.Type != 'Graded');
    
    aov_object <- aov(type_I_sgn_data[, feature] ~ Spike.Type, data = type_I_sgn_data);
    
    tukeyhsd_stats <- as.data.frame(TukeyHSD(aov_object)$Spike.Type);
    tukeyhsd_stats$pairwise <- row.names(tukeyhsd_stats);
    tukeyhsd_stats$feature <- feature;
    
    return(tukeyhsd_stats);
    }

## Hypothesis 3: Biophysical properties are associated w/ CDR
# By SGN Type
aov_cdr <- aov(CDR ~ Spike.Type, data = sgn_data);
TukeyHSD(aov_cdr); # Shot output

# Current clamp properties
lm_object <- lm(CDR ~ Current.Threshold, data = sgn_data);
summary(lm_object);

# Write function
sgn_cdr_lm_stats <- function(feature) {
    
    type_I_sgn_data <- subset(sgn_data, sgn_data$SGN.Type == 'Type I');
    aov_cdr <- aov(CDR ~ type_I_sgn_data[, feature], data = type_I_sgn_data);
    aov_summary <- summary(aov_cdr)
    
    stat_summary <- data.frame(
        feature = feature,
        aov.p = aov_summary[[1]]$`Pr(>F)`[1],
        aov.statistic = aov_summary[[1]]$`F value`[1]
        );
    return(stat_summary);
}

## Hypothesis 4: Biophysical properties are associated w/ age of animal
lm_object <- lm(AP.Latency ~ Age, data = sgn_data);
summary(lm_object);

# Model Comparison between AB and A model
# Format to omit NAs from data, update to impute medians
sgn_data_na_omit <- sgn_data[!(is.na(sgn_data$CDR) |
                                   is.na(sgn_data$AP.Latency.threshold)),];

AB.lm <- lm(CDR ~ Age + AP.Latency.threshold, data = sgn_data_na_omit);
A.lm  <- lm(CDR ~ Age + 1, data = sgn_data_na_omit);

# Tests whether the addition of effect provides additional information
lm_pval             <- as.data.frame(anova(AB.lm,A.lm))$`Pr(>F)`[2];
eta_squared         <- eta_squared(AB.lm)[1,2]; # Effect size
within.cluster.var  <- vcov(AB.lm)[2,2];

# Write function
sgn_cdr_lm_age_ancova <- function(feature) {
    
    sgn_data_na_omit <- sgn_data[!(is.na(sgn_data$CDR) |
                                       is.na(sgn_data[, feature])),];
    
    AB.lm <- lm(CDR ~ Age + sgn_data_na_omit[, feature], data = sgn_data_na_omit);
    A.lm  <- lm(CDR ~ Age + 1, data = sgn_data_na_omit);
    
    # Tests whether the addition of effect provides additional information
    lm_pval             <- as.data.frame(anova(AB.lm, A.lm))$`Pr(>F)`[2];
    eta_squared         <- eta_squared(AB.lm)[1,2]; # Effect size
    within_cluster_var  <- vcov(AB.lm)[2,2];
    
    stat_summary <- data.frame(
        feature = feature,
        lm_pval = lm_pval,
        eta_squared = eta_squared,
        within_cluster_var = within_cluster_var
        );
    return(stat_summary);
    }

### ITERATE STATISTICS #########################################################
cc_sgn_type_stats      <- NULL;
cc_spike_type_stats    <- NULL;
cc_sgn_cdr_stats       <- NULL;
cc_sgn_cdr_age_stats   <- NULL;

vc_sgn_type_stats      <- NULL;
vc_spike_type_stats    <- NULL;
vc_sgn_cdr_stats       <- NULL;
vc_sgn_cdr_age_stats   <- NULL;

# Gather CC stats
for (i in 1:length(sgn_cc_features)) {
    feature_sgn_type_stats      <- sgn_type_wilcox_stats(sgn_cc_features[i]);
    feature_spike_type_stats    <- spike_type_tukeyhsd_stats(sgn_cc_features[i]);
    feature_sgn_cdr_stats       <- sgn_cdr_lm_stats(sgn_cc_features[i]);
    feature_cdr_age_stats       <- sgn_cdr_lm_age_ancova(sgn_cc_features[i]);
    
    cc_sgn_type_stats      <- rbind.fill(cc_sgn_type_stats, feature_sgn_type_stats);
    cc_spike_type_stats    <- rbind.fill(cc_spike_type_stats, feature_spike_type_stats);
    cc_sgn_cdr_stats       <- rbind.fill(cc_sgn_cdr_stats, feature_sgn_cdr_stats);
    cc_sgn_cdr_age_stats   <- rbind.fill(cc_sgn_cdr_age_stats, feature_cdr_age_stats);
    }

# Gather VC stats
for (i in 1:length(sgn_vc_features)) {
    feature_sgn_type_stats      <- sgn_type_wilcox_stats(sgn_vc_features[i]);
    feature_spike_type_stats    <- spike_type_tukeyhsd_stats(sgn_vc_features[i]);
    feature_sgn_cdr_stats       <- sgn_cdr_lm_stats(sgn_vc_features[i]);
    feature_cdr_age_stats       <- sgn_cdr_lm_age_ancova(sgn_vc_features[i]);
    
    vc_sgn_type_stats      <- rbind.fill(vc_sgn_type_stats, feature_sgn_type_stats);
    vc_spike_type_stats    <- rbind.fill(vc_spike_type_stats, feature_spike_type_stats);
    vc_sgn_cdr_stats       <- rbind.fill(vc_sgn_cdr_stats, feature_sgn_cdr_stats);
    vc_sgn_cdr_age_stats   <- rbind.fill(vc_sgn_cdr_age_stats, feature_cdr_age_stats);
    }
### SAVE FILE ##################################################################
write.table(cc_sgn_type_stats, 'Derived/cc_sgn_type_stats.txt', sep = '\t');
write.table(cc_spike_type_stats, 'Derived/cc_spike_type_stats.txt', sep = '\t');
write.table(cc_sgn_cdr_stats, 'Derived/cc_sgn_cdr_stats.txt', sep = '\t');
write.table(cc_sgn_cdr_age_stats, 'Derived/cc_sgn_cdr_age_stats.txt', sep = '\t');

write.table(vc_sgn_type_stats, 'Derived/vc_sgn_type_stats.txt', sep = '\t');
write.table(vc_spike_type_stats, 'Derived/vc_spike_type_stats.txt', sep = '\t');
write.table(vc_sgn_cdr_stats, 'Derived/vc_sgn_cdr_stats.txt', sep = '\t');
write.table(vc_sgn_cdr_age_stats, 'Derived/vc_sgn_cdr_age_stats.txt', sep = '\t');
### END ########################################################################
