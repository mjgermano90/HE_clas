#########################################################################
################################################################################
# 1-1-1 Analisis de datos mediante classification of patients based on H&E TIL results ====
#########################################################################

#### Load library and functions ====
getwd()

library(DESeq2)
library(SummarizedExperiment)
library(org.Hs.eg.db) #Cambiar el simbolo de genes
library(ggrepel)
library(RColorBrewer)
library(survival)
library(survminer)
library(pheatmap)
library(tidyr)
library(stats)
library(dplyr)
library(rcompanion) # For posthoc.kruskal.dunn.test
library(mosaic) #for mosaic plot
library(chisq.posthoc.test) # to know which groups differ among them in proportion chisq test
library(TCGAbiolinks)
library(maftools)
library(stringr)
library(fgsea)
library(knitr) #to plot tables
library(readxl)
library(ggplot2)
library(EnhancedVolcano)
library(gridExtra)
library(grid)
library(patchwork) 
library(forcats)
library(ggpubr)
library(openxlsx)
library(broom)
library(coxphf)
library(forestmodel)
library(car)
library(DescTools)
library(forestplot)


pathLocalDb <- function(x) {
  completePath <- file.path(getwd(), "02_data", x)
  return(completePath)
}


pathLocalResults <- function(x) {
  completePath <- file.path(getwd(), "03_results", x)
  return(completePath)
}



#libray(survival)
# Function to perform multivariable Cox analysis
perform_multivariable_cox <- function(clinical_data, survival_object, col_indices) {
  # Create the formula string
  formula_str <- paste("surv_obj ~", paste(colnames(clinical_data)[col_indices], collapse = " + "))
  
  # Convert the string to a formula object
  formula_obj <- as.formula(formula_str)
  
  # Fit the Cox proportional hazards model
  multivar_cox <- coxph(formula_obj, data = clinical_data)
  
  # Return the fitted model
  return(multivar_cox)
}

create_grouped_violin_plot <- function(data, group, column_names, y_limit = NULL, annotations = NULL, y_text = NULL) {
  # Reshape the data into long format
  data_long <- gather(data, key = "variable", value = "value", all_of(column_names))
  
  # Mantener el orden de las variables tal como aparecen en 'column_names'
  data_long$variable <- factor(data_long$variable, levels = column_names)
  
  # Mantener el orden de los grupos como aparecen en el dataset
  group_order <- unique(data[[group]])
  data_long[[group]] <- factor(data_long[[group]], levels = group_order)
  
  # Create a color palette for fill
  pastel_palette <- brewer.pal(n = nlevels(data_long[[group]]), name = "Set2")
  
  # Create the grouped violin plot
  plot <- ggplot(data_long, aes(x = as.factor(.data[[group]]), y = value, fill = as.factor(.data[[group]]))) +
    geom_violin(scale = "width", trim = TRUE, adjust = 0.8, linewidth = 0.2) + 
    geom_boxplot(width = 0.08, fill = "white", size = 0.3, outlier.shape = NA, linewidth = 0.2, color = "black") + 
    labs(title = "Grouped Violin Plots",
         x = "Groups",
         y = y_text,
         fill = "Groups") +
    #theme_minimal() +
    theme_bw(base_size = 8) +
    theme(axis.text.x = element_blank()) +
    facet_wrap(~ variable, ncol = 4 ) + #ncol = length(column_names)
    scale_fill_manual(values = pastel_palette) +
    scale_color_manual(values = pastel_palette)
  
  # Add annotations if provided
  
  if (!is.null(annotations)) {
    for (ann in annotations) {
      plot <- plot + 
        geom_text(data = subset(data_long, variable == ann$variable),
                  aes(x = ann$x, y = ann$y),
                  label = strrep("*", ann$stars), 
                  color = "black", size = 6)
      
    }
  }
  
  # Apply y-axis limits if provided
  if (!is.null(y_limit)) {
    plot <- plot + scale_y_continuous(limits = c(0, y_limit))
  }
  
  print(plot)
}


create_grouped_boxplot <- function(data, group, column_names, y_limit = NULL, annotations = NULL, y_text = NULL) {
  # Reshape the data into long format
  data_long <- gather(data, key = "variable", value = "value", all_of(column_names))
  
  # Mantener el orden de las variables tal como aparecen en 'column_names'
  data_long$variable <- factor(data_long$variable, levels = column_names)
  
  # Mantener el orden de los grupos como aparecen en el dataset
  group_order <- unique(data[[group]])
  data_long[[group]] <- factor(data_long[[group]], levels = group_order)
  
  # Create a color palette for fill
  pastel_palette <- brewer.pal(n = nlevels(as.factor(data[[group]])), name = "Set2")
  
  # Create the grouped boxplot
  plot <- ggplot(data_long, aes(x = as.factor(.data[[group]]), y = value, fill = as.factor(.data[[group]]))) +
    geom_boxplot(outlier.shape = NA, linewidth = 0.3, width = 0.4, color = "black", position = position_dodge(width = 0.4)) +
    #geom_jitter(width = 0.12, size = 0.8, alpha = 0.6) +
    labs(title = group,
         x = "Groups",
         y = y_text,
         fill = "Groups") +
    #theme_minimal() +
    theme_bw(base_size = 8) +
    theme(
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      axis.text.y = element_text(size = 7),
      axis.title.y = element_text(size = 8),
      strip.text = element_text(size = 8, face = "bold"),
      legend.position = "none",
      panel.grid = element_blank()
    ) +
    facet_wrap(~ variable, ncol = 2) +  #ncol = length(column_names)
    scale_fill_manual(values = pastel_palette)
  
  # Add annotations if provided
  if (!is.null(annotations)) {
    for (ann in annotations) {
      plot <- plot + 
        geom_text(data = subset(data_long, variable == ann$variable),
                  aes(x = ann$x, y = ann$y),
                  label = strrep("*", ann$stars), 
                  color = "black", size = 6)
    }
  }
  
  # Apply y-axis limits if provided
  if (!is.null(y_limit)) {
    plot <- plot + scale_y_continuous(limits = c(0, y_limit))
  }
  
  print(plot)
}

create_grouped_barplot_median_sem <- function(data, group, column_names, y_limit = NULL, annotations = NULL, y_text = NULL) {
  # Reshape the data into long format
  data_long <- gather(data, key = "variable", value = "value", all_of(column_names))
  
  # Mantener el orden de las variables tal como aparecen en 'column_names'
  data_long$variable <- factor(data_long$variable, levels = column_names)
  
  # Mantener el orden de los grupos como aparecen en el dataset
  group_order <- unique(data[[group]])
  data_long[[group]] <- factor(data_long[[group]], levels = group_order)
  
  # Calcular mediana y SEM
  summary_data <- data_long %>%
    group_by(variable, .data[[group]]) %>%
    summarise(
      median = median(value, na.rm = TRUE),
      sem = sd(value, na.rm = TRUE), # / sqrt(n()),
      .groups = "drop"
    ) %>%
    rename(group_value = !!group)
  
  # Crear paleta de colores
  pastel_palette <- brewer.pal(n = nlevels(as.factor(data[[group]])), name = "Set2")
  
  # Crear gráfico
  plot <- ggplot(summary_data, aes(x = group_value, y = median, fill = group_value)) +
    geom_bar(stat = "identity", position = position_dodge(), width = 0.6) +
    geom_errorbar(aes(ymin = median - sem, ymax = median + sem), width = 0.2, position = position_dodge(0.6)) +
    facet_wrap(~ variable, ncol = length(column_names)) +
    labs(title = "Grouped Bar Plots (Median ± SEM)",
         x = "Groups",
         y = y_text,
         fill = "Groups") +
    theme_minimal() +
    theme(
      axis.text.x = element_blank(),
      strip.text = element_text(angle = 0, hjust = 1)
    ) +
    scale_fill_manual(values = pastel_palette)
  
  # Añadir anotaciones si están disponibles
  if (!is.null(annotations)) {
    for (ann in annotations) {
      plot <- plot + 
        geom_text(data = subset(summary_data, variable == ann$variable),
                  aes(x = ann$x, y = ann$y),
                  label = strrep("*", ann$stars), 
                  color = "black", size = 6)
    }
  }
  
  # Aplicar límites del eje Y si se proveen
  if (!is.null(y_limit)) {
    plot <- plot + scale_y_continuous(limits = c(0, y_limit))
  }
  
  print(plot)
}


perform_normality_tests <- function(data_frame, columns, group_column = NULL) {
  results <- list()  # Initialize an empty list to store results
  
  for (col in columns) {
    if (is.null(group_column)) {
      # Perform normality test without grouping
      normal_test <- shapiro.test(data_frame[[col]])
      results[[col]] <- list(
        overall = list(
          statistic = normal_test$statistic,
          p_value = normal_test$p.value
        )
      )
    } else {
      # Perform normality test within each group
      group_levels <- unique(data_frame[[group_column]])
      results[[col]] <- list()
      
      for (level in group_levels) {
        subset_data <- data_frame[data_frame[[group_column]] == level, col, drop = FALSE]
        normal_test <- shapiro.test(subset_data[[col]])
        results[[col]][[as.character(level)]] <- list(
          statistic = normal_test$statistic,
          p_value = normal_test$p.value
        )
      }
    }
  }
  
  return(results)
}


# wilcoxon_test_function <- function(data, column, condition_column, condition_value) {
#   # Extract values based on condition
#   condition_true <- data[data[[condition_column]] == condition_value, column]
#   condition_false <- data[data[[condition_column]] != condition_value, column]
#   
#   # Perform Wilcoxon test
#   wilcox_result <- wilcox.test(condition_true, condition_false)
#   
#   # Return the result
#   return(wilcox_result)
# }


wilcoxon_test_function <- function(data, column, condition_column, condition_value) {
  
  # Extract values
  condition_true <- data[data[[condition_column]] == condition_value, column]
  condition_false <- data[data[[condition_column]] != condition_value, column]
  
  # Remove NA
  condition_true <- condition_true[!is.na(condition_true)]
  condition_false <- condition_false[!is.na(condition_false)]
  
  # Sample sizes
  n1 <- length(condition_true)
  n2 <- length(condition_false)
  
  # Wilcoxon test
  wilcox_result <- wilcox.test(condition_true, condition_false)
  
  # U statistic (equivalente al W transformado)
  U <- wilcox_result$statistic
  
  # Rank-biserial correlation
  r_rb <- 1 - (2 * U) / (n1 * n2)
  
  return(list(
    test = wilcox_result,
    effect_size = r_rb,
    n1 = n1,
    n2 = n2
  ))
}

# # Call Wilcoxon funcion
# all_results <- list()
# for (col in normality_cols) {
#   # Perform Wilcoxon test for each column
#   wilcox_result <- wilcoxon_test_function(clinical, col, "immune_LumAB_clas", "low")
#   all_results[[col]] <- wilcox_result
#   
#   # Print results
#   cat("Wilcoxon test results for", col, ":\n")
#   print(paste0("p-value: ", all_results[[col]][["p.value"]]))
# }


#### Load and pre-process data ====
# Cargo rda
load(pathLocalDb("TCGA_BRCA_RNAseq_tumor_female_unique_immune.rda"))


#Elimino los que tienen tratamiento previo
TCGA_BRCA_RNAseq_tumor_female_unique_immune <- TCGA_BRCA_RNAseq_tumor_female_unique_immune[, colData(TCGA_BRCA_RNAseq_tumor_female_unique_immune)$prior_treatment != "Yes" ]


#TCGA_BRCA_RNAseq_tumor_female_unique_immune@colData$paper_BRCA_Subtype_PAM50 <- as.factor(TCGA_BRCA_RNAseq_tumor_female_unique_immune@colData$paper_BRCA_Subtype_PAM50)
## Filter out rda to keep only Luminal A and B
#TCGA_BRCA_RNAseq_tumor_female_unique_immune@colData@listData[["PAM50_Subtype"]] <- as.factor(TCGA_BRCA_RNAseq_tumor_female_unique_immune@colData@listData[["PAM50_Subtype"]])
TCGA_BRCA_LumAB_RNAseq <- TCGA_BRCA_RNAseq_tumor_female_unique_immune[, colData(TCGA_BRCA_RNAseq_tumor_female_unique_immune)$PAM50_Subtype == "LumA"| colData(TCGA_BRCA_RNAseq_tumor_female_unique_immune)$PAM50_Subtype == "LumB" ]
rm(TCGA_BRCA_RNAseq_tumor_female_unique_immune)
colData(TCGA_BRCA_LumAB_RNAseq)$PAM50_Subtype <- relevel(
  as.factor(colData(TCGA_BRCA_LumAB_RNAseq)$PAM50_Subtype),
  ref = "LumA"
)

clinical <- as.data.frame(TCGA_BRCA_LumAB_RNAseq@colData)
clinical[] <- lapply(clinical, function(x) {
  if (is.character(x) || is.factor(x)) {
    x <- as.character(x)
    x[x == "Indetermined"] <- NA
    return(x)
  }
  x
})



clinical$PAM50_Subtype <- as.factor(as.character(clinical$PAM50_Subtype))
levels(clinical$PAM50_Subtype)

#Los datos de porcentaje de TIL fueron descargados del paper  https://www.sciencedirect.com/science/article/pii/S2211124718304479#app2
# https://stonybrookmedicine.app.box.com/v/cellreportspaper/file/283628537637
# Tumor-Infiltrating Lymphocytes Maps from TCGA H&E Whole Slide Pathology Images

TCGA_BRCA_TIL_HE <-read.csv(pathLocalDb("brca_H&E.csv"))
#Using tidyverse to add a column with patient names
#library(dplyr)
#library(stringr)

TCGA_BRCA_TIL_HE <- TCGA_BRCA_TIL_HE %>%
  mutate(patient = str_extract(Slides, "^[^-]+-[^-]+-[^-]+"))

# Some patient has more than 1 sample of slide, select one for each sample based on the highest number of data points
TCGA_BRCA_TIL_HE<- TCGA_BRCA_TIL_HE %>%
  group_by(patient) %>%
  slice_max(order_by = number.of.data.points, n = 1) %>%
  ungroup()

# Select only "patient" and "til_percentage"
#Using library(dplyr)
library(dplyr)
TCGA_BRCA_TIL_HE  <- TCGA_BRCA_TIL_HE  %>% 
  select(patient, til_percentage)

# Merge with clinical data frame based on "patient"
clinical <- clinical %>%
  left_join(TCGA_BRCA_TIL_HE, by = "patient")


#################################################################################
#### Univariable Cox for HE TIL percentage, immune score and cell population  ====

### Survival Object 

censored_time <- pmin(clinical$overall_survival, 3650)
event_indicator <- clinical$vital_status_binary & (clinical$overall_survival <= 3650)

#creamos un objeto "surv" con la funci?n Surv()
surv_obj <- Surv(censored_time,
                 event_indicator)


# surv_obj <- Surv(clinical$censored_time,
#                  clinical$event_indicator)

#Obtain pvalue and coefficient for each cell 
immune_cell_columns <- colnames(clinical)[c(94:105,121)] #ojo 121 es para TIL si no he clasificado en high y low en funcion de Immune score
univ_cell_results <- list()
for (col in immune_cell_columns) {
  formula <- as.formula(paste("surv_obj ~", col))
  univ_cox <- coxph(formula, data = clinical)
  summary_univ_cox <- summary(univ_cox)
  p_value <- summary_univ_cox[["logtest"]][["pvalue"]]
  CI <- summary_univ_cox$conf.int[, c("lower .95", "upper .95")]
  HR <- summary_univ_cox$coefficients[1, "exp(coef)"]
  univ_cell_results[[col]] <- list(p_value = p_value, HR = HR, CI=CI)
}


# Convert the list to a data frame
cox_results_df <- do.call(rbind, lapply(names(univ_cell_results), function(col) {
  c(Cell_Type = col, 
    P_Value = univ_cell_results[[col]]$p_value, 
    HR = univ_cell_results[[col]]$HR,
    CI = univ_cell_results[[col]]$CI)
}))

# Convert columns to appropriate types
cox_results_df <- data.frame(cox_results_df, stringsAsFactors = FALSE)
cox_results_df$P_Value <- as.numeric(cox_results_df$P_Value)
cox_results_df$HR <- as.numeric(cox_results_df$HR)
cox_results_df$CI.lower..95 <- as.numeric(cox_results_df$CI.lower..95)
cox_results_df$CI.upper..95 <- as.numeric(cox_results_df$CI.upper..95)

# Sort the data frame by P_Value
cox_results_df <- cox_results_df[order(cox_results_df$P_Value), ]

# Print the sorted results
print(cox_results_df)


write.csv2(cox_results_df, file=pathLocalResults("Univariable Cox for immune score til percentages and quantiseq proportions patients without prior treatment.csv"), row.names = FALSE)


#########################################################################
#### Classification of patients based on til percentage ====

## Elimino los pacientes que no tienen dato para til_percentage
library(dplyr)
clinical <- clinical %>% filter(!is.na(til_percentage))

##### Comparar diferentes cut off

# Initialize a new column for immune classification
clinical$HE_clas <- NA

## Esto que aparece abajo, ser� ignorado al correr el codigo, es como comentarlo
if (FALSE) {
censored_time <- clinical$censored_time
clinical$event_indicator_binary <- ifelse(clinical$event_indicator == TRUE, 1,0)
event_indicator_binary <- clinical$event_indicator_binary

# Define cut off
surv_cut_point <- surv_cutpoint(
  clinical,
  time = "censored_time",
  event = "event_indicator_binary",
  colnames(clinical[122]), #122 es la columna con til_percentage
  minprop = 0.1,
  progressbar = TRUE
)

HE_clas_cutoff <-summary(surv_cut_point)$cutpoint 
HE_clas_cutoff#0.2338073 -> ES MUY BAJO, VER HISTOGRAMA Y CLASIFICAR POR LA MEDIANA
}


range(clinical$til_percentage, na.rm = TRUE)
hist(clinical$til_percentage, breaks = seq(0,35,1))
# "M1" = "#80b1d3",
# "M2" = "#b3de69", 
# "M3" = "#bc80bd", 
# "M4" = "#66c2a5",
# "M5" = "#fdb462"

p <- ggplot(clinical, aes(x = til_percentage)) +
  # Add histogram with proper styling
  geom_histogram(binwidth = 1, fill = "#66c2a5", color = "white", alpha = 1) +
  # Add smoothed density line (optional)
  geom_density(aes(y = ..density.. * 500), color = "black", size = 0.5) +
  # Use a clean theme
  theme_minimal(base_size = 14) +
  # Customize labelshttp://127.0.0.1:11019/graphics/81f4ede7-3f24-4861-af5a-dcdfa847eb91.png
  labs(title = "Distribution of TIL percentage", x = "Tumor-Infiltrating Lymphocyetes (%)", y = "Frequency") +
  # Refine theme for publication
  theme(
    plot.title = element_text(face = "bold"), panel.grid.minor = element_blank(), axis.line = element_line(color = "black")
  )

# 3. View the plot
print(p)

# 4. Save the plot for publication (High resolution)
ggsave(pathLocalResults("HE clas Results/FigureS1_histogram.pdf"), p, width = 8, height = 6, device = cairo_pdf)




summary(clinical$til_percentage)


sum(clinical$til_percentage > 10, na.rm = TRUE) #50

HE_clas_cutoff <- median(clinical$til_percentage, na.rm = TRUE)

#Columna HE_clas
clinical$HE_clas <- ifelse(
  clinical$til_percentage >= HE_clas_cutoff,
  "high", "low"
)
clinical$HE_clas <- as.factor(clinical$HE_clas)

# por cox

surv_obj <- Surv(clinical$censored_time,
                 clinical$event_indicator)

cox_fit <- coxph(Surv(censored_time, event_indicator) ~ HE_clas, data = clinical)
summary_cox <- summary(cox_fit)
summary_cox 


# Interpretation
# Call:
# coxph(formula = Surv(censored_time, event_indicator) ~ HE_clas, 
#     data = clinical)
# 
#   n= 669, number of events= 85 
# 
#              coef exp(coef) se(coef)   z Pr(>|z|)       exp(coef) is HR, se(coef) es la desviacion estandar del coef, z is Wald test statistic: higher absolute value means more significant.
# HE_claslow 0.8301    2.2934   0.2306 3.6 0.000319 ***
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
#            exp(coef) exp(-coef) lower .95 upper .95
# HE_claslow     2.293      0.436      1.46     3.604
# 
# Concordance= 0.6  (se = 0.03 )
# Likelihood ratio test= 13.76  on 1 df,   p=2e-04
# Wald test            = 12.96  on 1 df,   p=3e-04
# Score (logrank) test = 13.68  on 1 df,   p=2e-04
# 
# 
# Concordance (C-index) measures how well the model discriminates between patients with better vs worse survival.Range: 0.5 (no better than chance) to 1 (perfect prediction).
# Likelihood ratio test: Compares full vs. null model using log-likelihood. Significant improvement.
# Wald test: Tests if coef / se(coef) is different from 0. Confirms significance of effect.
# Score (Logrank) test: Based on rank data. Also confirms variable is significantly associated with surviva




# Evaluate different cut off for HE clas and generate the results data frame

# Define your cutoffs (fixed the sequence - removed 'by' parameter since you have specific values)
cutoffs <- c(0.5, 1, 1.074332, 2.0, 2.5, 3.0, 3.5, 4, 4.5, 5.0)
#cutoffs <- c(1.0, 1.074717, 2.0,3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10)

clinical$HE_clas <- relevel(clinical$HE_clas, ref = "low")

cox_cut_results <- lapply(cutoffs, function(cut) {
  # Create the dichotomized variable
  clinical$HE_clas <- ifelse(clinical$til_percentage >= cut, "high", "low") #use high as reference by default
  #clinical$HE_clas <- factor(clinical$HE_clas, levels = c("low", "high")) # si quiero cambiar el default, debo colocarlo dentro del loop
  
  # Get counts for each group
  n_high <- sum(clinical$HE_clas == "high")
  n_low <- sum(clinical$HE_clas == "low")
  
  # Fit Cox model
  cox_fit <- coxph(Surv(censored_time, event_indicator) ~ HE_clas, data = clinical)
  summary_cox <- summary(cox_fit)
  #print(summary_cox) # concondance index para cut off de 2, da menor (0.561), mientras que usando la media mejora un poco a 0.599
  
  # Return results with counts
  data.frame(
    cutoff = cut,
    n_high = n_high,
    n_low = n_low,
    HR = summary_cox$coefficients[2],
    lower_CI = summary_cox$conf.int[3],
    upper_CI = summary_cox$conf.int[4],
    p_value = summary_cox$coefficients[5],
    stringsAsFactors = FALSE
  )
})

# Combine results into one data frame
cox_cut_results_df <- do.call(rbind, cox_cut_results)

# View the results
cox_cut_results_df

write.csv2(cox_cut_results_df, file=pathLocalResults("HE clas Results/Univarable_Cox_cut_offs_without_prior_treatment.csv"), row.names = FALSE)


# Assuming you have a data frame with cutoff values, HRs, and p-values
library(forestplot)

# Create forest plot
# Create a text table
tabletext <- cbind(
  c("Cutoff", sprintf("%.2f", cox_cut_results_df$cutoff)),
  c("High (n)", as.character(cox_cut_results_df$n_high)),
  c("Low (n)", as.character(cox_cut_results_df$n_low)),
  c("HR (95% CI)", 
    sprintf("%.2f (%.2f–%.2f)", 
            cox_cut_results_df$HR, 
            cox_cut_results_df$lower_CI, 
            cox_cut_results_df$upper_CI)),
  c("p-value", 
    formatC(cox_cut_results_df$p_value, format = "e", digits = 2))
)

## Add til% continue variable to forest plot
cox_results_df_til <- cox_results_df[2,]
new_row <- matrix(c("TIL%", "-", "-", "0.94 (0.88-1.00)", "2.52e-02"), nrow = 1) #valores extraidos de cox_results_df_til 
# Combine: first header row, then new row, then rest
tabletext <- rbind(tabletext[1, , drop = FALSE], new_row, tabletext[-1, , drop = FALSE])

# Check result
tabletext

#Add a new row to cox_cut_results_df
# Get names from original df
col_names <- names(cox_cut_results_df)
# Create new row with same names and types (NA where needed)
new_row_df <- data.frame(
  cutoff = "TIL%",
  n_high = NA,
  n_low = NA,
  HR = 0.94,
  lower_CI = 0.88,
  upper_CI = 1.00,
  p_value = 2.52e-02,
  stringsAsFactors = FALSE
)

# Reorder to match original df
new_row_df <- new_row_df[, col_names]
# Combine
cox_cut_results_df <- rbind(new_row_df, cox_cut_results_df)

# Check result
cox_cut_results_df

dir.create(dirname(pathLocalResults("HE clas Results/Figure_1a.pdf")),
           recursive = TRUE,
           showWarnings = FALSE)

#Para figuras hechas cin forest (o que no sea con ggplot) definir el tamaño primero
cairo_pdf(pathLocalResults("HE clas Results/Figure_1a.pdf"),
          width = 9/2.54,
          height = 9/2.54)

tabletext[,1] <- stringr::str_wrap(tabletext[,1], width = 15)

forest <- forestplot(
  labeltext = tabletext,
  mean  = c(NA, cox_cut_results_df$HR),
  lower = c(NA, cox_cut_results_df$lower_CI),
  upper = c(NA, cox_cut_results_df$upper_CI),
  zero = 1,
  boxsize = 0.15,
  xlog = TRUE,
  colwidths = unit(c(4, 1.5, 2.5, 1.5), "cm"),   # <<< definir el tamaño de columnas
  txt_gp = fpTxtGp(                   # reduce el tamaño del texto
    label = gpar(cex = 0.5),
    ticks = gpar(cex = 0.5),
    xlab  = gpar(cex = 0.6),
    title = gpar(cex = 0.7)
  ),
  colgap = unit(2, "mm"), #reduce el tamaño entre columnas
  col = fpColors(box = "royalblue", line = "darkblue", summary = "royalblue"),
  title = "Univariable Cox HR by Cutoff",
  xlab = "Hazard Ratio (log scale)")



print(forest)
dev.off()


##### Eleccion del cut off y generacion de columna HE_clas
#Columna HE_clas

clinical$HE_clas <- NA
HE_clas_cutoff <- median(clinical$til_percentage, na.rm = TRUE)

clinical$HE_clas <- ifelse(
  clinical$til_percentage >= HE_clas_cutoff,
  "high", "low"
)
clinical$HE_clas <- as.factor(clinical$HE_clas)

#########################################################################
##### Curva Kaplan Meier entre high y low HE_clas (basado en TIL) ====
surv_obj <- Surv(clinical$censored_time,
                 clinical$event_indicator)

curva_immune<- survfit(surv_obj ~ HE_clas, data = clinical)
summary(curva_immune)

surv_results <- survdiff(surv_obj ~ HE_clas, data = clinical, rho=0) #rho= 0 indica que el test es Log-Rank test

#Creo una paleta de colores para el fill
pastel_palette <- brewer.pal(n = nlevels(clinical$HE_clas), name = "Set2")

# Curva immune
curva<- ggsurvplot(
           curva_immune, 
           data = clinical,
           size = 0.8,
           palette = pastel_palette,
           censor.shape = '|', censor.size = 2.5, #disminuyo el tamaño de la lina de censor
           conf.int = TRUE,
           pval = TRUE,
           pval.size = 6,
           pval.method = TRUE,
           risk.table = TRUE,
           risk.table.col = 'strata',
           legend.labs = list('0' = 'high', '1' = 'low'),
           risk.table.height = 0.25,
           xlim = c(0, 3650),
           break.time.by = 500,
           ggtheme = theme_bw(base_size = 7)) # base_size = 7  Eso reescala TODO.

curva$plot <- curva$plot +
  theme(
    legend.position = "top",
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 7),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 7)
  )

curva$table <- curva$table +
  theme(
    text = element_text(size = 6),
    axis.title = element_text(size = 7),
    axis.text = element_text(size = 6)
  )



dir.create(dirname(pathLocalResults("HE clas Results/Figure_1b.pdf")),
           recursive = TRUE,
           showWarnings = FALSE)
cairo_pdf(pathLocalResults("HE clas Results/Figure_1b.pdf"),
          width = 10/2.54,
          height = 9/2.54)
print(curva)
dev.off()


######################################################################
### Association between categorical variables =====

#### Comparison of stage proportions among node status ====

# Create a contingency table
contingency_table <- table(clinical$node_status, clinical$stage)
summary(contingency_table)


# Perform Chi-square test
chi_sq_test <- fisher.test(contingency_table) #ERROR NO SE PUEDE PORQUE NO HAY CASOS

#### Comparson of PAM50 subtypes proportions among node status ====

# Create a contingency table
contingency_table <- table(clinical$node_status, clinical$PAM50_Subtype)
summary(contingency_table)

# Perform Chi-square test
chi_sq_test <- fisher.test(contingency_table)

# View the results
print(chi_sq_test)


##### Plot PAM50 proportions
# Plot contingency table, adding shade, control the colors in function of residuals from contingency table in the graph 
mosaicplot(contingency_table,
           shade = c(1,2), 
           dir = "h",
           main="PAM50 proportions between node_satus",
           sub = paste("p-value by Chi-squared: ", round(chi_sq_test$p.value, 5)))


#### Comparson of PAM50 subtypes proportions among tumor size ====

# Create a contingency table
contingency_table <- table(clinical$tumor_size, clinical$PAM50_Subtype)
summary(contingency_table)

# Perform Chi-square test
chi_sq_test <- chisq.test(contingency_table)

# View the results
print(chi_sq_test)

#### Comparson of PAM50 subtypes proportions among stage ====

# Create a contingency table
contingency_table <- table(clinical$stage, clinical$PAM50_Subtype)
summary(contingency_table)

# Perform Chi-square test
chi_sq_test <- fisher.test(contingency_table)

# View the results
print(chi_sq_test)


##### Plot PAM50 proportions
# Plot contingency table, adding shade, control the colors in function of residuals from contingency table in the graph 
mosaicplot(contingency_table,
           shade = c(1,2), 
           dir = "h",
           main="PAM50 proportions among stages",
           sub = paste("p-value by Chi-squared: ", round(chi_sq_test$p.value, 5)))

#### Comparson of PAM50 subtypes proportions between high and low HE groups ====

# Create a contingency table
contingency_table <- table(clinical$HE_clas, clinical$PAM50_Subtype)
summary(contingency_table)

# Perform Chi-square test
chi_sq_test <- chisq.test(contingency_table)

# View the results
print(chi_sq_test)


##### Plot PAM50 proportions
# Plot contingency table, adding shade, control the colors in function of residuals from contingency table in the graph 
mosaicplot(contingency_table,
           shade = c(1,2), 
           dir = "h",
           main="PAM50 proportions between Immune groups",
           sub = paste("p-value by Chi-squared: ", round(chi_sq_test$p.value, 5)))


#### Comparson of node_status proportions between high and low HE groups ====

# Create a contingency table
contingency_table <- table(clinical$HE_clas, clinical$node_status)
summary(contingency_table)

# Perform Chi-square test
chi_sq_test <- chisq.test(contingency_table)

# View the results
print(chi_sq_test)


#### Comparson of stage proportions between high and low HE groups ====

# Create a contingency table
contingency_table <- table(clinical$HE_clas, clinical$stage)
summary(contingency_table)

# Perform Chi-square test
chi_sq_test <- chisq.test(contingency_table)

# View the results
print(chi_sq_test)

#### Comparson of tumor size proportions between high and low HE groups ====

# Create a contingency table
contingency_table <- table(clinical$HE_clas, clinical$tumor_size)
summary(contingency_table)

# Perform Chi-square test
chi_sq_test <- chisq.test(contingency_table)

# View the results
print(chi_sq_test)


##### Plot node_status proportions
# Plot contingency table, adding shade, control the colors in function of residuals from contingency table in the graph 
mosaicplot(contingency_table,
           shade = c(1,2), 
           dir = "h",
           main="Node Status proportions between Immune groups",
           sub = paste("p-value by Chi-squared: ", round(chi_sq_test$p.value, 3)))


# Association between node status and immune class with a direction
library(DescTools)

cochran_test <- CochranArmitageTest(contingency_table)
print(cochran_test)


#### Comparson of stage proportions between high and low HE groups ====

# Create a contingency table
contingency_table <- table(clinical$HE_clas, clinical$stage)
summary(contingency_table)

# Perform Chi-square test
chi_sq_test <- fisher.test(contingency_table)

# View the results
print(chi_sq_test)

# Plot contingency table, adding shade, control the colors in function of residuals from contingency table in the graph 
mosaicplot(contingency_table,
           shade = c(1,2), 
           dir = "h",
           main="Stage proportions between Immune groups",
           sub = paste("p-value by Fisher's test: ", round(chi_sq_test$p.value, 3)))

# Association between stage and immune class with a direction (higher infiltation -> higher stage? or viceversa)
cochran_test <- CochranArmitageTest(contingency_table)
print(cochran_test)


#### Comparson of histological type proportions between high and low HE groups ====

clinical$paper_BRCA_Pathology <- as.factor(as.character(clinical$paper_BRCA_Pathology))
levels(clinical$paper_BRCA_Pathology)
summary(clinical$paper_BRCA_Pathology)

# Create a contingency table

#Creo un nuevo data frame sin los "NA" -> OJO ESTAN COMO FACTOR
clinical_clean <- clinical[
  !is.na(clinical$paper_BRCA_Pathology) & 
    clinical$paper_BRCA_Pathology != "NA", 
]

clinical_clean$paper_BRCA_Pathology <- droplevels(clinical_clean$paper_BRCA_Pathology)

contingency_table <- table(clinical_clean$HE_clas, clinical_clean$paper_BRCA_Pathology)

summary(contingency_table)
contingency_table

# Perform Chi-square test
chi_sq_test <- chisq.test(contingency_table)

# View the results
print(chi_sq_test)

# Plot contingency table, adding shade, control the colors in function of residuals from contingency table in the graph 
mosaicplot(contingency_table,
           shade = c(1,2), 
           dir = "h",
           main="Pathology proportions between Immune groups",
           sub = paste("p-value by Chi-squared: ", round(chi_sq_test$p.value, 4)))

### ER, PR, HER2 status ====
hr_her2_status <- read.csv2(pathLocalDb("clinical_status.csv"), sep =",")

hr_her2_status <- hr_her2_status %>% 
  select(patient, er_status_by_ihc, pr_status_by_ihc, her2_status_ihc_fish)

clinical <- left_join(clinical, hr_her2_status, by= "patient")


clinical[c("er_status_by_ihc", "pr_status_by_ihc", "her2_status_ihc_fish")] <- lapply(clinical[c("er_status_by_ihc", "pr_status_by_ihc", "her2_status_ihc_fish")], factor)

clinical$hr_status <- ifelse (
  clinical$pr_status_by_ihc == "Positive"  | clinical$er_status_by_ihc == "Positive", "Positive",
  ifelse (clinical$pr_status_by_ihc == "Negative"  & clinical$er_status_by_ihc == "Negative" , "Negative",
          NA)
)

clinical$her2_status_ihc_fish_y <- ifelse (
  clinical$her2_status_ihc_fish  == "Positive"  | clinical$her2_status_ihc_fish  == "Negative", as.character(clinical$her2_status_ihc_fish),
  NA
)
clinical$her2_status_ihc_fish_y <- as.factor(clinical$her2_status_ihc_fish_y)

clinical$pr_status_by_ihc<- ifelse (
  clinical$pr_status_by_ihc  == "Positive"  | clinical$pr_status_by_ihc  == "Negative", as.character(clinical$pr_status_by_ihc),
  NA
)
clinical$pr_status_by_ihc <- as.factor(clinical$pr_status_by_ihc)

clinical$er_status_by_ihc<- ifelse (
  clinical$er_status_by_ihc  == "Positive"  | clinical$er_status_by_ihc  == "Negative", as.character(clinical$er_status_by_ihc),
  NA
)
clinical$er_status_by_ihc <- as.factor(clinical$er_status_by_ihc)

# Asociation between immune group and PR, ER, HER2 status
#ER
contingency_table <- table(
  droplevels(clinical$HE_clas[clinical$er_status_by_ihc != "[Not Evaluated]"]),
  droplevels(clinical$er_status_by_ihc[clinical$er_status_by_ihc != "[Not Evaluated]"])
)
contingency_table

summary(contingency_table)

# Perform Chi-square test
chi_sq_test <- fisher.test(contingency_table)

# View the results
print(chi_sq_test)

#PR
contingency_table <- table(
  droplevels(clinical$HE_clas[clinical$pr_status_by_ihc != "[Not Evaluated]" & clinical$pr_status_by_ihc != "Indeterminate"]),
  droplevels(clinical$pr_status_by_ihc[clinical$pr_status_by_ihc != "[Not Evaluated]" & clinical$pr_status_by_ihc != "Indeterminate"])
)
contingency_table

summary(contingency_table)

# Perform Chi-square test
chi_sq_test <- fisher.test(contingency_table)

# View the results
print(chi_sq_test)

#HER2
contingency_table <- table(
  droplevels(clinical$HE_clas[clinical$her2_status_ihc_fish == "Positive" | clinical$her2_status_ihc_fish == "Negative"]),
  droplevels(clinical$her2_status_ihc_fish[clinical$her2_status_ihc_fish == "Positive" | clinical$her2_status_ihc_fish == "Negative"])
)
contingency_table
summary(contingency_table)


######################################################################
### Analisis Cox Multivariado ====

# Debido a que "age_at_index" viola el supuesto de Proportional Hazard para el modelo de Cox los pacientes seran estratificados
# Crear categorias

# Define cut off
surv_cut_point <- surv_cutpoint(
  clinical,
  time = "censored_time",
  event = "event_indicator",
  colnames(clinical[48]), #48 edad
  minprop = 0.2,
  progressbar = TRUE
)

age_cutoff <-summary(surv_cut_point)$cutpoint 
age_cutoff#0.2338073 -> ES MUY BAJO, VER HISTOGRAMA Y CLASIFICAR POR LA MEDIANA

clinical$age_group <- as.factor(ifelse(clinical$age_at_index < 70, "< 70", "> 70"))
levels(clinical$age_group)
clinical$age_group <- relevel(clinical$age_group, ref = "< 70")


### Multivariable Cox model for  TIL (continuous) and HE clas ====

#### Multivariable Cox model for TIL coninuous ====

surv_obj <- Surv(clinical$censored_time,
                 clinical$event_indicator)
modelos_til <- list(
  MS1 = surv_obj ~ til_percentage + tumor_size,
  MS2 = surv_obj ~ til_percentage + stage,
  MS3 = surv_obj ~ til_percentage + PAM50_Subtype,
  MS4 = surv_obj ~ til_percentage + age_group + tumor_size,
  MS5 = surv_obj ~ til_percentage + age_group + stage,
  MS6 = surv_obj ~ til_percentage + age_group + PAM50_Subtype
)

# 1. Ejecutar los modelos
resultados_cox_til <- lapply(modelos_til , function(f) coxph(f, data = clinical))
resultados_cox_til 
# install.packages("broom")
library(broom)
library(dplyr)

# 2. Evaluación del Supuesto de Proporcionalidad (Schoenfeld)
# Guardamos los tests en una lista para revisarlos uno a uno
schoenfeld_tests_til <- lapply(resultados_cox_til, function(m) {
  test <- cox.zph(m)
  return(test)
})

# Imprimir los p-valores globales de cada modelo
cat("--- P-valores Globales del Test de Schoenfeld ---\n")
lapply(names(schoenfeld_tests_til), function(n) {
  cat(n, ": ", schoenfeld_tests_til[[n]]$table["GLOBAL", "p"], "\n")
})

# Imprimir la tabla detallada de cada modelo
for (nombre in names(schoenfeld_tests_til)) {
  cat("\n==========================================\n")
  cat("DETALLE TEST DE SCHOENFELD:", nombre, "\n")
  cat("==========================================\n")
  print(schoenfeld_tests_til[[nombre]])
}

# El modelo MS1 y el MS4 violan el supuesto de proporcionalidad para el modelo global pero no para sus variables
# Mientras que el modelo MS2 MS3 lo violan para til, 
# Graficar residuos para MS1
test_m1 <- cox.zph(resultados_cox_til[["MS1"]])
plot(test_m1[1]) # Residuos para til_percentage
abline(h = 0, col = "red")

# Graficar residuos para MS4
test_m5 <- cox.zph(resultados_cox_til[["MS4"]])
plot(test_m5[1]) # Residuos para til_percentage
abline(h = 0, col = "red")



# 3. Evaluación de Multicolinealidad (VIF)
library(car)
vif_resultados_til <- lapply(resultados_cox_til, function(m) {
  # Intentamos calcular VIF, capturando el error si la estratificación da problemas
  tryCatch(vif(m), error = function(e) return("Error en cálculo de VIF o Inf detectado"))
})

# Imprimir VIFs (enfocarse en GVIF^(1/(2*Df)))
cat("\n--- Resultados de VIF por Modelo ---\n")
print(vif_resultados_til)


# 4. Grafico de forest plot
library(forestmodel)
forest_model(resultados_cox_til[["MS1"]])
forest_model(resultados_cox_til[["MS2"]])
forest_model(resultados_cox_til[["MS3"]])
forest_model(resultados_cox_til[["MS4"]])
forest_model(resultados_cox_til[["MS5"]])
forest_model(resultados_cox_til[["MS6"]]) #este el el unico que no viola ningun supuesto

#Solo da significativo el modelo 3 (PAM50 age). Lo guardo para figura supplementaria
ggsave(pathLocalResults("HE clas Results/FigureS2_forest_til.pdf"), forest_model(resultados_cox_til[["MS3"]]),  width = 8, height = 6, device = cairo_pdf)

       
# 1. Unificar todos los modelos en un solo dataframe
library(dplyr)
library(ggplot2)
library(broom)

# 1. Extraer resultados y añadir las FILAS DE REFERENCIA manualmente
extraer_con_ref <- function(modelo, nombre) {
  res <- tidy(modelo, exponentiate = TRUE, conf.int = TRUE) %>%
    mutate(Modelo = nombre)
  
  # Creamos un dataframe con las referencias (HR = 1, IC = 1)
  referencias <- data.frame(
    term = c("Tumor Size T1", "Node Status 0", "Stage Stage I", "PAM50 LumA", "Age < 70"),
    estimate = 1, conf.low = 1, conf.high = 1, p.value = NA,
    Modelo = nombre
  )
  
  # Unimos y limpiamos nombres (ajustar según tus nombres de variables reales)
  res <- res %>%
    mutate(term = case_when(
      term == "til_percentage" ~ "TIL percentage",
      grepl("tumor_size", term) ~ gsub("tumor_size", "Tumor Size ", term),
      grepl("node_status", term) ~ gsub("node_status", "Node Status ", term),
      grepl("stage", term) ~ gsub("stage", "Stage ", term),
      grepl("PAM50", term) ~ gsub("PAM50_Subtype", "PAM50 ", term),
      grepl("age", term) ~ gsub("age_group", "Age ", term),
      TRUE ~ term
    )) %>%
    bind_rows(referencias) %>%
    # Solo nos quedamos con las variables que existen en ESTE modelo específico
    filter(term %in% c("TIL percentage",
                       "Tumor Size T1", "Tumor Size T2", "Tumor Size T3", "Tumor Size T4",
                       "Node Status 0", "Node Status 1 to 3", "Node Status more than 4",
                       "Stage Stage I", "Stage Stage II", "Stage Stage III", "Stage Stage IV",
                       "PAM50 LumA", "PAM50 LumB",
                       "Age < 70", "Age > 70"))
  return(res)
}

# 2. Unificar todos los modelos
df_grafico <- bind_rows(lapply(names(resultados_cox_til), function(x) extraer_con_ref(resultados_cox_til[[x]], x)))

# 3. DEFINIR EL ORDEN DEL EJE Y (De abajo hacia arriba para que Immune Group quede arriba)
orden_niveles <- rev(c(
  "TIL percentage",
  "Age < 70", "Age > 70",
  "Tumor Size T1", "Tumor Size T2", "Tumor Size T3", "Tumor Size T4",
  "Node Status 0", "Node Status 1 to 3", "Node Status more than 4",
  "Stage Stage I", "Stage Stage II", "Stage Stage III", "Stage Stage IV",
  "PAM50 LumA", "PAM50 LumB"
))

df_grafico$term <- factor(df_grafico$term, levels = orden_niveles)

# 4. Paleta de colores personalizada
mis_colores <- c(
  "MS1" = "#80b1d3",
  "MS2" = "#b3de69", 
  "MS3" = "#bc80bd", 
  "MS4" = "#66c2a5",
  "MS5" = "#fdb462",
  "MS6" = "#fdb462"
)

# 5. Generación del gráfico
grafico_til <- ggplot(df_grafico, aes(x = estimate, y = term)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey40") +
  #barras de error (CI)
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), 
                 height = 0.2, 
                 size = 0.3, # ACHIQUE EL TAMAÑO
                 show.legend = FALSE, 
                 color = "black") +
  geom_point(aes(color = Modelo), size = 3, shape = 15, show.legend = FALSE) +
  # Facetado por modelo para replicar la imagen b86dcd
  facet_grid(. ~ Modelo) +
  # Escala logaritmica constante para todos los modelos (0.8 a 25)
  scale_x_log10(limits = c(0.8, 25), breaks = c(1, 2, 5, 10, 20)) +
  # Aplicacion de LOS colores especificos
  scale_color_manual(values = mis_colores) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "grey95"),
    strip.text = element_text(face = "bold", size = 10),
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 9, color = "black"),
    panel.grid.minor = element_blank(),
    legend.position = "none" # Los colores ya identifican los paneles
  ) +
  labs(x = "Hazard Ratio (95% CI)")

# Mostrar y guardar para Illustrator
print(grafico_til)

# 4. Guardar para Illustrator
ggsave(pathLocalResults("HE clas Results/Figure_suppl_forest_comparativo_til_final.pdf"), grafico_til, width = 18, height = 9, units = "cm", device = "pdf")


## Obtener tabla para comparar modelos

get_cox_summary <- function(model, model_name = "Model") {
  coefs <- summary(model)$coefficients
  confs <- summary(model)$conf.int
  
  # Create clean formatted data frame
  df <- data.frame(
    Variable = rownames(coefs),
    HR_CI = sprintf("%.4f (%.4f-%.4f)", confs[, "exp(coef)"], confs[, "lower .95"], confs[, "upper .95"]),
    P_value = formatC(coefs[, "Pr(>|z|)"], format = "e", digits = 3),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
  
  colnames(df)[2:3] <- paste0(c("HR_CI_", "P_"), model_name)
  return(df)
}


Model_1_df <- get_cox_summary(resultados_cox_til[["MS1"]], "Model 1")
Model_2_df <- get_cox_summary(resultados_cox_til[["MS2"]], "Model 2")
Model_3_df <- get_cox_summary(resultados_cox_til[["MS3"]], "Model 3")
Model_4_df <- get_cox_summary(resultados_cox_til[["MS4"]], "Model 4")
Model_5_df <- get_cox_summary(resultados_cox_til[["MS5"]], "Model 5")
Model_6_df <- get_cox_summary(resultados_cox_til[["MS6"]], "Model 6")


# Join by Variable name (use full_join to preserve all rows)
library(dplyr)

summary_combined <- Model_1_df %>%
  full_join(Model_2_df, by = "Variable") %>%
  full_join(Model_3_df, by = "Variable") %>%
  full_join(Model_4_df, by = "Variable")  %>%
  full_join(Model_5_df, by = "Variable")



#### Multivariable Cox for HE clas ====

# Ver cuántos NAs hay por cada columna
colSums(is.na(clinical[, c("HE_clas", "tumor_size", "node_status", "age_group", "PAM50_Subtype", "stage")]))


modelos <- list(
  M1 = surv_obj ~ HE_clas + age_group + tumor_size + node_status,
  M2 = surv_obj ~ HE_clas + age_group + stage,
  M3 = surv_obj ~ HE_clas + age_group + PAM50_Subtype,
  M4 = surv_obj ~ HE_clas + age_group + PAM50_Subtype + stage,
  M5 = surv_obj ~ HE_clas + age_group + PAM50_Subtype + tumor_size + node_status
)

# 1. Ejecutar los modelos
resultados_cox <- lapply(modelos, function(f) coxph(f, data = clinical))

# install.packages("broom")
library(broom)
library(dplyr)

# 2. Evaluación del Supuesto de Proporcionalidad (Schoenfeld)
# Guardamos los tests en una lista para revisarlos uno a uno
schoenfeld_tests <- lapply(resultados_cox, function(m) {
  test <- cox.zph(m)
  return(test)
})

# Imprimir los p-valores globales de cada modelo
cat("--- P-valores Globales del Test de Schoenfeld ---\n")
lapply(names(schoenfeld_tests), function(n) {
  cat(n, ": ", schoenfeld_tests[[n]]$table["GLOBAL", "p"], "\n")
})

# Imprimir la tabla detallada de cada modelo
for (nombre in names(schoenfeld_tests)) {
  cat("\n==========================================\n")
  cat("DETALLE TEST DE SCHOENFELD:", nombre, "\n")
  cat("==========================================\n")
  print(schoenfeld_tests[[nombre]])
}

#Ninguno viola el supuesto ni global ni para cada variable!!!

# 3. Evaluación de Multicolinealidad (VIF)
library(car)
vif_resultados <- lapply(resultados_cox, function(m) {
  # Intentamos calcular VIF, capturando el error si la estratificación da problemas
  tryCatch(vif(m), error = function(e) return("Error en cálculo de VIF o Inf detectado"))
})

# Imprimir VIFs (enfocarse en GVIF^(1/(2*Df)))
cat("\n--- Resultados de VIF por Modelo ---\n")
print(vif_resultados)


# 4. Grafico de forest plot
library(forestmodel)
forest_model(resultados_cox[["M1"]])
forest_model(resultados_cox[["M2"]])
forest_model(resultados_cox[["M3"]])
forest_model(resultados_cox[["M4"]])
forest_model(resultados_cox[["M5"]])

# 1. Unificar todos los modelos en un solo dataframe
library(dplyr)
library(ggplot2)
library(broom)

# 1. Extraer resultados y añadir las FILAS DE REFERENCIA manualmente
extraer_con_ref <- function(modelo, nombre) {
  res <- tidy(modelo, exponentiate = TRUE, conf.int = TRUE) %>%
    mutate(Modelo = nombre)
  
  # Creamos un dataframe con las referencias (HR = 1, IC = 1)
  referencias <- data.frame(
    term = c("Immune Group High", "Tumor Size T1", "Node Status 0", "Stage Stage I", "PAM50 LumA", "Age < 70"),
    estimate = 1, conf.low = 1, conf.high = 1, p.value = NA,
    Modelo = nombre
  )
  
  # Unimos y limpiamos nombres (ajustar según tus nombres de variables reales)
  res <- res %>%
    mutate(term = case_when(
      term == "HE_claslow" ~ "Immune Group Low",
      grepl("tumor_size", term) ~ gsub("tumor_size", "Tumor Size ", term),
      grepl("node_status", term) ~ gsub("node_status", "Node Status ", term),
      grepl("stage", term) ~ gsub("stage", "Stage ", term),
      grepl("PAM50", term) ~ gsub("PAM50_Subtype", "PAM50 ", term),
      grepl("age", term) ~ gsub("age_group", "Age ", term),
      TRUE ~ term
    )) %>%
    bind_rows(referencias) %>%
    # Solo nos quedamos con las variables que existen en ESTE modelo específico
    filter(term %in% c("Immune Group High", "Immune Group Low", 
                       "Tumor Size T1", "Tumor Size T2", "Tumor Size T3", "Tumor Size T4",
                       "Node Status 0", "Node Status 1 to 3", "Node Status more than 4",
                       "Stage Stage I", "Stage Stage II", "Stage Stage III", "Stage Stage IV",
                       "PAM50 LumA", "PAM50 LumB",
                       "Age < 70", "Age > 70"))
  return(res)
}

# 2. Unificar todos los modelos
df_grafico <- bind_rows(lapply(names(resultados_cox), function(x) extraer_con_ref(resultados_cox[[x]], x)))

# 3. DEFINIR EL ORDEN DEL EJE Y (De abajo hacia arriba para que Immune Group quede arriba)
orden_niveles <- rev(c(
  "Immune Group High", "Immune Group Low",
  "Age < 70", "Age > 70",
  "Tumor Size T1", "Tumor Size T2", "Tumor Size T3", "Tumor Size T4",
  "Node Status 0", "Node Status 1 to 3", "Node Status more than 4",
  "Stage Stage I", "Stage Stage II", "Stage Stage III", "Stage Stage IV",
  "PAM50 LumA", "PAM50 LumB"
))

df_grafico$term <- factor(df_grafico$term, levels = orden_niveles)

# 4. Paleta de colores personalizada
mis_colores <- c(
  "M1" = "#80b1d3",
  "M2" = "#b3de69", 
  "M3" = "#bc80bd", 
  "M4" = "#66c2a5",
  "M5" = "#fdb462"
)

# 5. Generación del gráfico
grafico <- ggplot(df_grafico, aes(x = estimate, y = term)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey40") +
  #barras de error (CI)
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), 
                 height = 0.2, 
                 size = 0.3, # ACHIQUE EL TAMAnO
                 show.legend = FALSE, 
                 color = "black") +
  geom_point(aes(color = Modelo), size = 3, shape = 15, show.legend = FALSE) +
  # Facetado por modelo para replicar la imagen b86dcd
  facet_grid(. ~ Modelo) +
  # Escala logaritmica constante para todos los modelos (0.8 a 25)
  scale_x_log10(limits = c(0.6, 25), breaks = c(1, 2, 5, 10, 20)) +
  # Aplicación de LOS colores específicos
  scale_color_manual(values = mis_colores) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "grey95"),
    strip.text = element_text(face = "bold", size = 10),
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 9, color = "black"),
    panel.grid.minor = element_blank(),
    legend.position = "none" # Los colores ya identifican los paneles
  ) +
  labs(x = "Hazard Ratio (95% CI)")

# Mostrar y guardar para Illustrator
print(grafico)

# 4. Guardar para Illustrator
ggsave(pathLocalResults("HE clas Results/Figure_1c_forest_comparativo_final.pdf"), grafico, width = 18, height = 9, units = "cm", device = "pdf")


## Obtener tabla para comparar modelos

get_cox_summary <- function(model, model_name = "Model") {
  coefs <- summary(model)$coefficients
  confs <- summary(model)$conf.int
  
  # Create clean formatted data frame
  df <- data.frame(
    Variable = rownames(coefs),
    HR_CI = sprintf("%.4f (%.4f-%.4f)", confs[, "exp(coef)"], confs[, "lower .95"], confs[, "upper .95"]),
    P_value = formatC(coefs[, "Pr(>|z|)"], format = "e", digits = 3),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
  
  colnames(df)[2:3] <- paste0(c("HR_CI_", "P_"), model_name)
  return(df)
}


Model_1_df <- get_cox_summary(resultados_cox[["M1"]], "Model 1")
Model_2_df <- get_cox_summary(resultados_cox[["M2"]], "Model 2")
Model_3_df <- get_cox_summary(resultados_cox[["M3"]], "Model 3")
Model_4_df <- get_cox_summary(resultados_cox[["M4"]], "Model 4")
Model_5_df <- get_cox_summary(resultados_cox[["M5"]], "Model 5")


# Join by Variable name (use full_join to preserve all rows)
library(dplyr)

summary_combined <- Model_1_df %>%
  full_join(Model_2_df, by = "Variable") %>%
  full_join(Model_3_df, by = "Variable") %>%
  full_join(Model_4_df, by = "Variable")  %>%
  full_join(Model_5_df, by = "Variable")


write.csv(summary_combined, pathLocalResults("HE clas Results/multivariable_cox_models_age_table_without_prior_treatment.csv"), row.names = FALSE)
library(openxlsx)
write.xlsx(summary_combined, pathLocalResults("HE clas Results/multivariable_cox_models_age_table_excel_without_prior_treatment.xlsx"), rowNames = FALSE)


## Obtener p adj values para HE clas low-> esto permite confirmar que esta clasificacion es estable para todos los modelos
library(dplyr)
library(tidyr)

# 1. Identificar las columnas que contienen p-valores
columnas_p <- grep("^P_Model", names(summary_combined), value = TRUE)

# 2. Convertir a formato largo para procesar todos los p-valores juntos
# Esto facilita la aplicación de p.adjust a toda la bolsa de tests realizados
p_long <- summary_combined %>%
  select(Variable, all_of(columnas_p)) %>%
  pivot_longer(cols = all_of(columnas_p), 
               names_to = "Modelo", 
               values_to = "p_val_original") %>%
  # Quitamos los NA porque no se pueden ajustar
  filter(!is.na(p_val_original))

# 3. Calcular el P-Adjusted (FDR / Benjamini-Hochberg)
p_long <- p_long %>%
  mutate(p_val_adjusted = p.adjust(p_val_original, method = "fdr"))

# 4. Volver al formato ancho original para unirlo a tu tabla
p_adjusted_wide <- p_long %>%
  select(Variable, Modelo, p_val_adjusted) %>%
  mutate(Modelo = gsub("P_", "P_Adj_", Modelo)) %>% # Renombramos las columnas
  pivot_wider(names_from = Modelo, values_from = p_val_adjusted)

# 5. Unir con la tabla original
summary_final <- left_join(summary_combined, p_adjusted_wide, by = "Variable")

# Ver resultados de tu variable de interés
summary_final %>% 
  filter(Variable == "HE_claslow") %>% 
  select(Variable, starts_with("P_Model"), starts_with("P_Adj"))




#### Multivariable Cox for PAM50 ====

modelos <- list(
  M1 = surv_obj ~ PAM50_Subtype + age_group +  node_status,
  M2 = surv_obj ~ PAM50_Subtype + age_group + tumor_size,
  M3 = surv_obj ~ PAM50_Subtype + age_group + stage
)

# 1. Ejecutar los modelos
resultados_cox <- lapply(modelos, function(f) coxph(f, data = clinical))
resultados_cox
# install.packages("broom")
library(broom)
library(dplyr)

# 2. Evaluación del Supuesto de Proporcionalidad (Schoenfeld)
# Guardamos los tests en una lista para revisarlos uno a uno
schoenfeld_tests <- lapply(resultados_cox, function(m) {
  test <- cox.zph(m)
  return(test)
})

# Imprimir los p-valores globales de cada modelo
cat("--- P-valores Globales del Test de Schoenfeld ---\n")
lapply(names(schoenfeld_tests), function(n) {
  cat(n, ": ", schoenfeld_tests[[n]]$table["GLOBAL", "p"], "\n")
})

# Imprimir la tabla detallada de cada modelo
for (nombre in names(schoenfeld_tests)) {
  cat("\n==========================================\n")
  cat("DETALLE TEST DE SCHOENFELD:", nombre, "\n")
  cat("==========================================\n")
  print(schoenfeld_tests[[nombre]])
}


# 3. Evaluación de Multicolinealidad (VIF)
library(car)
vif_resultados <- lapply(resultados_cox, function(m) {
  # Intentamos calcular VIF, capturando el error si la estratificación da problemas
  tryCatch(vif(m), error = function(e) return("Error en cálculo de VIF o Inf detectado"))
})

# Imprimir VIFs (enfocarse en GVIF^(1/(2*Df)))
cat("\n--- Resultados de VIF por Modelo ---\n")
print(vif_resultados)


# 4. Grafico de forest plot
library(forestmodel)
pam50_stage <- forest_model(resultados_cox[["M3"]])
pam50_stage


# 4. Guardar para Illustrator
ggsave(pathLocalResults("HE clas Results/Figure_S3_forest_PAM50_stage.pdf"), pam50_stage , width = 18, height = 9, units = "cm", device = "pdf")


## Obtener tabla para comparar modelos

get_cox_summary <- function(model, model_name = "Model") {
  coefs <- summary(model)$coefficients
  confs <- summary(model)$conf.int
  
  # Create clean formatted data frame
  df <- data.frame(
    Variable = rownames(coefs),
    HR_CI = sprintf("%.4f (%.4f-%.4f)", confs[, "exp(coef)"], confs[, "lower .95"], confs[, "upper .95"]),
    P_value = formatC(coefs[, "Pr(>|z|)"], format = "e", digits = 3),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
  
  colnames(df)[2:3] <- paste0(c("HR_CI_", "P_"), model_name)
  return(df)
}


Model_1_df <- get_cox_summary(resultados_cox[["M1"]], "Model 1")
Model_2_df <- get_cox_summary(resultados_cox[["M2"]], "Model 2")
Model_3_df <- get_cox_summary(resultados_cox[["M3"]], "Model 3")
#Model_4_df <- get_cox_summary(resultados_cox[["M4"]], "Model 4")
#Model_5_df <- get_cox_summary(resultados_cox[["M5"]], "Model 5")


# Join by Variable name (use full_join to preserve all rows)
library(dplyr)

summary_combined <- Model_1_df %>%
  full_join(Model_2_df, by = "Variable") %>%
  full_join(Model_3_df, by = "Variable") 
summary_combined

write.csv(summary_combined, pathLocalResults("HE clas Results/multivariable_cox_models_PAM50_table_without_prior_treatment.csv"), row.names = FALSE)
library(openxlsx)
write.xlsx(summary_combined, pathLocalResults("HE clas Results/multivariable_cox_models_PAM50_table_excel_without_prior_treatment.xlsx"), rowNames = FALSE)









################################################################################
### TIL percentage between lumA and LumB ====

# Perform Normality test for total_perMB 
perform_normality_tests(clinical, "til_percentage", "PAM50_Subtype")

## Wilcoxon test compare TMB between Treg groups
wilcox_tmb_result <- wilcoxon_test_function(clinical, "til_percentage", "PAM50_Subtype", "LumA")
wilcox_tmb_result 

# Plot results TMB between high and low immune groups by violin plot
annotations <- list(
  list(x = 1.5, y = 2, stars = 1, variable = colnames(clinical)[122]) #122 es til percentage
)
create_grouped_boxplot(data = clinical, #a data frame with in columns should be "group_variable" and each "variable" to be plotted, in rows: samples
                       group = "PAM50_Subtype",
                       column_names = "til_percentage",
                       y_limit = 2.5,
                       annotations = annotations,
                       y_text = "TIL percentage")

annotations <- list(
  list(x = 1.5, y = 8, stars = 1, variable = colnames(clinical)[122]) #122 es til percentage
)
create_grouped_violin_plot(data = clinical, #a data frame with in columns should be "group_variable" and each "variable" to be plotted, in rows: samples
                       group = "PAM50_Subtype",
                       column_names = "til_percentage",
                       y_limit =8,
                       annotations = annotations,
                       y_text = "TIL percentage")

################################################################################
### MUTATION ANALYSIS ====
maf<- read.table(pathLocalDb("maf_BRCA.txt"))
MAF<- read.maf(maf)

#### TMB ====
tmb_patients <- tmb(MAF)

#library(stringr) # split the string and extract the desired part in column named

# Create the new "patient" column
tmb_patients$patient <- str_extract(tmb_patients$Tumor_Sample_Barcode, "^[^-]+-[^-]+-[^-]+")

# Display the head of the data frame to verify
head(tmb_patients)
summary(tmb_patients)

#Saco los hipermutadores (punto de corte tmb = 10)  -> No deberia eliminarlos en este caso

#tmb_patients <- tmb_patients[tmb_patients$total_perMB<11,]

tmb_patients <- tmb_patients[,c("total_perMB", "total_perMB_log", "patient", "Tumor_Sample_Barcode")]

#Agrego datos de TMB al clinical
clinical_tmb <- merge(clinical, tmb_patients, by = "patient")

# Perform Normality test for total_perMB 
clinical_tmb <- clinical_tmb[clinical_tmb$total_perMB!=0,]
perform_normality_tests(clinical_tmb, "total_perMB", "HE_clas")

## Wilcoxon test compare TMB between Treg groups
wilcox_tmb_result <- wilcoxon_test_function(clinical_tmb, "total_perMB", "HE_clas", "low")
wilcox_tmb_result 

# Wilcoxon rank sum test with continuity correction
# 
# data:  condition_true and condition_false
# W = 41027, p-value = 0.03353
# alternative hypothesis: true location shift is not equal to 0

# Plot results TMB between high and low immune groups by violin plot
annotations <- list(
  list(x = 1.5, y = 2, stars = 1, variable = colnames(clinical_tmb)[129]) #ojo con HER ER y PR status cambio el numero de columna
)
create_grouped_boxplot(data = clinical_tmb, #a data frame with in columns should be "group_variable" and each "variable" to be plotted, in rows: samples
                       group = "HE_clas",
                       column_names = "total_perMB",
                       y_limit = 2.5,
                       annotations = annotations,
                       y_text = "TMB per Mb")
library(ggplot2)

p <- ggplot(clinical_tmb, aes(x = HE_clas, y = total_perMB, fill = HE_clas)) +
  geom_boxplot(
    width = 0.6,
    outlier.shape = NA,   # ocultar outliers (ya están como puntos)
    alpha = 0.6,
    linewidth = 0.4
  ) +
  geom_jitter(
    width = 0.15,
    size = 1.2,
    alpha = 0.6
  ) +
  scale_fill_brewer(palette = "Set2") +
  labs(
    x = NULL,
    y = "Tumor Mutational Burden (mut/Mb)"
  ) +
  theme_bw(base_size = 8) +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 7),
    axis.title = element_text(size = 8),
    plot.title = element_text(size = 9, hjust = 0.5)
  )

p


## Univariable Cox for TMB 
censored_time <- pmin(clinical_tmb$overall_survival, 3650)
event_indicator <- clinical_tmb$vital_status_binary & (clinical_tmb$overall_survival <= 3650)

#creamos un objeto "surv" con la funci?n Surv()
surv_obj <- Surv(censored_time,
                 event_indicator)


formula <- as.formula(paste("surv_obj ~", "total_perMB"))
univ_tmb_cox <- coxph(formula, data = clinical_tmb)
univ_tmb_cox # No da significativo


#####################################################################
########## mutation frecuency analysis ###########
MAF_patients <- data.frame(
  Tumor_Sample_Barcode = getSampleSummary(MAF)$Tumor_Sample_Barcode,
  patient = str_extract(getSampleSummary(MAF)$Tumor_Sample_Barcode, "^[^-]+-[^-]+-[^-]+")
)

# Crear un objeto MAF que reconozca a todos los pacientes como datos clínicos
# incluso si no tienen mutaciones en la tabla principal.
filtered_MAF <- read.maf(
  maf = maf[maf$Tumor_Sample_Barcode %in% MAF_patients$Tumor_Sample_Barcode, ],
  clinicalData = MAF_patients # Esto "ancla" a todos los pacientes al objeto
)

# Ahora revisa el resumen de muestras
getSampleSummary(filtered_MAF)

# Debido a que hay pacientes sin mutaciones, entonces el MAF se va filtrando y se disminuye el n quedando como NA cuando deberia ser wt

# Identificar pacientes faltante
filtered_MAF_data <- as.data.frame(filtered_MAF@data)
pacientes_con_mutacion <- unique(filtered_MAF@data$Tumor_Sample_Barcode)
paciente_faltante <- setdiff(MAF_patients$Tumor_Sample_Barcode, pacientes_con_mutacion)
print(paciente_faltante)

# A. Obtener los top 20 genes
top20_genes <- getGeneSummary(filtered_MAF)$Hugo_Symbol[1:20]

# B. Extraer solo las mutaciones de esos genes
# Usamos el slot @data pero filtramos por nuestros genes de interés
filtered_MAF_data <- as.data.frame(filtered_MAF@data) %>% 
  filter(Hugo_Symbol %in% top20_genes)

# C. Crear una matriz base con TODOS los pacientes (986) y los 20 genes
# Inicializamos todo en 0 (Wild Type)
mut_status_matrix <- matrix(0, 
                            nrow = nrow(MAF_patients), 
                            ncol = length(top20_genes))

rownames(mut_status_matrix) <- MAF_patients$Tumor_Sample_Barcode
colnames(mut_status_matrix) <- top20_genes

# D. Llenar con 1 los casos donde sí hay mutación
# Recorremos el dataframe de mutaciones y marcamos la posición correspondiente
for(i in 1:nrow(filtered_MAF_data)) {
  sample <- as.character(filtered_MAF_data$Tumor_Sample_Barcode[i])
  gene <- as.character(filtered_MAF_data$Hugo_Symbol[i])
  
  if(sample %in% rownames(mut_status_matrix) & gene %in% colnames(mut_status_matrix)) {
    mut_status_matrix[sample, gene] <- 1
  }
}

# E. Convertir a DataFrame final
mutation_status<- as.data.frame(mut_status_matrix)
mutation_status$patient <- MAF_patients$patient[match(rownames(mutation_status), MAF_patients$Tumor_Sample_Barcode)]

# Merge with clinical 
clinical <- clinical %>% 
  left_join(mutation_status, by = "patient")

#convierto columns de genes de numerico a factor reemplazando 0 por wt y 1 por mutate
cols <- colnames(clinical[,129:148])

clinical[cols] <- lapply(clinical[cols], function(x) {
  factor(
    ifelse(x == 0, "wt",
           ifelse(x == 1, "mutated", NA)
    ),
    levels = c("wt", "mutated")
  )
})

## Análsisis estadisticos

#### Til_percentage between mutated and wt ====
# Perform Normality test 
perform_normality_tests <- function(df, outcome, col) {
  
  sub <- df[, c(outcome, col)]
  sub <- sub[!is.na(sub[[col]]) & !is.na(sub[[outcome]]), ]
  
  wt_vals  <- sub[sub[[col]] == "wt", outcome]
  mut_vals <- sub[sub[[col]] == "mutated", outcome]
  
  res <- list(
    gene = col,
    n_wt = length(wt_vals),
    n_mut = length(mut_vals),
    shapiro_wt = if (length(wt_vals) >= 3 && length(wt_vals) <= 5000)
      shapiro.test(wt_vals)$p.value else NA,
    shapiro_mut = if (length(mut_vals) >= 3 && length(mut_vals) <= 5000)
      shapiro.test(mut_vals)$p.value else NA
  )
  
  return(as.data.frame(res))
}

results <- do.call(rbind, lapply(top20_genes, function(i) {
  perform_normality_tests(clinical, "til_percentage", i)
}))
results 

## Wilcoxon test compare til_percentaje between mutated vs wt for each gene
### COPIAR FUNCION AL PRINCIPIO SI ES QUE FUNCIONA OK
wilcoxon_test_function <- function(data, column, condition_column, condition_value) {
  
  # Extract values
  condition_true <- data[data[[condition_column]] == condition_value, column]
  condition_false <- data[data[[condition_column]] != condition_value, column]
  
  # Remove NA
  condition_true <- condition_true[!is.na(condition_true)]
  condition_false <- condition_false[!is.na(condition_false)]
  
  # Sample sizes
  n1 <- length(condition_true)
  n2 <- length(condition_false)
  
  # Wilcoxon test
  wilcox_result <- wilcox.test(condition_true, condition_false)
  
  # U statistic (equivalente al W transformado)
  U <- wilcox_result$statistic
  
  # Rank-biserial correlation
  r_rb <- 1 - (2 * U) / (n1 * n2)
  
  return(list(
    test = wilcox_result,
    effect_size = r_rb,
    n1 = n1,
    n2 = n2
  ))
}


###########
all_results <- list()
for (col in top20_genes) { #top20_genes son los genes mas mutados en toda la cohorte
  # Perform Wilcoxon test
  wilcox_result <- wilcoxon_test_function(clinical, "til_percentage", col, "wt")
  
  # Calculate median TILs for wildtype and mutated groups
  wt_median <- median(clinical$til_percentage[clinical[[col]] == "wt"], na.rm = TRUE)
  mut_median <- median(clinical$til_percentage[clinical[[col]] != "wt"], na.rm = TRUE)
  
  # Store results in list
  all_results[[col]] <- list(
    p.value = wilcox_result$test$p.value,
    statistic = wilcox_result$test$statistic,
    effect_size = wilcox_result$effect_size,
    wt_median_til = wt_median,
    mut_median_til = mut_median
  )
  
  
  # Print results
  cat("Wilcoxon test results for", col, ":\n")
  print(paste0("p-value: ", wilcox_result$p.value))
}

# Extract components into vectors
genes <- names(all_results)
p_values <- sapply(all_results, function(x) x$p.value)
test_statistics <- sapply(all_results, function(x) x$statistic)
wt_til <- sapply(all_results, function(x) x$wt_median_til)
mut_til <- sapply(all_results, function(x) x$mut_median_til)
effect_sizes <- sapply(all_results, function(x) x$effect_size)

results_table <- data.frame(
  Gene = genes,
  P_Value = p_values,
  Test_Statistic = test_statistics,
  Effect_Size = effect_sizes,
  TIL_Median_WT = wt_til,
  TIL_Median_Mut = mut_til
)

# Add significance stars
results_table$Significance <- ifelse(results_table$P_Value < 0.05, "*", "")


result_table <- results_table[order(results_table$P_Value),]
result_table
dif_genes <- rownames(subset(result_table, P_Value <= 0.05))


#Save results
write.csv2(results_table, file=pathLocalResults("HE clas Results/Wilcoxon test of til percentages between wt and mutated genes without prior treatment.csv"), row.names = FALSE)

#Plot
library(ggplot2)
library(dplyr)

# Create a long-format data frame for plotting
plot_data <- clinical %>%
  filter(!is.na(til_percentage)) %>%
  select(patient, til_percentage, all_of(dif_genes)) %>%
  pivot_longer(
    cols = all_of(dif_genes),
    names_to = "Gene",
    values_to = "Mutation_Status"
  ) %>%
  filter(Mutation_Status %in% c("wt", "mutated"))  # adjust if labels are different

annotations <- data.frame(
  Gene = c("TP53", "HMCN1", "RYR2", "MAP3K1", "NCOR1"),
  label = c("*"),
  y = c(15, 15, 15, 15, 15),  # Y position of asterisk
  x = 1.5          # Center between 1 (wt) and 2 (mut)
)


til_mutated <- ggplot(plot_data, aes(x = Mutation_Status, y = til_percentage, fill = Mutation_Status)) +
  geom_boxplot(outlier.shape = NA, width = 0.7, linewidth = 0.2) +
  facet_wrap(~ Gene, nrow = 1) +
  scale_y_continuous(limits = c(0, 15)) +
  labs(
    x = "Mutation Status",
    y = "TIL Percentage",
    fill = "Mutation Status"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(face = "bold", size = 8),
    axis.text.x = element_text(size = 7), 
    axis.text.y = element_text(size = 7),
    axis.title.y = element_text(size = 8),
    # RELACIÓN DE ASPECTO: Controla que el carril no sea un tubo vertical
    aspect.ratio = 5, 
    legend.position = "bottom",
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 7),
    legend.key.size = unit(0.3, "cm"),
    panel.spacing = unit(0.8, "lines"), # Espacio entre genes
    panel.grid.major = element_line(linewidth = 0.1),
    panel.grid.minor = element_blank()
  ) +
  # ESTA LÍNEA aumenta el espacio a los lados de los boxplots dentro del carril
  scale_x_discrete(expand = expansion(mult = c(0.4, 0.4))) +
  scale_fill_manual(values = c("wt" = "#1f77b4", "mutated" = "#d62728")) +
  coord_cartesian(ylim = c(0, 15)) +
  geom_text(data = annotations, aes(x = x, y = y, label = label),
            inherit.aes = FALSE, size = 3)

# Mostrar y guardar para Illustrator
print(til_mutated)


# 4. Guardar para Illustrator
ggsave(pathLocalResults("HE clas Results/Figure_2b_til_beween_mutated_wt-without_treatment.pdf"), til_mutated, width = 9, height = 9, units = "cm", device = "pdf")


#### Mutation frequency between HE_clas by Fisher ====
# Chi cuadrado o fisher para HE_clas y genes mutados o wt

mutated_symbols <- top20_genes #lo habia hecho con dif genes, are those which TIL percentage changed between mutated and wt

## Loop de chi cuadrado para todos los genes mutados high /low vs mutated non mutated gene

groups <- c("high", "low")

results <- data.frame(
  Gene = character(),
  Score_High = numeric(),
  Score_Low = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

for (gene in mutated_symbols) {
  tbl <- clinical %>%
    select(HE_clas, !!sym(gene)) %>%
    mutate(mutated = ifelse(.data[[gene]] != "wt", 1, 0)) %>%
    group_by(HE_clas) %>%
    summarise(
      Mutated = sum(mutated, na.rm = TRUE),
      Total = n(),
      Score = Mutated / Total
    )
  
  
  # ensure both groups are present
  if (all(groups %in% tbl$HE_clas)) {
    score_high <- tbl$Score[tbl$HE_clas == "high"]
    score_low  <- tbl$Score[tbl$HE_clas == "low"]
    
    # test for association
    contingency <- table(clinical$HE_clas, clinical[[gene]])
    pval <- fisher.test(contingency)$p.value
    
    results <- rbind(results, data.frame(
      Gene = gene,
      Score_High = score_high,
      Score_Low = score_low,
      p_value = pval
    ))
  }
}

results <- results[order(results$p_value), ]
results_sign <- subset(results, p_value <= 0.05)

# Plot

library(ggplot2)
library(ggrepel)

mutation_frequency_plot <- ggplot(
  results, aes(x = Score_Low * 100, y = Score_High * 100, label = Gene)) +
  geom_abline(slope = 1, intercept = 0, color = "black") +
  geom_point(aes(color = p_value <= 0.05), size = 1.5) +
  # Use ggrepel with segments/arrows to connect labels to points
  geom_text_repel(
    #data = subset(results, p_value <= 0.05),
    aes(label = Gene),
    size = 2,
    segment.color = "grey50",      # Color of connecting lines
    segment.size = 0.2,            # Thickness of lines
    max.overlaps = 30,            # Increase to show more labels
    force = 3                    # Fuerza el alejamiento si hay mucho solapamiento
  ) +
  scale_color_manual(
    values = c("TRUE" = "#d62728", "FALSE" = "#1f77b4"),
    labels = c("p > 0.05","p ≤ 0.05"),
    name = "Significance"
  ) +
  scale_x_continuous(
    limits = c(0, 40),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    limits = c(0, 45),
    expand = c(0, 0)
  ) +
  labs(
    x = "Mutation Frequency (%) in Low immune group",
    y = "Mutation Frequency (%) in High immune group"
  ) +
  theme_minimal()+
  theme(
    axis.text.x = element_text(size = 7), 
    axis.text.y = element_text(size = 7),
    axis.title.y = element_text(size = 7),
    axis.title.x = element_text(size = 7),
    legend.position = "bottom",
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 7),
    legend.key.size = unit(0.3, "cm"),
  )

print(mutation_frequency_plot)
# 4. Guardar para Illustrator
ggsave(pathLocalResults("HE clas Results/Figure_2a_mutation_frequency_without_treatment.pdf"), mutation_frequency_plot, width = 9, height = 9, units = "cm", device = "pdf")


write.csv2(results, file =pathLocalResults("Chi-squared_results between mutated genes and immune group.csv"), row.names = FALSE)


# ### Otros analisis que pueden realizarse ###
# 
# 
# #### Exclusive/co-occurance event analysis on top 30 mutated genes ====
# all_interaction <- somaticInteractions(maf = filtered_MAF, top = 20, pvalue = c(0.05), countStats ="sig")
# sig_all_interaction <- all_interaction[all_interaction$pValue<=0.05,]
# 
# high_interaction <- somaticInteractions(maf = MAF_high, top = 20, pvalue = c(0.05), countStats ="sig")
# sig_high_interaction <- high_interaction[high_interaction$pValue<=0.05,]
# 
# low_interaction <- somaticInteractions(maf = MAF_low, top = 20, pvalue = c(0.05))
# sig_low_interaction <- low_interaction[low_interaction$pValue<=0.05,]
# 
# ## Analizar mas detalles en: https://academic.oup.com/bioinformatics/article/29/18/2238/240376
# 
# 
# #### Asociaciones entre mutaciones y características clínicas ====
# # DE ACUERDO AL PAQUETE MAFTOOLS https://bioconductor.org/packages/devel/bioc/vignettes/maftools/inst/doc/maftools.html#81_Reading_and_summarizing_gistic_output_files
# 
# 
# #### Plot genes cuya frecuencia de mutaciones varíe entre High vs low ====
# genes <- results_sign$Gene
# coOncoplot(m1 = MAF_high, m2 = MAF_low, m1Name = 'High', m2Name = 'Low', genes = genes, removeNonMutated = TRUE)
# coBarplot(m1 = MAF_high, m2 = MAF_low, m1Name = 'High', m2Name = 'Low', genes = genes)


####  Cox for each mutated gene ====
### Survival Object 

censored_time <- pmin(clinical$overall_survival, 3650)
event_indicator <- clinical$vital_status_binary & (clinical$overall_survival <= 3650)

#creamos un objeto "surv" con la funcion Surv()
surv_obj <- Surv(clinical$censored_time,
                 clinical$event_indicator)

#Obtain pvalue and coefficient for each cell 

mutated_symbols <- c("PIK3CA","CDH1","TTN","GATA3","TP53","MUC16","KMT2C", "MAP3K1","MAP2K4","MUC4","RUNX1","MUC5B","NEB","NCOR1","ARID1A","FLG","RYR2","HMCN1","USH2A", "ERBB2")

mutated_symbols <- top20_genes

univ_gene_results <- list()
for (col in mutated_symbols) {
  formula <- as.formula(paste("surv_obj ~", col))
  univ_cox <- coxph(formula, data = clinical)
  summary_univ_cox <- summary(univ_cox)
  p_value <- summary_univ_cox[["logtest"]][["pvalue"]]
  CI <- summary_univ_cox$conf.int[, c("lower .95", "upper .95")]
  HR <- summary_univ_cox$coefficients[1, "exp(coef)"]
  univ_gene_results[[col]] <- list(p_value = p_value, HR = HR, CI=CI)
}


# Convert the list to a data frame
cox_results_df <- do.call(rbind, lapply(names(univ_gene_results), function(col) {
  c(Gene = col, 
    P_Value = univ_gene_results[[col]]$p_value, 
    HR = univ_gene_results[[col]]$HR,
    CI_low = univ_gene_results[[col]][["CI"]][["lower .95"]],
    CI_upp = univ_gene_results[[col]][["CI"]][["upper .95"]])
}))

# Convert columns to appropriate types
cox_results_df <- data.frame(cox_results_df, stringsAsFactors = FALSE)
cox_results_df$P_Value <- as.numeric(cox_results_df$P_Value)
cox_results_df$HR <- as.numeric(cox_results_df$HR)
cox_results_df$CI_low <- as.numeric(cox_results_df$CI_low)
cox_results_df$CI_upp <- as.numeric(cox_results_df$CI_upp)

# Sort the data frame by P_Value
cox_results_df <- cox_results_df[order(cox_results_df$P_Value), ]

# Print the sorted results
print(cox_results_df)


# > print(cox_results_df)
#      Gene    P_Value        HR    CI_low     CI_upp
# 7   KMT2C 0.01388109 0.4546391 0.2559826  0.8074638
# 12  MUC5B 0.04871755 0.3480577 0.1403852  0.8629410
# 2    CDH1 0.06287353 1.9512986 0.9002415  4.2294942
# 13    NEB 0.08159456 3.9390048 0.5481765 28.3043109
# 5    TP53 0.17025070 0.6847209 0.4064577  1.1534847
# 18  HMCN1 0.22420205 2.7815412 0.3869061 19.9970284
# 15 ARID1A 0.30324294 1.9360828 0.4760052  7.8747395
# 17   RYR2 0.49883819 1.4548366 0.4594242  4.6069613
# 19  USH2A 0.53700221 0.7169448 0.2621887  1.9604573
# 1  PIK3CA 0.54587059 0.8715557 0.5594776  1.3577117
# 14  NCOR1 0.62892413 1.2697008 0.4650625  3.4665022
# 11  RUNX1 0.63224262 1.2669735 0.4639300  3.4600522
# 16    FLG 0.63469500 0.7971327 0.3224296  1.9707266
# 9  MAP2K4 0.65198488 1.2504545 0.4577803  3.4156914
# 6   MUC16 0.73902729 1.1229014 0.5619000  2.2440071
# 8  MAP3K1 0.77620274 0.8921538 0.4112192  1.9355572
# 10   MUC4 0.87186947 0.9080878 0.2862372  2.8809091
# 4   GATA3 0.96364580 1.0134399 0.5703065  1.8008920
# 3     TTN 0.97239438 0.9896092 0.5480886  1.7868031

library(coxphf)
clinical$censored_time <- pmin(clinical$overall_survival, 3650)
clinical$event_indicator <- clinical$vital_status_binary & (clinical$overall_survival <= 3650)


genes_to_test <- c("FBN3", "RYR1", "TP53")
penalized_results <- list()

for (gene in genes_to_test) {
  # Filtrar datos sin NA
  model_data <- clinical[, c("censored_time", "event_indicator", gene)]
  model_data <- na.omit(model_data)
  
  # Agregar columna 'surv_obj' dentro del dataframe para que coxphf la reconozca
  model_data$surv_obj <- with(model_data, Surv(censored_time, event_indicator))
  
  # Ajustar modelo penalizado
  formula <- as.formula(paste("surv_obj ~", gene))
  
  fit <- coxphf(formula, data = model_data)
  summary_fit <- summary(fit)
  
  penalized_results[[gene]] <- list(
    HR = summary_fit$coefficients["exp(coef)"],
    CI = c(summary_fit$conf.int["lower .95"], summary_fit$conf.int["upper .95"]),
    p_value = summary_fit$prob
  )
}


surv_obj <- Surv(clinical$censored_time,
                 clinical$event_indicator)

cox_fit <- coxph(Surv(censored_time, event_indicator) ~ TP53, data = clinical)
summary_cox <- summary(cox_fit)
summary_cox 

#Guardar usando write.csv2() (usa ; como separador de columnas y , como decimal)
#sto es ideal para abrir en Excel con configuración regional latinoamericana o europea:
write.csv2(cox_results_df, file =pathLocalResults("Cox results for each mutated genes.csv"), row.names = FALSE)


## Multivariado incorporando mutaciones
for (i in mutated_symbols) {
  clinical[[i]] <- relevel(
    as.factor(clinical[[i]]),
    ref = "wt"
  )
}

### Multivariable models incorporing TP53 as covariable ====
# Ver cuántos NAs hay por cada columna
colSums(is.na(clinical[, c("HE_clas", "tumor_size", "node_status", "age_group", "PAM50_Subtype", "stage", "TP53")]))

# Verificar N y eventos para el subgrupo TP53
clinical_tp53 <- clinical %>% filter(!is.na(TP53))
table(clinical_tp53$event_indicator) # ¿Cuántos eventos hay?

modelos <- list(
  M6 = surv_obj ~ HE_clas + age_group + TP53, #modelo de base genomica
  M7 = surv_obj ~ HE_clas + age_group + stage + TP53, #modelo clinico genomico (uso stage porque es lo mas potente y resume lo de TMN)
  M8 = surv_obj ~ HE_clas + age_group + PAM50_Subtype + TP53, #modelo de competencia genomica 
  M9 = surv_obj ~ HE_clas + age_group + stage + PAM50_Subtype + TP53 #Modelo compelto clinico genomico
)

# 1. Ejecutar los modelos
resultados_cox <- lapply(modelos, function(f) coxph(f, data = clinical))

# install.packages("broom")
library(broom)
library(dplyr)

# 2. Evaluación del Supuesto de Proporcionalidad (Schoenfeld)
# Guardamos los tests en una lista para revisarlos uno a uno
schoenfeld_tests <- lapply(resultados_cox, function(m) {
  test <- cox.zph(m)
  return(test)
})

# Imprimir los p-valores globales de cada modelo
cat("--- P-valores Globales del Test de Schoenfeld ---\n")
lapply(names(schoenfeld_tests), function(n) {
  cat(n, ": ", schoenfeld_tests[[n]]$table["GLOBAL", "p"], "\n")
})

# Imprimir la tabla detallada de cada modelo
for (nombre in names(schoenfeld_tests)) {
  cat("\n==========================================\n")
  cat("DETALLE TEST DE SCHOENFELD:", nombre, "\n")
  cat("==========================================\n")
  print(schoenfeld_tests[[nombre]])
}


# 3. Evaluación de Multicolinealidad (VIF)
library(car)
vif_resultados <- lapply(resultados_cox, function(m) {
  # Intentamos calcular VIF, capturando el error si la estratificación da problemas
  tryCatch(vif(m), error = function(e) return("Error en cálculo de VIF o Inf detectado"))
})

# Imprimir VIFs (enfocarse en GVIF^(1/(2*Df)))
cat("\n--- Resultados de VIF por Modelo ---\n")
print(vif_resultados)

# 1. Unificar todos los modelos en un solo dataframe
library(dplyr)
library(ggplot2)
library(broom)

# 1. Extraer resultados y añadir las FILAS DE REFERENCIA manualmente
extraer_con_ref <- function(modelo, nombre) {
  res <- tidy(modelo, exponentiate = TRUE, conf.int = TRUE) %>%
    mutate(Modelo = nombre)
  
  # Creamos un dataframe con las referencias (HR = 1, IC = 1)
  referencias <- data.frame(
    term = c("Immune Group High", "Stage Stage I", "PAM50 LumA", "Age < 70", "TP53 wt"),
    estimate = 1, conf.low = 1, conf.high = 1, p.value = NA,
    Modelo = nombre
  )
  
  # Unimos y limpiamos nombres (ajustar según tus nombres de variables reales)
  res <- res %>%
    mutate(term = case_when(
      term == "HE_claslow" ~ "Immune Group Low",
      grepl("stage", term) ~ gsub("stage", "Stage ", term),
      grepl("PAM50", term) ~ gsub("PAM50_Subtype", "PAM50 ", term),
      grepl("age", term) ~ gsub("age_group", "Age ", term),
      grepl("TP53", term) ~ gsub("TP53", "TP53 ", term),
      TRUE ~ term
    )) %>%
    bind_rows(referencias) %>%
    # Solo nos quedamos con las variables que existen en ESTE modelo específico
    filter(term %in% c("Immune Group High", "Immune Group Low", 
                       "Stage Stage I", "Stage Stage II", "Stage Stage III", "Stage Stage IV",
                       "PAM50 LumA", "PAM50 LumB",
                       "Age < 70", "Age > 70",
                       "TP53 wt", "TP53 mutated"))
  return(res)
}

# 2. Unificar todos los modelos
df_grafico <- bind_rows(lapply(names(resultados_cox), function(x) extraer_con_ref(resultados_cox[[x]], x)))

# 3. DEFINIR EL ORDEN DEL EJE Y (De abajo hacia arriba para que Immune Group quede arriba)
orden_niveles <- rev(c(
  "Immune Group High", "Immune Group Low",
  "Age < 70", "Age > 70",
  "TP53 wt", "TP53 mutated",
  "Stage Stage I", "Stage Stage II", "Stage Stage III", "Stage Stage IV",
  "PAM50 LumA", "PAM50 LumB"
))

df_grafico$term <- factor(df_grafico$term, levels = orden_niveles)

# 4. Paleta de colores personalizada
mis_colores <- c(
  "M6"  = "#b15928",
  "M7"  = "#ccebc5", #"#8dd3c7", "#ccebc5"
  "M8"  = "#fb8072",
  "M9"  = "#8dd3c7" #"#ffed6f" 
)


# 5. Generación del gráfico
grafico <- ggplot(df_grafico, aes(x = estimate, y = term)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey40") +
  #barras de error (CI)
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), 
                 height = 0.2, 
                 size = 0.3, # ACHIQUE EL TAMAÑO
                 show.legend = FALSE, 
                 color = "black") +
  geom_point(aes(color = Modelo), size = 3, shape = 15, show.legend = FALSE) +
  # Facetado por modelo para replicar la imagen b86dcd
  facet_grid(. ~ Modelo) +
  # Escala logarítmica constante para todos los modelos (0.8 a 25)
  scale_x_log10(limits = c(0.7, 25), breaks = c(1, 2, 5, 10, 20)) +
  # Aplicación de LOS colores específicos
  scale_color_manual(values = mis_colores) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "grey95"),
    strip.text = element_text(face = "bold", size = 10),
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 9, color = "black"),
    panel.grid.minor = element_blank(),
    legend.position = "none" # Los colores ya identifican los paneles
  ) +
  labs(x = "Hazard Ratio (95% CI)")

# Mostrar y guardar para Illustrator
print(grafico)

# 4. Guardar para Illustrator
ggsave(pathLocalResults("HE clas Results/FigureS3_forest_comparativo_TP53.pdf"), grafico, width = 18, height = 9, units = "cm", device = "pdf")


## Obtener tabla para comparar modelos

get_cox_summary <- function(model, model_name = "Model") {
  coefs <- summary(model)$coefficients
  confs <- summary(model)$conf.int
  
  # Extraer n y eventos (necesario para ver la representatividad)
  # R almacena esto en 'nevent' y 'n'
  n_total <- model$n
  eventos_totales <- model$nevent
  
  # Create clean formatted data frame
  df <- data.frame(
    Variable = rownames(coefs),
    HR_CI = sprintf("%.4f (%.4f-%.4f)", confs[, "exp(coef)"], confs[, "lower .95"], confs[, "upper .95"]),
    P_value = coefs[, "Pr(>|z|)"], # Lo dejamos numérico para ajustar FDR luego
    row.names = NULL,
    stringsAsFactors = FALSE
  )
  
  # Añadimos sufijo del modelo a las columnas
  colnames(df)[2:3] <- paste0(c("HR_CI_", "P_"), model_name)
  return(df)
}


Model_6_df <- get_cox_summary(resultados_cox[["M6"]], "Model 6")
Model_7_df <- get_cox_summary(resultados_cox[["M7"]], "Model 7")
Model_8_df <- get_cox_summary(resultados_cox[["M8"]], "Model 8")
Model_9_df <- get_cox_summary(resultados_cox[["M9"]], "Model 9")


# Join by Variable name (use full_join to preserve all rows)
library(dplyr)

summary_combined <- Model_6_df %>%
  full_join(Model_7_df, by = "Variable") %>%
  full_join(Model_8_df, by = "Variable") %>%
  full_join(Model_9_df, by = "Variable")


write.csv(summary_combined, pathLocalResults("HE clas Results/multivariable_cox_models_TP53_table_without_prior_treatment.csv"), row.names = FALSE)
library(openxlsx)
write.xlsx(summary_combined, pathLocalResults("HE clas Results/multivariable_cox_models_TP53_table_excel_without_prior_treatment.xlsx"), rowNames = FALSE)


## Obtener p adj values pra HE clas low-> esto permite confirmar que esta clasificacion es estable para todos los modelos
library(dplyr)
library(tidyr)

# 1. Identificar las columnas que contienen p-valores
columnas_p <- grep("^P_Model", names(summary_combined), value = TRUE)

# 2. Convertir a formato largo para procesar todos los p-valores juntos
# Esto facilita la aplicación de p.adjust a toda la bolsa de tests realizados
p_long <- summary_combined %>%
  select(Variable, all_of(columnas_p)) %>%
  pivot_longer(cols = all_of(columnas_p), 
               names_to = "Modelo", 
               values_to = "p_val_original") %>%
  # Quitamos los NA porque no se pueden ajustar
  filter(!is.na(p_val_original))

# 3. Calcular el P-Adjusted (FDR / Benjamini-Hochberg)
p_long <- p_long %>%
  mutate(p_val_adjusted = p.adjust(p_val_original, method = "fdr"))

# 4. Volver al formato ancho original para unirlo a tu tabla
p_adjusted_wide <- p_long %>%
  select(Variable, Modelo, p_val_adjusted) %>%
  mutate(Modelo = gsub("P_", "P_Adj_", Modelo)) %>% # Renombramos las columnas
  pivot_wider(names_from = Modelo, values_from = p_val_adjusted)

# 5. Unir con la tabla original
summary_final <- left_join(summary_combined, p_adjusted_wide, by = "Variable")

# Ver resultados de tu variable de interés
summary_final %>% 
  filter(Variable == "HE_claslow") %>% 
  select(Variable, starts_with("P_Model"), starts_with("P_Adj"))


########################################################################
### IMMUNE CELLS from HE comparison between HE immune clas ====
library(readxl)
immune_cell_HE <- read_excel(pathLocalDb("immune_cell_HE.xlsx"))
colnames(immune_cell_HE)
columns <- c("ParticipantBarcode","Leukocyte Fraction","Lymphocytes","Neutrophils","Eosinophils","Mast Cells",              
             "Dendritic Cells","Macrophages", "Global_Pattern", "Immune Subtype")

immune_cell_HE <- immune_cell_HE[,columns]


colnames(immune_cell_HE)[colnames(immune_cell_HE) == "ParticipantBarcode"] <- "patient"
colnames(immune_cell_HE)[colnames(immune_cell_HE) == "Lymphocytes"] <- "Lymphocytes_HE"
colnames(immune_cell_HE)[colnames(immune_cell_HE) == "Neutrophils"] <- "Neutrophils_HE"
colnames(immune_cell_HE)[colnames(immune_cell_HE) == "Eosinophils"] <- "Eosinophils_HE"
colnames(immune_cell_HE)[colnames(immune_cell_HE) == "Mast Cells"] <- "Mast_Cells_HE"
colnames(immune_cell_HE)[colnames(immune_cell_HE) == "Dendritic Cells"] <- "Dendritic_Cells_HE"
colnames(immune_cell_HE)[colnames(immune_cell_HE) == "Macrophages"] <- "Macrophages_HE"
colnames(immune_cell_HE)[colnames(immune_cell_HE) == "Leukocyte Fraction"] <- "Leukocyte_Fraction_HE"
colnames(immune_cell_HE)[colnames(immune_cell_HE) == "Immune Subtype"] <- "Immune_Subtype_HE"

library(dplyr)
clinical <- clinical %>%
  left_join(immune_cell_HE, by = "patient")

columns_to_convert <- colnames(clinical[,colnames(clinical) %in% colnames(immune_cell_HE) & 
                                          colnames(clinical)!= "patient" & 
                                          colnames(clinical)!= "Immune_Subtype_HE" & 
                                          colnames(clinical)!= "Global_Pattern"])

for (col_name in columns_to_convert) {
  clinical[[col_name]] <- as.numeric(clinical[[col_name]])
}

for (col_name in columns_to_convert){
  print(summary(clinical[[col_name]]))
  hist(clinical[[col_name]], main = col_name)
}

#Comparison of HE cells proportions between HE clas

HE_cell_wilcox_result <- list()

for (col_name in columns_to_convert) {
  result <- wilcoxon_test_function(clinical, col_name, "HE_clas", "high")
  HE_cell_wilcox_result[[col_name]] <- result
  print(HE_cell_wilcox_result) 
}

cell_types <- names(HE_cell_wilcox_result)
p_values <- sapply(HE_cell_wilcox_result, function(x) x[["test"]][["p.value"]])
test_statistics <- sapply(HE_cell_wilcox_result, function(x) x[["test"]][["statistic"]])
effect_sizes <- sapply(HE_cell_wilcox_result, function(x) x[["effect_size"]][["W"]])

HE_cell_results_table <- data.frame(
  cell_type = cell_types,
  p_value = p_values,
  test_statistic = test_statistics,
  Effect_Size = effect_sizes
)

HE_cell_results_table$Significance <- ifelse(HE_cell_results_table$p_value < 0.05, "*", "")
HE_cell_results_table <- HE_cell_results_table[order(HE_cell_results_table$p_value), ]
HE_cell_results_table

sign_cell_types <- c("Lymphocytes_HE", "Macrophages_HE", "Mast_Cells_HE")

annotations <- list(
  list(x = 1.5, y = 0.75, stars = 1, variable = c("Lymphocytes_HE", "Macrophages_HE", "Mast_Cells_HE"))
)
create_grouped_barplot_median_sem(data = clinical, #a data frame with in columns should be "group_variable" and each "variable" to be plotted, in rows: samples
                                  group = "HE_clas",
                                  column_names = sign_cell_types,
                                  #y_limit= 0.05,
                                  annotations = annotations,
                                  y_text = "Cells proportions")

create_grouped_boxplot(data = clinical, #a data frame with in columns should be "group_variable" and each "variable" to be plotted, in rows: samples
                       group = "HE_clas",
                       column_names = sign_cell_types,
                       #y_limit= 0.05,
                       annotations = annotations,
                       y_text = "Cells proportions")


## Barras apiladas

clinical_plot <- clinical[,c("HE_clas", "Lymphocytes_HE","Neutrophils_HE","Eosinophils_HE","Mast_Cells_HE",              
                             "Dendritic_Cells_HE","Macrophages_HE")]
head(clinical_plot)

# 1st transform data from wide to long format using pivot_longer() (from the tidyverse) or gather() (from tidyr).
#library(tidyverse)

# Transform data from wide to long format
clinical_long <- clinical_plot %>%
  pivot_longer(
    cols = c(Lymphocytes_HE, Neutrophils_HE, Eosinophils_HE, Mast_Cells_HE, Dendritic_Cells_HE, Macrophages_HE),
    names_to = "Cell_Type",
    values_to = "Proportion"
  ) %>%
  filter(!is.na(Proportion))
# Calculate mean proportion and percentage per cell type and HE_clas group
clinical_avg <- clinical_long %>%
  group_by(HE_clas, Cell_Type) %>%
  summarise(Percentage = mean(Proportion)* 100, .groups = "drop")
# Calculate label position (cumulative middle of each stacked bar)
clinical_avg <- clinical_avg %>%
  group_by(HE_clas) %>%
  arrange(HE_clas, desc(Percentage)) %>%
  mutate(pos = cumsum(Percentage) - (Percentage / 2)) %>%
  ungroup()
# Only keep labels for segments > 3%
clinical_avg <- clinical_avg %>%
  mutate(label = ifelse(Percentage > 2, sprintf("%.1f%%", Percentage), ""))

# Plot
HE_cell_plot <- ggplot(clinical_avg, aes(x = HE_clas, y = Percentage, fill = Cell_Type)) +
  geom_bar(stat = "identity", position = "stack", color = "white") +
  geom_text(aes(y = pos, label = label), color = "black", size = 3, na.rm = TRUE) +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  labs(x = "HE Classification", y = "Mean Cell Proportion (%)", fill = "Cell Type") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 7), 
    axis.text.y = element_text(size = 7),
    axis.title.y = element_text(size = 8),
    # RELACIÓN DE ASPECTO: Controla que el carril no sea un tubo vertical
    aspect.ratio = 1, 
    legend.position = "bottom",
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 7),
    legend.key.size = unit(0.3, "cm"))
HE_cell_plot

ggsave(pathLocalResults("HE clas Results/Figure_3a_HE_cells_without_treatment.pdf"), HE_cell_plot, width = 9, height = 9, units = "cm", device = "pdf")


########################################################################################
##### Correlation analysis between TIL HE and the sum of cells by Quantiseq ====

#Total TIL HE vs total quantiseq mononuclear
cor_quanti_til <-cor.test(rowSums(clinical[,c("B_cell", "NK_cell","T_cell_CD4_non_regulatory","T_cell_CD8","T_cell_regulatory_Tregs")]), clinical[,"til_percentage"])
plot(rowSums(clinical[,c("B_cell", "NK_cell","T_cell_CD4_non_regulatory","T_cell_CD8","T_cell_regulatory_Tregs")]), clinical[,"til_percentage"])
cor_quanti_til <-cor.test(1-(clinical[,"uncharacterized_cell"]), clinical[,"til_percentage"])
cor_quanti_til 
## Macrophages
hist(clinical$Macrophage_M1)
hist(clinical$Macrophage_M2)
hist(rowSums(clinical[,c("Macrophage_M1", "Macrophage_M2")]))
cor_quanti_Mo <-cor.test(rowSums(clinical[,c("Macrophage_M1", "Macrophage_M2")])/(1-(clinical[,"uncharacterized_cell"])), clinical[,"Macrophages_HE"])
cor_quanti_Mo

# Calculate the normalized macrophages-like proportion
clinical$estimated_macrophages <- rowSums(clinical[, c("Macrophage_M1", "Macrophage_M2")], na.rm = TRUE) / 
  (1 - clinical[, "uncharacterized_cell"])
# Plot with linear trend line
library(ggpubr)
ggplot(clinical, aes(x = estimated_macrophages, y = Macrophages_HE)) +
  geom_point(alpha = 0.6, color = "#66C2A5") +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed") +
  stat_cor(method = "spearman", label.x = 0.1, label.y = 0.9, size = 4) +  # adjust position if needed
  labs(
    x = "Estimated Macrophages Proportion (normalized)",
    y = "Macrophages (HE-stained slides)",
    title = "Correlation between Estimated and HE-derived Macrophages Proportions"
  ) +
  theme_minimal()

#Lymphocytes
hist(clinical$Lymphocytes_HE)
hist(rowSums(clinical[,c("B_cell","T_cell_CD4_non_regulatory","T_cell_CD8","T_cell_regulatory_Tregs", "NK_cell")])/(1-(clinical[,"uncharacterized_cell"])))
cor_quanti_Lym <-cor.test(rowSums(clinical[,c("B_cell","T_cell_CD4_non_regulatory","T_cell_CD8","T_cell_regulatory_Tregs", "NK_cell")])/(1-(clinical[,"uncharacterized_cell"])), clinical[,"Lymphocytes_HE"])

library(ggplot2)
# Calculate the normalized lymphocyte-like proportion
clinical$estimated_lymphocytes <- rowSums(clinical[, c("B_cell", "T_cell_CD4_non_regulatory", "T_cell_CD8", "T_cell_regulatory_Tregs", "NK_cell")], na.rm = TRUE) / 
  (1 - clinical[, "uncharacterized_cell"])
# Plot with linear trend line
library(ggpubr)
ggplot(clinical, aes(x = estimated_lymphocytes, y = Lymphocytes_HE)) +
  geom_point(alpha = 0.6, color = "#66C2A5") +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed") +
  stat_cor(method = "spearman", label.x = 0.1, label.y = 0.9, size = 4) +  # adjust position if needed
  labs(
    x = "Estimated Lymphocyte Proportion (normalized)",
    y = "Lymphocytes (HE-stained slides)",
    title = "Correlation between Estimated and HE-derived Lymphocyte Proportions"
  ) +
  theme_minimal()

#Dendritic cells
# Calculate the normalized lymphocyte-like proportion
clinical$estimated_DC <- clinical[, c("Myeloid_dendritic_cell")] / 
  (1 - clinical[, "uncharacterized_cell"])
# Plot with linear trend line
library(ggpubr)
ggplot(clinical, aes(x = estimated_DC, y = Dendritic_Cells_HE)) +
  geom_point(alpha = 0.6, color = "#66C2A5") +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed") +
  stat_cor(method = "spearman", label.x = 0.1, label.y = 0.9, size = 4) +  # adjust position if needed
  labs(
    x = "Estimated Dendritic cell Proportion (normalized)",
    y = "Dendritic cell (HE-stained slides)",
    title = "Correlation between Estimated and HE-derived Dentritic cell Proportions"
  ) +
  theme_minimal()


#### Comparison of Lymphocytes y Macrophages normalized to 1-Uncharacterized between HE Class ====

# Comparaci?n de proporcion de macr?fagos M1 y M2 normalizados entre High y Low 
clinical$estimated_M1 <- clinical[, c("Macrophage_M1")] / 
  (1 - clinical[, "uncharacterized_cell"])
clinical$estimated_M2 <- clinical[, c("Macrophage_M2")] / 
  (1 - clinical[, "uncharacterized_cell"])

Macrophages_columns <- c("estimated_M1", "estimated_M2")
Normalized_Macrophages_wilcox_result <- list()

for (col_name in Macrophages_columns) {
  result <- wilcoxon_test_function(clinical, col_name, "HE_clas", "high")
  Normalized_Macrophages_wilcox_result[[col_name]] <- result
  print(Normalized_Macrophages_wilcox_result) 
}

cell_types <- names(Normalized_Macrophages_wilcox_result)
p_values <- sapply(Normalized_Macrophages_wilcox_result, function(x) x[["test"]][["p.value"]])
test_statistics <- sapply(Normalized_Macrophages_wilcox_result, function(x) x[["test"]][["statistic"]])

Normalized_Macrophages_results_table <- data.frame(
  cell_type = cell_types,
  p_value = p_values,
  test_statistic = test_statistics
)

Normalized_Macrophages_results_table$Significance <- ifelse(Normalized_Macrophages_results_table$p_value < 0.05, "*", "")
Normalized_Macrophages_results_table <- Normalized_Macrophages_results_table[order(Normalized_Macrophages_results_table$p_value), ]
Normalized_Macrophages_results_table


annotations <- list(
  list(x = 1.5, y = 0.75, stars = 1, variable = c("estimated_M2"))
)

create_grouped_boxplot(data = clinical, #a data frame with in columns should be "group_variable" and each "variable" to be plotted, in rows: samples
                       group = "HE_clas",
                       column_names = Macrophages_columns,
                       #y_limit= 0.05,
                       annotations = annotations,
                       y_text = "Quantiseq Normalized to 1- Uncharacterized Cells proportions")

#Comparaci?n proporcion de TCD4, Treg, TCD4, B cells and NK cell normalizados
clinical$estimated_Bcells <- clinical[, c("B_cell")] / 
  (1 - clinical[, "uncharacterized_cell"])
clinical$estimated_TCD4 <- clinical[, c("T_cell_CD4_non_regulatory")] / 
  (1 - clinical[, "uncharacterized_cell"])
clinical$estimated_TCD8 <- clinical[, c("T_cell_CD8")] / 
  (1 - clinical[, "uncharacterized_cell"])
clinical$estimated_Tregs <- clinical[, c("T_cell_regulatory_Tregs")] / 
  (1 - clinical[, "uncharacterized_cell"])
clinical$estimated_NK <- clinical[, c("NK_cell")] / 
  (1 - clinical[, "uncharacterized_cell"])

#Comparison between HE clas
Lymphocytes_columns <- c("estimated_Bcells", "estimated_TCD4", "estimated_TCD8", "estimated_Tregs", "estimated_NK")
Normalized_Lymphocytes_wilcox_result <- list()

for (col_name in Lymphocytes_columns) {
  result <- wilcoxon_test_function(clinical, col_name, "HE_clas", "high")
  Normalized_Lymphocytes_wilcox_result[[col_name]] <- result
  print(Normalized_Lymphocytes_wilcox_result) 
}

cell_types <- names(Normalized_Lymphocytes_wilcox_result)
p_values <- sapply(Normalized_Lymphocytes_wilcox_result, function(x) x[["test"]][["p.value"]])
test_statistics <- sapply(Normalized_Lymphocytes_wilcox_result, function(x) x[["test"]][["statistic"]])

Normalized_Lymphocytes_results_table <- data.frame(
  cell_type = cell_types,
  p_value = p_values,
  test_statistic = test_statistics
)

Normalized_Lymphocytes_results_table$Significance <- ifelse(Normalized_Lymphocytes_results_table$p_value < 0.05, "*", "")
Normalized_Lymphocytes_results_table <- Normalized_Lymphocytes_results_table[order(Normalized_Lymphocytes_results_table$p_value), ]
Normalized_Lymphocytes_results_table


annotations <- list(
  list(x = 1.5, y = 0.25, stars = 1, variable = c("estimated_TCD8", "estimated_Tregs", "estimated_NK"))
)

create_grouped_boxplot(data = clinical, #a data frame with in columns should be "group_variable" and each "variable" to be plotted, in rows: samples
                       group = "HE_clas",
                       column_names = Lymphocytes_columns,
                       y_limit= 0.5,
                       annotations = annotations,
                       y_text = "Quantiseq Normalized to 1- Uncharacterized Cells proportions")

# General plots for Macrophages and Lympohocytes normalized to immune cell fraction
immune_cell_columns <- c("estimated_Bcells", "estimated_TCD4", "estimated_TCD8", "estimated_Tregs", "estimated_NK", "estimated_M1", "estimated_M2")

annotations <- list(
  list(x = 1.5, y = 0.5, stars = 1, variable = c("estimated_TCD8", "estimated_Tregs", "estimated_NK", "estimated_M2"))
)

create_grouped_boxplot(data = clinical,
                       group = "HE_clas",
                       column_names = immune_cell_columns,
                       y_limit= 0.5,
                       annotations = annotations,
                       y_text = "Quantiseq Normalized to 1- Uncharacterized Cells proportions")

quantiseq_violin_plot <- create_grouped_violin_plot(data = clinical,
                           group = "HE_clas",
                           column_names = immune_cell_columns,
                           y_limit= 0.8,
                           annotations = annotations,
                           y_text = "Quantiseq Normalized to 1- Uncharacterized Cells proportions")


ggsave(pathLocalResults("HE clas Results/Figure_3c_quantiseq_violin_4cols_without_treatment.pdf"), quantiseq_violin_plot, width = 18, height = 9, units = "cm", device = "pdf")


##### Comparisons of Quantiseq ratios between HE clas ====
# TCD8 to Treg
clinical$TCD8_Treg <- clinical$T_cell_CD8 / clinical$T_cell_regulatory_Tregs
#TCD8 to TCD4
clinical$TCD8_TCD4 <- clinical$T_cell_CD8 / clinical$T_cell_CD4_non_regulatory
clinical$TCD4_Treg <- clinical$T_cell_CD4_non_regulatory / clinical$T_cell_regulatory_Tregs

Lymphocytes_columns <- c("TCD8_Treg", "TCD8_TCD4")
Normalized_Lymphocytes_ratios_wilcox_result <- list()

for (col_name in Lymphocytes_columns) {
  result <- wilcoxon_test_function(clinical, col_name, "HE_clas", "high")
  Normalized_Lymphocytes_ratios_wilcox_result[[col_name]] <- result
  print(Normalized_Lymphocytes_ratios_wilcox_result) 
}

cell_types <- names(Normalized_Lymphocytes_ratios_wilcox_result)
p_values <- sapply(Normalized_Lymphocytes_ratios_wilcox_result, function(x) x[["test"]][["p.value"]])
test_statistics <- sapply(Normalized_Lymphocytes_ratios_wilcox_result, function(x) x[["test"]][["statistic"]])
effect_size <- sapply(Normalized_Lymphocytes_ratios_wilcox_result, function(x) x[["effect_size"]])

Normalized_Lymphocytes_ratios_results_table <- data.frame(
  cell_type = cell_types,
  p_value = p_values,
  test_statistic = test_statistics,
  effect_size =  effect_size
)

Normalized_Lymphocytes_ratios_results_table$Significance <- ifelse(Normalized_Lymphocytes_ratios_results_table$p_value < 0.05, "*", "")
Normalized_Lymphocytes_ratios_results_table <- Normalized_Lymphocytes_ratios_results_table[order(Normalized_Lymphocytes_ratios_results_table$p_value), ]
Normalized_Lymphocytes_ratios_results_table


annotations <- list(
  list(x = 1.5, y = 0.4, stars = 1, variable = c("TCD8_Treg", "TCD8_TCD4"))
)

T_ratio_plot <- create_grouped_boxplot(data = clinical, #a data frame with in columns should be "group_variable" and each "variable" to be plotted, in rows: samples
                       group = "HE_clas",
                       column_names = Lymphocytes_columns,
                       y_limit= 0.5,
                       annotations = annotations,
                       y_text = "T cells ratios")

create_grouped_violin_plot(data = clinical, #a data frame with in columns should be "group_variable" and each "variable" to be plotted, in rows: samples
                       group = "HE_clas",
                       column_names = Lymphocytes_columns,
                       y_limit= 2,
                       annotations = annotations,
                       y_text = "T cells ratios")

ggsave(pathLocalResults("HE clas Results/Figure_3b_T_ratio_without_treatment.pdf"), T_ratio_plot, width = 9, height = 9, units = "cm", device = "pdf")



######################################################

#### Cox univariado para cada subtipo celular HE ====
surv_obj <- Surv(clinical$censored_time,
                 clinical$event_indicator)

sign_cell_types <- c("Lymphocytes_HE", "Macrophages_HE", "Mast_Cells_HE", 
                     "til_percentage", "HE_clas","B_cell","Macrophage_M1","Macrophage_M2",
                     "NK_cell", "T_cell_CD4_non_regulatory","T_cell_CD8","T_cell_regulatory_Tregs",
                     "uncharacterized_cell")
univ_HE_cell_results <- list()
for (col in sign_cell_types) {
  formula <- as.formula(paste("surv_obj ~", col))
  univ_cox <- coxph(formula, data = clinical)
  summary_univ_cox <- summary(univ_cox)
  p_value <- summary_univ_cox[["logtest"]][["pvalue"]]
  CI <- summary_univ_cox$conf.int[, c("lower .95", "upper .95")]
  HR <- summary_univ_cox$coefficients[1, "exp(coef)"]
  univ_HE_cell_results[[col]] <- list(p_value = p_value, HR = HR, CI=CI)
}


# Convert the list to a data frame
cox_results_df <- do.call(rbind, lapply(names(univ_HE_cell_results), function(col) {
  c(Cell_Type = col, 
    P_Value = univ_HE_cell_results[[col]]$p_value, 
    HR = univ_HE_cell_results[[col]]$HR,
    CI = univ_HE_cell_results[[col]]$CI)
}))

# Convert columns to appropriate types
cox_results_df <- data.frame(cox_results_df, stringsAsFactors = FALSE)
cox_results_df$P_Value <- as.numeric(cox_results_df$P_Value)
cox_results_df$HR <- as.numeric(cox_results_df$HR)
cox_results_df$CI.lower..95 <- as.numeric(cox_results_df$CI.lower..95)
cox_results_df$CI.upper..95 <- as.numeric(cox_results_df$CI.upper..95)

# Sort the data frame by P_Value
cox_results_df <- cox_results_df[order(cox_results_df$P_Value), ]

# Print the sorted results
print(cox_results_df)


# print(cox_results_df)
#        Cell_Type     P_Value        HR CI.lower..95 CI.upper..95
# 1 Lymphocytes_HE 0.002184113 0.1011994   0.02299305     0.445409
# 3 Macrophages_HE 0.002503469 9.2955326   2.22265429    38.875558
# 2  Mast_Cells_HE 0.467840451 2.9934862   0.16576432    54.058433


##### Cox multivariado para immune clas, cell proportions, mutations, prognostic factors
# Modelo stage_node: HE, stage (90), age, TP53
estimated_M2_clas_cut_off <- mean(clinical$estimated_M2)
clinical$estimated_M2_clas <- ifelse(clinical$estimated_M2>=estimated_M2_clas_cut_off,
                                     "high", "low")
clinical$estimated_M2_clas <- as.factor(clinical$estimated_M2_clas)


estimated_M1_clas_cut_off <- median(clinical$estimated_M1)
clinical$estimated_M1_clas <- ifelse(clinical$estimated_M1>=estimated_M1_clas_cut_off,
                                     "high", "low")
clinical$estimated_M1_clas <- as.factor(clinical$estimated_M1_clas)

hist(clinical$estimated_Tregs)
estimated_Treg_clas_cut_off <- median(clinical$estimated_Tregs)
clinical$estimated_Treg_clas <- ifelse(clinical$estimated_Tregs>=estimated_Treg_clas_cut_off,
                                     "high", "low")
clinical$estimated_Treg_clas <- as.factor(clinical$estimated_Treg_clas)


surv_obj <- Surv(clinical$censored_time,
                 clinical$event_indicator)

multi_cox <- coxph(surv_obj ~ HE_clas +  estimated_M1 + stage + tt(age_at_index) + TP53, 
                           data = clinical,
                           tt = function(x, t, ...) x * log(t))

multi_cox <- coxph(surv_obj ~ HE_clas + estimated_M1 + node_status + tumor_size + tt(age_at_index) + TP53, 
                   data = clinical,
                   tt = function(x, t, ...) x * log(t))

multi_cox <- coxph(surv_obj ~ HE_clas + TCD8_Treg + node_status + tumor_size + tt(age_at_index) + TP53, 
                   data = clinical,
                   tt = function(x, t, ...) x * log(t))

# Print the summary of the model
summary(multi_cox)

# Forest Plot for multicox Results
ggforest(multi_cox, 
         data = clinical,
         main = "Multivariable cox analysis for HE clas",
         noDigits = 3)



################################################################################
###    DEG ANALYSIS ====

# Filtro rda para quedarme solo con los pacientes que tengo en clinical (con valores de TIL)
TCGA_BRCA_LumAB_RNAseq <- TCGA_BRCA_LumAB_RNAseq[, colData(TCGA_BRCA_LumAB_RNAseq)$patient %in% clinical$patient]
TCGA_BRCA_LumAB_RNAseq@colData@listData <- clinical

dds <- DESeqDataSet(TCGA_BRCA_LumAB_RNAseq, design = ~ HE_clas)


#Prefiltering
smallestGroupSize <- 3
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize # si no funciona Session -> restart R
dds <- dds[keep,]
dds <- DESeq(dds)

#Asigno Symbol a rownames
rownames(dds) <- gsub("\\..*","",rownames(dds))
rownames(dds) <- dds@rowRanges$gene_name
rownames(dds) <- dds@rowRanges@elementMetadata@listData$gene_name

results_names <- resultsNames(dds)
results_names
results_names <- results_names[-1]

##### DEG between high and low ====
res1 <- lfcShrink(dds, contrast = c('HE_clas', 'high', 'low'), type = 'ashr')

# Si es deseable se pueden ordenar los resultados de acuerdo al padj values (de menor a mayor)
res1_ordered <- res1[order(res1$padj),]
res1_ordered$symbol <- rownames(res1_ordered)

res1_ordered_df <- as.data.frame(res1_ordered@listData)

write.table(as.data.frame(res1_ordered), 
            file=pathLocalResults("HE clas Results/DESeq_high_vs_low_HE_clas_without_treatment.txt"))

##### DEG without shrunken para FGSEA ====
res_DEG_gsea <- results(dds, contrast = c('HE_clas', 'high', 'low'))

res_DEG_gsea <- res_DEG_gsea[order(res_DEG_gsea$padj),]
res_DEG_gsea$symbol <- rownames(res_DEG_gsea)

res_DEG_gsea<- as.data.frame(res_DEG_gsea@listData)

# library(org.Hs.eg.db)
# ens2symbol <- AnnotationDbi::select(org.Hs.eg.db,
#                                     key=rownames(res_DEG_gsea), 
#                                     columns="SYMBOL",
#                                     keytype="ENSEMBL")
# library(tidyverse)
# ens2symbol <- as_tibble(ens2symbol)
# ens2symbol
# 
# res_DEG_gsea <- as_tibble(res_DEG_gsea, rownames=NA)
# res_DEG_gsea$ENSEMBL <- rownames(res_DEG_gsea)
# res_DEG_gsea <- inner_join(res_DEG_gsea, ens2symbol, by="ENSEMBL")
# 
# res_DEG_gsea <- res_DEG_gsea %>% 
#   dplyr::select(SYMBOL, stat) %>% 
#   na.omit() %>% 
#   distinct() %>% 
#   group_by(SYMBOL) %>% 
#   summarize(stat=mean(stat))
# res_DEG_gsea
# 
# res_DEG_gsea <- res_DEG_gsea[order(res_DEG_gsea$stat), ]

write.table(as.data.frame(res_DEG_gsea), 
            file=pathLocalResults("HE clas Results/DESeq_gsea_high_vs_low_HE_clas_without_prior_treatment.txt"))

#Plot DEG (not gsea) Load DESeq results
res1_ordered <- read.csv(pathLocalResults("DESeq_high_vs_low_HE_clas_without_prior_treatment.txt"), sep="")

#### ...Volcano plots between high and low indicando genes inmunolgicos ====
#Importo genes inmunol?gicos
immune_genes_immport <- read.delim2(pathLocalDb("ImmuneGeneList.txt"), header = TRUE, sep="\t") # from Immport data base
immune_genes <- data.frame(
  symbol=  unique(immune_genes_immport$Symbol)
)
rm(immune_genes_immport)

##### ...Volcano plots para high vs low indicando genes inmunologicos ====
volcano_data <- data.frame(
  log2FoldChange = res1_ordered$log2FoldChange,
  pvalueadj = res1_ordered$padj,
  symbol = rownames(res1_ordered)
)

# Filter out rows where either log2FoldChange or padj is NA
volcano_data <- volcano_data[!is.na(volcano_data$log2FoldChange) | !is.na(volcano_data$pvalueadj), ]

volcano_data <- merge(volcano_data, immune_genes, by.x = "symbol", by.y = "symbol", all.x = TRUE)

log2cutoff <- 1
qvaluecutoff <- 0.05
library(dplyr)
input <- mutate(volcano_data, 
               immune = ifelse(volcano_data$symbol %in% immune_genes$symbol == TRUE 
                               & volcano_data$pvalueadj < qvaluecutoff 
                               & abs(log2FoldChange) > log2cutoff , 
                               "Immune DEG", 
                               ifelse(volcano_data$pvalueadj< qvaluecutoff & abs(log2FoldChange) > log2cutoff , "DEG", "Not Sig")))

keyvals.color <- ifelse(input$immune == "DEG", "#D3D3D3",
                        ifelse(input$immune == "Immune DEG", "#66c2a5", "#686868"))
keyvals.color[is.na(keyvals.color)] <- "#686868"

names(keyvals.color)[keyvals.color == "#66c2a5"] <- "Immune DEG"
names(keyvals.color)[keyvals.color == "#D3D3D3"] <- "Non-immune DEG"
names(keyvals.color)[keyvals.color == "#686868"] <- "Not Sig"


keyvals.shape <- ifelse(
  input$immune== "Immune DEG", 16, 1)
keyvals.shape[is.na(keyvals.shape)] <- 3

names(keyvals.shape)[keyvals.shape == 16] <- 'Immune DEG'
names(keyvals.shape)[keyvals.shape == 1] <- 'Non-immune'



# Plot without labels
library(EnhancedVolcano)
p <- EnhancedVolcano(input,
                lab = rep("", nrow(input)), # ithermanner is putting rownames(input), but a dash appear inside of the point
                labSize = 0,                 # hides labels completely
                x = 'log2FoldChange',
                y = 'pvalueadj',
                colCustom = keyvals.color,
                colAlpha = 1,                # solid points
                legendPosition = 'right',
                legendLabSize = 12,
                legendIconSize = 4.0,
                drawConnectors = FALSE,
                shapeCustom = keyvals.shape,
                xlim = c(-4,3),
                ylim = c(-1,40)
)

p
library(ggplot2)
p + geom_point(
  data = subset(input, immune == "Immune DEG"),
  aes(x = 'log2FoldChange', y = 'pvalueadj'),
  colour = "#4a90c2",
  size = 3
)


# #Plot usign ggplot
# genes_to_label <- c("ADORA2A", "LAG3", "HAVCR2", "PDCD1", "CD274", "PDCD1LG2", "CTLA4", "IDO1", "CD276", "VTCN1", "CD244", 
#                     "BTLA", "TIGIT", "CD80", "CD86", "VSIR", "CD28", "ICOS", "ICOSLG", "TNFRSF14", "CD160", "TNFSF14", 
#                     "TNFRSF9", "TNFSF9", "TNFRSF4", "CD70", "CD27", "CD40", "CD40LG", "LGALS9", "TNFSF18", "CEACAM1", 
#                     "CD47", "SIRPA", "DNAM1", "PVR", "CD244", "CD48", "TMIGD2", "HHLA2", "BTN2A1", "CD209", "BTN2A2", 
#                     "BTN3A1", "BTNL3", "BTNL9", "CD96", "TDO2", "CD200", "CD200R1", "GZMB","HMCN1", "TP53", "MAP3K1")
# genes_to_label <- immune_genes$symbol
# #genes_to_label <- c("CTLA4", "PDCD1", "CD274", "CD80", "CD86", "TIGIT", "KIR3DL2", "PRF1", "GZMB", "CD40")
# 
# library(ggplot2)
# volc = ggplot(input, aes(x=log2FoldChange, y=-log10(pvalueadj))) + #volcanoplot with log2Foldchange versus pvalue
#   geom_point(aes(col=sig)) + #add points colored by significance
#   
#   geom_point(data= input[!is.na(input$immune),], aes(x=log2FoldChange[!is.na(immune)], y=-log10(pvalueadj)[!is.na(immune)], col=immune), size = 3, shape=18) +
#   
#   
#   scale_colour_manual(name="",  
#                       values = c("DEG"="grey30", "Not Sig"="grey60", "Immune genes"=  #6699CC)) + 
# xlim(-5, 5) + 
#   ylim(0, 60) +
#   geom_hline(yintercept = -log(qvaluecutoff), linetype="dashed", 
#              color = "grey20", linewidth=0.5) + 
#   geom_vline(xintercept = -log2cutoff, linetype = "dashed", color = "grey20") +  # Add vertical line at log2FoldChange = -1
#   geom_vline(xintercept = log2cutoff, linetype = "dashed", color = "grey20") +   # Add vertical line at log2FoldChange = 1
#   labs(title = "DEG Analysis Volcano high vs low HE immune groups in luminal patients", x = "log2 FC", y = "-log10(pvalueadj)") 
# 
# # Use geom_label_repel to add labels in boxes with a white background
# + geom_label_repel(data=input %>% 
#                      filter((abs(log2FoldChange) > log2cutoff) & (pvalueadj < qvaluecutoff) 
#                             &(symbol %in% genes_to_label)
#                      ), 
#                    aes(log2FoldChange, -log10(pvalueadj) 
#                        ,label= ifelse(symbol %in% genes_to_label, symbol, "")
#                    ), 
#                    box.padding = 0.4,     # Adjust the padding inside the label
#                    label.padding = 0.2,   # Adjust padding around the label
#                    fill = "white",        # Set background color for the label
#                    color = "black",       # Label text color
#                    max.overlaps = 3000, 
#                    min.segment.length = 0.05, 
#                    size = 3, 
#                    segment.color = "grey50", 
#                    segment.size = 0.3)
# 
# 
# volc


# Extract genes from res_Immune with padj < 0.05 and Log2FoldChange > |1|
sign_genes <- res1_ordered$symbol[
  res1_ordered$padj <= qvaluecutoff &
    res1_ordered$log2FoldChange > log2cutoff |
    res1_ordered$log2FoldChange < -log2cutoff
]
sign_genes #208 L2FC 1, Padj 0.05


# Check if significant genes are present in Immune_genes_Immport
sign_genes_immune <- sign_genes[sign_genes %in% immune_genes$symbol]
sign_genes_immune #67 L2FC 1, Padj 0.05

up_sign_genes <- res1_ordered$symbol[
  res1_ordered$padj <= qvaluecutoff &
    res1_ordered$log2FoldChange > log2cutoff
]
up_sign_genes 
# 175 upregulados L2FC 1

sign_up_genes_immune <- up_sign_genes[up_sign_genes %in% immune_genes$symbol]
# 59 sign_up_genes_immune

sign_down_genes <- res1_ordered$symbol[
  res1_ordered$padj <= qvaluecutoff &
    res1_ordered$log2FoldChange < -log2cutoff
]
sign_down_genes
# 33 downregulados L2FC 1

sign_down_genes_immune <- sign_down_genes[sign_down_genes %in% immune_genes$symbol]
# 8 sign_down_genes_immune

#### DEG compared with normal ====

#### OJOOO MODIFICAR LA RUTA, AHORA ESTAN EN pathLocalDB  


TCGA_BRCA_normal_RNAseq <- load(pathLocalDb("TCGA_BRCA_RNAseq_normal_female.rda"))
TCGA_BRCA_normal_RNAseq <- TCGA_BRCA_RNAseq_tumor_female_unique_OS
rm(TCGA_BRCA_RNAseq_tumor_female_unique_OS)



# Replace "." with "_" in the column names of both datasets
colnames(TCGA_BRCA_LumAB_RNAseq@colData) <- gsub("\\.", "_", colnames(TCGA_BRCA_LumAB_RNAseq@colData))
colnames(TCGA_BRCA_LumAB_RNAseq@colData) <- gsub("\\-", "_", colnames(TCGA_BRCA_LumAB_RNAseq@colData))

colnames(TCGA_BRCA_normal_RNAseq@colData) <- gsub("\\.", "_", colnames(TCGA_BRCA_normal_RNAseq@colData))
colnames(TCGA_BRCA_normal_RNAseq@colData) <- gsub("\\ ", "_", colnames(TCGA_BRCA_normal_RNAseq@colData))
colnames(TCGA_BRCA_normal_RNAseq@colData) <- gsub("\\-", "_", colnames(TCGA_BRCA_normal_RNAseq@colData))

# Get column names from both datasets
colnames_tumor <- colnames(TCGA_BRCA_LumAB_RNAseq@colData)
colnames_normal <- colnames(TCGA_BRCA_normal_RNAseq@colData)

# Find columns missing in each dataset
missing_in_normal <- setdiff(colnames_tumor, colnames_normal)
missing_in_tumor <- setdiff(colnames_normal, colnames_tumor)

# Add missing columns with NA values to the normal dataset
for (col in missing_in_normal) {
  TCGA_BRCA_normal_RNAseq[[col]] <- NA
}

# Now both datasets have the same columns, perform rbind
TCGA_tumor_normal_RNAseq <- cbind(TCGA_BRCA_LumAB_RNAseq, TCGA_BRCA_normal_RNAseq)

# Add "normal" as a new level to the factor
colData(TCGA_tumor_normal_RNAseq)$HE_clas <- factor(colData(TCGA_tumor_normal_RNAseq)$HE_clas, levels = c(levels(colData(TCGA_tumor_normal_RNAseq)$HE_clas), "normal"))

# Replace NA values in the HE_clas column with "normal"
colData(TCGA_tumor_normal_RNAseq)$HE_clas[is.na(colData(TCGA_tumor_normal_RNAseq)$HE_clas)] <- "normal"
colData(TCGA_tumor_normal_RNAseq)$HE_clas <- relevel(
  as.factor(colData(TCGA_tumor_normal_RNAseq)$HE_clas),
  ref = "normal"
)

save(TCGA_tumor_normal_RNAseq, file=pathLocalDb("TCGA_tumor_normal_luminal_HE_RNAseq_without_treatment.rda"))

#TCGA_tumor_normal_RNAseq <-load(pathLocalDb("TCGA_tumor_normal_luminal_HE_RNAseq.rda"))

dds <- DESeqDataSet(TCGA_tumor_normal_RNAseq, design = ~ HE_clas)

#Prefiltering
smallestGroupSize <- 3
keep <- rowSums(counts(dds)>= 10) >= smallestGroupSize
dds <- dds[keep,]
dds <- DESeq(dds)

#Asigno Symbol a rownames
rownames(dds) <- gsub("\\..*","",rownames(dds))
rownames(dds) <- dds@rowRanges$gene_name

results_names <- resultsNames(dds)
results_names
results_names <- results_names[-1]

###### DEG high compared with normal ====

res1 <- lfcShrink(dds, contrast = c('HE_clas', 'high', "normal"), type = 'ashr')

# Si es deseable se pueden ordenar los resultados de acuerdo al padj values (de menor a mayor)
res1_ordered <- res1[order(res1$padj),]
res1_ordered$symbol <- rownames(res1_ordered)

write.table(as.data.frame(res1_ordered), 
            file=pathLocalResults("HE clas Results/DESeq_high_vs_normal_HE_clas_without_prior_treatment.txt"))

res_high_ordered <- res1_ordered 
# Plot
###### ...Volcano plots for high vs normal indicando genes inmunol?gicos ====

#res1_ordered <- read.csv(pathLocalResults("HE clas Results/DESeq_high_vs_normal_HE_clas.txt"), sep="")

#Importo genes inmunol?gicos
#immune_genes_immport <- read.delim2(pathLocalDb("ImmuneGeneList.txt"), header = TRUE, sep="\t")
#immune_genes <- data.frame(symbol=  unique(immune_genes_immport$Symbol))
#rm(immune_genes_immport)

# Create volcano data
volcano_data <- data.frame(
  log2FoldChange = res1_ordered$log2FoldChange,
  pvalueadj = res1_ordered$padj,
  symbol = rownames(res1_ordered)
)

# Filter out rows where either log2FoldChange or padj is NA
volcano_data <- volcano_data[!is.na(volcano_data$log2FoldChange) | !is.na(volcano_data$pvalueadj), ]

volcano_data <- merge(volcano_data, immune_genes, by.x = "symbol", by.y = "symbol", all.x = TRUE)


log2cutoff <- 1
qvaluecutoff <- 0.05
library(dplyr)
input<- mutate(volcano_data, 
               sig = ifelse(volcano_data$pvalueadj< qvaluecutoff & abs(log2FoldChange) > log2cutoff , "DEG", "Not Sig"),
               immune = ifelse(volcano_data$symbol %in% immune_genes$symbol == TRUE 
                               & volcano_data$pvalueadj < qvaluecutoff 
                               & abs(log2FoldChange) > log2cutoff , 
                               "Immune genes", NA))
genes_to_label <- c("ADORA2A", "LAG3", "HAVCR2", "PDCD1", "CD274", "PDCD1LG2", "CTLA4", "IDO1", "CD276", "VTCN1", "CD244", 
                    "BTLA", "TIGIT", "CD80", "CD86", "VSIR", "CD28", "ICOS", "ICOSLG", "TNFRSF14", "CD160", "TNFSF14", 
                    "TNFRSF9", "TNFSF9", "TNFRSF4", "CD70", "CD27", "CD40", "CD40LG", "LGALS9", "TNFSF18", "CEACAM1", 
                    "CD47", "SIRPA", "DNAM1", "PVR", "CD244", "CD48", "TMIGD2", "HHLA2", "BTN2A1", "CD209", "BTN2A2", 
                    "BTN3A1", "BTNL3", "BTNL9", "CD96", "TDO2", "CD200", "CD200R1","GZMB","HMCN1", "TP53", "MAP3K1", "MKI67")
library(ggplot2)
# Create the volcano plot with filtered labels
volc = ggplot(input, aes(x = log2FoldChange, y = -log10(pvalueadj))) + 
  geom_point(aes(col = sig)) + # Points colored by significance
  
  # Highlight immune genes
  geom_point(data = input[!is.na(input$immune),], 
             aes(x = log2FoldChange[!is.na(immune)], y = -log10(pvalueadj)[!is.na(immune)], 
                 col = immune), 
             size = 3, shape = 18) +
  
  # Set color scheme
  scale_colour_manual(name = "",  
                      values = c("DEG" = "grey30", "Not Sig" = "grey60", "Immune genes" = "#6699CC")) + 
  xlim(-8, 8) + 
  ylim(0, 360) +
  
  # Add significance lines
  geom_hline(yintercept = -log(qvaluecutoff), linetype = "dashed", color = "grey20", linewidth = 0.5) + 
  geom_vline(xintercept = -log2cutoff, linetype = "dashed", color = "grey20") +  
  geom_vline(xintercept = log2cutoff, linetype = "dashed", color = "grey20") +  
  
  # Labels and titles
  labs(title = "DEG Analysis Volcano high immune luminal group vs normal tissue", x = "log2 FC", y = "-log10(pvalueadj)") +
  
  # Use geom_label_repel to add labels in boxes with a white background
  geom_label_repel(data=input %>% 
                     filter((abs(log2FoldChange) > log2cutoff) & (pvalueadj < qvaluecutoff) & 
                              (symbol %in% genes_to_label)), 
                   aes(log2FoldChange, -log10(pvalueadj), 
                       label= ifelse(symbol %in% genes_to_label, symbol, "")), 
                   box.padding = 0.4,     # Adjust the padding inside the label
                   label.padding = 0.2,   # Adjust padding around the label
                   fill = "white",        # Set background color for the label
                   color = "black",       # Label text color
                   max.overlaps = 3000, 
                   min.segment.length = 0.05, 
                   size = 3, 
                   segment.color = "grey50", 
                   segment.size = 0.3)

# Display the plot
print(volc)

# Extract genes from res_Immune with padj < 0.05 and Log2FoldChange > |2|
sign_genes <- res1_ordered$symbol[
  res1_ordered$padj <= qvaluecutoff &
    res1_ordered$log2FoldChange > log2cutoff |
    res1_ordered$log2FoldChange < -log2cutoff
]

#sign_genes: 5072 desregulados

# Check if significant genes are present in Immune_genes_Immport
sign_genes_immune <- sign_genes[sign_genes %in% immune_genes$symbol]
#sing_genes_immune 515 desregulados

sign_up_genes <- res1_ordered$symbol[
  res1_ordered$padj <= qvaluecutoff &
    res1_ordered$log2FoldChange > log2cutoff
]
# 2506 upregulados

sign_up_genes_immune <- sign_up_genes[sign_up_genes %in% immune_genes$symbol]
#sing_genes_immune 227 upregulados

sign_down_genes <- res1_ordered$symbol[
  res1_ordered$padj <= qvaluecutoff &
    res1_ordered$log2FoldChange < -log2cutoff
]
# 2566 downregulados

sign_down_genes_immune <- sign_down_genes[sign_down_genes %in% immune_genes$symbol]
#288 down immuen genes


###### DEG low compared with normal ====
# Se usa el mismo dds (que incluye los normal) cambiando el contraste -> esto est? permitido cdo se usa ashr como m?todo de Shrink
res1 <- lfcShrink(dds, contrast = c('HE_clas', 'low', "normal"), type = 'ashr')

# Si es deseable se pueden ordenar los resultados de acuerdo al padj values (de menor a mayor)
res1_ordered <- res1[order(res1$padj),]
res1_ordered$symbol <- rownames(res1_ordered)

write.table(as.data.frame(res1_ordered), 
            file=pathLocalResults("HE clas Results/DESeq_low_vs_normal_HE_clas_without_treatment.txt"))

res_low_ordered <- res1_ordered
###### ...Volcano plots para low vs normal tissue indicating immune genes and labeling DEG immune checkpoint ====
volcano_data <- data.frame(
  log2FoldChange = res1_ordered$log2FoldChange,
  pvalueadj = res1_ordered$padj,
  symbol = rownames(res1_ordered)
)

# Filter out rows where either log2FoldChange or padj is NA
volcano_data <- volcano_data[!is.na(volcano_data$log2FoldChange) | !is.na(volcano_data$pvalueadj), ]

volcano_data <- merge(volcano_data, immune_genes, by.x = "symbol", by.y = "symbol", all.x = TRUE)


log2cutoff <- 1
qvaluecutoff <- 0.05

input<- mutate(volcano_data, 
               sig = ifelse(volcano_data$pvalueadj< qvaluecutoff & abs(log2FoldChange) > log2cutoff , "DEG", "Not Sig"),
               immune = ifelse(volcano_data$symbol %in% immune_genes$symbol == TRUE 
                               & volcano_data$pvalueadj < qvaluecutoff 
                               & abs(log2FoldChange) > log2cutoff , 
                               "Immune genes", NA))
# Create the volcano plot with filtered labels
volc = ggplot(input, aes(x = log2FoldChange, y = -log10(pvalueadj))) + 
  geom_point(aes(col = sig)) + # Points colored by significance
  
  # Highlight immune genes
  geom_point(data = input[!is.na(input$immune),], 
             aes(x = log2FoldChange[!is.na(immune)], y = -log10(pvalueadj)[!is.na(immune)], 
                 col = immune), 
             size = 3, shape = 18) +
  
  # Set color scheme
  scale_colour_manual(name = "",  
                      values = c("DEG" = "grey30", "Not Sig" = "grey60", "Immune genes" = "#6699CC")) + 
  xlim(-9, 9) + 
  ylim(0, 400) +
  
  # Add significance lines
  geom_hline(yintercept = -log(qvaluecutoff), linetype = "dashed", color = "grey20", linewidth = 0.5) + 
  geom_vline(xintercept = -log2cutoff, linetype = "dashed", color = "grey20") +  
  geom_vline(xintercept = log2cutoff, linetype = "dashed", color = "grey20") +  
  
  # Labels and titles
  labs(title = "DEG Analysis Volcano low vs normal immune groups in luminal patients", x = "log2 FC", y = "-log10(pvalueadj)") +
  
  # Add gene of interest label in white boxes  
  geom_label_repel(data=input %>% 
                     filter((abs(log2FoldChange) > log2cutoff) & (pvalueadj < qvaluecutoff) & 
                              (symbol %in% genes_to_label)), 
                   aes(log2FoldChange, -log10(pvalueadj), 
                       label= ifelse(symbol %in% genes_to_label, symbol, "")), 
                   box.padding = 0.4,     # Adjust the padding inside the label
                   label.padding = 0.2,   # Adjust padding around the label
                   fill = "white",        # Set background color for the label
                   color = "black",       # Label text color
                   max.overlaps = 3000, 
                   min.segment.length = 0.05, 
                   size = 3, 
                   segment.color = "grey50", 
                   segment.size = 0.3)

# Remove grey background
volc + theme_minimal(base_size = 15) + 
  theme(panel.background = element_rect(fill = "white")) # Set plot background to white

# Display the plot
volc

# Extract genes from res_Immune with padj < 0.05 and Log2FoldChange > |2|
sign_genes <- res1_ordered$symbol[
  res1_ordered$padj <= qvaluecutoff &
    res1_ordered$log2FoldChange > log2cutoff |
    res1_ordered$log2FoldChange < -log2cutoff
]

#sign_genes: 4864 desregulados

# Check if significant genes are present in Immune_genes_Immport
sign_genes_immune <- sign_genes[sign_genes %in% immune_genes$symbol]
#sing_genes_immune 485 desregulados

sign_up_genes <- res1_ordered$symbol[
  res1_ordered$padj <= qvaluecutoff &
    res1_ordered$log2FoldChange > log2cutoff
]
# 2309 upregulados

sign_up_genes_immune <- sign_up_genes[sign_up_genes %in% immune_genes$symbol]
#sing_genes_immune 177 upregulados

sign_down_genes <- res1_ordered$symbol[
  res1_ordered$padj <= qvaluecutoff &
    res1_ordered$log2FoldChange < -log2cutoff
]
# 2555 downregulados

sign_down_genes_immune <- sign_down_genes[sign_down_genes %in% immune_genes$symbol]
#308 down immune genes


### DEG and Heatmap for high vs Normal and Low vs Normal for immune function ====

## Heatmaps para High vs Normal y Low vs Normal for each Immune function
res_high_ordered <- read.csv(pathLocalResults("DESeq_high_vs_normal_HE_clas_without_prior_treatment.txt"), sep="")
res_low_ordered <- read.csv(pathLocalResults("DESeq_low_vs_normal_HE_clas_without_prior_treatment.txt"), sep="")
res1_ordered <- read.csv(pathLocalResults("DESeq_high_vs_low_HE_clas_without_prior_treatment.txt"), sep="")

res_high_ordered <- res_high_ordered[,c("log2FoldChange", "padj", "symbol")]
res_low_ordered <- res_low_ordered[,c("log2FoldChange", "padj", "symbol")]
res1_ordered <- res1_ordered[,c("log2FoldChange", "padj", "symbol")]
## Elijo que genes mostrar basado en el resultado de DEG entre high y low, manteniendo solo aquellos con Log2FC > 0.5 
res1_ordered <- res1_ordered %>%
  filter(abs(log2FoldChange) >= 0.5)

#Importo genes inmunol?gicos
immune_genes_immport <- read.delim2(pathLocalDb("ImmuneGeneList.txt"), header = TRUE, sep="\t")
library(readxl)
immune_checkpoint <- read_excel(pathLocalDb("Immune_checkpoint_genes_stimulatory_inhibitory_The_Immune_Landscape_Cancer_listo.xlsx"))

immune_genes_immport$Category <- as.factor(immune_genes_immport$Category)
immune_genes_immport$symbol <- immune_genes_immport$Symbol
immune_genes_immport <- immune_genes_immport[, c("symbol", "Category")]

immune_checkpoint$Category <- immune_checkpoint$Immune_Checkpoint
immune_checkpoint <- immune_checkpoint[, c("symbol", "Category")]
immune_checkpoint$Category <- as.factor(immune_checkpoint$Category)

immune_antigens <- immune_genes_immport[immune_genes_immport$Category=="Antigen_Processing_and_Presentation",]
immune_antigens$Category <- as.factor(as.character(immune_antigens$Category))

immune_chemokines <- immune_genes_immport[immune_genes_immport$Category=="Chemokines" | immune_genes_immport$Category=="Chemokine_Receptors",]
immune_chemokines$Category <- as.factor(as.character(immune_chemokines$Category))

#Combino DEG high vs Normal y Low vs Mormal y genero 3 data frames para checkpoint, antigens, chemokines categories

DEG_high_low_combined <- merge(res_high_ordered, res_low_ordered, by="symbol")

DEG_high_low_combined_checkpoint <- merge(DEG_high_low_combined, immune_checkpoint, by="symbol", all.x=FALSE, all.y = FALSE)
DEG_high_low_combined_checkpoint <- DEG_high_low_combined_checkpoint %>% distinct() #library(dplyr)
rownames(DEG_high_low_combined_checkpoint) <- DEG_high_low_combined_checkpoint$symbol

DEG_high_low_combined_antigens <- merge(DEG_high_low_combined, immune_antigens, by="symbol", all.x=FALSE, all.y = FALSE)
DEG_high_low_combined_antigens <- DEG_high_low_combined_antigens %>% distinct() #library(dplyr)
rownames(DEG_high_low_combined_antigens) <- DEG_high_low_combined_antigens$symbol

DEG_high_low_combined_chemokines <- merge(DEG_high_low_combined, immune_chemokines, by="symbol", all.x=FALSE, all.y = FALSE)
DEG_high_low_combined_chemokines<- DEG_high_low_combined_chemokines %>% distinct() #library(dplyr)
rownames(DEG_high_low_combined_chemokines) <- DEG_high_low_combined_chemokines$symbol

## Colors for annotations
brewer_palette <- brewer.pal(12, "Set3")

##### Checkpoint ====
res1_ordered_checkpoint <- merge(res1_ordered, immune_checkpoint, by="symbol", all.x=FALSE, all.y = FALSE)

#Filter genes based on those DEG between high and low in  res1_ordered_checkpoint
DEG_high_low_combined_checkpoint <- DEG_high_low_combined_checkpoint %>%
  filter(symbol %in% res1_ordered_checkpoint$symbol)

# Prepare matrix of log2FC with rownames
heatmap_matrix_checkpoint <- DEG_high_low_combined_checkpoint %>%
  select(log2FoldChange.x, log2FoldChange.y) %>%
  as.matrix()

# Rename columns for clarity
colnames(heatmap_matrix_checkpoint) <- c("High_vs_Normal", "Low_vs_Normal")
# Create annotation dataframe
row_ann_checkpoint <- DEG_high_low_combined_checkpoint %>%
  select(Category)
# Reorder genes by Category
gene_order <- row_ann_checkpoint %>%
  arrange(Category) %>%
  rownames()

heatmap_matrix_checkpoint <- heatmap_matrix_checkpoint[gene_order, ]
row_ann_checkpoint <- row_ann_checkpoint[gene_order, , drop = FALSE]

# Optional: set color for each category
ann_colors <- list(
  Category = c(
    Inhibitory = brewer_palette[1],
    Stimulatory = brewer_palette[2]
  )
)

checkpoint_heatmap <- pheatmap(
  heatmap_matrix_checkpoint,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = colorRampPalette(c("#2166ac","#6baed6", "#f7f7f7","#fb8072", "#d95f50"))(100),
  breaks = seq(-1.5, 1.5, length.out = 101),  # or use full range of your data
  annotation_row = row_ann_checkpoint,
  annotation_colors = ann_colors,
  show_rownames = TRUE,
  show_colnames = TRUE,
  border_color = NA,
  fontsize = 10,
  main = "Log2 Fold Change - Checkpoint Genes",
  cellwidth = 15,    # Fixed width for all columns
  cellheight = 10,   # Fixed height for all rows
  silent = TRUE
)


##### Antigen procesing ====
res1_ordered_antigens <- merge(res1_ordered, immune_antigens, by="symbol", all.x=FALSE, all.y = FALSE)

#Filter genes based on those DEG between high and low in  res1_ordered_antigens
DEG_high_low_combined_antigens <- DEG_high_low_combined_antigens %>%
  filter(symbol %in% res1_ordered_antigens$symbol)

# Prepare matrix of log2FC with rownames
heatmap_matrix_antigens <- DEG_high_low_combined_antigens %>%
  select(log2FoldChange.x, log2FoldChange.y) %>%
  as.matrix()

# Rename columns for clarity
colnames(heatmap_matrix_antigens) <- c("High_vs_Normal", "Low_vs_Normal")

# Create annotation dataframe
row_ann_antigens <- DEG_high_low_combined_antigens %>%
  select(Category)
# Reorder genes by Category
gene_order <- row_ann_antigens %>%
  arrange(Category) %>%
  rownames()

heatmap_matrix_antigens <- heatmap_matrix_antigens[gene_order, ]
row_ann_antigens <- row_ann_antigens[gene_order, , drop = FALSE]

# Optional: set color for each category
ann_colors <- list(
  Category = c(
    Antigen_Processing_and_Presentation = brewer_palette[3]
  )
)
# Plot
antigen_heatmap <-pheatmap(
  heatmap_matrix_antigens,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = colorRampPalette(c("#2166ac","#6baed6", "#f7f7f7","#fb8072", "#d95f50"))(100),
  breaks = seq(-1.5, 1.5, length.out = 101),  # or use full range of your data
  annotation_row = row_ann_antigens,
  annotation_colors = ann_colors,
  show_rownames = TRUE,
  show_colnames = TRUE,
  border_color = NA,
  fontsize = 10,
  main = "Log2 Fold Change - Antigen Genes",
  cellwidth = 15,    # Fixed width for all columns
  cellheight = 10,   # Fixed height for all rows
  silent = TRUE
)

##### Chemokines ====
res1_ordered_chemokines <- merge(res1_ordered, immune_chemokines, by="symbol", all.x=FALSE, all.y = FALSE)

#Filter genes based on those DEG between high and low in  res1_ordered_chemokines
DEG_high_low_combined_chemokines <- DEG_high_low_combined_chemokines %>%
  filter(symbol %in% res1_ordered_chemokines$symbol)

# Prepare matrix of log2FC with rownames
heatmap_matrix_chemokines <- DEG_high_low_combined_chemokines %>%
  select(log2FoldChange.x, log2FoldChange.y) %>%
  as.matrix()

# Rename columns for clarity
colnames(heatmap_matrix_chemokines) <- c("High_vs_Normal", "Low_vs_Normal")
# Create annotation dataframe
row_ann_chemokines <- DEG_high_low_combined_chemokines %>%
  select(Category)
# Reorder genes by Category
gene_order <- row_ann_chemokines %>%
  arrange(Category) %>%
  rownames()

heatmap_matrix_chemokines <- heatmap_matrix_chemokines[gene_order, ]
row_ann_chemokines <- row_ann_chemokines[gene_order, , drop = FALSE]

# Optional: set color for each category
ann_colors <- list(
  Category = c(
    Chemokine_Receptors = brewer_palette[4],
    Chemokines = brewer_palette[6]
  )
)

chemokine_heatmap <- pheatmap(
  heatmap_matrix_chemokines,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = colorRampPalette(c("#2166ac","#6baed6", "#f7f7f7","#fb8072", "#d95f50"))(100),
  breaks = seq(-1.5, 1.5, length.out = 101),  # or use full range of your data
  annotation_row = row_ann_chemokines,
  annotation_colors = ann_colors,
  show_rownames = TRUE,
  show_colnames = TRUE,
  border_color = NA,
  fontsize = 10,
  main = "Log2 Fold Change - Chemokine Genes",
  cellwidth = 15,    # Fixed width for all columns
  cellheight = 10,   # Fixed height for all rows
  silent = TRUE
)


# To combine multiple pheatmap() plots into a single figure, the most flexible and publication-quality approach is to use the {gridExtra} or {patchwork} package — 
# but since pheatmap() returns a grid object, the easiest method is to convert each pheatmap to a grob and then arrange them.


library(pheatmap)
library(gridExtra)
library(grid)

# 2. Extract the grobs
g1 <- checkpoint_heatmap[[4]]
g2 <- antigen_heatmap[[4]]
g3 <- chemokine_heatmap[[4]]

# 3. Arrange side by side or in a grid
grid.arrange(g1, g2, g3, ncol = 3)  # or nrow = 3 for vertical layout




###GSEA ANALYSIS  -> Ver Hallmark P53_PATHWAY

######################################################################
### GSEA: Gene Set Enrichment Analysis between Immune groups using DEG for GSEA ====
# Script based on https://biostatsquid.com/fgsea-tutorial-gsea/


# Usar DEG obtenidos desde result, no usar shrunken -> ranking hacerlo con stat en lugar de Log2FC y padj
# para eso usé la guia https://stephenturner.github.io/deseq-to-fgsea/ -> obtuve res_DEG_gsea

# Set relevant paths
list.files()
in_path <- "03_results" #Path where are DEG result
out_path <- "03_results" #Path where will be GSEA result
bg_path <- "02_data/Pathways" #Path where are gmt files (gene sets download from GSEA page)

# Functions ===================================================
## Function: Adjacency matrix to list -------------------------
matrix_to_list <- function(pws){
  pws.l <- list()
  for (pw in colnames(pws)) {
    pws.l[[pw]] <- rownames(pws)[as.logical(pws[, pw])]
  }
  return(pws.l)
}

## Function: prepare_gmt --------------------------------------
prepare_gmt <- function(gmt_file, genes_in_data, savefile = FALSE){
  # for debug
  #file <- gmt_files[1]
  #genes_in_data <- df$gene_symbol
  
  # Read in gmt file
  gmt <- gmtPathways(gmt_file)
  hidden <- unique(unlist(gmt))
  
  # Convert gmt file to a matrix with the genes as rows and for each go annotation (columns) the values are 0 or 1
  mat <- matrix(NA, dimnames = list(hidden, names(gmt)),
                nrow = length(hidden), ncol = length(gmt))
  for (i in 1:dim(mat)[2]){
    mat[,i] <- as.numeric(hidden %in% gmt[[i]])
  }
  
  #Subset to the genes that are present in our data to avoid bias
  hidden1 <- intersect(genes_in_data, hidden)
  mat <- mat[hidden1, colnames(mat)[which(colSums(mat[hidden1,])>5)]] # filter for gene sets with more than 5 genes annotated
  # And get the list again
  final_list <- matrix_to_list(mat) # for this we use the function we previously defined
  
  if(savefile){
    saveRDS(final_list, file = paste0(gsub('.gmt', '', gmt_file), '_subset_', format(Sys.time(), '%d%m'), '.RData'))
  }
  
  print('Wohoo! .gmt conversion successfull!:)')
  return(final_list)
}

# Analysis ====================================================
## High Immune vs Low Immune Comparison
## 1. Read in DEG Results-----------------------------------------------------------
list.files(in_path)
res_DEG_gsea <- read.csv(pathLocalResults("DESeq_gsea_high_vs_low_HE_clas_without_prior_treatment.txt"), sep="")
#categories for each gene set
hallmark_category <- read.csv(pathLocalDb("Pathways/hallmark_category.csv"), sep=";")

# Lista de KEGG gene sets con categorías
kegg_category <- c(
  # Cell Growth and Death
  "KEGG_APOPTOSIS" = "Cell Growth and Death",
  "KEGG_CELL_CYCLE" = "Cell Growth and Death",
  "KEGG_P53_SIGNALING_PATHWAY" = "Cell Growth and Death",
  
  # Replication and Repair
  "KEGG_DNA_REPLICATION" = "Replication and Repair",
  "KEGG_BASE_EXCISION_REPAIR" = "Replication and Repair",
  "KEGG_NUCLEOTIDE_EXCISION_REPAIR" = "Replication and Repair",
  "KEGG_MISMATCH_REPAIR" = "Replication and Repair",
  "KEGG_HOMOLOGOUS_RECOMBINATION" = "Replication and Repair",
  "KEGG_NON_HOMOLOGOUS_END_JOINING" = "Replication and Repair",
  
  # Immune System
  "KEGG_ANTIGEN_PROCESSING_AND_PRESENTATION" = "Immune System",
  "KEGG_T_CELL_RECEPTOR_SIGNALING_PATHWAY" = "Immune System",
  "KEGG_B_CELL_RECEPTOR_SIGNALING_PATHWAY" = "Immune System",
  "KEGG_NATURAL_KILLER_CELL_MEDIATED_CYTOTOXICITY" = "Immune System",
  "KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION" = "Immune System",
  "KEGG_CHEMOKINE_SIGNALING_PATHWAY" = "Immune System",
  "KEGG_LEUKOCYTE_TRANSENDOTHELIAL_MIGRATION" = "Immune System",
  "KEGG_HEMATOPOIETIC_CELL_LINEAGE" = "Immune System",
  "KEGG_COMPLEMENT_AND_COAGULATION_CASCADES" = "Immune System",
  "KEGG_FC_GAMMA_R_MEDIATED_PHAGOCYTOSIS" = "Immune System",
  "KEGG_FC_EPSILON_RI_SIGNALING_PATHWAY" = "Immune System",
  "KEGG_TOLL_LIKE_RECEPTOR_SIGNALING_PATHWAY" = "Immune System",
  "KEGG_JAK_STAT_SIGNALING_PATHWAY" = "Immune System",
  "KEGG_RIG_I_LIKE_RECEPTOR_SIGNALING_PATHWAY" = "Immune System",
  "KEGG_NOD_LIKE_RECEPTOR_SIGNALING_PATHWAY" = "Immune System",
  "KEGG_CYTOSOLIC_DNA_SENSING_PATHWAY" = "Immune System",
  "KEGG_INTESTINAL_IMMUNE_NETWORK_FOR_IGA_PRODUCTION" = "Immune System",
  
  # Cancer: Overview
  "KEGG_PATHWAYS_IN_CANCER" = "Cancer: Overview"
)

# Crear el data.frame
kegg_df <- data.frame(
  pathway = names(kegg_category),
  process_category = unname(kegg_category),
  stringsAsFactors = FALSE
)

# Ver
print(kegg_df)

# (Opcional) Guardar como CSV
write.csv(kegg_df, file = "02_data/Pathways/filtered_kegg_pathways_by_category.csv", row.names = FALSE)


## 2. Prepare background genes (gene sets)-----------------------------------------------

# Download gene sets .gmt files
#https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp

# For GSEA
# Filter out the gmt files for KEGG, Reactome and GOBP

list.files(bg_path) #Aca deben estar todos los set de genes que queramos analizar
gmt_files <- list.files(path = bg_path, pattern = '.gmt', full.names = TRUE)
gmt_files #nos va a devolver una lista de nobre de los archivos de todos los .gmt files que quiero analizar

# we need to use a set of background genes (gmt file) that contain ONLY the genes present in our dataframe (genes_in_data), to avoid bias.
genes_in_data <- res_DEG_gsea$symbol

#Prepare KEGG gmt
KEGG_genes <- prepare_gmt(gmt_files[1], genes_in_data, savefile = FALSE)

#Prepare GO BP gmt anfd filter GO terms
GO_BP_genes <- prepare_gmt(gmt_files[2], genes_in_data, savefile = FALSE)
# Filtrar GO:BP relacionados con linfocitos, T cells, B cells, macrófagos o mastocitos usando palabras claves
# Define immune-related terms
immune_keywords <- c("LYMPHOCYTE", "T_CELL", "B_CELL", "MACROPHAGE", "MAST_CELL")
# Build a regular expression pattern to match terms after "_" or at the start
# This will match e.g., "_T_CELL" or "T_CELL" at the beginning
pattern <- paste0("(^|_)(", paste(immune_keywords, collapse = "|"), ")($|_)")
# Filter the list by matching names
immune_GO_BP_genes <- GO_BP_genes[grepl(pattern, names(GO_BP_genes), ignore.case = TRUE)]
# Exclude those with "REGULATION_OF_"
immune_GO_BP_genes <- immune_GO_BP_genes[!grepl("REGULATION_OF_", names(immune_GO_BP_genes), ignore.case = TRUE)]


# # Obtener los ultimos valores de los nombres:
# # Get the names of the gene sets
# gene_set_names <- names(immune_GO_BP_genes)
# # Extract the last word after the last "_"
# last_words <- sub(".*_", "", gene_set_names)
# # Get unique values
# unique_last_words <- unique(last_words)
# # View them
# print(unique_last_words)


# Select interested suffixes
target_suffixes <- c(
  "PROLIFERATION", "ACTIVATION", "DIFFERENTIATION", "IMMUNE_RESPONSE",
  "CHEMOTAXIS", "MEDIATED_IMMUNITY", "PRODUCTION", "MIGRATION", "CYTOTOXICITY"
)

# Create regex pattern to match these at the end of the string
pattern_suffix <- paste0("(", paste(target_suffixes, collapse = "|"), ")$")

# Filter the list: keep only those whose names end in one of the target suffixes
immune_GO_BP_genes <- immune_GO_BP_genes[
  grepl(pattern_suffix, names(immune_GO_BP_genes), ignore.case = TRUE)
]


#Prepare Hallmarks gmt
Hallmark_genes <- prepare_gmt(gmt_files[6], genes_in_data, savefile = FALSE)


rankings <- res_DEG_gsea$stat
#rankings <- sign(res1_ordered$log2FoldChange)*(-log10(res1_ordered$padj)) # we will use the signed p values from spatial DGE as ranking
names(rankings) <- res_DEG_gsea$symbol# genes as names#

head(rankings)
rankings <- sort(rankings, decreasing = TRUE) # sort genes by ranking
plot(rankings)

#Check min and max ranking
max(rankings)
min(rankings)

# Some genes have such low p values that the signed pval is +- inf, we need to change it to the maximum * constant to avoid problems with fgsea
max_ranking <- max(rankings[is.finite(rankings)])
min_ranking <- min(rankings[is.finite(rankings)])
rankings <- replace(rankings, rankings > max_ranking, max_ranking * 10)
rankings <- replace(rankings, rankings < min_ranking, min_ranking * 10)
rankings <- sort(rankings, decreasing = TRUE) # sort genes by ranking

## 4. Run GSEA ---------------------------------------------------------------
# Function to perform GSEA and filter significant independent pathways
perform_gsea <- function(rankings, pathways, minSize = 10, maxSize = 500, nproc = 1) {
  gsea_res <- fgsea(pathways = pathways, stats = rankings, scoreType = 'std', minSize = minSize, maxSize = maxSize, nproc = nproc)
  collapsed_pathways <- collapsePathways(gsea_res[order(padj)][padj < 0.05], pathways, rankings)
  significant_pathways <- gsea_res[pathway %in% collapsed_pathways$mainPathways][order(-NES), ]
  return(significant_pathways)
}

## Un plot por cada set de genes
### Hallmarks
high_vs_low_Hallmark <- fgseaMultilevel(
  pathways = Hallmark_genes,
  stats = rankings,
  nPermSimple = 10000,
  eps = 0
)

### KEGG
high_vs_low_KEGG <- fgseaMultilevel(
  pathways = KEGG_genes,
  stats = rankings,
  nPermSimple = 10000,
  eps = 0
)

## GO_BP_genes
high_vs_low_GO_BP <- fgseaMultilevel(
  pathways = immune_GO_BP_genes, #I use only filter GO BP
  stats = rankings,
  nPermSimple = 10000,
  eps = 0
)

##  Multipanel Layout
# Load necessary libraries
library(patchwork) #Use patchwork (recommended) to layout the plots:
library(forcats)
library(RColorBrewer)

# Add "collection"  and "category" column to each -> category refers to biological function and it was downloaded from https://pmc.ncbi.nlm.nih.gov/articles/PMC4707969/

#Hallmark
hallmark <- high_vs_low_Hallmark %>%
  mutate(pathway = gsub("^HALLMARK_", "", pathway),
         collection = "Hallmark")
hallmark <- as.data.frame(merge(hallmark, hallmark_category, by ="pathway"))
#elimino la columna de tipo lista (contiene la lista de genes) ya que no se puede guardar como csv
hallmark<- hallmark[, !sapply(hallmark, is.list)]
write.csv2(hallmark, file=pathLocalResults("HE clas Results/GSEA_hallmarks_without_prior_treatment.csv"), row.names = FALSE)

#KEGG
high_vs_low_KEGG
# Merge filtrando sólo los pathways de interés (inner join)
KEGG <- merge(kegg_df, high_vs_low_KEGG, by = "pathway")
KEGG<- KEGG %>%
  mutate(pathway = gsub("^KEGG_", "", pathway),
         collection = "KEGG")
KEGG <- KEGG[, !sapply(KEGG, is.list)]
write.csv2(KEGG, file=pathLocalResults("HE clas Results/GSEA_KEGG_without_prior_treatment.csv"), row.names = FALSE)

#GO BP
gobp <- high_vs_low_GO_BP %>%
  mutate(pathway = gsub("^GOBP_", "", pathway),
         collection = "GO_BP")
gobp <- gobp[, !sapply(gobp, is.list)]
write.csv2(gobp, file=pathLocalResults("HE clas Results/GSEA_GO_BP_without_prior_treatment.csv"), row.names = FALSE)


# Combine all
all_pathways <- bind_rows(hallmark, KEGG) %>%
  mutate(
    neglog10padj = -log10(padj),
    collection = factor(collection, levels = c("KEGG", "Hallmark")) # set facet order
  )

# Define consistent palette
colors_set2 <- brewer.pal(8, "Set2")
green_up <- colors_set2[1]     # Green (upregulated)
orange_down <- colors_set2[2]  # Orange (downregulated)

# Helper function to plot each group separately
plot_pathways <- function(data, collection_name) {
  data %>%
    filter(collection == collection_name) %>%
    mutate(
      neglog10padj = -log10(padj),
      pathway = fct_reorder(pathway, NES)
    ) %>%
    ggplot(aes(x = NES, y = pathway, size = neglog10padj, color = NES)) +
    geom_point(alpha = 0.9) +
    scale_color_gradient2(
      low = orange_down, mid = "white", high = green_up, midpoint = 0,
      name = "NES (high vs low)"
    ) +
    scale_size(range = c(3, 8), name = "-log10(padj)") +
    theme_minimal() +
    theme(
      axis.text.y = element_text(size = 8),
      axis.text.x = element_text(size = 10),
      axis.title = element_text(size = 12),
      plot.title = element_text(size = 14, face = "bold"),
      legend.position = "right"
    ) +
    labs(
      title = paste(collection_name, "Pathway Enrichment"),
      x = "Normalized Enrichment Score (NES)", y = NULL
    )
}

# If you already have all_pathways as a single dataframe:
# Make sure "category" has values: KEGG, Hallmark, GO_BP
# and NES, padj columns are present

# Create individual plots
kegg_plot <- plot_pathways(all_pathways, "KEGG")
hallmark_plot <- plot_pathways(all_pathways, "Hallmark")
gobp_plot <- plot_pathways(all_pathways, "GO_BP")

# Save plots (SVG recommended for Illustrator)
ggsave("KEGG_bubbleplot.pdf", kegg_plot, width = 9, height = 9, path = out_path)
ggsave("Hallmark_bubbleplot.pdf", hallmark_plot, width = 9, height = 9, path = out_path)
ggsave("GO_BP_bubbleplot.pdf", gobp_plot, width = 9, height = 9, path = out_path)
