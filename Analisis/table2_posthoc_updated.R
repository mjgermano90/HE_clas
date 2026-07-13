# Define your 10 clinical factor columns
clinical_factors <- c("stage", "tumor_size", "node_status", "PAM50_Subtype", 
                      "paper_BRCA_Pathology", "er_status_by_ihc", 
                      "pr_status_by_ihc", "hr_status", "her2_status_ihc_fish_y")


# 2. Part A: Global Statistical Tests
# NOTA: Ya NO se ajusta el p_value global entre factores (se elimina el BH aqui).
# El p_value que sale de aca es el crudo (t-test / ANOVA / Wilcoxon / Kruskal-Wallis).
global_results <- map_df(clinical_factors, function(factor_col) {
  df_clean <- clinical_df %>% 
    filter(!is.na(.data[[factor_col]]), !is.na(til_percentage)) %>%
    mutate(!!factor_col := as.factor(.data[[factor_col]]))
  
  levels_count <- n_distinct(df_clean[[factor_col]])
  if (levels_count < 2) return(NULL)
  
  formula_str <- as.formula(paste("til_percentage ~", factor_col))
  model <- lm(formula_str, data = df_clean, na.action = na.exclude)
  shapiro_p <- shapiro.test(residuals(model))$p.value
  is_normal <- shapiro_p >= 0.05
  
  if (is_normal) {
    if (levels_count == 2) {
      p_val <- t.test(formula_str, data = df_clean)$p.value
      method_used <- "Independent t-test"
    } else {
      p_val <- summary(aov(formula_str, data = df_clean))[[1]]["Pr(>F)"][1, 1]
      method_used <- "ANOVA"
    }
  } else {
    if (levels_count == 2) {
      p_val <- wilcox.test(formula_str, data = df_clean)$p.value
      method_used <- "Wilcoxon rank-sum"
    } else {
      p_val <- kruskal.test(formula_str, data = df_clean)$p.value
      method_used <- "Kruskal-Wallis"
    }
  }
  
  tibble(Factor = factor_col, Levels = levels_count, Method = method_used, p_value = p_val)
})
# <-- Se elimino: %>% mutate(adj_p_value = p.adjust(p_value, method = "BH"))


# 3. Part B: Conditional Post-Hoc Testing
# Condicion actualizada: se dispara el post-hoc con el p_value CRUDO (no adj_p_value)
# siempre que Levels > 2 y p_value <= 0.05.
# El ajuste BH ahora vive DENTRO del post-hoc (pairwise.wilcox.test ya lo aplica
# a nivel de comparaciones por pares dentro de cada factor).
post_hoc_summaries <- map_chr(1:nrow(global_results), function(i) {
  row <- global_results[i, ]
  factor_col <- row$Factor
  
  # Condicion: correr post-hoc solo si hay mas de 2 niveles Y el p_value crudo es <= 0.05
  if (row$Levels > 2 && !is.na(row$p_value) && row$p_value <= 0.05) {
    
    df_clean <- clinical_df %>% 
      filter(!is.na(.data[[factor_col]]), !is.na(til_percentage))
    
    # Run pairwise Wilcoxon con ajuste BH para las comparaciones entre grupos
    # (este es el adj_p_value del post-hoc que queda reflejado en cada par)
    pw_test <- pairwise.wilcox.test(df_clean$til_percentage, df_clean[[factor_col]], 
                                    p.adjust.method = "BH", exact = FALSE)
    p_matrix <- pw_test$p.value
    
    # Parse matrix to extract significantly different pairs (adj_p < 0.05)
    sig_pairs <- c()
    pairs_idx <- which(!is.na(p_matrix), arr.ind = TRUE)
    
    if (nrow(pairs_idx) > 0) {
      for (j in 1:nrow(pairs_idx)) {
        r_idx <- pairs_idx[j, 1]
        c_idx <- pairs_idx[j, 2]
        p_pair <- p_matrix[r_idx, c_idx]  # ya es adj_p_value (BH) del post-hoc
        
        if (p_pair < 0.05) {
          group1 <- colnames(p_matrix)[c_idx]
          group2 <- rownames(p_matrix)[r_idx]
          p_str <- if(p_pair < 0.001) "<0.001" else sprintf("%.3f", p_pair)
          sig_pairs <- c(sig_pairs, paste0(group1, " vs ", group2, " (adj.p=", p_str, ")"))
        }
      }
    }
    
    if (length(sig_pairs) == 0) return("No significant pairwise differences")
    return(paste(sig_pairs, collapse = "; "))
    
  } else if (row$Levels == 2 && !is.na(row$p_value) && row$p_value <= 0.05) {
    return("N/A (Only 2 levels)")
  } else {
    return("-") # No significativo globalmente (p crudo > 0.05), no se corre post-hoc
  }
})

# Bind post-hoc results back to global results
global_results$Post_Hoc_Analysis <- post_hoc_summaries


# 4. Part C: Calculate Descriptive Stats per Level
descriptive_stats <- map_df(clinical_factors, function(factor_col) {
  clinical_df %>%
    filter(!is.na(.data[[factor_col]]), !is.na(til_percentage)) %>%
    group_by(Level = as.character(.data[[factor_col]])) %>%
    summarise(
      N = n(),
      Summary_Stat = sprintf("%.1f%% (%.1f-%.1f)", 
                             median(til_percentage), 
                             quantile(til_percentage, 0.25), 
                             quantile(til_percentage, 0.75)),
      .groups = "drop"
    ) %>%
    mutate(Factor = factor_col)
})


# 5. Part D: Merge and Clean for Display
# NOTA: ya no existe adj_p_value a nivel global -> se muestra el p_value crudo del test global.
# El p ajustado (BH) que importa por comparaciones multiples queda dentro de
# Post_Hoc_Analysis (cada par ya trae su "adj.p=").
table_data <- descriptive_stats %>%
  left_join(global_results, by = "Factor") %>%
  mutate(
    p_value = if_else(p_value < 0.001, "< 0.001", sprintf("%.3f", p_value)),
    Factor = gsub("_", " ", Factor),
    Factor = tools::toTitleCase(Factor)
  ) %>%
  select(Factor, Level, N, Summary_Stat, Method, p_value, Post_Hoc_Analysis)


# 6. Part E: Generate the Updated Table 2
table_2_updated <- table_data %>%
  gt(groupname_col = "Factor") %>%
  tab_header(
    title = "Table 2: Association of TIL Percentage with Patient Clinical Characteristics",
    subtitle = "Post-hoc pairwise comparisons (BH-adjusted) run when global test p <= 0.05"
  ) %>%
  cols_label(
    Level = "Characteristic Level",
    N = "N",
    Summary_Stat = "Median TIL % (IQR)",
    Method = "Global Test",
    p_value = "Global p-value",
    Post_Hoc_Analysis = "Significant Pairwise Differences (Post-Hoc, BH-adj.)"
  ) %>%
  cols_align(align = "left", columns = c(Level, Post_Hoc_Analysis)) %>%
  cols_align(align = "center", columns = c(N, Summary_Stat, Method, p_value)) %>%
  tab_options(
    row_group.font.weight = "bold",
    heading.title.font.weight = "bold",
    table.border.top.color = "black",
    table.border.bottom.color = "black",
    table_body.border.bottom.color = "black",
    column_labels.border.bottom.color = "black",
    column_labels.border.top.color = "black"
  )

# Preview Table
print(table_2_updated)
