## Fitting time course linear models with limma and splines 

fit_limma_poly <- function(type, .x = 2, wide_imp = wide_imp, wide_dat = wide_dat, top_prots = 18, pval_cutoff = 0.01){
          
          # load required packages
          
          library(splines)
          library(limma)
          library(dplyr)
          
          # reformating ----
          mat_imp <- wide_imp %>% 
                    tibble::column_to_rownames("ID") %>%
                    as.matrix()
          
          # Extract time as a continuous variable  -----
          
          levs <- c("06h", "12h","18h", 
                    "24h", "48h", "72h", "96h")
          
          thevars <- colnames(wide_imp)[-1]
          
          timesf <- str_split_fixed(thevars, pattern = "\\_", n = 2)[,1]
          
          repsf <- str_split_fixed(thevars, pattern = "\\_", n = 2)[,2]
          
          factors <- factor(x = timesf, levels = levs)
          
          time <- str_remove(factors, "h") %>% as.numeric()
          
          # create design matrix ----
          
          if (type == "spline"){
                    
                    cubic <- ns(time, df = .x)
                    
                    design <- model.matrix(~cubic)
                    
                    row.names(design) <- thevars
                    
                    design
                    
                    coefs <- 2:(1+.x)
                    
          } else if (type == "polynome"){
                    
                    polyn <- poly(time, degree = .x)
                    
                    design <- model.matrix(~polyn)
                    
                    row.names(design) <- thevars
                    
                    design
                    
                    coefs <- 2:(1+.x)
                    
          } else if (type == "linear"){
                    
                    design <- model.matrix(~time)
                    
                    row.names(design) <- thevars
                    
                    design
                    
                    coefs <- 2
          }
          
          # Fit linear model (`limma`)  ----
          
          fit <- lmFit(mat_imp, 
                       design = design, 
                       method = "ls")  
          
          fit <- eBayes(fit)
          
          toptable <- topTable(fit, coef = coefs, number = Inf) %>% 
                    tibble::rownames_to_column("ID")
          
          fited_vals <- fitted.MArrayLM(fit) %>%
                    as.data.frame() %>%
                    tibble::rownames_to_column("ID") %>%
                    pivot_longer(cols = 2:ncol(.),
                                 names_sep = "\\_",
                                 names_to = c("Time","Rep"),
                                 values_to = "Fitted_Abundance") %>% 
                    mutate(Run = paste0(Time,"_",Rep),
                           protein = ID)
          
          limma_out <- list(limma_fit = fit,
                            toptable = toptable,
                            fitted_values = fited_vals)
          
          return(limma_out)
          
}

# Visualizations ----

## histogram ----

histogram_toptable <- function(toptable, type, .x){
          
          
          np005 <- toptable %>% 
                    filter(adj.P.Val <= 0.05) %>% 
                    pull(ID) %>% 
                    length()
          
          np001 <- toptable %>% 
                    filter(adj.P.Val <= 0.01) %>% 
                    pull(ID) %>% 
                    length()
          
          histo_pvals <- ggplot(toptable) + 
                    geom_histogram(aes(x = adj.P.Val), 
                                   binwidth = 0.005,
                                   fill = "#2AB7CA", 
                                   color = "#e9ecef", 
                                   alpha=0.9) +
                    geom_histogram(data = toptable %>% filter(adj.P.Val <= 0.05),
                                   aes(x = adj.P.Val),
                                   binwidth = 0.005,
                                   fill = "#FE8585", 
                                   color = "#e9ecef", 
                                   alpha=0.9) +
                    geom_histogram(data = toptable %>% filter(adj.P.Val <= 0.01),
                                   aes(x = adj.P.Val),
                                   binwidth = 0.005,
                                   fill = "#FE4A49", 
                                   color = "#e9ecef", 
                                   alpha=0.9) + 
                    geom_vline(xintercept = 0.05,
                               color = "blue",
                               size = 1, linetype = "dashed") +
                    geom_vline(xintercept = 0.01,
                               color = "red",
                               size = 1, linetype = "dashed") +
                    ggtitle(paste0("Distribution of adjusted p-values: ",
                                   " proteins <= 0.05 = ",np005," ;proteins <= 0.01 = ", np001),
                            subtitle = paste0("Fitting = ", type,
                                              " ; Degrees of Freedom = ", .x)) + 
                    theme(axis.text.x = element_text(hjust = 0.5, vjust = 0.1, size = 10),
                          panel.background = element_blank(),
                          panel.grid.major = element_line(color = "grey"),
                          panel.border = element_rect(colour = "black", fill=NA, size=1.5),
                          axis.title=element_text(size=12,face="bold"))  
          
          return(histo_pvals)
          
}

## profile plots ----

vis_profs <- function(data, 
                      toptable, 
                      pval_cutoff, 
                      title1 = NULL, 
                      subtitle1 = NULL, 
                      top_nr = 18, 
                      fited_values, 
                      method = NULL,
                      interesting){
          
          
          ori_long <- data %>% 
                    pivot_longer(cols = 2:ncol(.),
                                 names_sep = "\\_",
                                 names_to = c("Time","Rep"),
                                 values_to = "Abundance") %>% 
                    mutate(Run = paste0(Time,"_",Rep),
                           protein = ID)
          
          
          filt_top <- toptable %>% 
                    dplyr::filter(ID %in% interesting) %>%
                    filter(adj.P.Val <= pval_cutoff) %>%
                    slice_min(order_by = adj.P.Val,
                              n = top_nr)
          
          proteins <- pull(filt_top, ID)
          
          fitednreal <- left_join(ori_long, fited_values) %>% 
                    suppressMessages() %>%
                    suppressWarnings()
          
          toprofplot <- fitednreal %>%
                    dplyr::filter(ID %in% proteins)
          
          bitrout <- clusterProfiler::bitr(toprofplot$ID,
                                               fromType = "UNIPROT",
                                               toType = "SYMBOL", 
                                               OrgDb = org.Hs.eg.db) %>% 
                    dplyr::rename(ID = UNIPROT) %>% 
                    suppressMessages() %>%
                    suppressWarnings()
          
          toprofplot <- left_join(toprofplot, bitrout, by = "ID") %>% 
                    suppressMessages() %>%
                    suppressWarnings() %>% 
                    mutate(SYMBOL = ifelse(is.na(SYMBOL),
                                           yes = ID,
                                           no = SYMBOL))
          
          
          profile_plot1 <- ggplot(data = toprofplot, 
                                  aes(x = Time, y = Abundance, group=ID, color = Time)) +
                    geom_line(stat = "summary") +
                    geom_point() +
                    geom_smooth(data = toprofplot,
                                mapping = aes(x = Time, y = Fitted_Abundance),
                                method = method, color = "red") +
                    facet_wrap(~SYMBOL, ncol = 3) +
                    ggtitle(label = title1,
                            subtitle = subtitle1) + 
                    theme(axis.text.x = element_text(hjust = 0.5, vjust = 0.1, size =10),
                          panel.background = element_blank(),
                          panel.grid.major = element_line(color = "grey"),
                          panel.border = element_rect(colour = "black", fill=NA, size=1),
                          axis.title=element_text(size=12,face="bold"))
          
          return(profile_plot1)
}

## function to prep list for intersection analysis 

fit_limma_tocomp <- function(.x, type, wide_imp = wide_imp, wide_dat = wide_dat, pval_cutoff = 0.01){
          
          # load required packages
          
          library(splines)
          library(limma)
          library(dplyr)
          
          # reformating ----
          mat_imp <- wide_imp %>% 
                    tibble::column_to_rownames("ID") %>%
                    as.matrix()
          
          # Extract time as a continuous variable  -----
          
          levs <- c("06h", "12h","18h", 
                    "24h", "48h", "72h", "96h")
          
          thevars <- colnames(wide_imp)[-1]
          
          timesf <- str_split_fixed(thevars, pattern = "\\_", n = 2)[,1]
          
          repsf <- str_split_fixed(thevars, pattern = "\\_", n = 2)[,2]
          
          factors <- factor(x = timesf, levels = levs)
          
          time <- str_remove(factors, "h") %>% as.numeric()
          
          # create design matrix ----
          
          if (type == "spline"){
                    
                    cubic <- ns(time, df = .x)
                    
                    design <- model.matrix(~cubic)
                    
                    row.names(design) <- thevars
                    
                    design
                    
                    coefs <- 2:(1+.x)
                    
          } else if (type == "polynome"){
                    
                    polyn <- poly(time, degree = .x)
                    
                    design <- model.matrix(~polyn)
                    
                    row.names(design) <- thevars
                    
                    design
                    
                    coefs <- 2:(1+.x)
                    
          } else if (type == "linear"){
                    
                    design <- model.matrix(~time)
                    
                    row.names(design) <- thevars
                    
                    design
                    
                    coefs <- 2
          }
          
          # Fit linear model (`limma`)  ----
          
          fit <- lmFit(mat_imp, 
                       design = design, 
                       method = "ls")  
          
          fit <- eBayes(fit)
          
          toptable <- topTable(fit, coef = coefs, number = Inf) %>% 
                    tibble::rownames_to_column("ID")
          
          ## proteins per condition ----
          
          significan_proteins <- toptable %>% 
                    filter(adj.P.Val <= pval_cutoff) %>% 
                    pull(ID) 
          
          return(significan_proteins)
          
          
          
}

#prof_plots <- vis_profs(data = wide_dat, 
#                        toptable = toptable, 
#                        pval_cutoff = pval_cutoff, 
#                        title1 = paste0("Profile plots of top ",top_prots," with lowest adjusted p-values"),
#                        subtitle1 = paste0("Fitting = ", type," ; Degrees of Freedom = ", .x),
#                        top_nr = top_prots)
#
#
#print(design)
#print(histo_pvals)
#print(prof_plots)