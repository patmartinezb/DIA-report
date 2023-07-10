


all_analysis <- function(df, meta, comp, thres, org){
  
  # This function performs the whole analysis, and generates the figures and tables as lists, so they can be used within the app to generate the excel and pptx
  # The enrich_dotplot() and gen_pptx() can be found in /R

# QC plots ----------------------------------------------------------------
  
  # Prep the data to be plotted
  
  df_plot <- df %>%
    select(Protein.Ids,
           matches(meta$BioReplicate)) %>%
    pivot_longer(!Protein.Ids, names_to = "BioReplicate", values_to = "vals") %>%
    left_join(meta, by = "BioReplicate") %>%
    select(
      Protein.Ids,
      BioReplicate,
      vals,
      Group) %>%
    filter(!is.na(vals))
  
  # Number of proteins per sample
  
  fig.1 <- df_plot %>% 
    group_by(BioReplicate, Group) %>%
    count() %>%
    ggplot(aes(reorder(BioReplicate, n), n)) +
    geom_bar(aes(fill = Group), stat = "identity", color = "black") +
    xlab("Sample") +
    ylab("Count") +
    ggtitle("Number of proteins per sample") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90)) +
    scale_fill_brewer(palette = "Set2") +
    geom_text(aes(label = n), vjust = 0.5, hjust = 1.5, angle = 90)
  
  # Number of proteins per Group
  
  fig.2 <- df_plot %>% 
    group_by(Group) %>%
    summarise(n = n_distinct(Protein.Ids)) %>%
    ggplot(aes(reorder(Group, n), n)) +
    geom_bar(aes(fill = Group), stat = "identity", color = "black") +
    xlab("Sample") +
    ylab("Count") +
    ggtitle("Number of proteins per Group") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90),
          legend.position = "none") +
    scale_fill_brewer(palette = "Set2") +
    geom_text(aes(label = n), vjust = 0.5, hjust = 1.5, angle = 90)
  
  # Plot NA
  
  df.isna <- df %>%
    select_if(is.numeric) %>%
    gather(key = "BioReplicate", value = "val") %>%
    mutate(isna = is.na(val)) %>%
    group_by(BioReplicate) %>%
    mutate(total = n()) %>%
    group_by(BioReplicate, total, isna) %>%
    summarise(num.isna = n()) %>%
    mutate(pct = num.isna / total * 100) %>%
    left_join(meta, by = "BioReplicate") %>%
    select(BioReplicate, total, isna, num.isna, pct, Group)
  
  levels <- (df.isna %>% filter(isna == T) %>% arrange(desc(pct)))$BioReplicate
  
  fig.3 <- df.isna %>%
    ggplot() +
    geom_bar(aes(x = reorder(BioReplicate, -pct), y = pct, fill = isna),
             stat = "identity", alpha = 0.8
    ) +
    scale_fill_manual(
      name = "",
      values = c("gray43", "mediumpurple1"),
      labels = c("Present", "Missing")
    ) +
    coord_flip() +
    labs(x = "Samples", y = "% of missing values") +
    theme_minimal() +
    scale_x_discrete(labels = levels, limits = levels) +
    geom_text(
      aes(x = reorder(BioReplicate, -pct), y = pct, label = round(pct, 2)), 
      hjust = .3, size = 2.5,
      position = position_dodge(width = 1),
      inherit.aes = TRUE
    )
  
  # NA heatmap
  
  plot_na <- df %>%
    select_if(is.numeric) %>%
    mutate_all(list(~ ifelse(is.na(.), 0, 1)))



  my_sample_col <- data.frame(Group = meta$Group)



  if (identical(meta$BioReplicate, colnames(plot_na))){

    rownames(my_sample_col) <- meta$BioReplicate

  }


  fig.4 <- as.ggplot(
    pheatmap(plot_na,
           color = c("skyblue", "grey"),
           show_rownames = FALSE,
           annotation_col = my_sample_col)
    )
  
  # Boxplot
  
  fig.5 <- df %>%
    select(Protein.Ids, which(plyr::colwise(is.numeric)(.) == TRUE)) %>%
    pivot_longer(!Protein.Ids, names_to = "BioReplicate", values_to = "counts") %>%
    left_join(meta, by = "BioReplicate") %>%
    ggplot(aes(x = BioReplicate, y = counts, fill = Group)) +
    geom_boxplot() +
    labs(y = "Log2 of abundance", x = "Samples") +
    scale_fill_brewer(palette = "Set2") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 12))
  
  # Density plot
  
  fig.6 <- df %>%
    select(Protein.Ids, which(plyr::colwise(is.numeric)(.) == TRUE)) %>%
    pivot_longer(!Protein.Ids, names_to = "BioReplicate", values_to = "counts") %>%
    left_join(meta, by = "BioReplicate") %>%
    ggplot(aes(x = counts, color = Group)) +
    geom_density() +
    labs(y = "Log2 of abundance", x = "") +
    scale_fill_brewer(palette = "Set2") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 60, vjust = 0.5, size = 12))
  
  # CV plot
  
  cv.mean <- df %>% 
    select_if(is.numeric) %>% 
    mutate(mean = rowMeans(., na.rm = TRUE)) %>% 
    select(mean)
  
  cv <- apply(select_if(df, is.numeric), 1, function(x){
    
    sd <- sd(x, na.rm = TRUE)
    m <- mean(x, na.rm = TRUE)
    
    cv <- sd/m*100
    
    return(cv)
    
  })
  
  cv.df <- data.frame(Protein.Ids = df$Protein.Ids,
                      logmean = cv.mean$mean,
                      cv = cv)
  
  fig.7 <- cv.df %>% 
    ggplot(aes(logmean, cv)) +
    geom_point() +
    theme_classic() +
    xlab("Log2 protein abundance") + ylab("% coefficient of variation") +
    geom_smooth(method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)
  
  
  # Venn diagram
  
  df2 <- df %>% 
    select(Protein.Ids,
           matches(meta$BioReplicate)) %>% 
    pivot_longer(!Protein.Ids, names_to = "BioReplicate", values_to = "vals") %>% 
    na.omit() %>% 
    inner_join(meta, by = "BioReplicate")
  
  
  # Loop because it's not know before hand the number of groups involved in the analysis
  
  list.group <- list()
  
  for (i in 1:length(levels(as.factor(df2$Group)))){
    
    list.group[[i]] <- unique(df2$Protein.Ids[df2$Group == levels(as.factor(df2$Group))[i]])
    
    names(list.group)[i] <- levels(as.factor(df2$Group))[i]
    
  }
  
  
  fig.13 <- as.ggplot(venn::venn(list.group, zcolor = "style", ggplot = TRUE))

  
  

# Filtering ---------------------------------------------------------------

  # Based on the number of NAs - if the ratio of NAs is above or equals a certain pre-definied value, it will be removed
  
  df.filt <- df %>%
    mutate(filt = case_when(
      rowSums(!is.na(across(where(is.numeric)))) / ncol(across(where(is.numeric))) >= thres ~ TRUE,
      TRUE ~ FALSE
    )) %>%
    filter(filt == TRUE) %>%
    select(-filt)
  

# Normalization -----------------------------------------------------------

  # Median normalization
  
  prots <- df.filt %>%
    select_if(is.character)
  
  normData <- select_if(df.filt, is.numeric)
  
  medianNorm <- function(logDat) {
    # Find medians of each sample
    # Divide by median
    # Multiply by mean of medians
    sampleMed <- apply(logDat, 2, median, na.rm=TRUE)
    meanMed <- mean(sampleMed, na.rm=TRUE)
    out <- t(t(logDat) / sampleMed)
    out <- out * meanMed
    return(as.matrix(out))
  }
  
  df.norm <- as.data.frame(medianNorm(normData))
  
  df.norm <- cbind(prots, df.norm)
  
  # Plotting the before and after normalization
  
  p1 <- df.filt %>%
    select(Protein.Ids, which(plyr::colwise(is.numeric)(.) == TRUE)) %>%
    pivot_longer(!Protein.Ids, names_to = "BioReplicate", values_to = "counts") %>%
    left_join(meta, by = "BioReplicate") %>%
    ggplot(aes(x = BioReplicate, y = counts, fill = Group)) +
    geom_boxplot() +
    labs(y = "Log2 of abundance", x = "Samples") +
    scale_fill_brewer(palette = "Set2") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 8))
  
  p2 <- df.norm %>%
    select(Protein.Ids, which(plyr::colwise(is.numeric)(.) == TRUE)) %>%
    pivot_longer(!Protein.Ids, names_to = "BioReplicate", values_to = "counts") %>%
    left_join(meta, by = "BioReplicate") %>%
    ggplot(aes(x = BioReplicate, y = counts, fill = Group)) +
    geom_boxplot() +
    labs(y = "Log2 of abundance", x = "Samples") +
    scale_fill_brewer(palette = "Set2") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 8))
  
  fig.8 <- ggpubr::ggarrange(p1, p2, labels = c("Pre-normalization", "Post-normalization"),
                           common.legend = TRUE, legend = "bottom")
  

# Imputation --------------------------------------------------------------

  # Imputation is using QRILC
  
  df.imp <- select_if(df.norm, is.numeric)
  
  set.seed(42)
  
  df.imp <- as.data.frame(imputeLCMD::impute.QRILC(df.imp)[[1]])
  
  
  df.imp <- cbind(prots, df.imp)
  
  
  # Plotting the before and after imputation
  
  h3 <- df.norm %>%
    select(Protein.Ids, which(plyr::colwise(is.numeric)(.) == TRUE)) %>%
    pivot_longer(!Protein.Ids, names_to = "BioReplicate", values_to = "counts") %>%
    left_join(meta, by = "BioReplicate") %>%
    ggplot(aes(x = counts, color = Group)) +
    geom_density() +
    labs(y = "Log2 of abundance", x = "") +
    scale_fill_brewer(palette = "Set2") +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 60, vjust = 0.5, size = 12),
      legend.position = "bottom"
    )
  
  h4 <- df.imp %>%
    select(Protein.Ids, which(plyr::colwise(is.numeric)(.) == TRUE)) %>%
    pivot_longer(!Protein.Ids, names_to = "BioReplicate", values_to = "counts") %>%
    left_join(meta, by = "BioReplicate") %>%
    ggplot(aes(x = counts, color = Group)) +
    geom_density() +
    labs(y = "Log2 of abundance", x = "") +
    scale_fill_brewer(palette = "Set2") +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 60, vjust = 0.5, size = 12),
      legend.position = "bottom"
    )
  
  fig.9 <- ggpubr::ggarrange(h3, h4, labels = c("Pre-imputation", "Post-imputation"),
                                common.legend = TRUE)

# PCA ---------------------------------------------------------------------
  
  pca_prot <- function(df, meta, var_p, scale = FALSE, labels = "point"){
    
    df <- na.exclude(df)
    
    tt <- as.data.frame(t(select_if(df, is.numeric)))
    
    tt <- tt %>%
      mutate(BioReplicate = rownames(.)) %>%
      inner_join(meta, by = "BioReplicate")
    
    rownames(tt) <- colnames(select_if(df, is.numeric))
    
    res.pca <- prcomp(select_if(tt, is.numeric), scale = FALSE)
    
    fviz_pca_ind(res.pca,
                 geom.ind = labels,
                 habillage = tt[[var_p]],
                 legend.title = "", title = "",
                 invisible = "quali",
                 pointsize = 4,
                 pointshape = 19
    ) +
      theme_classic() +
      theme(legend.title = element_blank()) +
      scale_fill_brewer(palette = "Set2")
    
  }
  
  fig.10 <- pca_prot(df.imp, meta, var_p = "Group")  

# DE analysis -------------------------------------------------------------

  # Careful with the number of rows of the comparison file, if it's just one the make.names() function will mess with it, no returning the necessary df 
  
  if (nrow(comp) == 1){
    
    comp <- as.data.frame(t(data.frame(apply(comp, 2, make.names))))
    
    rownames(comp) <- NULL
    
  } else {
    
    comp <- as_tibble(apply(comp, 2, make.names))
    
  }
  

  comp.v <- comp %>%
    mutate(comp = str_c(Condition, Control, sep = "-")) %>%
    select(comp) %>%
    pull()

  rownames(df.imp) <- df.imp$Protein.Ids

  anno <- rownames(df.imp)

  # Model
  mod <- model.matrix(~0 + as.factor(Group), data = meta)
  rownames(mod) <- colnames(select_if(df.imp, is.numeric))

  colnames(mod) <- sort(unique(meta$Group))


  # We use limma::lmFit to fit the model Fit the model
  fit <- lmFit(select_if(df.imp, is.numeric), mod)

  # Comparisons
  contrast.matrix <- makeContrasts(contrasts = comp.v,
                                   levels = mod)


  # Fit the comparisons
  fit2 <- contrasts.fit(fit, contrast.matrix)

  # Test these coefficients being != 0 with limma::eBayes
  fit2 <- eBayes(fit2)


  # Get results:
  
  # Everything
  
  list.res.all <- list()
  
  for (i in 1:length(comp.v)){
    
    df.all <- topTable(fit2, coef = i, number = Inf, genelist = as.data.frame(anno))
    
    list.res.all[[i]] <- df.all
    
    names(list.res.all)[i] <- comp.v[i]
    
  }
  
  
  # Only DE proteins per comparison
  res <- as.data.frame(summary(decideTests(fit2)))
  
  res.sig <- res %>% 
    filter(Var1 == "NotSig") %>% 
    mutate(sig = case_when(Freq == dim(df.imp)[1] ~ "nosig",
                           TRUE ~ "sig")) %>% 
    filter(sig == "sig") %>% 
    pull(Var2) %>% 
    as.character()
  
  
  # Obtain results dfs depending on the number of comparisons - coef argument 
  
  list.res.sig <- list()
  
  for (i in 1:length(res.sig)){
    
    df.all <- topTable(fit2, coef = match(res.sig[i], comp.v), number = Inf, genelist = as.data.frame(anno))
    
    list.res.sig[[i]] <- df.all
    
    names(list.res.sig)[i] <- res.sig[i]
    
  }
  
  # Volcano plots
  
  list.volcanos <- list()
  
  for (i in 1:length(res.sig)){
    
    df.all <- list.res.sig[[i]]
    
    fig <- as.ggplot(
      EnhancedVolcano(df.all, lab = rownames(df.all), x = "logFC", 
                      y = "adj.P.Val", xlim = c(-2,2), ylim = c(0, 2), 
                      title = res.sig[i], subtitle = "", 
                      xlab = bquote(~Log[2] ~ "fold change"), 
                      ylab = bquote(~-Log[10] ~ adjusted ~ italic(P)), pCutoff = 0.05, 
                      colAlpha = 1, legendPosition = "bottom", legendLabSize = 10, legendIconSize = 3,
                      caption = "")
    )
    
    list.volcanos[[i]] <- fig
    
  }
  
  # Venn diagram
  
  if (length(res.sig) > 1){
  
  venn.list <- list()
  
  for (i in 1:length(res.sig)){
    
    prots <- list.res.sig[[i]] %>% 
      filter(adj.P.Val < .05) %>% 
      pull(anno)
    
    venn.list[[i]] <- prots
    
    names(venn.list)[i] <- gsub("-.*","", res.sig[i])
    
  }
  
  
  fig.11 <- as.ggplot(venn::venn(venn.list, zcolor = "style", ggplot = TRUE)) # It has to be transformed to a ggplot object
  
  } else {
    
    fig.11 <- NULL
    
  }
  
  # Enrichment --------------------------------------------------------------
  
  # Enrichment analysis depending on the results
  
  if (org == "Human"){
    
    org <- org.Hs.eg.db
    
  } else {
    
    org <- org.Mm.eg.db
    
  }
  
  enrich.list <- list()
  
  for (i in seq_along(res.sig)){
    
    nlist <- list()
    
    df.de <- list.res.sig[[i]] %>% 
      filter(adj.P.Val < .05)
    
    go.cc <- enrichGO(gene = df.de$anno,
                      universe = df.imp$Protein.Ids,
                      OrgDb         = org,
                      keyType       = 'UNIPROT',
                      ont           = "CC",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.01,
                      qvalueCutoff  = 0.05)
    
    if (!is.null(go.cc)){
      
      go.cc <- go.cc@result
      
      nlist[[1]] <- go.cc
      names(nlist)[1] <- paste0(res.sig[i], "_GO_cc")

    }
    
    
    go.bp <- enrichGO(gene = df.de$anno,
                      universe = df.imp$Protein.Ids,
                      OrgDb         = org,
                      keyType       = 'UNIPROT',
                      ont           = "BP",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.01,
                      qvalueCutoff  = 0.05)
    
    if (!is.null(go.bp)){

      go.bp <- go.bp@result
      nlist[[2]] <- go.bp
      names(nlist)[2] <- paste0(res.sig[i], "_GO_bp")
      
    }
    
    go.mf <- enrichGO(gene = df.de$anno,
                      universe = df.imp$Protein.Ids,
                      OrgDb         = org,
                      keyType       = 'UNIPROT',
                      ont           = "MF",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.01,
                      qvalueCutoff  = 0.05)
    
    if (!is.null(go.mf)){

      go.mf <- go.mf@result
      nlist[[3]] <- go.mf
      names(nlist)[3] <- paste0(res.sig[i], "_GO_mf")
      
    }
    
    if (length(nlist) > 0){
      
      enrich.list[[i]] <- nlist
      
    }
  }
  
  
  enrich.list <- unlist(enrich.list, recursive = FALSE)
  

# Enrichment figures ------------------------------------------------------
  
  enrich.figs <- list()
  
  if (length(enrich.list) > 0){
  
  for (i in 1:length(enrich.list)){
    
    fig <- enrich_dotplot(enrich.list[[i]], title = names(enrich.list)[i])
    
    enrich.figs[[i]] <- fig
    
    }
  }


# Create list of tables and figures to download ---------------------------
  
  # Pass tables to list
  
  listy <- list(filtered = df.filt,
                norm = df.norm,
                norm_imputed = df.imp)
  
  if (length(list.res.all) > 0){
    
    listy <- append(listy, list.res.all)
    
  }
  
  # if (length(list.res.sig) > 0){
  #   
  #   listy <- append(listy, list.res.sig)
  #   
  # }
  
  if (length(enrich.list) > 0){
    
    listy <- append(listy, enrich.list)
    
  }
  
  listy <- listy %>% discard(is.null) # NULL objects within lists throw errors when compiling
  
  
  # Pass figures to list
  
  listy2 <- list(fig1 = fig.1,
                 fig2 = fig.2,
                 fig3 = fig.3,
                 fig4 = fig.4,
                 fig5 = fig.5,
                 fig6 = fig.6,
                 fig7 = fig.7,
                 fig13 = fig.13,
                 fig8 = fig.8,
                 fig9 = fig.9,
                 fig10 = fig.10,
                 fig11 = fig.11)
  
  if (length(list.volcanos) > 0){
    
    listy2 <- append(listy2, list.volcanos)
    
  }
  
  
  if (length(enrich.figs) > 0){
    
    listy2 <- append(listy2, enrich.figs)
    
  }
  
  listy2 <- listy2 %>% discard(is.null) # NULL objects within lists throw errors when compiling
  
  
  # Create list of lists
  
  list_all <- list(tables = listy,
                   figs = listy2)
  
  return(list_all)
  
}
