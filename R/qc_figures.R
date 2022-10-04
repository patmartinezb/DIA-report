

plot_n_prots_sample <- function(df, metadata){
  
  df_plot <- df %>%
    select(Protein.Ids,
           matches(metadata$BioReplicate)) %>%
    pivot_longer(!Protein.Ids, names_to = "BioReplicate", values_to = "vals") %>%
    left_join(metadata, by = "BioReplicate") %>%
    select(
      BioReplicate,
      vals,
      Group) %>%
    filter(!is.na(vals)) %>% 
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
  
  return(df_plot)
  
}

plot_na <- function(df, metadata){
  
  df.isna <- df %>%
    select_if(is.numeric) %>%
    gather(key = "BioReplicate", value = "val") %>%
    mutate(isna = is.na(val)) %>%
    group_by(BioReplicate) %>%
    mutate(total = n()) %>%
    group_by(BioReplicate, total, isna) %>%
    summarise(num.isna = n()) %>%
    mutate(pct = num.isna / total * 100) %>%
    left_join(metadata, by = "BioReplicate") %>%
    select(BioReplicate, total, isna, num.isna, pct, Group)
  
  levels <- (df.isna %>% filter(isna == T) %>% arrange(desc(pct)))$BioReplicate
  
  df.plot <- df.isna %>%
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
  
  return(df.plot)
  
}
