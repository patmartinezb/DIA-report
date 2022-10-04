
# Enrichment dotplot ------------------------------------------------------

enrich_dotplot <- function(df, title){
  
  # Creates dotplot figure, similar to ClusterProfiler dotplot, for every GO that has a p.adjust value < 0.05
  
  if (nrow(filter(df, p.adjust < .05)) > 0){
    
    df$Description <- toTitleCase(df$Description)
    
    fig <- df %>% 
      filter(p.adjust < .05) %>% 
      slice(1:30) %>% 
      separate(GeneRatio, c("v1", "v2"), convert = TRUE) %>% 
      mutate(GeneRatio = v1/v2) %>% 
      select(-v1,
             -v2) %>% 
      ggplot(aes(GeneRatio, reorder(Description, GeneRatio))) +
      geom_point(aes(size = Count, color = p.adjust)) +
      theme_bw() +
      ylab("") +
      scale_color_gradient(high = "blue", low = "red") +
      ggtitle(title)
    
    return(fig)
    
  }
  
}
