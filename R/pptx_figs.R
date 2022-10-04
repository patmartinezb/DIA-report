
# Generate .pptx  ---------------------------------------------------------

gen_pptx <- function(figs, file) {
  
  # Generates .pptx looping inside the list of figures, and creating a new slide for every single one of them
  
  pres <- read_pptx()
  
  for (i in 1:length(figs)){
    
    pres <- add_slide(pres, layout = "Blank", master = "Office Theme")
    
    pres <- ph_with(pres, value = figs[[i]], location = ph_location_fullsize())
 
  }
  
  print(pres, target = file)
  
}
