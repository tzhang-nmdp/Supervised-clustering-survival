# The function to save heatmap  
         
save_pheatmap_pdf <- function(x, filename, width=20, height=8) 
   {
   stopifnot(!missing(x))
   stopifnot(!missing(filename))
   pdf(filename, width=width, height=height)
   grid::grid.newpage()
   grid::grid.draw(x$gtable)
   dev.off()
   }  
