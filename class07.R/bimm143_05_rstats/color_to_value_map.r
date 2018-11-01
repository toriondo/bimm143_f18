map.colors2 <- function (x, high.low, palette) {
 
   #determine percent values of the "high.low" range
  percent <- ((x -high.low[1]) / (high.low[2] - high.low[1]))
  
  #find corresponding index position in the color "palette" (plus one so that it isnt zero)
  index <- round ((length (palette)-1) * percent) + 1
  
  return (palette[index])
}

