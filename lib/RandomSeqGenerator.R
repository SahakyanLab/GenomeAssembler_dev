###########################################################################################
###### This function generates a random sequenced composed of the four DNA bases 	   ######
###### A, T, G and C. 																                               ######
###########################################################################################

RandomSeqGenerator <- function(len = len, G.cont = G.cont, C.cont = C.cont, A.cont = A.cont){
  # multiply base content by len of sequence
  G.cont <- round(G.cont*len)
  C.cont <- round(C.cont*len)
  A.cont <- round(A.cont*len)
  T.cont <- len-(G.cont+C.cont+A.cont)
  
  if(!identical(len, (G.cont+C.cont+A.cont+T.cont))){
    stop("RandomSeqGenerator: sum of base content not equal to len of genome!")
  }
  
  seq <- c(rep("G", G.cont), rep("C", C.cont),
           rep("A", A.cont), rep("T", T.cont))
  
  # random sampling of base positions in the sequence
  seq <- sample(x = seq, size = len, replace = FALSE)
  
  return(seq)
}