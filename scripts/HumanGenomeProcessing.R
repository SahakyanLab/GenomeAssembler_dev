setwd("/Volumes/Paddy_Backup/ProjectBoard_Patrick/02-Proof_of_principle/scripts/")
suppressPackageStartupMessages(library(Biostrings))

print("Unpacking the toplevel combined fasta file...", quote=F)
system(paste0("gunzip ../data/GRCh38_latest_genomic.fna.gz"))

# Initialising vectors to store values
A.cont      <- numeric()
C.cont      <- numeric()
G.cont      <- numeric()
T.cont      <- numeric()
vec.length  <- numeric()

# More memory efficient to request a small subset of sequences per iteration
fai <- fasta.index(paste0("../data/GRCh38_latest_genomic.fna"))
pb  <- txtProgressBar(min = 1, max = dim(fai)[1], style=3) 

for(seq in 1:dim(fai)[1]){
  seq.string <- readDNAStringSet(fai[seq, ],format="fasta")
  all.letters <- letterFrequency(seq.string, letters="ACGT", OR=0)
  if(seq==1){
    A.cont      <- as.numeric(all.letters[,"A"])
    C.cont      <- as.numeric(all.letters[,"C"])
    G.cont      <- as.numeric(all.letters[,"G"])
    T.cont      <- as.numeric(all.letters[,"T"])
    vec.length  <- as.numeric(sum(unname(all.letters[1,1:4])))
  } else {
    A.cont      <- A.cont+as.numeric(all.letters[,"A"])
    C.cont      <- C.cont+as.numeric(all.letters[,"C"])
    G.cont      <- G.cont+as.numeric(all.letters[,"G"])
    T.cont      <- T.cont+as.numeric(all.letters[,"T"])
    vec.length  <- vec.length+as.numeric(sum(unname(all.letters[1,1:4])))
  }
  setTxtProgressBar(pb, seq)
}
close(pb)

df <- data.frame(A.cont = A.cont/vec.length,
                 C.cont = C.cont/vec.length,
                 G.cont = G.cont/vec.length,
                 T.cont = T.cont/vec.length)

write.csv(df, file = "../data/HumanGenomeATGC.csv", row.names = FALSE)