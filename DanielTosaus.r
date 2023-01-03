library("Biostrings")
v <- DNAString("ATTGTC")
w <- DNAString("ATCGTG")

HammingDist <- function(v, w){
  d <- 0
  for(i in 1:(length(v))){
    if(v[i]!=w[i]){
      d <- d +1
    }
  }
  return(d)
}

HammingDist(v,w)

v <- DNAString("GCA")
DNA <- DNAStringSet(c("ATCCAGCT", "GGGCAACT", "ATGGATCT", "AAGCAACC", "TTGGAACT"))

TotalDist <- function(v, DNA){
  accmindist <- 0
  for(i in 1:(length(DNA))){
    min <- 10
    for(j in 1:(length(DNA[[i]])-length(v))){
      
      #For each DNA sequence we will compute the Hamming distance to v, starting from index j
      #We will advance j until j+length(v) == length of the sequence and we will keep the minimum distance
      #seq <- as.character(DNA[[i]])
      seq <- DNA[[i]]
      seq <- seq[j:(j+length(v))]
      d <- HammingDist(v, seq)
      if(d < min){
        min <- d
      }
    }
    accmindist <- accmindist + min

  }
  return(accmindist)
}

a <- TotalDist(v,DNA)




PatternSearch <- function(DNA, k){
  
  # We will call TotalDist will all the posible combination of v[1:k]
  chars <- c("A","T","C","G")
  # We obtain all possible permutations with length k made with the elements in chars
  combs <- expand.grid(rep(list(chars), 3)) 
  min <- 20
  row <- 0
  for(s in 1:(nrow(combs))){
    seq <- paste(as.vector(unlist(combs[s,])), collapse = "")
    seq <- DNAString(seq)
    d <- TotalDist(seq, DNA)
    # If the total distance of the seq with the set is lower that the previous we update it
    if(d<min){
      min <- d
      row <- s
    }
    
  }
  # We return the sequence with minimum distance
  return(paste(as.vector(unlist(combs[row,])), collapse = ""))
}

PatternSearch(DNA, 3)

