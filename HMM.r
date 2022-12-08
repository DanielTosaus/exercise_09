

HMM <- function(N, M, A, B, pi){
  # We create the probability matrix that will store the 
  # emission probabilities for each char in each state
  P <- matrix(0, length(M), length(N))
  for(i in (1:length(M))){
    for(j in (1:length(N))){
      # In the first case we calculate the probabilities only
      # with emission probability and transition probability
      # from the start to each possible state
      if(i == 1){
        P[i,j] <- A[j,i] + B[j,i]
      }
      else{
        # We have to loop for every state to form the 
        # max{} expression
        # The current state has a different expression
        maxvec <- c((B[j,i-1]+A[j,j]))
        for(k in (1:length(N))){
          if(k!=j){
            maxvec <- append(maxvec, B[k,i-1] + A[k,j])
          }
        }
        # We compute the probability
        P[i,j] <- B[j,i] + max(maxvec)
      }
    }
  }
  # We could have saved the max component of maxvec for every
  # state so we know de direction of the paths buf for full
  # interconected states we can simply compute the max of every
  # column of the matrix P
  
  gen_states <- c()
  for(i in (1:length(M))){
    gen_state <- which.max(P[i,1:(length(N))])
    gen_states <- append(gen_states, gen_state)
  }
  print(paste("Sequence: ", toString(M)))
  print(paste("Generator state: ", toString(gen_states)))
  print(paste("Log Probabilities: ", toString(apply(P,1,max))))
  
  return(P)
}
N <- HMM1$N
M <- HMM1$M
A <- HMM1$A
B <- HMM1$B
pi <- HMM1$pi

P <- HMM(N,M,A,B,pi)


N <- HMM2$N
M <- HMM2$M
A <- HMM2$A
B <- HMM2$B
pi <- HMM2$pi

P <- HMM(N,M,A,B,pi)
