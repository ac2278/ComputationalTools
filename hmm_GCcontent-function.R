#' @title HMM GC content
#'
#' @description The function uses a two-state Hidden Markov model (HMM) to identify regions of high GC content and 
#'              could be considered a crude gene-predictor since genes tend to have higher GC content than the surrounding noncoding DNA. 
#'              The state-transition parameters of this model are defined by one parameter: mu, the probability of switching states. 
#'              The probabilities of the self-transitions are 1-mu. Let xi represent the hidden state at position i of the sequence; 
#'              xi can take one of two values: H (high GC content) or L (low GC content). Let yi represent 
#'              the observed nucleotide at position i of the sequence (i.e., A, C, G, or T).
#'              The function implements the Viterbi algorithm to parse the sequence y into GC-rich and -poor regions (i.e., it finds the most probable state path)
#'              and outputs a list of nonoverlapping intervals (u,v) for the predicted GC-rich regions, 
#'              where each u indicates the first base of a GC-rich interval and each v indicates the last base. 
#'              The function also implements the forward and backward algorithms to compute the marginal posterior probability of 
#'              state H at every position in the sequence. The function includes an option to output a plot showing 
#'              these probabilities along the sequence, with position i on the horizontal axis and P( xi=H | y={y1,...yn} ) on the verticle axis.
#' 
#' @param            y (character vector) specifies the bases of the sequence. Element i of y specifies the ith nucleotide of the sequence.
#' @param      initial (named numeric vector; length 2) specifies the initial probabilities of each states, L and H.
#'                     Example:                 L   H 
#'                                            0.5 0.5
#' @param     emission (named numeric matrix) specifies the transition probabilities of the HMM; a 2x4 matrix, where rows represent hidden 
#'                     states (i.e. X = {L,H}) and columns represent emitted symbols (i.e. Y = {A,C,G,T})
#'                     Example:      A    C    G    T
#'                              L 0.32 0.18 0.18 0.32
#'                              H 0.13 0.37 0.37 0.13
#' @param           mu (numeric scalar) specifies the probability of switching states; takes values in [0,1]
#' @param         plot (logical) specifies whether or not to generate a plot showing the marginal posterior probability of 
#'                     state H at every position in the sequence
#' @param       pdf.fn (character) specifies the file name of the pdf to be generated; must be specified if ${plot} parameter is set to TRUE
#' 
#' @return The function _outputs_ a list with five elements:   
#'                (1)     viterbiPATH: (character vector) the most probable state path
#'                (2) GCrichintervals: (character matrix) a list of nonoverlapping intervals (u,v) for the predicted GC-rich regions, 
#                                      where each u indicates the first base of a GC-rich interval and each v indicates the last base. 
#'                (3)         F.final: (numeric) the final value of P(y) as computed by the termination step of the forward algorithm. 
#'                (4)         B.final: (numeric) the final value of P(y) as computed by the termination step of the backward algorithm. 
#'                (5)           gamma: (numeric matrix) the marginal posterior probability of state L anf H at every position in the sequence 
#'                                     (i.e. P( xi | y={y1,...yn} ); dimensions: 2 x n, where n is the length of the sequence
#'         
#' @author Ariel W Chan, \email{ac2278@@cornell.edu}
#' 
#' @export
################################################################################
hmmgc = function(y, initial, emission, mu, plot, pdf.fn){
  if(length(dimnames(emission))==0){ print("Error: The emission matrix must have named rows (hidden states of the HMM) and columns (observation symbols A, C, G, and T, listed in that (alphabetal) order).") } else
    if(length(names(initial))==0){ print("Error: Elements of the initial vector must be named (hidden states of the HMM: L and H).") } else
      n = length(y)
    mode(y) = "character"
    # The hidden states of the HMM: 
    X = rownames(emission)
    
    # The observation symbols of the HMM, listed in alphabetical order: 
    Y = colnames(emission)
    
    # Define the transition matrix of the HMM: 
    # The state-transition parameters of this mdeo are defined by one parameter: mu, the probability of switching states.
    # The probabilities of the self-transitions are 1-mu, as needed for proper conditional probability distributions. 
    transition = matrix(c(1-mu, mu, 
                          mu, 1-mu), byrow=T, nrow=length(X), ncol=length(X), dimnames=list(X, X))
    
    # Use the log transformation of all probabilities to avoid numerical problems (i.e. underflow). 
    # It is more efficient to take the log of all model parameters before running the algorithms, to avoid calling the 
    # logarithm function repeatedly during the dynamic programming iteration. 
    tilde.transition = log(transition)
    tilde.initial = log(initial)
    tilde.emission = log(emission)
    
    # VITERBI ALGORITHM: 
    # Implement the Viterbi algorithm using log probabilities to find the most probable state path.
    # Parse the sequence y into GC-rich and -poor regions. 
    # (1) Initialization: 
    V = POINTER = matrix(nrow=length(X), ncol=n, dimnames=list(X, paste("v_k(i=", 1:n, ")", sep="")))
    for (state in X){
      V[state,1] = tilde.initial[state] + tilde.emission[state,y[1]]
    }
    
    # (2) Recursion: 
    for (t in 2:n){
      for (state in X){
        fromk = c("fromL" = V["L",(t-1)] + tilde.transition["L",state], # F(1,1) + log(P( transition from L to L )) = F(1,2)
                  "fromH" = V["H",(t-1)] + tilde.transition["H",state]) # F(2,1) + log(P( transition from H to L )) = F(1,2)
        V[state,t] = tilde.emission[state, y[t]] + sort(fromk, dec=T)[1] # Maximization step
        POINTER[state,t] = names(sort(fromk, dec=T)[1])
      }
    }
    
    POINTER[,1] = paste("from", rownames(POINTER), sep="")
    
    # (3) Termination: 
    # Find the most probable state path, denoted {x1,...,xn}*.  
    # {x1,...,xn}* = argmax_{x1,...,xn} log(P({y1,...,yn}, {x1,...,xn}))
    # Return the log(P({y1,...,yn}, {x1,...,xn}*)): 
    sort(V[,n], dec=T)[1]
    content = c()
    
    # Parse the sequence y into GC-rich and -poor regions by implementing the traceback procedure: 
    coordinate = c(names(sort(V[,n], dec=T)[1]), n) # current cell 
    while(as.numeric(coordinate[2])>=1){
      if(POINTER[coordinate[1], as.numeric(coordinate[2])]
         == "fromH"){ content = c(content, "rich");
         coordinate = c("H", as.numeric(coordinate[2])-1) } else
           
           if(POINTER[coordinate[1], as.numeric(coordinate[2])] 
              == "fromL"){ content = c(content, "poor");
              coordinate = c("L", as.numeric(coordinate[2])-1) }
    }
    
    
    content = rev(content)
    names(content) = 1:n
    
    intervals = paste(content[1:(n-1)], content[2:n])
    if(content[n]=="rich"){ intervals = c(intervals, "rich poor");
    names(intervals) = paste(1:n, 2:(n+1), sep=":")
    } else names(intervals) = paste(1:(n-1), 2:n, sep=":")
    
    u = gsub('[0-9]+:', '', names(grep("poor rich", intervals, value=T)))
    v = gsub(':[0-9]+', '', names(grep("rich poor", intervals, value=T)))
    intervals = cbind(u, v)
    
    
    
    # FORWARD ALGORITHM: 
    # Implement the forward algorithm in log space. 
    # (1) Initialization: 
    F = matrix(nrow=length(X), ncol=n, dimnames=list(X, paste("f_k(i=", 1:n, ")", sep="")))
    for (state in X){
      F[state,1] = tilde.initial[state] + tilde.emission[state,y[1]]
    }
    
    # (2) Recursion: 
    for (t in 2:n){
      for (state in X){
        fromk = c("fromL" = F["L",(t-1)] + tilde.transition["L",state], # F(1,1) + log(P( transition from L to L )) = F(1,2)
                  "fromH" = F["H",(t-1)] + tilde.transition["H",state]) # F(2,1) + log(P( transition from H to L )) = F(1,2)
        # Summation Step (contrast this to the Viterbi Maximization Step)
        tilde.a = tilde.emission[state,y[t]]+fromk["fromL"]
        tilde.b = tilde.emission[state,y[t]]+fromk["fromH"]
        F[state,t] = tilde.a + log(1+exp(tilde.b-tilde.a)) # notice that F[state,t] does include emission probability for time t
      }
    }
    
    # (3) Final value of P(y) as computed by the termination step of the forward algorithm:
    tilde.a = F["L",n]
    tilde.b = F["H",n]
    F.finalvalue = tilde.a + log(1+exp(tilde.b-tilde.a)) # Compute the log probability of the sequence y, log(P(y)) 
    
    # BACKWARD ALGORITHM: 
    # Implement the backward algorithm in log space. 
    # (1) Initialization: 
    B = matrix(nrow=length(X), ncol=n, dimnames=list(X, paste("b_k(i=", 1:n, ")", sep="")))
    for (state in X){
      B[state,n] = 1
    }
    
    # (2) Recursion: 
    for (t in rev(1:(n-1))){
      for (state in X){
        fromk = c("fromL" = tilde.emission[state,y[(t+1)]] + # P( emit observation y(i+1) at time i+1 given state l at time i+1 ) 
                    B["L",(t+1)] +                           # recursion relation; summation across all possible paths leading to state L at time i+1 
                    tilde.transition["L",state],             # P( transition FROM state L at time i+1 TO state k at time i ) 
                  
                  "fromH" = tilde.emission[state,y[(t+1)]] + # P( emit observation y(i+1) at time i+1 given state l at time i+1 ) 
                    B["H",(t+1)] +                           # recursion relation; summation across all possible paths leading to state H at time i+1 
                    tilde.transition["H",state])             # P( transition FROM state H at time i+1 TO state k at time i ) 
        
        # Summation Step (contrast this to the Viterbi Maximization Step)
        tilde.a = fromk["fromL"]
        tilde.b = fromk["fromH"]
        B[state,t] = tilde.a + log(1+exp(tilde.b-tilde.a))   # Notice that B[state,t] does not include emission probability for time t
      }
    }
    
    # (3) Final value of P(y) as computed by the termination step of the backward algorithm:
    tilde.a = B["L",1] + tilde.initial["L"] + tilde.emission["L",y[1]]
    tilde.b = B["H",1] + tilde.initial["H"] + tilde.emission["H",y[1]]
    B.finalvalue = tilde.a + log(1+exp(tilde.b-tilde.a)) # Compute the log probability of the sequence y, log(P(y)) 
    
    # Compute the marginal posterior probability of state H at every position in the sequence 
    # (i.e. compute P( xi=H | y={y1,...yn} ).
    # See Equation (27) of Rabiner for details. 
    # f_k(i) accounts for the partial observation sequence y1,...yi and state k at time i
    # b_k(i) accounts for the remainder of the observation sequence y(i+1),y(i+2)...y(n) given 
    # state k at time i
    # Use the local estimate of P(y) as the denominator in computing the posterior probabilities for 
    # position i to ensure that P( xi=H | y,model ) + P( xi=L | y,model ) = 1
    # (i.e. the probability of being in state H at time i, given the observation
    # sequence y={y1,...yn} and the model plus the probability of being in state L at time i, 
    # given the observation sequence y={y1,...yn} sum to one). 
    state.L = F["L",]+B["L",] 
    state.H = F["H",]+B["H",] 
    normalization = state.L + log(1+exp(state.H-state.L)) # Compute the local estimate of P(y) and use this as the normalization factor. 
    tilde.gammaH = state.H-normalization
    tilde.gammaL = state.L-normalization
    gamma.H = exp(tilde.gammaH)
    gamma.L = exp(tilde.gammaL)
    gamma = rbind(gamma.L, gamma.H)
    colnames(gamma) = 1:n
    rownames(gamma) = c("L", "H")
    
    
    # Plot the marginal posterior probability of state H at every position i in the sequence. 
    if ((plot==T)==TRUE){ pdf(pdf.fn)
      plot(1:n, gamma.H, cex=0.50, ylab="P( xi=H | y={y1,...yn} )", xlab="position i",
           sub="Red data points indicate the GC-rich regions predicted by the Viterbi algorithm.",
           main="The marginal posterior probability of state H \nat every position i in the sequence.")
      pnt = as.numeric(names(content[which(content=="rich")]))
      points(c(1:n)[pnt], gamma.H[pnt], col = "red", cex=0.50)
      dev.off()
      return(list("viterbiPATH"=content, "GCrichintervals"=intervals, "F.final"=F.finalvalue, "B.final"=B.finalvalue, "gamma"=gamma)) } else
        
        return(list("viterbiPATH"=content, "GCrichintervals"=intervals, "F.final"=F.finalvalue, "B.final"=B.finalvalue, "gamma"=gamma))
}

################################################################################
# Example: parse the sequence y into GC-rich and -poor regions.
# (1) Define the sequence y (i.e. import the sequence y using the funtion importFASTA): 
importFASTA = function(FASTAfilename){
  N = as.numeric(system(paste("grep -c '>'", FASTAfilename, sep=" "), intern=T))
  print(N)
  species = system(paste("grep -o '>[[:alpha:]]*'", FASTAfilename, "| sed -E 's/>//g'", sep=" "), intern=T)
  sequences = system(paste("grep '^[[:alpha:]]'", FASTAfilename, sep=" "), intern=T)
  FASTA = vector(mode="list")
  for (i in 1:N){
    FASTA[species[i]] = strsplit(sequences[i], "")
    names(FASTA[[i]]) = 1:length(FASTA[[i]])
  }
  return(FASTA) 
}
y <- importFASTA("hmmgc_example.fasta")[[1]]  # a DNA sequence of length 1000

# (2) Define the emission probabilities of the HMM: 
# First define the hidden states of the HMM: 
X = c("L", "H")
# Then define the observation symbols of the HMM. List the symbols in alphabetical order: 
Y = c("A", "C", "G", "T")
emission <- matrix(c(0.32, 0.18, 0.18, 0.32,
                     0.13, 0.37, 0.37, 0.13), 
                   byrow=T, nrow=length(X), ncol=length(Y), dimnames=list(X, Y)) 

# (3) Define the state-transition parameters of this model (the transition matrix is defined by one parameter: 
mu <- 0.01

# (4) Define the initial probabilities of the HMM: 
initial <- c("L"=0.50, "H"=0.50) 


example <- hmmgc(y, initial, emission, mu, pdf.fn="hmmgc_example.pdf")

# view detected GC rich intervals
example$GCrichintervals

# Compare P(y) as computed by the forward and backward algorithms to check for correctness.
# The two quantities may differ slightly because of numerical error but should be close.
example$F.final
exampleB.final
