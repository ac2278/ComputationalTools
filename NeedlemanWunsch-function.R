# This function reports the number of sequences found in the FASTA file and imports the FASTA sequences into R:
# Specify the name of the FASTA file [character], FASTAfilename. 
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

# This function implements the Needleman-Wunsch global alignment algorithm. 
# The function returns a list of two elements: 
#   (1) best.globalscore 
#   (2) best.globalalignment
# Specify a value for the gap-open penalty, d [numeric]. 
# Specify a for the gap-extension penalty, e [numeric]. 
# The object sequences must be a list of length 2 [character].
# The object scorematrix must be a matrix consisting of four named rows and four named columns [numeric].
align = function(sequences, d, e, scorematrix){
  if(!length(sequences)==2){ print("Error: The object sequences must be a list of length 2."); break } else
    if(e==0){ print("Implement the Needleman-Wunsch global alignment algorithm -- linear gap penalty.");  
      x = sequences[[1]]; n = length(x); names(x) = 1:n
      y = sequences[[2]]; m = length(y); names(y) = 1:m
      
      # Define the F and POINTER matrix; matrix F is a lookup table that stores the 
      # results of solving sub-problems: 
      F = matrix(nrow=m+1, ncol=n+1, dimnames=list(0:m, 0:n))
      POINTER = matrix(nrow=m+1, ncol=n+1, dimnames=list(0:m, 0:n))  
      
      # Define the boundary conditions: 
      F["0",] = as.numeric(colnames(F))*(-d)
      F[,"0"] = as.numeric(rownames(F))*(-d)
      
      # Precompute the score function for matches and mismatches; 
      # store values in matrix S: 
      S = matrix(nrow=m, ncol=n)
      for (i in 1:n){
        S[,i] = scorematrix[x[as.character(i)], y]  
      }
      mode(S) = "numeric"
      S = cbind(NA, S); S = rbind(NA, S)
      dimnames(S) = dimnames(F)
      
      # Compute the elements of the F matrix (from left to right and top to bottom), 
      # using the recurrence relation given in equation 2.8 of Durbin et al. Each element 
      # of the matrix has THREE precursors. 
      # Maintain back pointers (required to implement traceback method for identifying 
      # the best alignment): 
      for (i in 1:n){
        for (j in 1:m){
          scores = c()
          scores = c(scores, diagonal = F[as.character(j-1), as.character(i-1)] + S[as.character(j), as.character(i)])
          scores = c(scores, vertical = F[as.character(j-1), as.character(i)] - d)
          scores = c(scores, horizontal = F[as.character(j), as.character(i-1)] - d)
          maxscore = sort(scores, dec=T)[1]
          F[as.character(j), as.character(i)] = maxscore
          POINTER[as.character(j), as.character(i)] = names(maxscore)
        }
      }
      
      POINTER[,1] = "vertical"
      POINTER[1,] = "horizontal"
      POINTER[1,1] = NA
      
      # The value in the final cell of the F matrix is by definition the best score for an alignment 
      # of x to y (i.e. the score of the best GLOBAL alignment of x to y). 
      best.score = F[as.character(m), as.character(n)]
      
      # Implement traceback procedure to find the best global alignment: 
      # Build the alignment in reverse, starting from the final cell of the POINTER matrix. 
      coordinate = c(m, n) # current cell
      x.alignment = c()
      y.alignment = c()
      
      while(sum(coordinate)>=2){
        if(POINTER[as.character(coordinate[1]), as.character(coordinate[2])] 
           == "diagonal"){ x.alignment = c(x.alignment, x[as.character(coordinate[2])]);
           y.alignment = c(y.alignment, y[as.character(coordinate[1])]);
           coordinate = c(coordinate[1]-1, coordinate[2]-1) } else
             
             if(POINTER[as.character(coordinate[1]), as.character(coordinate[2])] 
                == "horizontal"){ x.alignment = c(x.alignment, x[as.character(coordinate[2])]);
                y.alignment = c(y.alignment, "-");
                coordinate = c(coordinate[1], coordinate[2]-1) } else
                  
                  if(POINTER[as.character(coordinate[1]), as.character(coordinate[2])] 
                     == "vertical"){ x.alignment = c(x.alignment, "-");
                     y.alignment = c(y.alignment, y[as.character(coordinate[1])]);
                     coordinate = c(coordinate[1]-1, coordinate[2]) }
      }
      
      best.alignment = rbind(rev(x.alignment), rev(y.alignment))
      dimnames(best.alignment) = list(names(sequences), 1:ncol(best.alignment))
      return(list("best.globalscore" = best.score, "best.globalalignment" = best.alignment))
    } else print("Implement the Needleman-Wunsch global alignment algorithm -- affine gap penalty.");  
  x = sequences[[1]]; n = length(x); names(x) = 1:n
  y = sequences[[2]]; m = length(y); names(y) = 1:m
  
  # Define the M, Ix, Iy, and POINTER matrix: 
  M = Ix = Iy = POINTER = matrix(nrow=m+1, ncol=n+1, dimnames=list(0:m, 0:n))
  
  # Precompute the score function for matches and mismatches; 
  # store values in matrix S: 
  S = matrix(nrow=m, ncol=n)
  for (i in 1:n){
    S[,i] = scorematrix[x[as.character(i)], y]  
  }
  mode(S) = "numeric"
  S = cbind(NA, S); S = rbind(NA, S)
  dimnames(S) = dimnames(POINTER)
  
  # Define the boundary conditions: 
  M["0",] = -Inf
  M[,"0"] = -Inf
  M["0","0"] = 0
  
  Ix["0",] = -d-(as.numeric(colnames(Ix))-1)*e
  Ix[,"0"] = -Inf
  
  Iy[,"0"] = -d-(as.numeric(rownames(Iy))-1)*e
  Iy["0",] = -Inf
  
  #######################################
  #######################################
  # Compute the elements of the three matrices, 
  # using the recurrence relation given in equation 2.16 of Durbin et al. 
  # Maintain back pointers (required to implement traceback method for identifying 
  # the best alignment): 
  #######################################
  #######################################
  for (i in 1:n){
    for (j in 1:m){
      M[as.character(j),as.character(i)] = 
        sort(c(diagonal = M[as.character(j-1),as.character(i-1)] + S[as.character(j),as.character(i)], 
               horizontal = Ix[as.character(j-1),as.character(i-1)] + S[as.character(j),as.character(i)], 
               vertical = Iy[as.character(j-1),as.character(i-1)] + S[as.character(j),as.character(i)]), dec=T)[1] 
      
      Ix[as.character(j),as.character(i)] = 
        sort(c(diagonal = M[as.character(j),as.character(i-1)] - d,
               horizontal = Ix[as.character(j),as.character(i-1)] - e), dec=T)[1]
      
      Iy[as.character(j),as.character(i)] = 
        sort(c(diagonal = M[as.character(j-1),as.character(i)] - d,   
               vertical = Iy[as.character(j-1),as.character(i)] - e), dec=T)[1]
      
      argmax = sort(c(diagonal = M[as.character(j),as.character(i)], horizontal = Ix[as.character(j),as.character(i)], vertical = Iy[as.character(j),as.character(i)]), dec=T)[1]
      POINTER[as.character(j),as.character(i)] = names(argmax)
      
    }
  }
  
  POINTER[,1] = "vertical"
  POINTER[1,] = "horizontal"
  POINTER[1,1] = NA
  # Find the score of the best GLOBAL alignment of x to y. 
  best.score = sort(c(Mmatrix = M[as.character(m), as.character(n)], 
                      Ixmatrix = Ix[as.character(m), as.character(n)],
                      Iymatrix = Iy[as.character(m), as.character(n)]), dec=T)[1]
  
  # Implement traceback procedure to find the best global alignment: 
  # Build the alignment in reverse, starting from the final cell of the POINTER matrix. 
  coordinate = c(m, n) # current cell
  x.alignment = c()
  y.alignment = c()
  
  while(sum(coordinate)>=2){
    if(POINTER[as.character(coordinate[1]), as.character(coordinate[2])] 
       == "diagonal"){ x.alignment = c(x.alignment, x[as.character(coordinate[2])]);
       y.alignment = c(y.alignment, y[as.character(coordinate[1])]);
       coordinate = c(coordinate[1]-1, coordinate[2]-1) } else
         
         if(POINTER[as.character(coordinate[1]), as.character(coordinate[2])] 
            == "horizontal"){ x.alignment = c(x.alignment, x[as.character(coordinate[2])]);
            y.alignment = c(y.alignment, "-");
            coordinate = c(coordinate[1], coordinate[2]-1) } else
              
              if(POINTER[as.character(coordinate[1]), as.character(coordinate[2])] 
                 == "vertical"){ x.alignment = c(x.alignment, "-");
                 y.alignment = c(y.alignment, y[as.character(coordinate[1])]);
                 coordinate = c(coordinate[1]-1, coordinate[2]) }
  }
  
  best.alignment = rbind(rev(x.alignment), rev(y.alignment))
  dimnames(best.alignment) = list(names(sequences), 1:ncol(best.alignment))
  return(list("best.globalscore" = best.score, "best.globalalignment" = best.alignment))
  
} 

###########
# optional example
sequences = importFASTA("example.fasta")         # read the example FASTA sequences into R

scorematrix <- matrix(c(91, -114, -31, -123,     # define the score matrix; note that the score matrix defined here is symmetric
                        -114, 100, -125, -31, 
                        -31, -125, 100, -114, 
                        -123, -31, -114, 91), 
                      byrow=T, nrow=4, ncol=4,
                      dimnames=list(c("A", "C", "G", "T"), c("A", "C", "G", "T")))

# linear penalty alignment
alignment.linear <- align(sequences, d=100, e=0, scorematrix)

# gap penalty alignmnet
alignment.affine <- align(sequences, d=430, e=30, scorematrix)

# NOTE 1: The alignments are different. A natrual question is 'which alignment is better?'. The answer
#         to this question will vary depending on how one measures 'better' (i.e. 'better' with respect to what?)
#         The linear penalty alignment is 'better' than the gap penalty alignmnet 
#         with respect to alignment scores. The affine gap penalty alignment (alignment.affine) 
#         has an alignment score (3077) lower than that of the linear penalty alignment (4225). 
#         One could argue that comparing the alignment scores of these two alignments
#         is not an adequate approach to selecting the 'better' alignment. 
#         The linear penalty alignment assumes the simplest gap model, in which the gap score 
#         is a simple multiple of the length. This type of scoring scheme is not ideal for 
#         biological sequences since it penalizes additional gap steps as much as the first. 
#         In reality, when gaps do occur, they are often longer than one residue. 
#         One could therefore reason that the affine gap penalty alignment is the better alignment
#         since it captures aspects of biological reality.

# NOTE 2: If one has access to a dataset with alignments that are “correct”,
#         one could decide on the “right” score matrix and gap penalties. One could treat the score 
#         matrix and gap penalties as unknown parameters and use the data to estimate parameter values 
#         (i.e. model training).
