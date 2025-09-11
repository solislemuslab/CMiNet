# Bayesian estimation function(based on smart PCA)
#SimMap: a function for estimating the precision matrix using the similarities between features
#Recall: Baysian estimation of precision matrix
#Suppose Y = (y_1,...,y_n)^T ~ N(0,thete^-1)   Multivariate normal
#A = YY^T ~ W( n , thete^-1 ) , E(A|thete) = n*thete^-1
#Prior dist : theta ~ W(nue, G )  ; nue > p-1 , G is a p*p positive-definite matrix
# G = (nue*omega)^-1 ,  omega is a p*p positive-determined matrix  
# E(theta) = omega ^ -1  ; a predefined structure for omega gives the "structural information" about the mean of theta. 
#Posterior dist : theta ~ W(nue.star , (nue.star*omega.star)^-1)
# nue.star = n+nue  ; omega.star = (n/(n+nue))*S + (nue/(nue+n)*omega) ; S is covariance matrix of samples
#Baysian Estamtion for theta : MAP of posterior dist
# MAP = (nue.star - p - 1)*(nue.star*omega.star)
SimMap <- function(Y,nue,epsilon1,epsilon2,C) { 
  
  p<-ncol(Y)   #number of features
  n<-nrow(Y)   #number of samples 
  #nue<- (2*p)+1    #nue > p-1 (or 2p)
  
  #Kernel function
  #D2<- distance            #distance matrix includes the distance between functions

  
  #predefining structure for omega according to smartPCA 
  S1<-cov(Y)               #S1 is covariance matrix of samples
  V <- sqrt(diag(S1))      # Vector of standard deviations
  V_diag <- diag(V)        # Convert V into a diagonal matrix
  omega.hat <- V_diag %*% C %*% V_diag
  omega.hat<-omega.hat+diag(epsilon1,p)   #add epsilon for being sure the matrix is positive-definite   
  omega.star<-((n/(n+nue))*S1)+((nue/(nue+n))*omega.hat)
  nue.star<-nue+n
  #Baysian estimation for precision matrix
  IB<-as.matrix((nue.star)*omega.star)+epsilon2   #add epsilon for being sure the matrix is positive-definite
  B <-solve(IB)
  mode.star<-(nue.star-p-1)*B
  
}
 
################################################################################
 #sparse function
 sparse_quantile <- function(Y, Isigma_star, quantile_level ) {
   p <- ncol(Y)
   
   # Extract the upper triangle of Isigma_star
   upper_tri_values <- abs(Isigma_star[upper.tri(Isigma_star)])
   
   # Determine the quantile threshold based on the specified quantile level
   ths <- quantile(upper_tri_values, quantile_level)
   
   # Apply thresholding to make Isigma sparse
   Isigma_sparse <- ifelse(abs(Isigma_star) < ths, 0, Isigma_star)
   
   # Create binary matrix for theta
   theta <- ifelse(Isigma_sparse != 0, 1, 0)
   
   # Return both the sparse precision matrix and the binary adjacency matrix
   results <- list(Isigma_sparse = Isigma_sparse, theta = theta, threshold = ths)
   return(results)
 }
     
################################################################################
#Function Evaluation
 
evaluate<-function(number_of_genes,AD,TrueNet){
          
#####Diagnostic Measures #####
  
  matches <- matrix(0,nrow(TrueNet),number_of_genes)
  for ( i in 1:nrow(TrueNet)) {
    for ( j in i:number_of_genes) {
      if (AD[i,j] == 1 && TrueNet[i,j]==1) {matches [i,j]=1} else {
        if (AD[i,j] == 0 && TrueNet[i,j]==0) {matches [i,j]=0} else {
          if (AD[i,j] == 1 && TrueNet[i,j]==0) {matches [i,j]=2} else {
            matches [i,j]=3}}}
    }}
  
  matches_UP <-  matches[upper.tri(matches,diag=FALSE)]
  
  TP <- length (which(matches_UP == 1))
  TN <- length (which(matches_UP == 0))
  FP <- length (which(matches_UP == 2))
  FN <- length (which(matches_UP == 3))
  
  Diagnostic_measures_undir <- cbind(TP,TN,FP,FN)
  
  ###################
  Recall      <- TP / (TP+FN)
  Precision   <- TP / (TP+FP)
  Specificity <- TN / (TN+FP)
  Accuracy    <- (TP+TN) / (TP+FP+TN+FN)
  F_score     <- 2*(Recall*Precision)/(Recall+Precision)
  
  Accuracy_measures_undir <- cbind(Recall,Specificity,Precision,Accuracy,F_score)
  Diag_results <- list(Diagnostic_measures_undir,Accuracy_measures_undir)
  
  return(Diag_results)
}
       
################################################################################
##evaluate 
#function to set TFs target in estimated precision matrix
  fly_evaluate <- function(precision){
    evalute_matrix <- c()
    for ( i in goldTF[,2]){
      evalute_matrix <- rbind(evalute_matrix , precision[i,] )
    
  }
  
  rownames(evalute_matrix) <- goldTF[,1]
  return(evalute_matrix)
}

  ################################################################################  
  
# Function to perform bootstrapping with truncated normal replacement
  generate_bootstrap_data_with_truncnorm <- function(data, per, n_bootstrap) {
    bootstrap_data <- list()
    
    for (i in 1:n_bootstrap) {
      modified_data <- data
      for (col in 1:ncol(data)) {
        col_values <- data[, col]
        n_replace <- floor(per * length(col_values))
        replace_indices <- sample(1:length(col_values), n_replace, replace = FALSE)
        col_mean <- mean(col_values)
        col_var <- var(col_values)
        col_sd <- sqrt(col_var)
        
        # Generate truncated normal values
        trunc_values <- rtruncnorm(n = n_replace, a = min(col_values), b = max(col_values), mean = col_mean, sd = col_sd)
        
        # Replace selected values
        modified_data[replace_indices, col] <- col_values[replace_indices] + trunc_values
      }
      bootstrap_data[[i]] <- modified_data
    }
    
    return(bootstrap_data)
  }
  # Function to perform bootstrapping and generate bootstrap datasets
  generate_bootstrap_data <- function(data, n_bootstrap) {
    bootstrap_data <- list()
    for (i in 1:n_bootstrap) {
      # Sample rows with replacement
      sampled_indices <- sample(1:nrow(data), size = nrow(data), replace = TRUE)
      bootstrap_data[[i]] <- data[sampled_indices, , drop = FALSE]
    }
    return(bootstrap_data)
  }         
         