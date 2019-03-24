#Functions for performing Direct correlation method 

#build corpcor
buildCorpcor <- function(data,lam,write_output=1){
  corpcor_pcor=NULL
  if(lam!=0){
    corpcor_pcor=pcor.shrink(data,lambda=lam)
  }else{
    corpcor_pcor=pcor.shrink(data)
  }
  #space strips the colnames off, so replace those
  colnames(corpcor_pcor) <- colnames(data)
  rownames(corpcor_pcor) <- rownames(data)
  Result_Corpcor<<-corpcor_pcor
  if(write_output==1){
    Result_Corpcor_uppertri=Result_Corpcor[upper.tri(Result_Corpcor)]
    write_file(Result_Corpcor,paste0(dir_direct_bootstrap,"Corpcor_Result.csv"))
    write_file(Result_Corpcor_uppertri,paste0(dir_direct_bootstrap_uppertri,"Corpcor_Result.csv"))
  }
  return(Result_Corpcor)
}

#build space
buildSpace <- function(data,l1,l2,w,iter,write_output=1){
  
  result <- space.joint(data, lam1=l1, lam2=l2,weight=w, iter=iter)$ParCor 
  
  #space strips the colnames off, so replace those
  colnames(result) <- colnames(data)
  rownames(result) <- rownames(data)
  Result_Space<<- as.matrix(result)
  if(write_output==1){
    Result_Space_uppertri=Result_Space[upper.tri(Result_Space)]
    write_file(Result_Space,paste0(dir_direct_bootstrap,"Space_Result.csv"))
    write_file(Result_Space_uppertri,paste0(dir_direct_bootstrap_uppertri,"Space_Result.csv"))
  }
  return(Result_Space)
}

#build CorND
buildCorND <- function (data,corND_alpha,corND_beta,write_output=1){
  
  # Compute Spearman's correlation coefficient
  mat <- cor(data, method="spearman")
  
  # Set input paraemters
  alpha <- corND_alpha
  beta <- corND_beta
  
  
  if (min(mat) != max(mat)){
    # Preprocessing the input matrix - Linearly mapping the input matrix to be between -1 and 1
    abs_mat <- abs(mat)
    scale_mat <- (abs_mat-min(abs_mat))/(max(abs_mat)-min(abs_mat))
    # Thresholding the input matrix
    y <- quantile(scale_mat,1-alpha)
    mat_th <- scale_mat*(scale_mat>=y)
    mat_th <- sign(mat) * mat_th
  }else{
    print("The input matrix is a constant matrix")
  }
  
  # Making the matrix symetric if already not
  mat_th <- (mat_th+t(mat_th))/2
  
  # Eigen decomposition
  r <- eigen(mat_th)
  D <- diag(r$values)
  U <- r$vectors
  
  lam_n <- abs(min(D,0))
  lam_p <- abs(max(D,0))
  m1 <- lam_p*(1-beta)/beta
  m2 <- lam_n*(1+beta)/beta
  m <- max(m1,m2)
  
  # Network deconvolution
  for(i in 1:ncol(D)){
    D[i,i] <- (D[i,i]/(m+D[i,i]))
  }
  mat_new1 <- U%*%D%*%solve(U)
  
  
  # Linearly mapping the deconvolved matrix to be between -1 and 1
  result <- norm_mat(mat_new1)
  
  colnames(result)=colnames(data)
  rownames(result)=colnames(data)
  Result_CorND<<-data.matrix(result)
  if(write_output==1){
    Result_CorND_uppertri=Result_CorND[upper.tri(Result_CorND)]
    write_file(Result_CorND,paste0(dir_direct_bootstrap,"CorND_Result.csv"))
    write_file(Result_CorND_uppertri,paste0(dir_direct_bootstrap_uppertri,"CorND_Result.csv"))
  }
  return(Result_CorND)
}

#build IDA
buildIDA <- function(data, alpha,write_output=1){
  # Infer causal structure
  suffStat <- list(C=cor(data), n=nrow(data))
  pc.fit <- pc_stable(suffStat, indepTest=gaussCItest, p=ncol(data), alpha=alpha, skel.method = c("stable"))
  
  #Get a set of possible causal effects and infer a bound on causal effects based on obervation data
  n <- ncol(data)
  cov_data <- cov(data)
  IDA_result <- matrix(0,n,n)
  
  for (i in 1:n) {
    # infer causal effects
    ce <- idaFast(i, 1:n, cov_data, pc.fit@graph)
    
    # Select minimum of absolute values among possible causal effects.
    ce_min <- matrix(0,1,n)
    for (j in 1:n){
      index <- which(abs(ce) == min(abs(ce[j,])), arr.ind = TRUE)
      ce_min[,j] <- ce[j,index[1,2]]
    }
    IDA_result[,i] <- ce_min
  }
  # Linearly mapping the IDA result to be between -1 and 1
  result <- norm_mat(IDA_result)
  
  colnames(result)=colnames(data)
  rownames(result)=colnames(data)
  Result_IDA <- data.matrix(result)
  if(write_output==1){
    Result_IDA_uppertri=Result_IDA[upper.tri(Result_IDA)]
    write_file(Result_IDA,paste0(dir_direct_bootstrap,"IDA_Result.csv"))
    write_file(Result_IDA_uppertri,paste0(dir_direct_bootstrap_uppertri,"IDA_Result.csv"))
  }
  return(result)
}

#build Silencing 
#build space
buildSilencing <- function(data,Norm,Maximum,Minimum,write_output=1){
  
  G <- cor(data, method="spearman")
  
  # Step 1 - Preprocessing 
  #print("Preprocessing. Checking singularity of G")
  N <- nrow(G)
  R <- qr(G)$rank
  Go <- G - diag(1,N)
  #R <- N 
  if (R<N) {# G matrix is singular ==> Renormalizing off-diagonal terms
    #print("G matrix is singular - renormalizing off-diagonal terms")
    Go <- Go * Norm
  }
  
  
  # Step 2 - Silencing
  
  MaxValue <- 2
  n <- 1
  alpha <- 1 # alpha = the final normalizing factor Go ---> alpha * Go
  
  while(MaxValue <= Minimum || MaxValue >= Maximum){
    G1 = Go + diag(1,N)
    
    D <- Go %*% G1
    i <- 1
    while(i >=1 && i<=N){
      j <- i+1
      while(j >=i+1 && j<=N){
        D[i,j] = 0
        D[j,i] = 0
        j=j+1
      }
      i=i+1
    }
    
    MaxDiag <- max(max(abs(D)))
    #print("Maximum of Diagonal is Dmax = ")
    #print(MaxDiag)
    S <- (Go +D) %*% solve(G1, tol=1e-25)
    #print("Obtaining spectrum of S.")
    
    LS <- eigen(S)$values
    LS <- diag(LS)
    MaxValue <- max(abs(LS))
    #print("Maximum eigen value of S is Ls = ") 
    #print(MaxValue)
    
    if (MaxValue > Minimum && MaxValue < Maximum){
      print("Silencing done. Off-diagonal terms of G were renormalized by alpha = ") 
      print(alpha)
    } else{
      # We renormalize off-diagonal terms. The greater is the difference
      # between MaxValue and the desired range, the smaller is Norm
      Norm = 1 / (sqrt(MaxValue / (0.5 * (Maximum + Minimum))))
      #print("Renormalizing G for next iteration. Norm = ")
      #print(Norm)
      Go = Go * Norm;
      alpha = alpha * Norm;
    }
    
    n <- n+1
  }
  # Step 3 - POSTPROCESSING
  # Linearly mapping the IDA result to be between -1 and 1
  result <- norm_mat(S)
  colnames(result)=colnames(data)
  rownames(result)=colnames(data)
  Result_Silencing <- data.matrix(result)
  if(write_output==1){
    Result_Silencing_uppertri=Result_Silencing[upper.tri(Result_Silencing)]
    write_file(Result_Silencing,paste0(dir_direct_bootstrap,"Silencing_Result.csv"))
    write_file(Result_Silencing_uppertri,paste0(dir_direct_bootstrap_uppertri,"Silencing_Result.csv"))
  }
  return(result)
}

#Bootstrapping of Direct Correlation Results
#corpcor
corpcor_bootstrap<-function (data,bootstrap_type,iterations,sample.percentage,topk,params){
  boot<- bootstrap(data,"buildCorpcor",bootstrap_type, sample.percentage, iterations,topk,params)
  Result <- uptri2mat(colnames(data), boot)
  Result_uppertri=Result[upper.tri(Result)]
  write_file(Result,paste0(dir_direct_bootstrap,"Corpcor_",bootstrap_type,"_Result.csv"))
  write_file(Result_uppertri,paste0(dir_direct_bootstrap_uppertri,"Corpcor_",bootstrap_type,"_Result.csv"))
  return(Result)
}

#space
space_bootstrap<-function (data,bootstrap_type,iterations,sample.percentage,topk,params){
  boot<- bootstrap(data,"buildSpace" ,bootstrap_type, sample.percentage, iterations,topk,params)
  Result <- uptri2mat(colnames(data), boot)
  Result_uppertri=Result[upper.tri(Result)]
  write_file(Result,paste0(dir_direct_bootstrap,"Space_",bootstrap_type,"_Result.csv"))
  write_file(Result_uppertri,paste0(dir_direct_bootstrap_uppertri,"Space_",bootstrap_type,"_Result.csv"))
  return(Result)
}
#corND
corND_bootstrap<-function (data,bootstrap_type,iterations,sample.percentage,topk,params){
  boot<- bootstrap(data,"buildCorND" ,bootstrap_type, sample.percentage, iterations,topk,params)
  Result <- uptri2mat(colnames(data), boot)
  Result_uppertri=Result[upper.tri(Result)]
  write_file(Result,paste0(dir_direct_bootstrap,"CorND_",bootstrap_type,"_Result.csv"))
  write_file(Result_uppertri,paste0(dir_direct_bootstrap_uppertri,"CorND_",bootstrap_type,"_Result.csv"))
  return(Result)
}
#IDA
IDA_bootstrap<-function (data,bootstrap_type,iterations,sample.percentage,topk,params){
  boot<- bootstrap(data,"buildIDA" ,bootstrap_type, sample.percentage, iterations,topk,params)
  Result <- uptri2mat(colnames(data), boot)
  Result_uppertri=Result[upper.tri(Result)]
  write_file(Result,paste0(dir_direct_bootstrap,"IDA_",bootstrap_type,"_Result.csv"))
  write_file(Result_uppertri,paste0(dir_direct_bootstrap_uppertri,"IDA_",bootstrap_type,"_Result.csv"))
  return(Result)
}
#Silencing
Silencing_bootstrap<-function (data,bootstrap_type,iterations,sample.percentage,topk,params){
  boot<- bootstrap(data,"buildSilencing" ,bootstrap_type, sample.percentage, iterations,topk,params)
  Result <- uptri2mat(colnames(data), boot)
  Result_uppertri=Result[upper.tri(Result)]
  write_file(Result,paste0(dir_direct_bootstrap,"Silencing_",bootstrap_type,"_Result.csv"))
  write_file(Result_uppertri,paste0(dir_direct_bootstrap_uppertri,"Silencing_",bootstrap_type,"_Result.csv"))
  return(Result)
}

