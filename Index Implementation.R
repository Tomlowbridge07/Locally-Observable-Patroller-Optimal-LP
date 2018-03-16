source("Optimal solution by LP.R")

FindVMax<-function(Lambda,b,x)
{
  #First calculate B and the remainder
  B=ceiling(x)
  R=B-x
  
  CurrentMax=-1
  for(v in 0:(b+1))
  {
    if(v <= (Lambda * (1-R)))
    {
      CurrentMax=v
    }
    else
    {
      return(CurrentMax)
    }
  }
  return(CurrentMax)
}

CreateVMaxVector<-function(n,LambdaVec,bVec,xVec)
{
  VMaxVec=vector(length=n)
  for(i in 1:n)
  {
    VMaxVec[i]=FindVMax(LambdaVec[i],bVec[i],xVec[i])
  }
  return(VMaxVec)
}

Delta<-function(tilde=FALSE,CostAtNode,Lambda,b,x,v,vMax)
{
  #First calculate B and the remainder
  B=ceiling(x)
  R=B-x
  
  if(tilde==FALSE)
  {
    #Calculate Sum
    Sum=Lambda * R * (B+1) + v * (B+1-TruncPoissonHazard(Lambda,b,v))
    if(v>0)
    {
    for(i in 0:(v-1))
    { 
      Sum=Sum- i * TruncPoissonPMF(Lambda,b,i)
    }
    }
    Sum=CostAtNode * Sum
    
    return (Sum)
  }
  if(tilde==TRUE)
  {
    #Now calculate sum
    Sum= Lambda * (B+1 - R * TruncPoissonHazard(Lambda,b,vMax))
    if(vMax>0)
    {
    for(i in 0:(vMax-1))
    {
      Sum=Sum- i * TruncPoissonPMF(Lambda,b,i)
    }
    }
    Sum=CostAtNode * Sum
    return(Sum)
  }
}

PlainIndexForNode<-function(s,v,Cost,Lambda,b,x,vMax)
{
  #First calculate B 
  B=ceiling(x)
  
  if(s < B)
  {
    return(0)
  }
  else if(s==B && v < vMax)
  {
    return(Delta(tilde = FALSE,Cost,Lambda,b,x,v,vMax))
  }
  else
  {
    return(Delta(tilde = TRUE,Cost,Lambda,b,x,v,vMax))
  }
}

EqualStepIndexForNode<-function(s,v,Cost,Lambda,b,x,vMax)
{
  #First calculate B 
  B=ceiling(x)
  
  
  if(s <= B && v < vMax)
  {
    
    return(Delta(tilde = FALSE,Cost,Lambda,b,x,v,vMax)*(s/B))
  }
  else if(s <= B && v >=vMax)
  {
    
    return(Delta(tilde = TRUE,Cost,Lambda,b,x,v,vMax)*(s/B))
  }
  else
  {
   
    return(Delta(tilde = TRUE,Cost,Lambda,b,x,v,vMax))
  }
}

IncreasingStepIndexForNode<-function(s,v,Cost,Lambda,b,x,vMax)
{
  #First calculate B 
  B=ceiling(x)
  
  if(s <= B && v < vMax)
  {
    return(Delta(tilde = FALSE,Cost,Lambda,b,x,v,vMax)/(B+1-s))
  }
  else if(s <= B && v >=vMax)
  {
    return(Delta(tilde = TRUE,Cost,Lambda,b,x,v,vMax)/(B+1-s))
  }
  else
  {
    return(Delta(tilde = TRUE,Cost,Lambda,b,x,v,vMax))
  }
}

RandomIndexForNode<-function(s,v,Cost,Lambda,b,x,vMax)
{
  return(runif(1,min=0,max=1))
}

WaitingIndexForNode<-function(s,v,Cost,Lambda,b,x,vMax)
{
  if(s==1)
  {
    return(1)
  }
  else
  {
    return(0)
  }
}

IndicesForNodes<-function(n,IndexForNodeFunction,sVec,vVec,CostVec,LambdaVec,bVec,xVec,vMaxVec=NULL)
{
  if(is.null(vMaxVec))
  {
    #Create vMaxVec
    vMaxVec=CreateVMaxVector(n,LambdaVec,bVec,xVec)
  }
  
  
  #Now for each node we find its index
  Indices=vector(length=n)
  for(i in 1:n)
  { 
    Indices[i]=IndexForNodeFunction(sVec[i],vVec[i],CostVec[i],LambdaVec[i],bVec[i],xVec[i],vMaxVec[i])
  }
  
  return(Indices)
}