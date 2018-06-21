source("Optimal solution by LP.R")

FindVMax<-function(Lambda,b,x)
{
  #First calculate B and the remainder
  B=ceiling(x)
  R=B-x
  
  CurrentMax=-1
  for(v in 0:b)
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
    if(v==-1)
    {
      return(0)
    }
    #Calculate Sum
    Sum=Lambda * R * B

    Sum=Sum + (v * (B+1-TruncPoissonHazard(Lambda,b,v)))

    if(v>0)
    {
     for(i in 0:(v-1))
     { 
       Sum=Sum- (i * TruncPoissonPMF(Lambda,b,i))
     }
    }
    Sum=CostAtNode * Sum
    return (Sum)
  }
  if(tilde==TRUE)
  {
    #Now calculate sum
    Sum= Lambda * (B + 1 - R + (R-1)*TruncPoissonHazard(Lambda,b,vMax+1))
    for(i in 0:vMax)
    {
      Sum=Sum- (i * TruncPoissonPMF(Lambda,b,i))
    }
    Sum=CostAtNode * Sum
    return(Sum)
  }
}

#This index is for use with the no local-observations
PlainIndexForDetCostNode<-function(s,v,Cost,Lambda,b,x,vMax)
{
  B=ceiling(x)
  if(v!=0)
  {
    print("Error")
  }
  
  if(s<B)
  {
    return(0)
  }
  else if(s==B)
  {
    return(Cost*Lambda*(B-x))
  }
  else if(s==B+1)
  {
    return(Cost*Lambda)
  }
}

#This index is for use with the no local-observations
EqualPenaltyIndexForDetCostNode<-function(s,v,Cost,Lambda,b,x,vMax)
{
  B=ceiling(x)
  if(v!=0)
  {
    print("Error")
  }
  
  if(s <= B)
  {
    return((Cost*Lambda*(B-x))/B)
  }
  else if(s==B+1)
  {
    return(Cost*Lambda)
  }
}

#This index is for use with the no local-observations
EqualBenefitIndexForDetCostNode<-function(s,v,Cost,Lambda,b,x,vMax)
{
  B=ceiling(x)
  if(v!=0)
  {
    print("Error")
  }
  
  if(s <= B)
  {
    return((s*Cost*Lambda*(B-x))/B)
  }
  else if(s==B+1)
  {
    return(Cost*Lambda)
  }
}

#This index is for use with the no local-observations
UnequalPenaltyForDetCostNode<-function(s,v,Cost,Lambda,b,x,vMax)
{
  B=ceiling(x)
  if(v!=0)
  {
    print("Error")
  }
  
  if(s <= B)
  {
    return((2*s*Cost*Lambda*(B-x))/(B+1))
  }
  else if(s==B+1)
  {
    return(Cost*Lambda)
  }
}

#This index is for use with the no local-observations
UnequalBenefitIndexForDetCostNode<-function(s,v,Cost,Lambda,b,x,vMax)
{
  B=ceiling(x)
  if(v!=0)
  {
    print("Error")
  }
  
  if(s <= B)
  {
    #We must work out the sum of all lesser and equal s
    Sum=0
    for(i in 1:s)
    {
      Sum=Sum+((2*i*B)/(B+1)) 
    }
    return((Sum*Cost*Lambda*(B-x))/B)
  }
  else if(s==B+1)
  {
    return(Cost*Lambda)
  }
}

#For use with LOCAL-OBSERVATIONS :)
OldPlainIndexForNode<-function(s,v,Cost,Lambda,b,x,vMax)
{
  #First calculate B 
  B=ceiling(x)
  
  
  stopifnot(s <= (B+1))
  stopifnot(v <= (b+1))
  # print(s)
  # print(v)
  # print(vMax)
  
  if(s < B)
  {
    return(0)
  }
  else if(s==B && v < (vMax+1))
  {
    return(Delta(tilde = FALSE,Cost,Lambda,b,x,v,vMax))
  }
  else if(s==B && v >= (vMax+1))
  {
    return(Delta(tilde = TRUE,Cost,Lambda,b,x,v,vMax))
  }
  else if(s==(B+1) && v < (vMax))
  {
    return(Delta(tilde = TRUE,Cost,Lambda,b,x,v+1,vMax))
  }
  else if(s==(B+1) && v >= (vMax))
  {
    return(Delta(tilde = TRUE,Cost,Lambda,b,x,v+1,vMax))
  }
  else
  {
    print("Error")
  }
}  

PlainIndexForNode<-function(s,v,Cost,Lambda,b,x,vMax)
{
  #First calculate B 
  B=ceiling(x)

  
 stopifnot(s <= (B+1))
 stopifnot(v <= (b+1))
 # print(s)
 # print(v)
 # print(vMax)
  
  if(s < B)
  {
    return(0)
  }
  else if(s==B && v < (vMax+1))
  {
    return(Delta(tilde = FALSE,Cost,Lambda,b,x,v,vMax))
  }
  else if(s==B && v >= (vMax+1))
  {
    return(Delta(tilde = FALSE,Cost,Lambda,b,x,v,vMax))
  }
  else if(s==(B+1) && v < (vMax))
  {
    return(Delta(tilde = FALSE,Cost,Lambda,b,x,v+1,vMax))
  }
  else if(s==(B+1) && v >= (vMax))
  {
    return(Delta(tilde = TRUE,Cost,Lambda,b,x,v+1,vMax))
  }
  else
  {
    print("Error")
  }
}

EqualStepIndexForNode<-function(s,v,Cost,Lambda,b,x,vMax)
{
  #First calculate B 
  B=ceiling(x)

  
  if(s <= B && v < (vMax+1))
  {
    return(Delta(tilde = FALSE,Cost,Lambda,b,x,v,vMax)*(s/B))
  }
  else if(s <= B && v >=(vMax+1))
  {
    return(Delta(tilde = FALSE,Cost,Lambda,b,x,v,vMax)*(s/B))
  }
  else if(s==(B+1) && v < vMax)
  {
    return(Delta(tilde = FALSE,Cost,Lambda,b,x,v+1,vMax))
  }
  else if(s==(B+1) && v >= vMax)
  {
    return(Delta(tilde = TRUE,Cost,Lambda,b,x,v+1,vMax))
  }
  else
  {
    print("Error")
  }
}

FlatIndexForNode<-function(s,v,Cost,Lambda,b,x,vMax)
{
  #First calculate B 
  B=ceiling(x)
  
  
  if(s <= B && v < (vMax+1))
  {
    return(Delta(tilde = FALSE,Cost,Lambda,b,x,v,vMax)*(1/B))
  }
  else if(s <= B && v >=(vMax+1))
  {
    return(Delta(tilde = FALSE,Cost,Lambda,b,x,v,vMax)*(1/B))
  }
  else if(s==(B+1) && v < vMax)
  {
    return(Delta(tilde = FALSE,Cost,Lambda,b,x,v+1,vMax))
  }
  else if(s==(B+1) && v >= vMax)
  {
    return(Delta(tilde = TRUE,Cost,Lambda,b,x,v+1,vMax))
  }
  else
  {
    print("Error")
  }
}

IncreasingStepIndexForNode<-function(s,v,Cost,Lambda,b,x,vMax)
{
  #First calculate B 
  B=ceiling(x)
  
  if(s <= B && v <= vMax)
  {
    return(Delta(tilde = FALSE,Cost,Lambda,b,x,v,vMax)/(B+1-s))
  }
  else if(s <= B && v >= (vMax+1))
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
    #print(Indices[i])
  }
  
  return(Indices)
}

# #This function creates a matrix of s on row and v on column for a particular node
# CreateIndexTable<-function(IndexForNodeFunction,CostVec,LambdaVec,bVec,xVec,vMaxVec=NULL)
# {
# 
# }