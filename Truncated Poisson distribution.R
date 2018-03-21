#Create truncated poisson pmf/cdf/hazard
TruncPoissonPMF<-function(lambda,b,i)
{
  if(i > b)
  {
    return(0)
  }
  else if(i==0 && b==0)
  {
    return(1)
  }
  else
  {
    if(i!=b)
    {
      return(((lambda^i)/factorial(i))*exp(-lambda))
    } 
    if(i==b)
    { 
      Prob=1
      for(j in 0:(b-1))
      {
        Prob=Prob-((lambda^j)/factorial(j))*exp(-lambda)
      }
      return(Prob)
    }
  }
}

TruncPoissionCDF<-function(lambda,b,i)
{
  if(i < b)
  {
    return(0)
  }
  else if(i!=b)
  { 
    Sum=0
    for(j in 0:i)
    {
      Sum=Sum+TruncPoissonPMF(lambda,b,j)
    }
    return(Sum)
  }
  else if(i==b)
  {
    return(1)
  }  
}

TruncPoissonHazard<-function(lambda,b,i)
{
  if(i > b)
  {
    return(0)
  }
  else if(i!=b)
  {
    Sum=0
    for(j in i:b)
    {
      Sum=Sum+TruncPoissonPMF(lambda,b,j)
    }
    return(Sum)
  }
  else if(i==b)
  {
    return(TruncPoissonPMF(lambda,b,i))
  }  
}

TruncPoissionMean<-function(lambda,b)
{
  Mean=0
  for(i in 0:b)
  {
    Mean=Mean+i*TruncPoissonPMF(lambda,b,i)
  }
  return(Mean)
}