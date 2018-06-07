source("Testing Heuristics.R")

library(ggplot2)
library(reshape)
library(tikzDevice)

#This R file contains functions which take in the Test Matrix and return values
#Test Matrices come in the order MinError,BestHeuristics,AdjMatrix,AttackTimes,Capacity,Lambda,Costs,Whole Errors

#This function returns the matrix of errors for one scenario
DecodeTestMatrixErrors<-function(ScenarioNumber,TestMatrix)
{
  print(TestMatrix)
  print(TestMatrix[ScenarioNumber,8])
 return(TestMatrix[ScenarioNumber,8][[1]])
}

#This function just returns the vector of errors for a scenario
DecodeTestVectorErrors<-function(ScenarioNumber,TestMatrix)
{
  Matrix=DecodeTestMatrixErrors(ScenarioNumber,TestMatrix)
  return(Matrix[,4])
}

#This function collects the test matrix error data for all heuristics per row
DecodeErrorDataPerHeuristic<-function(TestMatrix)
{
  #First we form a matrix of errors rows=heuristicnumber,col=scenarionumber
  NumberOfScenarios=nrow(TestMatrix)
  NumberOfHeuristics=length(DecodeTestVectorErrors(1,TestMatrix))
  MatrixOfErrors=matrix(0,nrow=NumberOfHeuristics,ncol=NumberOfScenarios)
  for(ScenarioNumber in 1:NumberOfScenarios)
  {
    VectorOfErrors=DecodeTestVectorErrors(ScenarioNumber,TestMatrix)
    MatrixOfErrors[,ScenarioNumber]=VectorOfErrors
  }

 return(MatrixOfErrors)  
}

#This function takes a heuristic vector and categorizes it
CategorizeHeuristicErrorData<-function(HeuristicErrors,NumberOfCategories,MaxErrorCategory)
{
  #Work out categories
  #Note we use the number of categories+1 to allow for the maximum category to be picked
  CategoryBoundaries=seq(0,MaxErrorCategory,length.out = NumberOfCategories+1)
  MidPointOfCategories=((CategoryBoundaries[2]-CategoryBoundaries[1])/2) + CategoryBoundaries
  NumberInEachCategory=vector(length=NumberOfCategories+1)
  
  #Now find out how many are in each category
  for(CategoryNum in 1:(NumberOfCategories+1))
  {
    #Set upper and lower bounds
    if(CategoryNum==0)
    {
     LowerBound=-1
    }
    else
    {
     LowerBound=CategoryBoundaries[CategoryNum]
    }

    if(CategoryNum==(NumberOfCategories+1))
    {
     UpperBound=Inf
    }
    else
    {
     UpperBound=CategoryBoundaries[CategoryNum+1]
    }
    NumberInEachCategory[CategoryNum]=length(HeuristicErrors[HeuristicErrors>LowerBound & HeuristicErrors<=UpperBound])
  }
  return(list(CategoryMidPoints=MidPointOfCategories,NumberInEachCategory=NumberInEachCategory))
}

#Plot a single heuristic
PlotHeuristicErrorData<-function(HeuristicErrors,NumberOfCategories,MaxErrorCategory)
{
  #We run the categorize to get the data to plot
  Categorized=CategorizeHeuristicErrorData(HeuristicErrors,NumberOfCategories,MaxErrorCategory)
  XCoordinates=Categorized$CategoryMidPoints
  YCoordinates=Categorized$NumberInEachCategory
  
  DataFrame=data.frame(XCoordinates,YCoordinates)
  Plot<-ggplot(DataFrame,show.legend='True') + geom_point(aes(x = XCoordinates, y = YCoordinates)) +
    geom_line(aes(x = XCoordinates, y = YCoordinates))
  print(Plot)
  return(Plot)
}
  
#Plot a group of heuristics
#Note the heuristic erros should be provided in a matrix of rows for each heuristic
PlotMultipleHeuristicErrorData<-function(HeuristicErrors,NumberOfCategories,MaxErrorCategory,HeuristicNames=NULL,SaveTexImage=F,FileName=NULL,Size=c(3,2))
{
  #First we look at how many graphs we are going to plot
  NumHeuristics=nrow(HeuristicErrors)

  XCoordinates=matrix(0,ncol=length(CategorizeHeuristicErrorData(HeuristicErrors[1,],NumberOfCategories,MaxErrorCategory)$NumberInEachCategory),nrow=NumHeuristics)
  YCoordinates=matrix(0,ncol=length(CategorizeHeuristicErrorData(HeuristicErrors[1,],NumberOfCategories,MaxErrorCategory)$NumberInEachCategory),nrow=NumHeuristics)
  
  #Now categorize the data
  for(i in 1:NumHeuristics)
  { 
    if(i==1)
    {
      #Initialize
    Categorized=CategorizeHeuristicErrorData(HeuristicErrors[i,],NumberOfCategories,MaxErrorCategory)
    YCoordinates=matrix(0,ncol=length(Categorized$NumberInEachCategory),nrow=NumHeuristics)
    XCoordinates=Categorized$CategoryMidPoints
    YCoordinates[i,]=Categorized$NumberInEachCategory
    YCoordinates[i,]=YCoordinates[i,]/sum(YCoordinates[i,])
    Ydataframecol<-data.frame(YCoordinates[i,])
    
    if(is.null(HeuristicNames))
    {
      names(Ydataframecol)=paste("YCoordinates",toString(i))
    }
    else
    {
      names(Ydataframecol)=toString(HeuristicNames[i])
    }    
                                          
                                          
    XYCoordinates=data.frame(XCoordinates)
    XYCoordinates<-cbind(XYCoordinates,Ydataframecol)
    
    }
    else
    {
      Categorized=CategorizeHeuristicErrorData(HeuristicErrors[i,],NumberOfCategories,MaxErrorCategory)
      YCoordinates[i,]=Categorized$NumberInEachCategory
      YCoordinates[i,]=YCoordinates[i,]/sum(YCoordinates[i,])
      Ydataframecol<-data.frame(YCoordinates[i,])
      if(is.null(HeuristicNames))
      {
        names(Ydataframecol)=paste("YCoordinates",toString(i))
      }
      else
      {
        names(Ydataframecol)=toString(HeuristicNames[i])
      }
      
      XYCoordinates<-cbind(XYCoordinates,Ydataframecol)
    }
 print(XYCoordinates)
  }
  
  #Now we plot the data
  DataFrame<-XYCoordinates
  print(DataFrame)
  MeltedDataFrame<-melt(DataFrame,id="XCoordinates")
  print(MeltedDataFrame)
  
  if(SaveTexImage)
  {
    tikz(file=paste("/local/pmxtol/Dropbox/",toString(FileName),".tex",sep=""),width=Size[1],height=Size[2])
  }
  
  Plot<-ggplot(MeltedDataFrame,aes(x=XCoordinates,y=value,color=variable),show.legend='True') + #geom_point()+
       geom_line() +xlab("Percentage Error")+ylab("Frequency Density") + labs(color="Heuristic Type")
  if(SaveTexImage)
  {
    print(Plot)
    dev.off()  
  }
  
  print(Plot)
  return(Plot)
  
}








