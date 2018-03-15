source("Heuristics.R")
library(stats)
library(ggplot2)
library(reshape2)

SimulateVStateEvolution<-function(NodeToEvolve,NewsVec,vVec,xVec,LambdaVec,bVec)
{
  BVec=ceiling(xVec)
  NewvVec=(NewVState(vVec,NewsVec,NodeToEvolve,BVec,bVec,LambdaVec)$State)[1,]
  NewvVec[NodeToEvolve]=min(bVec[NodeToEvolve],rpois(1,LambdaVec[NodeToEvolve]))
  return(NewvVec)
}

SimulationForEvolution<-function(NumberOfRunSteps,HeuristicDepth,HeuristicFunction,n,AdjacencyMatrix,IndexForNodeFunction,CostVec,LambdaVec,bVec,xVec,vMaxVec=NULL)
{
  BVec=ceiling(xVec)
  if(is.null(vMaxVec))
  {
    #Create vMaxVec
    vMaxVec=CreateVMaxVector(n,LambdaVec,bVec,xVec)
  }
  
  RunCost=vector(length=NumberOfRunSteps)
  
  for(run in 0:NumberOfRunSteps)
  {
    if(run==0)
    {
      #Set up initial States
      StartNode=StartingNodeHeuristic(n,IndexForNodeFunction,CostVec,LambdaVec,bVec,xVec,vMaxVec=NULL)
      
      #print(paste("Starting at ",toString(StartNode)))
      
      sVec=BVec+1
      sVec[StartNode]=1
      
      vVec=vector(length=n)
      for(i in 1:n)
      {
        vVec[i]=TruncPoissionMean(LambdaVec[i],bVec[i])
      }
    }
    else
    {
      #Perform a run
      OldsVec=sVec
      OldvVec=vVec
      
      #decide where to move to using heuristic
      MoveToNode=HeuristicFunction(HeuristicDepth,n,AdjacencyMatrix,IndexForNodeFunction,OldsVec,OldvVec,CostVec,LambdaVec,bVec,xVec,vMaxVec)
      #print(MoveToNode)
      
      #Calculate cost of doing such an action
      RunCost[run]=CostOfAction(c(OldsVec,OldvVec),MoveToNode,n,CostVec,xVec,LambdaVec)
      
      #Evolve System
      sVec=NewSState(OldsVec,MoveToNode,BVec)
      vVec=SimulateVStateEvolution(MoveToNode,sVec,OldvVec,xVec,LambdaVec,bVec)
      
      #print(paste("Moved to ",toString(MoveToNode),"Costing ",toString(RunCost[run])))
    }
    
  }
  AverageCost=sum(RunCost)/NumberOfRunSteps
  print(paste("Average cost is ",toString(AverageCost)))
  return(list(Average=AverageCost,CostForStep=RunCost))
}

SimulationExperiment<-function(NumberOfTrials,NumberOfRunSteps,HeuristicDepth,HeuristicFunction,n,AdjacencyMatrix,IndexForNodeFunction,CostVec,LambdaVec,bVec,xVec,vMaxVec=NULL)
{
 #We repeat the simulation
  RunningCostMatrix=matrix(0,nrow = 0,ncol=NumberOfRunSteps)
  AveragecostVec=vector(length=NumberOfTrials)
  
  for(trial in 1:NumberOfTrials)
  {
    Simulation=SimulationForEvolution(NumberOfRunSteps,HeuristicDepth,HeuristicFunction,n,AdjacencyMatrix,IndexForNodeFunction,CostVec,LambdaVec,bVec,xVec)
   RunningCostMatrix=rbind(RunningCostMatrix,Simulation$CostForStep)
   AveragecostVec[trial]=Simulation$Average
   #print(RunningCostMatrix)
  }
  AverageAmongSimulations=sum(AveragecostVec)/NumberOfTrials
  return(AverageAmongSimulations)
}














