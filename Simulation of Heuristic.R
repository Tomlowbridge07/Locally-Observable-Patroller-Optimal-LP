source("Heuristics.R")
library(stats)

SimulateVStateEvolution<-function(NodeToEvolve,vVec,LambdaVec,bVec)
{
  NewvVec=vVec
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
      
      print(paste("Starting at ",toString(StartNode)))
      
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
      vVec=SimulateVStateEvolution(MoveToNode,OldvVec,LambdaVec,bVec)
      
      print(paste("Moved to ",toString(MoveToNode),"Costing ",toString(RunCost[run])))
    }
    
  }
  AverageCost=sum(RunCost)/NumberOfRunSteps
  print(paste("Average cost is ",toString(AverageCost)))
}

















