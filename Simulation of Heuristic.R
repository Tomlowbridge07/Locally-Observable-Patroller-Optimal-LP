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
    print(paste("On run ",toString(run)))
    if(run==0)
    {
      #Set up initial States
      StartNode=StartingNodeHeuristic(n,IndexForNodeFunction,CostVec,LambdaVec,bVec,xVec,vMaxVec)
      
      print(paste("Starting at ",toString(StartNode)))
      
      sVec=BVec+1
      sVec[StartNode]=1
      print("Current S is ")
      print(sVec)
      
      vVec=vector(length=n)
      for(i in 1:n)
      {
        if(i==StartNode)
        {
         vVec[i]=TruncPoissionMean(LambdaVec[i],bVec[i])
         print(paste("I have just set the V to",toString(vVec[i]),"Using lam=",toString(LambdaVec[i])," and b=",toString(bVec[i])))
        }
        else
        {
          vVec[i]=0
        }

      }
      print("Current v is")
      print(vVec)
    }
    else
    {
      #Perform a run
      OldsVec=sVec
      OldvVec=vVec
      print("Current S is ")
      print(OldsVec)
      print("Current V is ")
      print(OldvVec)
      
      
      #decide where to move to using heuristic
      MoveToNode=HeuristicFunction(HeuristicDepth,n,AdjacencyMatrix,IndexForNodeFunction,OldsVec,OldvVec,CostVec,LambdaVec,bVec,xVec,vMaxVec)
      print(paste("Choosing to move to",toString(MoveToNode)))
      
      #Calculate cost of doing such an action
      RunCost[run]=CostOfAction(c(OldsVec,OldvVec),MoveToNode,n,CostVec,xVec,LambdaVec)
      print(paste("Cost of action is",toString(RunCost[run])))
      
      #Evolve System
      sVec=NewSState(OldsVec,MoveToNode,BVec)
      print("Evolved S is ")
      print(sVec)
      vVec=SimulateVStateEvolution(MoveToNode,sVec,OldvVec,xVec,LambdaVec,bVec)
      #vVec=NewMeanVState(OldvVec,sVec,MoveToNode,BVec,bVec,LambdaVec)
      print("Evolved v is")
      print(vVec)
      
      print(paste("Moved to ",toString(MoveToNode),"Costing ",toString(RunCost[run])))
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














