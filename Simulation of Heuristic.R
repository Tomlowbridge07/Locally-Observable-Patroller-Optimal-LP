source("Heuristics.R")
library(stats)
library(ggplot2)
library(reshape2)

#This function returns a vector for the simulation to use
CreateSimulationScenario<-function(NoSteps,n,LambdaVec)
{
  Scenario=matrix(0,nrow=n,ncol=NoSteps)
  
  for(i in 1:n)
  {
   for(j in 1:NoSteps)
   {
     Scenario[i,j]=rpois(1,LambdaVec[i])
   }
  }
  return(Scenario)
}

# SimulateVStateEvolution<-function(NodeToEvolve,NewsVec,vVec,xVec,LambdaVec,bVec)
# {
#   BVec=ceiling(xVec)
#   #NewvVec=(NewVState(vVec,NewsVec,NodeToEvolve,BVec,bVec,LambdaVec)$State)[1,]
#   NewvVec=vVec
#   NewvVec[NodeToEvolve]=min(bVec[NodeToEvolve],rpois(1,LambdaVec[NodeToEvolve]))
#   # NewvVec=NewMeanVState(vVec,NewsVec,NodeToEvolve,BVec,bVec,LambdaVec)
#   return(NewvVec)
# }

SimulationForEvolution<-function(NumberOfRunSteps,HeuristicDepth,HeuristicFunction,n,AdjacencyMatrix,IndexForNodeFunction,CostVec,LambdaVec,bVec,xVec,SimulationScenario=NULL,vMaxVec=NULL)
{
  #Set up simulation matrix and tracking
  if(is.null(SimulationScenario))
  {
    SimulationScenario=CreateSimulationScenario(NumberOfRunSteps,n,LambdaVec)
  }
  
  #This vector tracks which column to currently use in each row
  TrackingScenario=rep(1,n)
  
  
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
         #vVec[i]=TruncPoissionMean(LambdaVec[i],bVec[i])
         vVec[i]=0
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
      #print(sVec)
      
      #Creating new vVec using scenario
      vVec=OldvVec
      vVec[MoveToNode]=min(SimulationScenario[MoveToNode,TrackingScenario[MoveToNode]],bVec[MoveToNode])
      TrackingScenario[MoveToNode]=TrackingScenario[MoveToNode]+1
      #vVec=NewMeanVState(OldvVec,sVec,MoveToNode,BVec,bVec,LambdaVec)
      print("Evolved v is")
      #print(vVec)
      
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














