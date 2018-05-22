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

SimulationForEvolution<-function(NumberOfRunSteps,HeuristicDepth,HeuristicFunction,n,AdjacencyMatrix,IndexForNodeFunction,CostVec,LambdaVec,bVec,xVec,BurnOut=0,SimulationScenario=NULL,vMaxVec=NULL)
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
  SeperatedRunCost=matrix(0,nrow=NumberOfRunSteps,ncol=2)
  CumulativeSeperatedRunCost=matrix(0,nrow=NumberOfRunSteps,ncol=2)
  AllInfo=matrix(list(),nrow=NumberOfRunSteps,ncol=4) #All info stores; state, decision , cost from arrival , cost from observed
  
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
      SeperatedRunCost[run,1]=SeperatedCostOfAction(c(OldsVec,OldvVec),MoveToNode,n,CostVec,xVec,LambdaVec)$CostDueToArrivals
      SeperatedRunCost[run,2]=SeperatedCostOfAction(c(OldsVec,OldvVec),MoveToNode,n,CostVec,xVec,LambdaVec)$CostDueToObserved
      if(run==1)
      {
        CumulativeSeperatedRunCost[run,1]=SeperatedRunCost[run,1]
        CumulativeSeperatedRunCost[run,2]=SeperatedRunCost[run,2]        
      }
      else
      {
        CumulativeSeperatedRunCost[run,1]=CumulativeSeperatedRunCost[run-1,1]+SeperatedRunCost[run,1]
        CumulativeSeperatedRunCost[run,2]=CumulativeSeperatedRunCost[run-1,2]+SeperatedRunCost[run,2]
      }
      
      
      print(paste("Cost of action is",toString(RunCost[run])))
      
      #Evolve System
      sVec=NewSState(OldsVec,MoveToNode,BVec)
      #print("Evolved S is ")
      #print(sVec)
      
      #Creating new vVec using scenario
      vVec=OldvVec
      vVec[MoveToNode]=min(SimulationScenario[MoveToNode,TrackingScenario[MoveToNode]],bVec[MoveToNode])
      TrackingScenario[MoveToNode]=TrackingScenario[MoveToNode]+1
      #vVec=NewMeanVState(OldvVec,sVec,MoveToNode,BVec,bVec,LambdaVec)
      #print("Evolved v is")
      #print(vVec)
      
      AllInfo[[run,1]]=c(OldsVec,OldvVec)
      AllInfo[[run,2]]=MoveToNode
      AllInfo[[run,3]]=SeperatedRunCost[run,1]
      AllInfo[[run,4]]=SeperatedRunCost[run,2]
      
      print(paste("Moved to ",toString(MoveToNode),"Costing ",toString(RunCost[run])))
    }
    
  }
  
  #We apply a burn out to remove the first lot of results (to allow it to reach some equilibrium)
  Keep=(BurnOut+1):NumberOfRunSteps
  RunCost=RunCost[Keep]
  SeperatedRunCost=SeperatedRunCost[Keep,]
  CumulativeSeperatedRunCost=CumulativeSeperatedRunCost[Keep,]
  AllInfo=AllInfo[Keep,]
  
  AverageCost=sum(RunCost)/NumberOfRunSteps
  print(paste("Average cost is ",toString(AverageCost)))
  return(list(Average=AverageCost,CostForStep=RunCost,SeperatedCostForStep=SeperatedRunCost,CumulativeSeperatedCostForStep=CumulativeSeperatedRunCost,FullInfoMatrix=AllInfo))
}

#This function is desinged to take in the FullInfoMatrix and compare its decisions to actions from some policy
CompareSimulationInfoToPolicy<-function(FullInfo,ActionPolicy,StateSpace,BVec)
{
  AgreeAtStep=vector(length=nrow(FullInfo))
  for(i in 1:nrow(FullInfo))
  {
    #For each step we now compare the decision to the policy
    CurrentState=FullInfo[[i,1]]
    Decision=FullInfo[[i,2]]
    l=length(CurrentState)
    n=l/2
    #Note here we retain the information for the observed when in state B+1, we will now remove this
    for(j in 1:n)
    {
      if(CurrentState[j]==(BVec[j]+1))
      {
        CurrentState[j+n]=0
      }
    }
    
    CurrentStateID=IdenityRow(CurrentState,StateSpace)
    PolicyDecision=ActionPolicy[CurrentStateID]

    
    print(paste("Our Decision at state ",toString(CurrentState)))
    print(paste("Proposed by simulation is: ",toString(Decision)," and policy suggests: ",toString(PolicyDecision)))
    
    if(Decision==PolicyDecision)
    {
      AgreeAtStep[i]=1
    }
    else
    {
      AgreeAtStep[i]=0
    }
    
  }
  
  #Disagree stores which we disagree at
  DisagreeAt=which(AgreeAtStep==0)
  print("The Policys disagree at the following")
  if(length(DisagreeAt)!=0)
  {
   for(i in 1:length(DisagreeAt))
   {
     print(paste("Step: ",toString(DisagreeAt[i]),", State: ",toString(FullInfo[[DisagreeAt[i],1]])))
     print(paste("Decisions made where, Index proposes: ",toString(FullInfo[[DisagreeAt[i],2]])," but policy suggests: ",toString(ActionPolicy[IdenityRow(FullInfo[[DisagreeAt[i],1]],StateSpace)])))
   }
  }
  return(AgreeAtStep)
}

#Similar to Simulation but runs using a given polciy (note it will be of list form, we will use the first choice)
#Note the index function is used to determine the starting node of the scenario
RunPolicyForScenario<-function(NumberOfRunSteps,Policy,n,AdjacencyMatrix,IndexForNodeFunction,CostVec,LambdaVec,bVec,xVec,SimulationScenario,StateSpace=NULL,vMaxVec=NULL)
{
  #Set up simulation matrix and tracking
  if(is.null(SimulationScenario))
  {
    SimulationScenario=CreateSimulationScenario(NumberOfRunSteps,n,LambdaVec)
  }
  #This vector tracks which column to currently use in each row
  TrackingScenario=rep(1,n)
  BVec=ceiling(xVec)
  if(is.null(StateSpace))
  {
    StateSpace=CreateSVStates(n,BVec,bVec)
  }
  
  
  
  if(is.null(vMaxVec))
  {
    #Create vMaxVec
    vMaxVec=CreateVMaxVector(n,LambdaVec,bVec,xVec)
  }
  
  RunCost=vector(length=NumberOfRunSteps)
  SeperatedRunCost=matrix(0,nrow=NumberOfRunSteps,ncol=2)
  CumulativeSeperatedRunCost=matrix(0,nrow=NumberOfRunSteps,ncol=2)
  AllInfo=matrix(list(),nrow=NumberOfRunSteps,ncol=4) #All info stores; state, decision , cost from arrival , cost from observed
  
  for(run in 0:NumberOfRunSteps)
  {
    print(paste("On run ",toString(run)))
    if(run==0)
    {
      #Set up initial States
      #This is set up as before even though we are running a policy
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
      MoveToNode=Policy[[IdenityRow(c(OldsVec,OldvVec),StateSpace)]]
      print(paste("Choosing to move to",toString(MoveToNode)))
      
      #Calculate cost of doing such an action
      RunCost[run]=CostOfAction(c(OldsVec,OldvVec),MoveToNode,n,CostVec,xVec,LambdaVec)
      SeperatedRunCost[run,1]=SeperatedCostOfAction(c(OldsVec,OldvVec),MoveToNode,n,CostVec,xVec,LambdaVec)$CostDueToArrivals
      SeperatedRunCost[run,2]=SeperatedCostOfAction(c(OldsVec,OldvVec),MoveToNode,n,CostVec,xVec,LambdaVec)$CostDueToObserved
      if(run==1)
      {
        CumulativeSeperatedRunCost[run,1]=SeperatedRunCost[run,1]
        CumulativeSeperatedRunCost[run,2]=SeperatedRunCost[run,2]        
      }
      else
      {
        CumulativeSeperatedRunCost[run,1]=CumulativeSeperatedRunCost[run-1,1]+SeperatedRunCost[run,1]
        CumulativeSeperatedRunCost[run,2]=CumulativeSeperatedRunCost[run-1,2]+SeperatedRunCost[run,2]
      }
      
      
      print(paste("Cost of action is",toString(RunCost[run])))
      
      #Evolve System
      sVec=NewSState(OldsVec,MoveToNode,BVec)
      #print("Evolved S is ")
      #print(sVec)
      
      #Creating new vVec using scenario
      vVec=OldvVec
      vVec[MoveToNode]=min(SimulationScenario[MoveToNode,TrackingScenario[MoveToNode]],bVec[MoveToNode])
      TrackingScenario[MoveToNode]=TrackingScenario[MoveToNode]+1
      #vVec=NewMeanVState(OldvVec,sVec,MoveToNode,BVec,bVec,LambdaVec)
      #print("Evolved v is")
      #print(vVec)
      
      AllInfo[[run,1]]=c(OldsVec,OldvVec)
      AllInfo[[run,2]]=MoveToNode
      AllInfo[[run,3]]=SeperatedRunCost[run,1]
      AllInfo[[run,4]]=SeperatedRunCost[run,2]
      
      print(paste("Moved to ",toString(MoveToNode),"Costing ",toString(RunCost[run])))
    }
    
  }
  AverageCost=sum(RunCost)/NumberOfRunSteps
  print(paste("Average cost is ",toString(AverageCost)))
  return(list(Average=AverageCost,CostForStep=RunCost,SeperatedCostForStep=SeperatedRunCost,CumulativeSeperatedCostForStep=CumulativeSeperatedRunCost,FullInfoMatrix=AllInfo))
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














