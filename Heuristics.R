source("Index Implementation.R")

#One step Benefit Heuristic
#This function decides given the current state to look at the available actions and pick the action that produces the highest index
OneStepBenefitHeuristic<-function(n,AdjacencyMatrix,IndexForNodeFunction,sVec,vVec,CostVec,LambdaVec,bVec,xVec,vMaxVec=NULL)
{
  if(is.null(vMaxVec))
  {
    #Create vMaxVec
    vMaxVec=CreateVMaxVector(n,LambdaVec,bVec,xVec)
  }
  
  #Get indices
  NodeIndexes=IndicesForNodes(n,IndexForNodeFunction,sVec,vVec,CostVec,LambdaVec,bVec,xVec,vMaxVec)
  #print(NodeIndexes)
  
  #Identify which state we are at by
  CurrentNode=which.min(sVec)
  
  #Actions available
  Actions=AdjacencyMatrix[CurrentNode,]
  
  IndexForActions=Actions * NodeIndexes
  
  #Find best action
  BestAction=which.max(IndexForActions)
  
  return(BestAction)
}

OneStepPenaltyHeuristic<-function(n,AdjacencyMatrix,IndexForNodeFunction,sVec,vVec,CostVec,LambdaVec,bVec,xVec,vMaxVec=NULL)
{
  if(is.null(vMaxVec))
  {
    #Create vMaxVec
    vMaxVec=CreateVMaxVector(n,LambdaVec,bVec,xVec)
  }
  
  #Get indices
  NodeIndexes=IndicesForNodes(n,IndexForNodeFunction,sVec,vVec,CostVec,LambdaVec,bVec,xVec,vMaxVec)
  
  #Identify which state we are at by
  CurrentNode=which.min(sVec)
  
  #Actions available
  Actions=AdjacencyMatrix[CurrentNode,]
  
  IndexForActions=Actions * NodeIndexes
  
  PenaltyForActions=rep(sum(IndexForActions),length(Actions))-IndexForActions
  #print(PenaltyForActions)
  
  #Find best action
  BestAction=which.min(PenaltyForActions)
  
  return(BestAction)
}

DeterministicCostEvaluationOfPath<-function(Path,n,sVec,vVec,CostVec,LambdaVec,bVec,xVec,vMaxVec=NULL)
{
  BVec=ceiling(xVec)
  #We follow the path for as many steps as it has
  i=1
  SumOfCost=0
  while(Path[i]!=0 && i<=length(Path))
  {
    #Add actions cost
    SumOfCost=SumOfCost+CostOfAction(c(sVec,vVec),Path[i],n,CostVec,xVec,LambdaVec)
    
    #Evolve DETERMINISTICALLY
    sVec=NewSState(sVec,Path[i],BVec)
    vVec=NewMeanVState(vVec,sVec,Path[i],BVec,bVec,LambdaVec)

    i=i+1
  }
  
  return(list(Average=SumOfCost/(i-1),Overall=SumOfCost))
}

#Here we sum the indices up over the number of steps for all paths of the Number of steps
MultiStepBenefitHeuristic<-function(NoSteps,n,AdjacencyMatrix,IndexForNodeFunction,sVec,vVec,CostVec,LambdaVec,bVec,xVec,vMaxVec=NULL)
{
  if(is.null(vMaxVec))
  {
    #Create vMaxVec
    vMaxVec=CreateVMaxVector(n,LambdaVec,bVec,xVec)
  }
  
  BVec=ceiling(xVec)
  
  
  Paths=matrix(0,nrow=0,ncol=NoSteps)
  #Note. We will be using mean v Evolution for use with the index
  EvolvedStates=matrix(0,nrow=0,ncol=2*n)
  BenefitForPath=vector(length=0)
  BestPathforStep=matrix(0,nrow=NoSteps,ncol=NoSteps)
  
  for(Step in 1:NoSteps)
  {
    if(Step==1)
    {
      #We run the initial set up
      CurrentNode=which.min(sVec)
      
      #Actions
      Actions=AdjacencyMatrix[CurrentNode,]
      

      NodeIndexes=IndicesForNodes(n,IndexForNodeFunction,sVec,vVec,CostVec,LambdaVec,bVec,xVec,vMaxVec)
      print(NodeIndexes)
      
      BenefitForAction=Actions * NodeIndexes
      
      #Form Paths
      for(action in 1:n)
      {
        if(Actions[action]==1)
        {
          Paths=rbind(Paths,c(action,rep(0,NoSteps-1)))
          NewsVec=NewSState(sVec,action,BVec)
          EvolvedStates=rbind(EvolvedStates,c(NewsVec,NewMeanVState(vVec,NewsVec,action,BVec,bVec,LambdaVec)))
          BenefitForPath=c(BenefitForPath,BenefitForAction[action])

        }
      }
      print(paste("I am about to compare all paths of length ",toString(Step)))
      print(Paths)
      print("They have a benefit of")
      print(BenefitForPath)
      BestPath=Paths[which.max(BenefitForPath),]
      BestPathforStep[Step,]=BestPath
      print("I have chosen the path")
      print(BestPath)
    }
    else
    {
      
      #We need to copy and expand each
      OldPaths=Paths
      OldEvolvedStates=EvolvedStates
      OldBenefitForPath=BenefitForPath
      
      Paths=matrix(0,nrow=0,ncol=NoSteps)
      EvolvedStates=matrix(0,nrow=0,ncol=2*n)
      BenefitForPath=vector(length=0)
      
      #We now look at expanding each
      for(row in 1:nrow(OldPaths))
      {
        #for each row we expand it to allow all possible actions
        Actions=AdjacencyMatrix[OldPaths[row,Step-1],]
        

        NodeIndexes=IndicesForNodes(n,IndexForNodeFunction,OldEvolvedStates[row,1:n],OldEvolvedStates[row,(n+1):(2*n)],CostVec,LambdaVec,bVec,xVec,vMaxVec)

        
        BenefitForAction=Actions * NodeIndexes

        
        #Form Paths
        for(action in 1:n)
        {
          if(Actions[action]==1)
          {
            Paths=rbind(Paths,c(OldPaths[row,1:(Step-1)],action,rep(0,NoSteps-Step)))
            NewsVec=NewSState(OldEvolvedStates[row,1:n],action,BVec)
            EvolvedStates=rbind(EvolvedStates,c(NewsVec,NewMeanVState(OldEvolvedStates[row,(n+1):(2*n)],NewsVec,action,BVec,bVec,LambdaVec)))
            BenefitForPath=c(BenefitForPath,OldBenefitForPath[row]+BenefitForAction[action])
          }
        }
        
      }
      print(paste("I am about to compare all paths of length ",toString(Step)))
      print(Paths)
      print("They have a benefit of")
      print(BenefitForPath)
      BestPath=Paths[which.max(BenefitForPath),]
      BestPathforStep[Step,]=BestPath
      print("I have chosen the path")
      print(BestPath)

    }
    
  }
  
  #For each look ahead step we have a path
  #print(BestPathforStep)
  AverageCostforPath=vector(length=NoSteps)
  #We now need to see how good they perform
  for(i in 1:NoSteps)
  {
    #We compute the average cost of following such a strategy to decide which paths to pick
    #We use determinsitic evolution to the mean state in v
    AverageCostforPath[i]=DeterministicCostEvaluationOfPath(BestPathforStep[i,],n,sVec,vVec,CostVec,LambdaVec,bVec,xVec,vMaxVec)$Average
  }
  print("about to print paths and determinisitic cost of paths")
  print(BestPathforStep)
  print(AverageCostforPath)
  OverallBestPath=BestPathforStep[which.min(AverageCostforPath),]
  return(OverallBestPath[1])
}






#Here we sum the indices up over the number of steps for all paths of the Number of steps
MultiStepPenaltyHeuristic<-function(NoSteps,n,AdjacencyMatrix,IndexForNodeFunction,sVec,vVec,CostVec,LambdaVec,bVec,xVec,vMaxVec=NULL)
{
  if(is.null(vMaxVec))
  {
    #Create vMaxVec
    vMaxVec=CreateVMaxVector(n,LambdaVec,bVec,xVec)
  }
  
  BVec=ceiling(xVec)
  
  
  Paths=matrix(0,nrow=0,ncol=NoSteps)
  #Note. We will be using mean v Evolution for use with the index
  EvolvedStates=matrix(0,nrow=0,ncol=2*n)
  PenaltyForPath=vector(length=0)
  BestPathforStep=matrix(0,nrow=NoSteps,ncol=NoSteps)
  
  for(Step in 1:NoSteps)
  {
    if(Step==1)
    {
      #We run the initial set up
      CurrentNode=which.min(sVec)
      
      #Actions
      Actions=AdjacencyMatrix[CurrentNode,]
      
      NodeIndexes=IndicesForNodes(n,IndexForNodeFunction,sVec,vVec,CostVec,LambdaVec,bVec,xVec,vMaxVec)
      
      IndexForActions= Actions * NodeIndexes
      
      PenaltyForAction=rep(sum(IndexForActions),length(Actions))-IndexForActions
      
      #Form Paths
      for(action in 1:n)
      {
        if(Actions[action]==1)
        {
          Paths=rbind(Paths,c(action,rep(0,NoSteps-1)))
          NewsVec=NewSState(sVec,action,BVec)
          EvolvedStates=rbind(EvolvedStates,c(NewsVec,NewMeanVState(vVec,NewsVec,action,BVec,bVec,LambdaVec)))
          PenaltyForPath=c(PenaltyForPath,PenaltyForAction[action])
        }
      }
      BestPath=Paths[which.min(PenaltyForPath),]
      BestPathforStep[Step,]=BestPath
    }
    else
    {
      #We need to copy and expand each
      OldPaths=Paths
      OldEvolvedStates=EvolvedStates
      OldPenaltyForPath=PenaltyForPath
      
      Paths=matrix(0,nrow=0,ncol=NoSteps)
      EvolvedStates=matrix(0,nrow=0,ncol=2*n)
      PenaltyForPath=vector(length=0)
      
      #We now look at expanding each
      for(row in 1:nrow(OldPaths))
      {
        #for each row we expand it to allow all possible actions
        Actions=AdjacencyMatrix[OldPaths[row,Step-1],]
        
        NodeIndexes=IndicesForNodes(n,IndexForNodeFunction,OldEvolvedStates[row,1:n],OldEvolvedStates[row,(n+1):(2*n)],CostVec,LambdaVec,bVec,xVec,vMaxVec)
        
        IndexForActions= Actions * NodeIndexes
        
        PenaltyForAction=rep(sum(IndexForActions),length(Actions))-IndexForActions
        
        #Form Paths
        for(action in 1:n)
        {
          if(Actions[action]==1)
          {
            Paths=rbind(Paths,c(OldPaths[row,1:(Step-1)],action,rep(0,NoSteps-Step)))
            NewsVec=NewSState(OldEvolvedStates[row,1:n],action,BVec)
            EvolvedStates=rbind(EvolvedStates,c(NewsVec,NewMeanVState(OldEvolvedStates[row,(n+1):(2*n)],NewsVec,action,BVec,bVec,LambdaVec)))
            PenaltyForPath=c(PenaltyForPath,OldPenaltyForPath[row]+PenaltyForAction[action])
          }
        }
        
      }
      BestPath=Paths[which.min(PenaltyForPath),]
      BestPathforStep[Step,]=BestPath
    }
    
  }
  
  #For each look ahead step we have a path
  #print(BestPathforStep)
  AverageCostforPath=vector(length=NoSteps)
  #We now need to see how good they perform
  for(i in 1:NoSteps)
  {
    #We compute the average cost of following such a strategy to decide which paths to pick
    #We use determinsitic evolution to the mean state in v
    AverageCostforPath[i]=DeterministicCostEvaluationOfPath(BestPathforStep[i,],n,sVec,vVec,CostVec,LambdaVec,bVec,xVec,vMaxVec)$Average
  }
  #print(AverageCostforPath)
  OverallBestPath=BestPathforStep[which.min(AverageCostforPath),]
  return(OverallBestPath[1])
}


StartingNodeHeuristic<-function(n,IndexForNodeFunction,CostVec,LambdaVec,bVec,xVec,vMaxVec=NULL)
{
  if(is.null(vMaxVec))
  {
    #Create vMaxVec
    vMaxVec=CreateVMaxVector(n,LambdaVec,bVec,xVec)
  }
  
  #Elapsed state, we will consider the indices when we have not visited in a very long time.
  StartingSVec=ceiling(bVec)+1
  #and we will assume that all v's are in their mean state.
  StartingVVec=vector(length=n)
  for(i in 1:n)
  {
    StartingVVec[i]=TruncPoissionMean(LambdaVec[i],bVec[i])
  }
  
  #To decide a starting node we work out the index for all nodes (when S is maximum) and pick the biggest
  NodeIndexes=IndicesForNodes(n,IndexForNodeFunction,StartingSVec,StartingVVec,CostVec,LambdaVec,bVec,xVec,vMaxVec)
  
  BestStart=which.max(NodeIndexes)
  return(BestStart)
}