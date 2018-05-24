source("Simulation of Heuristic.R")

#This function runs our optimality and heuristic policy to compare the answers
#It runs the optimal policy to find the optimal answer then runs the policy in value iteration
RunTest<-function(AdjacencyMatrix,xVec,bVec,CostVec,LambdaVec,HeuristicFunction,HeuristicDepth,IndexForNodeFunction,MaxStepsForIteration)
{
  n=nrow(AdjacencyMatrix)
  #We first solve for optimality using the dual
  print("We are going to solve the dual problem first")
  DualSolved=SolveDualLP(AdjacencyMatrix,n,xVec,bVec,CostVec,LambdaVec)
  DualObjectiveValue=DualSolved$Value
  StateSpace=DualSolved$StateSpace
  print(paste("Dual has been solved for:",toString(DualObjectiveValue)))
  
  #We now create the Heuristic policy
  print("We are creating the Heuristic Policy")
  PolicyByHeuristic=HeuristicPolicy(HeuristicDepth,HeuristicFunction,n,AdjacencyMatrix,IndexForNodeFunction,CostVec,LambdaVec,bVec,xVec,StateSpace)
  print("Policy Has been created")
  
  #Run the Heurstic policy
  #We will select the tolerance to depend on the answer above, because we care about the size of the error we will work to a minimum error of 0.01% so that is we care about a tenthousandth error
  ToleranceForIt=10^floor(log10(DualObjectiveValue)-4)

  
  ValueItByHeuristic=ValueIterationForPolicy(MaxStepsForIteration,ToleranceForIt,StateSpace,AdjacencyMatrix,xVec,bVec,CostVec,LambdaVec,PolicyByHeuristic)
  
  #ValueFuncByHeuristic=ValueFunctionForPolicy(100,StateSpace,AdjacencyMatrix,xVec,bVec,CostVec,LambdaVec,PolicyByHeuristic)
  ValueFuncByHeuristic=ValueItByHeuristic$ValueFunction
  ValueFunSteps=ValueItByHeuristic$StepsRun
  AverageByFunc=mean(ValueFuncByHeuristic)/ValueFunSteps
  
  #We now work out the level of error and return it
  AbsError=ValueItByHeuristic$UpperBound - DualObjectiveValue
  PercentageError=(AbsError/DualObjectiveValue) *100
  
  AltAbsError=AverageByFunc-DualObjectiveValue
  AltPercentageError=(AltAbsError/DualObjectiveValue) *100
  
  print(paste("Percentage Error by Iteration is:",toString(PercentageError)))
  print(paste("Percentage Error by Function is:",toString(AltPercentageError)))
  return(PercentageError)
}

#The aim of this function is to the run the test (on a scenario) for multiple heuristics
RunTestForMultipleHeuristics<-function(AdjacencyMatrix,xVec,bVec,CostVec,LambdaVec,ListOfHeuristicFunctions,ListOfHeuristicDepths,ListOfIndexForNodeFunctions,MaxStepsForIteration)
{
  NumberOfHeuristicFuncs=length(ListOfHeuristicFunctions)
  NumberOfHeuristicsDepths=length(ListOfHeuristicDepths)
  NumberOfIndexFuncs=length(ListOfIndexForNodeFunctions)
  
  #We solve the dual problem once
  n=nrow(AdjacencyMatrix)
  #We first solve for optimality using the dual
  print("We are going to solve the dual problem first")
  DualSolved=SolveDualLP(AdjacencyMatrix,n,xVec,bVec,CostVec,LambdaVec)
  DualObjectiveValue=DualSolved$Value
  StateSpace=DualSolved$StateSpace
  print(paste("Dual has been solved for:",toString(DualObjectiveValue)))
  
  ToleranceForIt=10^floor(log10(DualObjectiveValue)-4)
  
  #Create Storage for errors
  Errors=matrix(0,ncol=4,nrow=NumberOfHeuristicFuncs*NumberOfHeuristicsDepths*NumberOfIndexFuncs)
  counter=1
  

  for(HeuristicFuncNum in 1:NumberOfHeuristicFuncs)
  {
    for(HeuristicDepthNum in 1:NumberOfHeuristicsDepths)
    {
      for(IndexFuncNum in 1:NumberOfIndexFuncs)
      {
        #We now form the policy for the heuristic, at the depth , using the index
        
        print("We are creating the Heuristic Policy")
        PolicyByHeuristic=HeuristicPolicy(ListOfHeuristicDepths[HeuristicDepthNum],ListOfHeuristicFunctions[[HeuristicFuncNum]],
                                          n,AdjacencyMatrix,ListOfIndexForNodeFunctions[[IndexFuncNum]],CostVec,LambdaVec,bVec,xVec,StateSpace)
        print("Policy Has been created")
        
        #Run the heuristic
        ValueItByHeuristic=ValueIterationForPolicy(MaxStepsForIteration,ToleranceForIt,StateSpace,AdjacencyMatrix,xVec,bVec,CostVec,LambdaVec,PolicyByHeuristic)
        
        ValueFuncByHeuristic=ValueItByHeuristic$ValueFunction
        ValueFunSteps=ValueItByHeuristic$StepsRun
        AverageByFunc=mean(ValueFuncByHeuristic)/ValueFunSteps
        
        #We now work out the level of error and return it
        AbsError=ValueItByHeuristic$UpperBound - DualObjectiveValue
        PercentageError=(AbsError/DualObjectiveValue) *100
        
        AltAbsError=AverageByFunc-DualObjectiveValue
        AltPercentageError=(AltAbsError/DualObjectiveValue) *100
        
        print(paste("Percentage Error by Iteration is:",toString(PercentageError)))
        print(paste("Percentage Error by Function is:",toString(AltPercentageError)))
        
        Errors[counter,1]=HeuristicFuncNum
        Errors[counter,2]=HeuristicDepthNum
        Errors[counter,3]=IndexFuncNum
        Errors[counter,4]=PercentageError
        counter=counter+1
        
      }
    }
  }
  
  #It is worth noting that by the ordering the structure is BigBlocks (with heuristic func), small blocks (with heuristic depth) and elements (with index func)
  
  #Now identify the best Heuristic
  MinError=min(Errors[,4])
  IDBestHeurisitic=which(Errors[,4]==MinError)
  print(IDBestHeurisitic)
  
  BestHeurisitcs=Errors[IDBestHeurisitic,1:3]
    
  return(list(Errors=Errors,MinError=MinError,BestHeurisitcs=BestHeurisitcs))
}

#This functions generates a random adjacency matrix
GenerateAdjConnectedMatrix<-function(NumNodes,NumEdges)
{
  stopifnot(NumEdges>(NumNodes-1))
  stopifnot(NumEdges<(NumNodes+1)*(NumNodes/2))
  
  AdjacencyMatrix=matrix(0,nrow=NumNodes,ncol=NumNodes)
  
  S=sample.int(NumNodes,size=1)
  NS=seq(1,NumNodes,1)
  NS=NS[any(NS!=S)]
  #We now construct it by connecting
  while(length(NS)!=0)
  {
    #Pick a node at random to be include
    NodeToBeAdded=sample(NS,size=1)
    
    #Now pick a random node to connect it to
    ConnectionToS=sample(S,size=1)
    
    #Now add this edge into the graph
    AdjacencyMatrix[NodeToBeAdded,ConnectionToS]=1
    AdjacencyMatrix[ConnectionToS,NodeToBeAdded]=1
    
    #We now removed this node from not selected and add it to selected
    S=c(S,NodeToBeAdded)
    NS=NS[NS!=NodeToBeAdded]
  }
  
  #Now we have a spanning tree, so
  NodesNeeded=NumEdges-(NumNodes-1)
  
  Consider=lower.tri(AdjacencyMatrix,diag=FALSE)
  Consider=Consider & (AdjacencyMatrix==0)
  #Now we add at random these edges into this random spanning tree
  ZerosInMatrix=which(Consider==TRUE)

  #Create a matrix which we will later transpose and add on
  AddOn=matrix(0,nrow=NumNodes,ncol=NumNodes)
  
  if(NodesNeeded>0)
  {
   EdgesToAdd=sample(ZerosInMatrix,size=NodesNeeded)
   AddOn[EdgesToAdd]=1
  }

  #Add on the edges symmetrically
   AdjacencyMatrix=AdjacencyMatrix +AddOn + t(AddOn)

  #We now add diagonal ones
  for(i in 1:NumNodes)
  {
    AdjacencyMatrix[i,i]=1
  }
  
  return(AdjacencyMatrix)
}

#Generate Scenarios- This function generates a collection of adjacency matrices, xVec, bVec, LambdaVec and CostVec within a given range.
#Note. The matrix will be connected.
GenerateTestScenarios<-function(MinNumNodes,MaxNumNodes,MinAttackTime,MaxAttackTime,MinObservedSize,MaxObservedSize,MinArrivalRate,MaxArrivalRate,MinCost,MaxCost)
{
  #Generating the Matrix
  NumberOfNodes
}
  
#This function is 
RunTestForMultipleScenarios