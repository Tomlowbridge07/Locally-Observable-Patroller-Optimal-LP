source("Simulation of Heuristic.R")
source("Value Iteration approach.R")

#This function runs our optimality and heuristic policy to compare the answers
#It runs the optimal policy to find the optimal answer then runs the policy in value iteration
RunTest<-function(AdjacencyMatrix,xVec,bVec,CostVec,LambdaVec,HeuristicFunction,HeuristicDepth,IndexForNodeFunction,MaxStepsForIteration,
                  UseValueItForOptimal=FALSE,ValueItOptMaxSteps=500,ValueItOptTolerance=0.001,PrintOutput=FALSE)
{
  n=nrow(AdjacencyMatrix)
  
  
  #We first solve for optimality using the dual
  if(UseValueItForOptimal==FALSE)
  {
   if(PrintOutput)
   {
     print("We are going to solve the dual problem first")
   }
   DualSolved=SolveDualLP(AdjacencyMatrix,n,xVec,bVec,CostVec,LambdaVec)
   DualObjectiveValue=DualSolved$Value
   StateSpace=DualSolved$StateSpace
   if(PrintOutput)
   {
     print(paste("Dual has been solved for:",toString(DualObjectiveValue)))
   }
  }
  else if(UseValueItForOptimal==TRUE)
  {
    #Instead we use value iteration to get the 'optimal' answer
    if(PrintOutput)
    {
      print("We are going to solve value iteration to near optimum")
    }
    DualSolved=ValueIterationForGame(ValueItOptMaxSteps,ValueItOptTolerance,AdjacencyMatrix,xVec,bVec,CostVec,LambdaVec)
    DualObjectiveValue=DualSolved$LowerBound
    StateSpace=DualSolved$StateSpace
    if(PrintOutput)
    {
     print(paste("Value iteration solved for:",toString(DualObjectiveValue))) 
    }
    
  }

  #We now create the Heuristic policy
  if(PrintOutput)
  {
   print("We are creating the Heuristic Policy")
  }
  PolicyByHeuristic=HeuristicPolicy(HeuristicDepth,HeuristicFunction,n,AdjacencyMatrix,IndexForNodeFunction,CostVec,LambdaVec,bVec,xVec,StateSpace)
  if(PrintOutput)
  {
   print("Policy Has been created")
  }
  
  
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
  print(AbsError)
  PercentageError=(AbsError/DualObjectiveValue) *100
  
  AltAbsError=AverageByFunc-DualObjectiveValue
  AltPercentageError=(AltAbsError/DualObjectiveValue) *100
  
  if(PrintOutput)
  {
   print(paste("Percentage Error by Iteration is:",toString(PercentageError)))
   print(paste("Percentage Error by Function is:",toString(AltPercentageError)))
  }

  return(PercentageError)
}

#The aim of this function is to the run the test (on a scenario) for multiple heuristics
RunTestForMultipleHeuristics<-function(AdjacencyMatrix,xVec,bVec,CostVec,LambdaVec,ListOfHeuristicFunctions,ListOfHeuristicDepths,ListOfIndexForNodeFunctions,
                                       MaxStepsForIteration,PrintOutput=TRUE)
{
  NumberOfHeuristicFuncs=length(ListOfHeuristicFunctions)
  NumberOfHeuristicsDepths=length(ListOfHeuristicDepths)
  NumberOfIndexFuncs=length(ListOfIndexForNodeFunctions)
  
  #We solve the dual problem once
  n=nrow(AdjacencyMatrix)
  #We first solve for optimality using the dual
  if(PrintOutput)
  {
   print("We are going to solve the dual problem first")
  }
  
  DualSolved=SolveDualLP(AdjacencyMatrix,n,xVec,bVec,CostVec,LambdaVec)
  DualObjectiveValue=DualSolved$Value
  StateSpace=DualSolved$StateSpace
  if(PrintOutput)
  {
    print(paste("Dual has been solved for:",toString(DualObjectiveValue)))
  }
  
  
  if(DualObjectiveValue>0)
  {
    ToleranceForIt=10^floor(log10(DualObjectiveValue)-6)
  }
  else
  {
    ToleranceForIt=10^floor(-1000)
  }
  
  
  #Create Storage for errors
  Errors=matrix(0,ncol=4,nrow=NumberOfHeuristicFuncs*NumberOfHeuristicsDepths*NumberOfIndexFuncs)
  PolicyList=list() #Note this list is a list of lists
  counter=1
  

  for(HeuristicFuncNum in 1:NumberOfHeuristicFuncs)
  {
    for(HeuristicDepthNum in 1:NumberOfHeuristicsDepths)
    {
      for(IndexFuncNum in 1:NumberOfIndexFuncs)
      {
        #We now form the policy for the heuristic, at the depth , using the index
        
        if(PrintOutput)
        {
          print(paste("We are in heurisitic func:",toString(HeuristicFuncNum),
                      " depth:",toString(HeuristicDepthNum)," index type:",toString(IndexFuncNum)))
          
          print("We are creating the Heuristic Policy")
        }
        
        PolicyByHeuristic=HeuristicPolicy(ListOfHeuristicDepths[HeuristicDepthNum],ListOfHeuristicFunctions[[HeuristicFuncNum]],
                                          n,AdjacencyMatrix,ListOfIndexForNodeFunctions[[IndexFuncNum]],CostVec,LambdaVec,bVec,xVec,StateSpace)
        if(PrintOutput)
        {
         print("Policy Has been created")
        }
        
        
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
        
        if(PrintOutput)
        {
         print(paste("Percentage Error by Iteration is:",toString(PercentageError)))
         print(paste("Percentage Error by Function is:",toString(AltPercentageError)))
        }

        
        Errors[counter,1]=HeuristicFuncNum
        Errors[counter,2]=HeuristicDepthNum
        Errors[counter,3]=IndexFuncNum
        Errors[counter,4]=PercentageError
        PolicyList[[counter]]=PolicyByHeuristic
        #print(PolicyList)
        counter=counter+1
        
      }
    }
  }
  
  #It is worth noting that by the ordering the structure is BigBlocks (with heuristic func), small blocks (with heuristic depth) and elements (with index func)
  
  #Now identify the best Heuristic
  MinError=min(Errors[,4])
  IDBestHeurisitic=which(Errors[,4]==MinError)
  print(IDBestHeurisitic)
  
  BestHeuristics=Errors[IDBestHeurisitic,1:3]
    
  return(list(Errors=Errors,MinError=MinError,BestHeuristics=BestHeuristics))
}

#This functions generates a random adjacency matrix
GenerateAdjConnectedMatrix<-function(NumNodes,NumEdges)
{
  stopifnot(NumEdges>=(NumNodes-1))
  stopifnot(NumEdges<=(NumNodes+1)*(NumNodes/2))
  
  AdjacencyMatrix=matrix(0,nrow=NumNodes,ncol=NumNodes)
  
  S=sample.int(NumNodes,size=1)
  NS=seq(1,NumNodes,1)
  NS=NS[NS!=S]
  #We now construct it by connecting
  while(length(NS)!=0)
  {

    #Pick a node at random to be include
    if(length(NS)==1)
    {
      NodeToBeAdded=NS
    }
    else
    {
     NodeToBeAdded=sample(NS,size=1)
    }
    
    
    #Now pick a random node to connect it to
    if(length(S)==1)
    {
      ConnectionToS=S
    }
    else
    {
      ConnectionToS=sample(S,size=1)
    }
    
    
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
    if(length(ZerosInMatrix)==1)
    {
      EdgesToAdd=ZerosInMatrix
    }
    else
    {
      EdgesToAdd=sample(ZerosInMatrix,size=NodesNeeded)
    }
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
  #Note becase of the ceiling we will be using the min-1 instead of min
  
  #Generating the Matrix
  NumberOfNodes=ceiling(runif(1,min=MinNumNodes-1,max=MaxNumNodes))
  NumberOfEdges=ceiling(runif(1,min=(NumberOfNodes-1)-1,max=((NumberOfNodes-1)*NumberOfNodes/2)))
  AdjacencyMatrix=GenerateAdjConnectedMatrix(NumberOfNodes,NumberOfEdges)
  
  n=NumberOfNodes
  
  #Generate the xvec,bvec,lambdavec,costvec
  xVec=runif(n,min=MinAttackTime,max=MaxAttackTime)
  bVec=ceiling(runif(n,min=MinObservedSize-1,max=MaxObservedSize))
  LambdaVec=runif(n,min=MinArrivalRate,max=MaxArrivalRate)
  CostVec=runif(n,min=MinCost,max=MaxCost)
  
  return(list(AdjacencyMatrix=AdjacencyMatrix,xVec=xVec,bVec=bVec,LambdaVec=LambdaVec,CostVec=CostVec))
}

#This generates a scenario for a fixed
GenerateTestScenariosVec<-function(NumNodes,MinAttackTimeVec,MaxAttackTimeVec,MinObservedSizeVec,MaxObservedSizeVec,MinArrivalRateVec,
                                   MaxArrivalRateVec,MinCostVec,MaxCostVec)
{
  #Note becase of the ceiling we will be using the min-1 instead of min
  
  #Generating the Matrix
  NumberOfNodes=NumNodes
  NumberOfEdges=ceiling(runif(1,min=(NumberOfNodes-1)-1,max=((NumberOfNodes-1)*NumberOfNodes/2)))
  AdjacencyMatrix=GenerateAdjConnectedMatrix(NumberOfNodes,NumberOfEdges)
  
  n=NumberOfNodes
  # print(MinAttackTimeVec)
  #       print(MaxAttackTimeVec)
  #       print(MinObservedSizeVec)
  #       print(MaxObservedSizeVec)
  #       print(MinArrivalRateVec)
  #       print(MaxArrivalRateVec)
  #       print(MinCostVec)
  #       print(MaxCostVec)
  #Generate the xvec,bvec,lambdavec,costvec
        xVec=vector(length=n)
        bVec=vector(length=n)
        LambdaVec=vector(length=n)
        CostVec=vector(length=n)
  for(i in 1:n)
  {
   xVec[i]=runif(1,min=MinAttackTimeVec[i],max=MaxAttackTimeVec[i])
   bVec[i]=ceiling(runif(1,min=MinObservedSizeVec[i]-1,max=MaxObservedSizeVec[i]))
   LambdaVec[i]=runif(1,min=MinArrivalRateVec[i],max=MaxArrivalRateVec[i])
   CostVec[i]=runif(1,min=MinCostVec[i],max=MaxCostVec[i])    
  }

  
  return(list(AdjacencyMatrix=AdjacencyMatrix,xVec=xVec,bVec=bVec,LambdaVec=LambdaVec,CostVec=CostVec))
  
}
  
#This function is going to run the test for multiple scenarios
RunTestForMultipleScenarios<-function(NumberOfScenarios,ListOfHeuristicFunctions,ListOfHeuristicDepths,ListOfIndexForNodeFunctions,MaxStepsForIteration,
                                      MinNumNodes,MaxNumNodes,MinAttackTime,MaxAttackTime,MinObservedSize,MaxObservedSize,MinArrivalRate,MaxArrivalRate,MinCost,MaxCost)
{
  #This  matrix  stores the minerror,the best heuristic and the scenario
  ScenarioRecording=matrix(list(),nrow=NumberOfScenarios,ncol=8)
  for(ScenarioNumber in 1:NumberOfScenarios)
  {
    #For each scenario we generate  the  scenario
    Scenario=GenerateTestScenarios(MinNumNodes,MaxNumNodes,MinAttackTime,MaxAttackTime,MinObservedSize,MaxObservedSize,MinArrivalRate,MaxArrivalRate,MinCost,MaxCost)
    print("Generated Scenario")
    AdjacencyMatrix=Scenario$AdjacencyMatrix
    xVec=Scenario$xVec
    bVec=Scenario$bVec
    LambdaVec=Scenario$LambdaVec
    CostVec=Scenario$CostVec
    
    ScenarioTest=RunTestForMultipleHeuristics(AdjacencyMatrix,xVec,bVec,CostVec,LambdaVec,ListOfHeuristicFunctions,ListOfHeuristicDepths,ListOfIndexForNodeFunctions,MaxStepsForIteration)
    BestHeuristics=ScenarioTest$BestHeuristics
    MinError=ScenarioTest$MinError
    Errors=ScenarioTest$Errors
    print(Errors)
    
    ScenarioRecording[[ScenarioNumber,1]]=MinError
    ScenarioRecording[[ScenarioNumber,2]]=BestHeuristics
    ScenarioRecording[[ScenarioNumber,3]]=AdjacencyMatrix
    ScenarioRecording[[ScenarioNumber,4]]=xVec
    ScenarioRecording[[ScenarioNumber,5]]=bVec
    ScenarioRecording[[ScenarioNumber,6]]=LambdaVec
    ScenarioRecording[[ScenarioNumber,7]]=CostVec
    ScenarioRecording[[ScenarioNumber,8]]=Errors
    #print(Scenario)
    
  }
  return(ScenarioRecording)
}

#This function is going to run the test for multiple scenarios- With a fixed complete graph
RunTestForMultipleScenariosCompleteGraphs<-function(NumberOfScenarios,ListOfHeuristicFunctions,ListOfHeuristicDepths,ListOfIndexForNodeFunctions,MaxStepsForIteration,
                                      SizeOfCompleteGraph,MinAttackTimeVec,MaxAttackTimeVec,MinObservedSizeVec,MaxObservedSizeVec,MinArrivalRateVec,MaxArrivalRateVec,MinCostVec,MaxCostVec)
{
  #This  matrix  stores the minerror,the best heuristic and the scenario
  ScenarioRecording=matrix(list(),nrow=NumberOfScenarios,ncol=8)
  AdjacencyMatrix=matrix(rep(1,SizeOfCompleteGraph*SizeOfCompleteGraph),nrow=SizeOfCompleteGraph)
  for(ScenarioNumber in 1:NumberOfScenarios)
  {
    #For each scenario we generate  the  scenario
    Scenario=GenerateTestScenariosVec(SizeOfCompleteGraph,MinAttackTimeVec,MaxAttackTimeVec,MinObservedSizeVec,MaxObservedSizeVec,MinArrivalRateVec,
                                   MaxArrivalRateVec,MinCostVec,MaxCostVec)
    print("Generated Scenario")
    xVec=Scenario$xVec
    bVec=Scenario$bVec
    LambdaVec=Scenario$LambdaVec
    CostVec=Scenario$CostVec
    
    ScenarioTest=RunTestForMultipleHeuristics(AdjacencyMatrix,xVec,bVec,CostVec,LambdaVec,ListOfHeuristicFunctions,ListOfHeuristicDepths,ListOfIndexForNodeFunctions,MaxStepsForIteration)
    BestHeuristics=ScenarioTest$BestHeuristics
    MinError=ScenarioTest$MinError
    Errors=ScenarioTest$Errors
    print(Errors)
    
    ScenarioRecording[[ScenarioNumber,1]]=MinError
    ScenarioRecording[[ScenarioNumber,2]]=BestHeuristics
    ScenarioRecording[[ScenarioNumber,3]]=AdjacencyMatrix
    ScenarioRecording[[ScenarioNumber,4]]=xVec
    ScenarioRecording[[ScenarioNumber,5]]=bVec
    ScenarioRecording[[ScenarioNumber,6]]=LambdaVec
    ScenarioRecording[[ScenarioNumber,7]]=CostVec
    ScenarioRecording[[ScenarioNumber,8]]=Errors
    #print(Scenario)
    
  }
  return(ScenarioRecording)
}

#This function creates the adjacency matrix for a line graph
CreateLineGraph<-function(NumNodes)
{
  AdjacencyMatrix=matrix(0,nrow=NumNodes,ncol=NumNodes)
  
  for(i in 1:NumNodes)
  {
    if(i==1)
    {
      AdjacencyMatrix[i,i]=1
      AdjacencyMatrix[i,i+1]=1
    }
    else if(i<NumNodes)
    {
      AdjacencyMatrix[i,i-1]=1
      AdjacencyMatrix[i,i]=1
      AdjacencyMatrix[i,i+1]=1
    }
    else if(i==NumNodes)
    {
      AdjacencyMatrix[i,i-1]=1
      AdjacencyMatrix[i,i]=1 
    }
  }
  print(AdjacencyMatrix)
  return(AdjacencyMatrix)
}

#This function is going to run the test for multiple scenarios- With a fixed complete graph
RunTestForMultipleScenariosLineGraphs<-function(NumberOfScenarios,ListOfHeuristicFunctions,ListOfHeuristicDepths,ListOfIndexForNodeFunctions,MaxStepsForIteration,
                                                    SizeOfLineGraph,MinAttackTimeVec,MaxAttackTimeVec,MinObservedSizeVec,MaxObservedSizeVec,MinArrivalRateVec,MaxArrivalRateVec,MinCostVec,MaxCostVec)
{
  #This  matrix  stores the minerror,the best heuristic and the scenario
  ScenarioRecording=matrix(list(),nrow=NumberOfScenarios,ncol=8)
  AdjacencyMatrix=CreateLineGraph(SizeOfLineGraph)
  for(ScenarioNumber in 1:NumberOfScenarios)
  {
    #For each scenario we generate  the  scenario
    Scenario=GenerateTestScenariosVec(SizeOfLineGraph,MinAttackTimeVec,MaxAttackTimeVec,MinObservedSizeVec,MaxObservedSizeVec,MinArrivalRateVec,MaxArrivalRateVec,MinCostVec,MaxCostVec)
    print("Generated Scenario")
    xVec=Scenario$xVec
    bVec=Scenario$bVec
    LambdaVec=Scenario$LambdaVec
    CostVec=Scenario$CostVec
    
    ScenarioTest=RunTestForMultipleHeuristics(AdjacencyMatrix,xVec,bVec,CostVec,LambdaVec,ListOfHeuristicFunctions,ListOfHeuristicDepths,ListOfIndexForNodeFunctions,MaxStepsForIteration)
    BestHeuristics=ScenarioTest$BestHeuristics
    MinError=ScenarioTest$MinError
    Errors=ScenarioTest$Errors
    print(Errors)
    
    ScenarioRecording[[ScenarioNumber,1]]=MinError
    ScenarioRecording[[ScenarioNumber,2]]=BestHeuristics
    ScenarioRecording[[ScenarioNumber,3]]=AdjacencyMatrix
    ScenarioRecording[[ScenarioNumber,4]]=xVec
    ScenarioRecording[[ScenarioNumber,5]]=bVec
    ScenarioRecording[[ScenarioNumber,6]]=LambdaVec
    ScenarioRecording[[ScenarioNumber,7]]=CostVec
    ScenarioRecording[[ScenarioNumber,8]]=Errors
    #print(Scenario)
    
  }
  return(ScenarioRecording)
}

