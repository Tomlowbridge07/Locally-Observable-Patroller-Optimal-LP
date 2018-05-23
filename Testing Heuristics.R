source("Simulation of Heuristic.R")

#This function runs our optimality and heuristic policy to compare the answers
#It runs the optimal policy to find the optimal answer then runs the policy in value iteration
RunTest<-function(AdjacencyMatrix,xVec,bVec,CostVec,LambdaVec,HeuristicFunction,HeuristicDepth,IndexForNodeFunction,MaxStepsForIteration)
{
  n=nrow(AdjacencyMatrix)
  #We first solve for optimality using the dual
  DualSolved=SolveDualLP(AdjacencyMatrix,n,xVec,bVec,CostVec,LambdaVec)
  DualObjectiveValue=DualSolved$Solution
  
  #We now create the Heuristic policy
  PolicyByHeuristic=HeuristicPolicy(HeuristicDepth,n,HeuristicFunction,AdjacencyMatrix,IndexForNodeFunction,CostVec,LambdaVec,bVec,xVec)
  
  #Run the Heurstic policy
  #We will select the tolerance to depend on the answer above, because we care about the size of the error we will work to a minimum error of 0.01% so that is we care about a tenthousandth error
  ToleranceForIt=10^floor(log10(DualObjectiveValue))
  StateSpace=DualSolved$StateSpace
  ValueItByHeuristic=ValueIterationForPolicy(MaxStepsForIteration,ToleranceForIt,StateSpace,AdjacencyMatrix,xVec,bVec,CostVec,LambdaVec,PolicyByHeuristic)
  
  #We now work out the level of error and return it
  AbsError=ValueItByHeuristic$UpperBound - DualObjectiveValue
  PercentageError=(AbsError/DualObjectiveValue) *100
  
  return(PercentageError)
}