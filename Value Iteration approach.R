source("Optimal solution by LP.R")

#Function to work out the value for a particular number of steps
#Expects states to be passed as a matrix (with rows being a state with s_1 , s_2,...,s_n,v_1,...,v_2)
ValueFunction<-function(Steps,StateSpace,AdjacencyMatrix,xVec,bVec,CostVec,LambdaVec,PriorValueFunction=NULL)
{
  StateSpaceSize=nrow(StateSpace)
  n=nrow(AdjacencyMatrix)
  BVec=ceiling(xVec)
  
  ValueVector=rep(0,StateSpaceSize)
  
  if(Steps==0) #Base case
  {
    return(ValueVector)
  }
  else if(Steps!=0) 
  {
    #Work out previous step
    if(is.null(PriorValueFunction))
    {
      PriorValueFunction=ValueFunction(Steps-1,StateSpace,AdjacencyMatrix,xVec,bVec,CostVec,LambdaVec)
    }
    
    
    #Form a vector from which to take the minimum for all states
    for(state in 1:StateSpaceSize)
    {
      #Current state is
      CurrentState=StateSpace[state,]
      
      #Current node is
      CurrentNode=which.min(CurrentState[1:n])
    
      #For each state we calculate all the values
      OptionsVector=vector(length=n)
      

      #for each option of node to move to calculate the cost and add the previous cost
      for(option in 1:n)
      {
        
        if(AdjacencyMatrix[CurrentNode,option]==1)
        {
          OptionsVector[option]=CostOfAction(CurrentState,option,n,CostVec,xVec,LambdaVec) #Cost of action
          
          #We now add the expected future cost, this means a proportion for each possible v state we can transistion to by taking this action
          #We can retrive the probability and states from a function
          
          EvolvedStates=NewSVState(CurrentState,option,BVec,bVec,LambdaVec)$State
          EvolvedStatesProb=NewSVState(CurrentState,option,BVec,bVec,LambdaVec)$Prob
          
          #We now need to identify the states and add the cost
          

          for(i in 1:nrow(EvolvedStates))
          {
            IDEvolvedState=IdenityRow(EvolvedStates[i,],StateSpace)
            
            #Add its prior costs * probability
            OptionsVector[option]=OptionsVector[option]+PriorValueFunction[IDEvolvedState]* EvolvedStatesProb[i]
          }
        }
        else
        {
          OptionsVector[option]=NaN
        }
        
        
        
      }
      
      
      
      #Set the Value Vector for that state 
      ValueVector[state]=min(OptionsVector)
    }
  }
  
  #Return all values
  return(ValueVector)
}

ValueIteration<-function(MaxNoSteps,Tolerance,StateSpace,AdjacencyMatrix,xVec,bVec,CostVec,LambdaVec)
{
  step=1
  BoundWidthError=Tolerance+1
  PriorValueFunction=ValueFunction(0,StateSpace,AdjacencyMatrix,xVec,bVec,CostVec,LambdaVec)
  while(step<=MaxNoSteps && BoundWidthError>=Tolerance)
  {
    print(paste("On Step ",toString(step)))
    #Work out the value vector for this number of steps
    NewValueFunction=ValueFunction(step,StateSpace,AdjacencyMatrix,xVec,bVec,CostVec,LambdaVec,PriorValueFunction)

    
    #For this state calcluate min and max for all states
    CostBetweenSteps=NewValueFunction-PriorValueFunction
    MaxForStates=max(CostBetweenSteps)
    MinForStates=min(CostBetweenSteps)
    
    BoundWidth=MaxForStates-MinForStates
    BoundWidthError=BoundWidth/MinForStates
    step=step+1
    
    PriorValueFunction=NewValueFunction
    
    print(MinForStates)
    print(MaxForStates)
  }
  if(BoundWidthError<Tolerance)
  {
    print("Returning due to tolerance reached")
    return(list(LowerBound=MinForStates,UpperBound=MaxForStates))
  }
  else
  {
    print("Returning due to time out")
    return(list(LowerBound=MinForStates,UpperBound=MaxForStates))
  }
  
}