source("Optimal solution by LP.R")

#Function to work out the value for a particular number of steps
#Expects states to be passed as a matrix (with rows being a state with s_1 , s_2,...,s_n,v_1,...,v_2)
ValueFunction<-function(Steps,StateSpace,AdjacencyMatrix,xVec,bVec,CostVec,LambdaVec,PriorValueFunction=NULL,PriorActionsMatrix=NULL)
{
  StateSpaceSize=nrow(StateSpace)
  n=nrow(AdjacencyMatrix)
  BVec=ceiling(xVec)
  
  #Stores the value of this iteration
  ValueVector=rep(0,StateSpaceSize)
  
  #Store a list of actions for this iteration
  AddOnActionsMatrix=matrix(list(),nrow=1,ncol=StateSpaceSize)
  
  if(Steps==0) #Base case
  {
    return(list(Values=ValueVector,Actions=AddOnActionsMatrix))
  }
  else if(Steps!=0) 
  {
    #Work out previous step
    if(is.null(PriorValueFunction) || is.null(PriorActionsMatrix))
    {
      Prior=ValueFunction(Steps-1,StateSpace,AdjacencyMatrix,xVec,bVec,CostVec,LambdaVec)
      PriorValueFunction=Prior$Values
      PriorActionsMatrix=Prior$Actions
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
      #We will also store the action used to achieve this.
      AddOnActionsMatrix[[1,state]]=which(OptionsVector==ValueVector[state])
    }
    ActionsMatrix=rbind(PriorActionsMatrix,AddOnActionsMatrix)

  }
  
  #Return all values
  return(list(Values=ValueVector,Actions=ActionsMatrix))
}

ValueIteration<-function(MaxNoSteps,Tolerance,StateSpace,AdjacencyMatrix,xVec,bVec,CostVec,LambdaVec)
{
  step=1
  BoundWidthError=Tolerance+1
  PriorValueFunction=ValueFunction(0,StateSpace,AdjacencyMatrix,xVec,bVec,CostVec,LambdaVec)$Values
  PriorActionsMatrix=ValueFunction(0,StateSpace,AdjacencyMatrix,xVec,bVec,CostVec,LambdaVec)$Actions
  while(step<=MaxNoSteps && BoundWidthError>=Tolerance)
  {
    print(paste("On Step ",toString(step)))
    
    #Work out the value vector for this number of steps
    New=ValueFunction(step,StateSpace,AdjacencyMatrix,xVec,bVec,CostVec,LambdaVec,PriorValueFunction,PriorActionsMatrix)
    NewValueFunction=New$Values
    NewActionsMatrix=New$Actions

    
    #For this state calcluate min and max for all states
    CostBetweenSteps=NewValueFunction-PriorValueFunction
    MaxForStates=max(CostBetweenSteps)
    MinForStates=min(CostBetweenSteps)
    
    BoundWidth=MaxForStates-MinForStates
    BoundWidthError=BoundWidth/MinForStates
    step=step+1
    
    PriorValueFunction=NewValueFunction
    PriorActionsMatrix=NewActionsMatrix
    
    print(MinForStates)
    print(MaxForStates)
  }
  #We need to remove the top row of the NewActionsMatrix
  
  TrueActionsMatrix=NewActionsMatrix[-1,]
  
  if(BoundWidthError<Tolerance)
  {
    print("Returning due to tolerance reached")
    return(list(LowerBound=MinForStates,UpperBound=MaxForStates,Actions=TrueActionsMatrix))
  }
  else
  {
    print("Returning due to time out")
    return(list(LowerBound=MinForStates,UpperBound=MaxForStates,Actions=TrueActionsMatrix))
  }
}

ValueIterationForGame<-function(MaxNoSteps,Tolerance,AdjacencyMatrix,xVec,bVec,CostVec,LambdaVec)
{
  #Set up games statespace
  BVec=ceiling(xVec)
  n=nrow(AdjacencyMatrix)
  StateSpace=CreateSVStates(n,BVec,bVec)
  
  #Solve the iteration
  ValueIt=ValueIteration(MaxNoSteps,Tolerance,StateSpace,AdjacencyMatrix,xVec,bVec,CostVec,LambdaVec)
  LowerBound=ValueIt$LowerBound
  UpperBound=ValueIt$UpperBound
  ActionsMatrix=ValueIt$Actions
  
  #This is assumed to be the actions to be performed in the limit (Theoretical justification may be needed)
  EndActions=ActionsMatrix[nrow(ActionsMatrix),]

 return(list(LowerBound=LowerBound,UpperBound=UpperBound,Actions=ActionsMatrix,EndActions=EndActions))
}