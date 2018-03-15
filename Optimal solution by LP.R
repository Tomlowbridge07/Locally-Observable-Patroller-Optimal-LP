library(gtools)
library(utils)
library(lpSolve)

source("Truncated Poisson distribution.R")

#Function to work out the size of the state space, we will assume that B>n
#Note. B is definded as being the smallest integer such that we are less than it, so B+1 is the last state space
SizeOfOMEGA<-function(n,B,b)
{
  SumOfOmegas=0
  for(i in 0:(n-1))
  {
   SumOfOmegas=SumOfOmegas+(choose(n,1)*choose(n-1,i)*choose(B-1,n-1-i)*factorial(n-1-i))
  }
  SumOfOmegav=b^n
  SumOfOmega=SumOfOmegas*SumOfOmegav
  return(SumOfOmega)
}


#Creating a matrix with all stats for s
CreateSStatesForConstantB<-function(n,B)
{
  StateMatrix=matrix(nrow=0,ncol=n)
  for(ChoiceToBeOne in 1:n) #choose which element is 1
  {
    #print(paste("Choosing element ",toString(ChoiceToBeOne)," to be set to 1"))
    OneRow=vector(length = n)
    #We set on of the elements to be 1
    OneRow[ChoiceToBeOne]=1
    
    #Create a vector of values not yet assigned
    NotOneInRow=1:n
    NotOneInRow=NotOneInRow[NotOneInRow!=ChoiceToBeOne]
    
    #Choice 0 to (n-1) to be B+1 (Note. The choice of 0 may be reduced if n-1 > B-1 ,as we will not have enough unique to pick)
    for(ChoiceToBeFull in max(0,n-B):(n-1))
    {
      #print(paste("Our current Choice to be full is ",toString(ChoiceToBeFull)))
      if(ChoiceToBeFull==0)
      {
        #We do not need to assign any to be Full i.e B+2 , so no changes
        NotAssignedInRow=NotOneInRow
        Row=OneRow
        
        if(length(NotAssignedInRow)==0)
        {
          #If we have non left to asign we can just store the Row
          StateMatrix=rbind(StateMatrix,Row)
          #print("Saving the Row in the State Matrix")
          #print("Row is ")
          #print(Row)
        }
        else
        {
        #For the remaining Not yet assinged elements we need a vector of elements that they can uniquely pick from
        UniqueStatesLeft=2:B
        
        #Now we want all combinations where we pick the remaining elements
        UniqueCombinationsMatrix<-combinations(length(UniqueStatesLeft),length(NotAssignedInRow),UniqueStatesLeft)
        
        #for each combination we now need to permutate them and assign them
        for(UniqueCombRow in 1:nrow(UniqueCombinationsMatrix))
        {
          #Get the chosen values to be assigned
          ChosenVector=UniqueCombinationsMatrix[UniqueCombRow,]
          
          #Create the perumatutation matrix
          ChosenPermutationMatrix<-permutations(length(ChosenVector),length(ChosenVector),ChosenVector)
          
          #for each permutation we now want to assign them to the elements we have not yet assigned
          for(ChosenPermutationRow in 1:nrow(ChosenPermutationMatrix))
          {
            #Assign values to the row
            Row[NotAssignedInRow]=ChosenPermutationMatrix[ChosenPermutationRow,]
            
            #As we have no completed making one state we will store it in a matrix
            StateMatrix=rbind(StateMatrix,Row)
            #print("Saving the Row in the State Matrix")
            #print("Row is ")
            #print(Row)
          }
        }
        }
      }
      else
      {
      #we now generate vectors to know which elements will be set to B+1
      FullCombinationsMatrix<-combinations(length(NotOneInRow),ChoiceToBeFull,NotOneInRow)
      #print("Full combinations matrix is")
      #print(FullCombinationsMatrix)
      
      for(FullCombRow in 1:nrow(FullCombinationsMatrix))
      {
        #for each combination set these to be B+2
        Row=OneRow

        Row[FullCombinationsMatrix[FullCombRow,]]=B+1
        #print("Setting selected elements to be Full")
        #print(Row)
        
        #Set that these combinations have now been assigned
        NotAssignedInRow=NotOneInRow
        #print(NotAssignedInRow)
        
        for(ValueToRemove in FullCombinationsMatrix[FullCombRow,])
        {
          #Remove this from the assigned vector
          NotAssignedInRow=NotAssignedInRow[NotAssignedInRow!=ValueToRemove]
        }
        #print("We have not yet assigned")
        #print(NotAssignedInRow)
        
        if(length(NotAssignedInRow)==0)
        {
          StateMatrix=rbind(StateMatrix,Row)
          #print("Saving the Row in the State Matrix")
          #print("Row is ")
          #print(Row)
        }
        else
        {
        
        #For the remaining Not yet assinged elements we need a vector of elements that they can uniquely pick from
        UniqueStatesLeft=2:B
        
        #Now we want all combinations where we pick the remaining elements
        UniqueCombinationsMatrix<-combinations(length(UniqueStatesLeft),length(NotAssignedInRow),UniqueStatesLeft)
        
        #for each combination we now need to permutate them and assign them
        for(UniqueCombRow in 1:nrow(UniqueCombinationsMatrix))
        {
          #Get the chosen values to be assigned
          ChosenVector=UniqueCombinationsMatrix[UniqueCombRow,]
          
          #Create the perumatutation matrix
          ChosenPermutationMatrix<-permutations(length(ChosenVector),length(ChosenVector),ChosenVector)
          
          #for each permutation we now want to assign them to the elements we have not yet assigned
          for(ChosenPermutationRow in 1:nrow(ChosenPermutationMatrix))
          {
            #Assign values to the row
            Row[NotAssignedInRow]=ChosenPermutationMatrix[ChosenPermutationRow,]
            
            #As we have no completed making one state we will store it in a matrix
            StateMatrix=rbind(StateMatrix,Row)
            #print(StateMatrix)
            #print("Saving the Row in the State Matrix")
            #print("Row is ")
            #print(Row)
          }
          
        }
        }
        
        
        
      }
      }
    }
  }
  return(StateMatrix)
}

#Using a B vector and not just a flat value
CreateSStates<-function(n,BVec)
{
  #First create for the highest B
  BMax=max(BVec)
  StateSpace=CreateSStatesForConstantB(n,BMax)
  
  #Then we have to go through the state space and reduce the columns to the maximum allowed for the node
  for(i in 1:nrow(StateSpace))
  {
    for(j in 1:n)
    {
      if(StateSpace[i,j] > BVec[j]+1)
      {
        StateSpace[i,j]=BVec[j]+1
      }
    }
  }
  #However now we may have repeated states, so we will go through and remove duplicate rows
  StateSpace=unique(StateSpace)
  
  return(StateSpace)
}

CreateVStatesForConstantb<-function(n,b)
{
  #Setting initial sequences
  StateMatrix=matrix(nrow=(b+1),ncol=n)
  for(i in 1:(b+1))
  {
    StateMatrix[i,1]=i-1
  }

  for(j in 2:n)
  {
   TempMatrix=matrix(nrow=0,ncol=n)
   
   #Copy prior sequence
   for(row in 1:nrow(StateMatrix))
   {
    BlockMatrix=matrix(rep(StateMatrix[row,],b+1),nrow=(b+1),ncol=n,byrow = TRUE)
    TempMatrix=rbind(TempMatrix,BlockMatrix)
   }

   #now for each row in the Temp Matrix we add 0,..,b
   for(row in 1:nrow(TempMatrix))
   {
     TempMatrix[row,j]=(row-1) %% (b+1)
   }
   
   StateMatrix=TempMatrix
   #print(StateMatrix)
  }
  return(StateMatrix)
}

CreateVStates<-function(n,bVec)
{
  #Setting initial sequences
  StateMatrix=matrix(nrow=(bVec[1]+1),ncol=n)
  for(i in 1:(bVec[1]+1))
  {
    StateMatrix[i,1]=i-1
  }
  
  for(j in 2:n)
  {
    TempMatrix=matrix(nrow=0,ncol=n)
    
    #Copy prior sequence
    for(row in 1:nrow(StateMatrix))
    {
      BlockMatrix=matrix(rep(StateMatrix[row,],bVec[j]+1),nrow=(bVec[j]+1),ncol=n,byrow = TRUE)
      TempMatrix=rbind(TempMatrix,BlockMatrix)
    }
    
    #now for each row in the Temp Matrix we add 0,..,b
    for(row in 1:nrow(TempMatrix))
    {
      TempMatrix[row,j]=(row-1) %% (bVec[j]+1)
    }
    
    StateMatrix=TempMatrix
    #print(StateMatrix)
  }
  return(StateMatrix)
}

CreateSVStatesForConstantBb<-function(n,B,b)
{
  SStates=CreateSStatesForConstantB(n,B)
  VStates=CreateVStatesForConstantb(n,b)
  
  #Then for each row in SStates we all all rows in V states to be added
  
  FullStateSpace=matrix(nrow=(nrow(SStates)*nrow(VStates)),ncol=2*n)
  
  for(i in 1:nrow(FullStateSpace))
  {
     #split i into its Block and remainder
     Block=ceiling(i/nrow(VStates))
     Remainder=((i-1) %% nrow(VStates)) +1
     FullStateSpace[i,1:n]=SStates[Block,]
     
     #Before we can insert the VState, we need to check that if s=B+2, then we set v=0
     WorkingVState=VStates[Remainder,]
     for(j in 1:n)
     {
       if(FullStateSpace[i,j]==B+1)
       {
         WorkingVState[j]=0
       }
     }
     
     FullStateSpace[i,(n+1):(2*n)]=WorkingVState
  }
  FullStateSpace=unique(FullStateSpace)
  return(FullStateSpace)
}

CreateSVStates<-function(n,BVec,bVec)
{
  SStates=CreateSStates(n,BVec)
  VStates=CreateVStates(n,bVec)
  
  #Then for each row in SStates we all all rows in V states to be added
  
  FullStateSpace=matrix(nrow=(nrow(SStates)*nrow(VStates)),ncol=2*n)
  
  for(i in 1:nrow(FullStateSpace))
  {
    #split i into its Block and remainder
    Block=ceiling(i/nrow(VStates))
    Remainder=((i-1) %% nrow(VStates)) +1
    FullStateSpace[i,1:n]=SStates[Block,]
    
    #Before we can insert the VState, we need to check that if s=B+2, then we set v=0
    WorkingVState=VStates[Remainder,]
    for(j in 1:n)
    {
      if(FullStateSpace[i,j]==BVec[j]+1)
      {
        WorkingVState[j]=0
      }
    }
    
    FullStateSpace[i,(n+1):(2*n)]=WorkingVState
  }
  FullStateSpace=unique(FullStateSpace)
  return(FullStateSpace) 
}

# #create function for expectation of evolved state
# ListOfEvolvedState<-function(StateVector,NodeMovedTo,n,B,b)
# {
#   #for each state and action we evolve to b+1 possible states (0,...,b)
#   NewStateVector<-StateVector
#   
#   #Evolve the s state part
#   for(i in 1:n)
#   {
#     if(i==NodeMovedto)
#     {
#       NewStateVector[i]=1
#     }
#     else
#     {
#       NewStateVector[i]=min(NewStateVector[i]+1,B+1)
#     }
#   }
#   
#   NewStatesMatrix=matrix(nrow=b+1,ncol=2*n)
#   for(i in 1:(b+1))
#   {
#     NewStatesMatrix[i,]=NewStateVector
#     NewStatesMatrix[i,NodeMovedto]=i-1
#   }
#   
#   return(NewStatesMatrix)
# }
# 
# ListOfEvolvedStateIDs<-function(StateVector,NodeMovedTo,n,B,b)
# {
#   NewStatesMatrix=ListOfEvolvedState(StateVector,NodeMovedTo,n,B,b)
#   EvolvedIDs=vector(length=b+1)
#   for(i in 1:(b+1))
#   {
#     EvolvedIDs[i]=StateNumberID(NewStatesMatrix[i,])
#   }
#   
#   return(EvolvedIDs)
# }

#Evolution of S states function
NewSState<-function(CurrentSState,NodeMovedTo,BVec)
{
  NewSState=vector(length=length(CurrentSState))
  for(i in 1:length(CurrentSState))
  {
    if(i==NodeMovedTo)
    {
      NewSState[i]=1
    }
    else
    {
      NewSState[i]=min(CurrentSState[i]+1,BVec[i]+1)
    }
  }
  return(NewSState)
}

#Evolution of V States function
#note. We need to know the current S state , as if one value is B+2, then v is set to 0
NewVState<-function(CurrentVState,NewSState,NodeMovedTo,BVec,bVec,lambdaVec)
{
  #We aim to store the New V states and the probability of ending up in one.
  NewVState=matrix(nrow=(bVec[NodeMovedTo]+1),ncol=length(CurrentVState))
  NewVStateProb=vector(length=(bVec[NodeMovedTo]+1))
  
  WorkingVState=CurrentVState
  #We now need to set any v=0 if there s=B+1
  for(i in 1:length(NewSState))
  {
    if(NewSState[i]==BVec[i]+1)
    {
      WorkingVState[i]=0
    }
  }
  
  for(i in 1:(bVec[NodeMovedTo]+1))
  {
    #copy the current state
    NewVState[i,]=WorkingVState
    #Store the new evolved value
    NewVState[i,NodeMovedTo]=(i-1)
   
    #Also store the probability in a vector
    NewVStateProb[i]=TruncPoissonPMF(lambdaVec[NodeMovedTo],bVec[NodeMovedTo],i-1)
  }

  return(list(State=NewVState,Prob=NewVStateProb))
}

#NOTE. This function does not return an actual possible state, but a mean state
NewMeanVState<-function(CurrentVState,NodeMovedTo,bVec,lambdaVec)
{
  NewVState=CurrentVState
  NewVState[NodeMovedTo]=TruncPoissionMean(lambdaVec[NodeMovedTo],bVec[NodeMovedTo])
  return(NewVState)
}

#Evolution of SV States
NewSVState<-function(CurrentSVState,NodeMovedTo,BVec,bVec,lambdaVec)
{
  #Get seperate information from prior functions
  n=length(CurrentSVState)/2
  NewSState=NewSState(CurrentSVState[1:n],NodeMovedTo,BVec)
  NewVStates=NewVState(CurrentSVState[(n+1):(2*n)],NewSState,NodeMovedTo,BVec,bVec,lambdaVec)$State
  NewVStatesProbs=NewVState(CurrentSVState[(n+1):(2*n)],NewSState,NodeMovedTo,BVec,bVec,lambdaVec)$Prob
  
  #Rejoin the S and V together
  NewSVState=matrix(nrow=(bVec[NodeMovedTo]+1),ncol=(2*n))
  
  for(i in 1:(bVec[NodeMovedTo]+1))
  {
    NewSVState[i,1:n]=NewSState
    NewSVState[i,(n+1):(2*n)]=NewVStates[i,]
  }
  
  return(list(State=NewSVState,Prob=NewVStatesProbs))
}

#Identify the row which a vector comes from
IdenityRow<-function(Vec,Mat)
{
  for(i in 1:nrow(Mat))
  {
    if(identical(Vec,Mat[i,]))
    {
      return(i)
    }
  }
  print("Error has occured, the vector is not in the matrix ")
  return(-1)
}

#Action Cost C_j
CostOfActionOnNode<-function(Node,StateVector,NodeMovedTo,n,CostVec,xVec,LambdaVec)
{
  stopifnot(length(StateVector)==2*n)
  
  BVec=ceiling(xVec)
  
  if(StateVector[Node]<BVec[Node])
  {
    return(0)
  }
  else if(StateVector[Node]==BVec[Node])
  {
    if(NodeMovedTo==Node)
    {
      if(xVec[Node]>1)
      {
        return(CostVec[Node] * (LambdaVec[Node] * (BVec[Node]-xVec[Node]-1) + StateVector[n+Node]))
      }
      else if(xVec[Node] <= 1)
      {
        return(CostVec[Node] * (LambdaVec[Node] * (BVec[Node]-2*xVec[Node]) + StateVector[n+Node]))
      }
    }
    else
    {
      return(CostVec[Node] * (LambdaVec[Node] * (BVec[Node]-xVec[Node]) + StateVector[n+Node]))
    }
  }
  else if(StateVector[Node]==BVec[Node]+1)
  {
    if(NodeMovedTo==Node)
    {
      if(xVec[Node]>1)
      {
        return(0)
      }
      else if(xVec[Node] <= 1)
      {
        return(CostVec[Node] * LambdaVec[Node] * (1-xVec[Node]))
      }
    }
    else
    {
      return(CostVec[Node] * LambdaVec[Node])
    }
  }
  else
  {
    print("An error has occured in the cost of a node function")
  }
}

#Total action cost
CostOfAction<-function(StateVector,NodeMovedTo,n,CostVec,xVec,LambdaVec)
{
  Sum=0
  for(j in 1:n)
  {

    Sum=Sum+CostOfActionOnNode(j,StateVector,NodeMovedTo,n,CostVec,xVec,LambdaVec)

  }

  return(Sum)
}
  
#create for a state the constraint matrix parts for a particular state
CreateConstraintMatrixForState<-function(StateVector,AdjMatrix,n,xVec,bVec,CostVec,LambdaVec,SVStateSpace=NULL)
{
  if(is.null(SVStateSpace))
  {
    SVStateSpace=CreateSVStates(n,BVec,bVec)
  }
  
  BVec=ceiling(xVec)

  #To find out out the list of nodes we can choose to move to we need to look which s_i=1 in our state
  CurrentNode=min(StateVector[1:n]) #just look in the S state space
  CurrentNodeRow=AdjMatrix[CurrentNode,]
  
  #Now create a vector of nodes we can move to
  Actions=vector(length=0)
  for(Node in 1:n)
  {
    if(CurrentNodeRow[Node]==1)
    {
      Actions=c(Actions,Node)
    }
  }  

  #Number of variables is 1 plus all the state spaces
  NoOfVariables=1+nrow(SVStateSpace)
  
  #creating the LHS matrix part A
  ALHS=matrix(nrow=0,ncol=NoOfVariables)
  
  #storage for the RHS vector part b
  bRHS=vector(length=0)
  
  #Create the matrix row by row
  for(NodeToMoveTo in Actions)
  {
   #We will now create the constraint for taking that action
   ConstraintRow=vector(length=NoOfVariables)
   
   #new SV state and probability
   NewSVStateList=NewSVState(StateVector,NodeToMoveTo,BVec,bVec,LambdaVec)$State
   NewSVStateProbList=NewSVState(StateVector,NodeToMoveTo,BVec,bVec,LambdaVec)$Prob
   
   #Now we need the identity of such a state
   IDNewSVStateVector=vector(length=nrow(NewSVStateList))
   for(i in 1:nrow(NewSVStateList))
   {
     IDNewSVStateVector[i]=IdenityRow(NewSVStateList[i,],SVStateSpace)+1 #Note we add one as we have the variable g
   }
   
   #now we can create the Constraint row
   ConstraintRow[1]=1 #as the g is always part of it
   ConstraintRow[IdenityRow(StateVector,SVStateSpace)+1]=1 #as we have one lot of h(current state)
   
   for(i in 1:length(IDNewSVStateVector))
   {
     ConstraintRow[IDNewSVStateVector[i]]=ConstraintRow[IDNewSVStateVector[i]]-NewSVStateProbList[i]
   }
   
   ALHS=rbind(ALHS,ConstraintRow)
   bRHS=c(bRHS,CostOfAction(StateVector,NodeToMoveTo,n,CostVec,xVec,LambdaVec))
  }
  return(list(LHS=ALHS,RHS=bRHS))
}

CreateConstraintMatrix<-function(AdjMatrix,n,xVec,bVec,CostVec,LambdaVec)
{
  BVec=ceiling(xVec)
  #We now repeat for each possible state the construction and bind together
  print("Constructing SV State Space")
  SVStateSpace=CreateSVStates(n,BVec,bVec)
  print("Constructed")
  
  NoOfVariables=1+nrow(SVStateSpace)
  
  FullConstraintMatrix=matrix(nrow=0,ncol=NoOfVariables)
  Fullbounds=vector(length=0)
  
  for(i in 1:nrow(SVStateSpace))
  {
    print("On state")
    print(SVStateSpace[i,])

    
    #For each state run prior construction and retreive
    FullConstraintMatrix=rbind(FullConstraintMatrix,CreateConstraintMatrixForState(SVStateSpace[i,],AdjMatrix,n,xVec,bVec,CostVec,LambdaVec,SVStateSpace)$LHS)
    Fullbounds=c(Fullbounds,CreateConstraintMatrixForState(SVStateSpace[i,],AdjMatrix,n,xVec,bVec,CostVec,LambdaVec,SVStateSpace)$RHS)
    
  }
  return(list(MatrixConstraints=FullConstraintMatrix,VectorBounds=Fullbounds))
}

SolveLP<-function(AdjMatrix,n,xVec,bVec,CostVec,LambdaVec)
{
  CreatedAb=CreateConstraintMatrix(AdjMatrix,n,xVec,bVec,CostVec,LambdaVec)
  A=CreatedAb$MatrixConstraints
  b=CreatedAb$VectorBounds
  
  
  Objdir="max"
  Objective=c(1,rep(0,(ncol(A)-1)))
  Constdir=rep("<=",nrow(A))
  
  Solved=lp(Objdir,Objective,A,Constdir,b)
  return(list(Value=Solved ,Solution=Solved$solution))
}

