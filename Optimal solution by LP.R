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
       if(FullStateSpace[i,j]==(B+1))
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
      if(FullStateSpace[i,j]==(BVec[j]+1))
      {
        WorkingVState[j]=0
      }
    }
    
    FullStateSpace[i,(n+1):(2*n)]=WorkingVState
  }
  FullStateSpace=unique(FullStateSpace)
  return(FullStateSpace) 
}

#This function is a variation on the above with the possibility of having an unkown transitionary state.
CreateTransSVStates<-function(n,BVec,bVec)
{
  
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
#note. We need to know the current S state , as if one value is B+1, then v is set to 0
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
NewMeanVState<-function(CurrentVState,NewSState,NodeMovedTo,BVec,bVec,lambdaVec)
{
  NewV=(NewVState(CurrentVState,NewSState,NodeMovedTo,BVec,bVec,lambdaVec)$State)[1,]
  NewV[NodeMovedTo]=TruncPoissionMean(lambdaVec[NodeMovedTo],bVec[NodeMovedTo])
  return(NewV)
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

#New state if current state is neglected
NewSVStateIfNeglect<-function(CurrentSVState,BVec,bVec,lambdaVec)
{
  #Note. We cannot simply use the prior functions
  n=length(CurrentSVState)/2
  sVec=CurrentSVState[1:n]
  vVec=CurrentSVState[(n+1):(2*n)]
  
  #Now we evolve sVec , by increasing it to cap, if its capped we removed the observed info
  for(i in 1:n)
  {
    if(sVec[i] < BVec[i])
    {
      sVec[i]=sVec[i]+1
    }
    else if(sVec[i]==BVec[i])
    {
      sVec[i]=sVec[i]+1
      vVec[i]=0
    }
    else
    {
      vVec[i]=0
    }
  }
  return(c(sVec,vVec))
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
  RVec=BVec-xVec
  
  if(StateVector[Node]<BVec[Node])
  {
    return(0)
  }
  else if(StateVector[Node]==BVec[Node])
  {
    if(NodeMovedTo==Node)
    {
      # if(xVec[Node]>1)
      # {
      #   return(CostVec[Node] * (LambdaVec[Node] * (BVec[Node]-xVec[Node]) + StateVector[n+Node]))
      # }
      # else if(xVec[Node] <= 1)
      # {
      #   return(CostVec[Node] * StateVector[n+Node])
      # }
      return(0)
    }
    else
    {
      return(CostVec[Node] * ((LambdaVec[Node] * RVec[Node]) + StateVector[n+Node]))
    }
  }
  else if(StateVector[Node]==BVec[Node]+1)
  {
    if(NodeMovedTo==Node)
    {
      # if(xVec[Node]>1)
      # {
      #   return(CostVec[Node] * LambdaVec[Node])
      # }
      # else if(xVec[Node] <= 1)
      # {
      #   return(CostVec[Node] * LambdaVec[Node] * xVec[Node])
      # }
      return(0)
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

#Cost if the state is neglected
CostToNeglect<-function(StateVector,n,CostVec,xVec,LambdaVec)
{
  
  #We can use prior formula but with no node every visited
  return(CostOfAction(StateVector,-1,n,CostVec,xVec,LambdaVec))
}


SeperatedCostOfActionOnNode<-function(Node,StateVector,NodeMovedTo,n,CostVec,xVec,LambdaVec)
{
  stopifnot(length(StateVector)==2*n)
  
  BVec=ceiling(xVec)
  RVec=BVec-xVec
  
  if(StateVector[Node]<BVec[Node])
  {
    return(list(CostDueToArrivals=0,CostDueToObserved=0))
  }
  else if(StateVector[Node]==BVec[Node])
  {
    if(NodeMovedTo==Node)
    {
      return(list(CostDueToArrivals=0,CostDueToObserved=0))
    }
    else
    {
      return(list(CostDueToArrivals=CostVec[Node] * LambdaVec[Node] * RVec[Node],CostDueToObserved=CostVec[Node] *StateVector[n+Node]))
    }
  }
  else if(StateVector[Node]==BVec[Node]+1)
  {
    if(NodeMovedTo==Node)
    {
      return(list(CostDueToArrivals=0,CostDueToObserved=0))
    }
    else
    {
      return(list(CostDueToArrivals=CostVec[Node] * LambdaVec[Node],CostDueToObserved=0))
    }
  }
  else
  {
    print("An error has occured in the cost of a node function")
  }
}
 
SeperatedCostOfAction<-function(StateVector,NodeMovedTo,n,CostVec,xVec,LambdaVec)
{
  SumDueToArrivals=0
  SumDueToObserved=0
  for(j in 1:n)
  {
    SumDueToArrivals=SumDueToArrivals+SeperatedCostOfActionOnNode(j,StateVector,NodeMovedTo,n,CostVec,xVec,LambdaVec)$CostDueToArrivals
    SumDueToObserved=SumDueToObserved+SeperatedCostOfActionOnNode(j,StateVector,NodeMovedTo,n,CostVec,xVec,LambdaVec)$CostDueToObserved
  }
  
  return(list(CostDueToArrivals=SumDueToArrivals,CostDueToObserved=SumDueToObserved))
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
  CurrentNode=which.min(StateVector[1:n]) #just look in the S state space
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
  print(SVStateSpace)
  
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
  return(list(MatrixConstraints=FullConstraintMatrix,VectorBounds=Fullbounds,SVStateSpace=SVStateSpace))
}

SolveLP<-function(AdjMatrix,n,xVec,bVec,CostVec,LambdaVec)
{
  CreatedAb=CreateConstraintMatrix(AdjMatrix,n,xVec,bVec,CostVec,LambdaVec)
  A=CreatedAb$MatrixConstraints
  b=CreatedAb$VectorBounds
  
  print(A)
  print(b)
  
  Objdir="max"
  Objective=c(1,rep(0,(ncol(A)-1)))
  Constdir=rep("<=",nrow(A))

  
  print("Starting to solve")
  Solved=lp(Objdir,Objective,A,Constdir,b)
  return(list(Value=Solved ,Solution=Solved$solution))
}

ActiveConstraints<-function(AdjMatrix,n,xVec,bVec,CostVec,LambdaVec)
{
  #Retrive the solution to the LP, Solve for the solution
  
  CreatedAb=CreateConstraintMatrix(AdjMatrix,n,xVec,bVec,CostVec,LambdaVec)
  A=CreatedAb$MatrixConstraints
  b=CreatedAb$VectorBounds
  SVStateSpace=CreatedAb$SVStateSpace
  
  print(A)
  print(b)
  
  Objdir="max"
  Objective=c(1,rep(0,(ncol(A)-1)))
  Constdir=rep("<=",nrow(A))
  
  
  print("Starting to solve")
  Solved=lp(Objdir,Objective,A,Constdir,b)
  Solution=Solved$solution
  print("Solution is")
  print(Solution)
  
  #Create the solutions b values
  Solutionsb=as.vector(A %*% Solution)
  
  print("calculated solution")
  print(Solutionsb)
  print("b")
  print(b)
  
  #Compare values to see which b are active
  ActiveConstraints=vector(length=0)
  StateBestAction=matrix(data=list(),nrow=nrow(SVStateSpace),ncol=2)
  # counter=1
  # for(i in 1:length(b))
  # {
  #   if(abs(Solutionsb[i]-b[i])<0.001)
  #   {
  #     print(paste("An active constraint exists at constraint",toString(i)))
  #     #We record that the constraint is active
  #     ActiveConstraints=c(ActiveConstraints,i)
  #     #We want to find the corresponding state and action for this active constraint
  #     #Note. Each state has n actions available, so we use blocksize n to get the state
  #     ActionIndex=i %% n
  #     if(ActionIndex==0)
  #     {
  #       Action=n
  #     }
  #     else
  #     {
  #       Action=ActionIndex
  #     }
  #     
  #     StateBlock=1+(i-Action)/n
  #     State=SVStateSpace[StateBlock,]
  #     print(State)
  #     
  #     #Store the state and action in the data frame
  #     StateBestAction[[counter,1]]=toString(State)
  #     StateBestAction[[counter,2]]=Action
  #     counter=counter+1
  #   }
  # }
  
  #for each state block we will look at which actions make the constraints active
  counter=1
  for(statenumber in 1:nrow(SVStateSpace))
  {
    #First from statenumber we need to identify the current node (for actions)
    CurrentState=SVStateSpace[statenumber,]
    CurrentNode=min(CurrentState[1:n])
    StateBestAction[[statenumber,1]]=CurrentState
    #For each possible action see if the constraint is active and record if so
    ActionsActive=vector(length=0)
    for(action in 1:ncol(AdjMatrix))
    {
      
      if(AdjMatrix[CurrentNode,action]==1)
      {
        #Now check if this action constraint is active
        if(abs(Solutionsb[counter]-b[counter])<0.001)
        {
          ActionsActive=c(ActionsActive,action)
        }
        #Increment the counter if we checked to see if it was active or not
        counter=counter+1
      }
    }
    StateBestAction[[statenumber,2]]=ActionsActive
    
  }

  print(SVStateSpace)
  return(StateBestAction)
}

#We solve the LP , then insert the value for g and then solve to maximize the sum of the h's
ActiveConstraintsExperimental<-function(AdjMatrix,n,xVec,bVec,CostVec,LambdaVec)
{
  #Retrive the solution to the LP, Solve for the solution
  
  CreatedAb=CreateConstraintMatrix(AdjMatrix,n,xVec,bVec,CostVec,LambdaVec)
  A=CreatedAb$MatrixConstraints
  b=CreatedAb$VectorBounds
  SVStateSpace=CreatedAb$SVStateSpace
  
  print(A)
  print(b)
  
  Objdir="max"
  Objective=c(1,rep(0,(ncol(A)-1)))
  Constdir=rep("<=",nrow(A))
  
  
  print("Starting to solve")
  Solved=lp(Objdir,Objective,A,Constdir,b)
  InitialSolution=Solved$solution
  print("Inital Solution is")
  print(InitialSolution)
  
  #Create the solutions b values
  Solutionsb=as.vector(A %*% InitialSolution)
  
  print("calculated solution")
  print(Solutionsb)
  print("b")
  print(b)
  
  #Compare values to see which b are active
  ActiveConstraints=vector(length=0)
  StateBestAction=matrix(data=list(),nrow=nrow(SVStateSpace),ncol=2)
  
  #for each state block we will look at which actions make the constraints active
  counter=1
  NoActions=vector(length=0)
  for(statenumber in 1:nrow(SVStateSpace))
  {
    #First from statenumber we need to identify the current node (for actions)
    CurrentState=SVStateSpace[statenumber,]
    CurrentNode=min(CurrentState[1:n])
    StateBestAction[[statenumber,1]]=CurrentState
    #For each possible action see if the constraint is active and record if so
    ActionsActive=vector(length=0)
    for(action in 1:ncol(AdjMatrix))
    {
      
      if(AdjMatrix[CurrentNode,action]==1)
      {
        #Now check if this action constraint is active
        if(abs(Solutionsb[counter]-b[counter])<0.001)
        {
          ActionsActive=c(ActionsActive,action)
        }
        #Increment the counter if we checked to see if it was active or not
        counter=counter+1
      }
    }
    if(length(ActionsActive)==0)
    {
      NoActions=c(NoActions,statenumber)
    }
    StateBestAction[[statenumber,2]]=ActionsActive
    
  }
  print(StateBestAction)
  
  #We now seek to get Active constrains for the current states which are not assigned an action
  HasActions=1:nrow(SVStateSpace)
  HasActions=HasActions[-NoActions]
  
  #Now we extract g
  g=InitialSolution[1]
  
  #Now we form an LP with this value in the constraints and try to maximize the sum of h's
  NewA=A[,-1] #remove the first column of A to remove g
  NumberofOriginalConstraints=nrow(NewA)
  Newb=b-g #Subtract g from the LHS
  for(i in HasActions)
  {
    #For each currently assigned value, we force the value to take the current solutions value
    NewA=rbind(NewA,c(rep(0,i-1),1,rep(0,ncol(NewA)-i)))
    Newb=c(Newb,InitialSolution[i+1])
  }
  
  print("New A is")
  print(NewA)

  print("New b is")
  print(Newb)
  
  NewObjdir="max"
  NewObjective=rep(1,ncol(NewA))
  #NewObjective=c(0,0,10000,rep(0,ncol(NewA)-3))
  NewConstdir=c(rep("<=",NumberofOriginalConstraints),rep("=",length(HasActions)))
  
  print("Starting to solve new problem for h's")
  HSolved=lp(NewObjdir,NewObjective,NewA,NewConstdir,Newb)
  HSolution=HSolved$solution
  print("Solution for h's is")
  print(HSolution)
  
  #Now combine for full solution
  Solution=c(g,HSolution)
  
  
  
  #Create the solutions b values
  Solutionsb=as.vector(A %*% Solution)
  
  print("calculated solution")
  print(Solutionsb)
  print("b")
  print(b)
  
  #Compare values to see which b are active
  ActiveConstraints=vector(length=0)
  StateBestAction=matrix(data=list(),nrow=nrow(SVStateSpace),ncol=2)
  # counter=1
  # for(i in 1:length(b))
  # {
  #   if(abs(Solutionsb[i]-b[i])<0.001)
  #   {
  #     print(paste("An active constraint exists at constraint",toString(i)))
  #     #We record that the constraint is active
  #     ActiveConstraints=c(ActiveConstraints,i)
  #     #We want to find the corresponding state and action for this active constraint
  #     #Note. Each state has n actions available, so we use blocksize n to get the state
  #     ActionIndex=i %% n
  #     if(ActionIndex==0)
  #     {
  #       Action=n
  #     }
  #     else
  #     {
  #       Action=ActionIndex
  #     }
  #     
  #     StateBlock=1+(i-Action)/n
  #     State=SVStateSpace[StateBlock,]
  #     print(State)
  #     
  #     #Store the state and action in the data frame
  #     StateBestAction[[counter,1]]=toString(State)
  #     StateBestAction[[counter,2]]=Action
  #     counter=counter+1
  #   }
  # }
  
  #for each state block we will look at which actions make the constraints active
  counter=1
  for(statenumber in 1:nrow(SVStateSpace))
  {
    #First from statenumber we need to identify the current node (for actions)
    CurrentState=SVStateSpace[statenumber,]
    CurrentNode=min(CurrentState[1:n])
    StateBestAction[[statenumber,1]]=CurrentState
    #For each possible action see if the constraint is active and record if so
    ActionsActive=vector(length=0)
    for(action in 1:ncol(AdjMatrix))
    {
      
      if(AdjMatrix[CurrentNode,action]==1)
      {
        #Now check if this action constraint is active
        if(abs(Solutionsb[counter]-b[counter])<0.001)
        {
          ActionsActive=c(ActionsActive,action)
        }
        #Increment the counter if we checked to see if it was active or not
        counter=counter+1
      }
    }
    StateBestAction[[statenumber,2]]=ActionsActive
    
  }
  
  #print(SVStateSpace)
  return(StateBestAction)
}


#We Note our list of x,y are in the form of blocks; x->state->action
CreateDualSetup<-function(AdjMatrix,n,xVec,bVec,CostVec,LambdaVec,SVStateSpace=NULL,AlphaVec=NULL)
{
  
  BVec=ceiling(xVec)
  if(is.null(SVStateSpace))
  {
    SVStateSpace=CreateSVStates(n,BVec,bVec)
  }
  if(is.null(AlphaVec))
  {
    AlphaVec=vector(length=nrow(SVStateSpace))
    AlphaVec=rep(1/length(AlphaVec),length(AlphaVec))
  }
  
  #We now work out the number of x's and y's needed (i.e the size of (s,a))
  NumberOfXVariables=0
  NumberOfActionsFromState=vector(length=0)
  for(state in 1:nrow(SVStateSpace))
  {
    #For each state we need to work out the number of actions
    #First we work out the current node (and then how many actions can be taken)
    StateVector=SVStateSpace[state,]
    CurrentNode=which.min(StateVector[1:n])
    #print(paste("Current Node is ",toString(CurrentNode)))
    
    AdjRow=AdjMatrix[CurrentNode,]
  
    NumberOfXVariables=NumberOfXVariables+length(AdjRow[AdjRow==1])
    NumberOfActionsFromState=c(NumberOfActionsFromState,length(AdjRow[AdjRow==1]))
    #print(NumberOfActionsFromState)
  }

  NumberOfYVariables=NumberOfXVariables
  NumberOfVariables=NumberOfXVariables+NumberOfYVariables
  
  
  
  #creating the LHS matrix part A
  ALHS=matrix(nrow=0,ncol=NumberOfVariables)
  
  #storage for the RHS vector part b
  bRHS=vector(length=0)
  
  #FIRST DUAL CONSTRAINT
  for(StateNumber in 1:nrow(SVStateSpace))
  {
    
    #We will now create the constraint for taking that action
    ConstraintRow=vector(length=NumberOfVariables)
    
    #We store our state we are working on the constraint row for
    State=SVStateSpace[StateNumber,]
    CurrentNode=which.min(State[1:n])
    CurrentNodeRow=AdjMatrix[CurrentNode,]
    
    #Actions available from this state
    Actions=vector(length=0)
    for(Node in 1:n)
    {
      if(CurrentNodeRow[Node]==1)
      {
        Actions=c(Actions,Node)
      }
    }  
    
    #We place 1's in the all the actions possible for this state x
    for(i in 1:NumberOfActionsFromState[StateNumber])
    {
      if(StateNumber==1)
      {
        ConstraintRow[i]=ConstraintRow[i]+1
      }
      else
      {
        ConstraintRow[sum(NumberOfActionsFromState[1:(StateNumber-1)])+i]=
          ConstraintRow[sum(NumberOfActionsFromState[1:(StateNumber-1)])+i]+1
      }
    }
    
    #We now subtract prob(moving to state summed over all states and actions) for x
    for(oldstatenumber in 1:nrow(SVStateSpace))
    {
      OldState=SVStateSpace[oldstatenumber,]
      OldNode=which.min(OldState[1:n])
      OldNodeRow=AdjMatrix[OldNode,]
      #Actions available from this state
      OldActions=vector(length=0)
      for(Node in 1:n)
      {
        if(OldNodeRow[Node]==1)
        {
          OldActions=c(OldActions,Node)
        }
      }  
      
      for(actionnumber in 1:NumberOfActionsFromState[oldstatenumber])
      {
        #Using this oldstate and action , is it possible that the new state is the current working state
        New=NewSVState(OldState,OldActions[actionnumber],BVec,bVec,LambdaVec)
        NewState=New$State
        NewStateProb=New$Prob
        
        for(NewStateNumber in 1:nrow(NewState))
        {
          if(all(NewState[NewStateNumber,]==SVStateSpace[StateNumber,]))
          {
            #If we can possibly move to the state we are working on we subtract the probability
            if(oldstatenumber==1)
            {
             ConstraintRow[actionnumber]=ConstraintRow[actionnumber]-NewStateProb[NewStateNumber]
            }
            else
            {
             ConstraintRow[sum(NumberOfActionsFromState[1:(oldstatenumber-1)])+actionnumber]=
               ConstraintRow[sum(NumberOfActionsFromState[1:(oldstatenumber-1)])+actionnumber]-NewStateProb[NewStateNumber]
            }
          }
        }
        
        
        

      }
    }
    
    ALHS=rbind(ALHS,ConstraintRow)
    bRHS=c(bRHS,0)
  }
  
  #SECOND DUAL CONSTRAINT
  for(StateNumber in 1:nrow(SVStateSpace))
  {
    #We will now create the constraint for taking that action
    ConstraintRow=vector(length=NumberOfVariables)
    
    #We store our state we are working on the constraint row for
    State=SVStateSpace[StateNumber,]
    CurrentNode=which.min(State[1:n])
    CurrentNodeRow=AdjMatrix[CurrentNode,]
    
    #Actions available from this state
    Actions=vector(length=0)
    for(Node in 1:n)
    {
      if(CurrentNodeRow[Node]==1)
      {
        Actions=c(Actions,Node)
      }
    }  
    
    #We place 1's in the all the actions possible for this state for x
    for(i in 1:NumberOfActionsFromState[StateNumber])
    {
      if(StateNumber==1)
      {
        ConstraintRow[i]=ConstraintRow[i]+1
      }
      else
      {
        ConstraintRow[sum(NumberOfActionsFromState[1:(StateNumber-1)])+i]=
          ConstraintRow[sum(NumberOfActionsFromState[1:(StateNumber-1)])+i]+1
      }
    }
    
    #We place 1's in the all the actions possible for this state for y
    for(i in 1:NumberOfActionsFromState[StateNumber])
    {
      if(StateNumber==1)
      {
        ConstraintRow[NumberOfXVariables+i]=ConstraintRow[NumberOfXVariables+i]+1
      }
      else
      {
        ConstraintRow[NumberOfXVariables+sum(NumberOfActionsFromState[1:(StateNumber-1)])+i]=
          ConstraintRow[NumberOfXVariables+sum(NumberOfActionsFromState[1:(StateNumber-1)])+i]+1
      }
    }
    
    #We now subtract prob(moving to state summed over all states and actions) for y
    for(oldstatenumber in 1:nrow(SVStateSpace))
    {
      OldState=SVStateSpace[oldstatenumber,]
      OldNode=which.min(OldState[1:n])
      OldNodeRow=AdjMatrix[OldNode,]
      #Actions available from this state
      OldActions=vector(length=0)
      for(Node in 1:n)
      {
        if(OldNodeRow[Node]==1)
        {
          OldActions=c(OldActions,Node)
        }
      }  
      for(actionnumber in 1:NumberOfActionsFromState[oldstatenumber])
      {
        #Using this oldstate and action , is it possible that the new state is the current working state
        
        New=NewSVState(OldState,OldActions[actionnumber],BVec,bVec,LambdaVec)
        NewState=New$State
        NewStateProb=New$Prob
        
        for(NewStateNumber in 1:nrow(NewState))
        {
          if(all(NewState[NewStateNumber,]==SVStateSpace[StateNumber,]))
          {
            #If we can possibly move to the state we are working on we subtract the probability
            if(oldstatenumber==1)
            {
              ConstraintRow[NumberOfXVariables+actionnumber]=ConstraintRow[NumberOfXVariables+actionnumber]-NewStateProb[NewStateNumber]
            }
            else
            {
              ConstraintRow[NumberOfXVariables+sum(NumberOfActionsFromState[1:(oldstatenumber-1)])+actionnumber]=
                ConstraintRow[NumberOfXVariables+sum(NumberOfActionsFromState[1:(oldstatenumber-1)])+actionnumber]-NewStateProb[NewStateNumber]
            }
          }
        }
        
        
        
        
      }
    }
    
    
    
    
    ALHS=rbind(ALHS,ConstraintRow)
    bRHS=c(bRHS,AlphaVec[StateNumber])
  }
  
  Objective=vector(length=NumberOfVariables)
  #We now create the Objective
  for(StateNumber in 1:nrow(SVStateSpace))
  {
    #We store our state we are working on the constraint row for
    State=SVStateSpace[StateNumber,]
    CurrentNode=which.min(State[1:n])
    CurrentNodeRow=AdjMatrix[CurrentNode,]
    
    #Actions available from this state
    Actions=vector(length=0)
    for(Node in 1:n)
    {
      if(CurrentNodeRow[Node]==1)
      {
        Actions=c(Actions,Node)
      }
    }  
    
    for(ActionNumber in 1:NumberOfActionsFromState[StateNumber])
    {
      #For this state and this action store the cost in the objective function
      if(StateNumber==1)
      {
        Objective[ActionNumber]=CostOfAction(State,Actions[ActionNumber],n,CostVec,xVec,LambdaVec)
      }
      else
      {
        Objective[sum(NumberOfActionsFromState[1:StateNumber-1])+ActionNumber]=CostOfAction(State,Actions[ActionNumber],n,CostVec,xVec,LambdaVec)
      }
    }
  }
  
  return(list(Objective=Objective,MatrixConstraints=ALHS,VectorBounds=bRHS,StateSpace=SVStateSpace,NumberOfActionsFromState=NumberOfActionsFromState,AlphaVec=AlphaVec))
}
  
SolveDualLP<-function(AdjMatrix,n,xVec,bVec,CostVec,LambdaVec)
{
  print("Starting to Set up the dual problem")
  CreatedDual=CreateDualSetup(AdjMatrix,n,xVec,bVec,CostVec,LambdaVec)
  A=CreatedDual$MatrixConstraints
  b=CreatedDual$VectorBounds
  StateSpace=CreatedDual$StateSpace
  
  print(A)
  print(b)
  
  Objdir="min"
  Objective=CreatedDual$Objective
  
  
  Constdir=rep("=",nrow(A))
  
  
  print("Starting to solve")
  Solved=lp(Objdir,Objective,A,Constdir,b)
  return(list(Value=Solved$objval ,Solution=Solved$solution,StateSpace=StateSpace))
}

OptimalDualDesicionPolicy<-function(AdjMatrix,n,xVec,bVec,CostVec,LambdaVec)
{
  #Solve the Dual LP
  CreatedDual=CreateDualSetup(AdjMatrix,n,xVec,bVec,CostVec,LambdaVec)
  A=CreatedDual$MatrixConstraints
  b=CreatedDual$VectorBounds
  SVStateSpace=CreatedDual$StateSpace
  NumberOfActionsFromState=CreatedDual$NumberOfActionsFromState
  
  print(A)
  print(b)
  
  Objdir="min"
  Objective=CreatedDual$Objective
  
  
  Constdir=rep("=",nrow(A))
  
  
  print("Starting to solve")
  Solved=lp(Objdir,Objective,A,Constdir,b)
  Solution=Solved$solution
  Value=Solved$objval
  print("Solution found")
  print(Solution)
  print("For objective value")
  print(Value)
  
  #We now split the solution into x's and y's
  NumberOfVariables=length(Solution)/2
  OptimalX=Solution[1:NumberOfVariables]
  print("Optimal x's are")
  print(OptimalX)
  OptimalY=Solution[(NumberOfVariables+1):length(Solution)]
  print("Optimal y's are")
  print(OptimalY)
  
  #From the optimal x/y's we create a decision rule
  
  #First identify Recurrent states and transient states
  #Record 1 if recurrent, 0 if transient
  StateType=vector(length=nrow(SVStateSpace))
  #OptimalDecision=list(length=nrow(SVStateSpace))
  OptimalDecision=list()
  for(StateNumber in 1:nrow(SVStateSpace))
  {
    State=SVStateSpace[StateNumber,]
    CurrentNode=which.min(State[1:n])
    CurrentActions=AdjMatrix[CurrentNode,]
    
    NodesCanMoveTo=which(CurrentActions==1)
    
    #Work out the sum of the x's
    if(StateNumber==1)
    {
      SumOfX=sum(OptimalX[1:NumberOfActionsFromState[1]])
    }
    else
    {
      SumOfX=sum(OptimalX[(sum(NumberOfActionsFromState[1:(StateNumber-1)])+1):sum(NumberOfActionsFromState[1:StateNumber])])
    }
    print(SumOfX)
    if(SumOfX>0)
    {
      #print("We have a recurrent State")
      #If the sum is postive then this state is recurrent and put into S_{x}
      StateType[StateNumber]=1
    
      StatesOptimalDecisions=vector(length=0)
      #We then run through and record all x>0
      for(ActionNumber in 1:NumberOfActionsFromState[StateNumber])
      {
       if(StateNumber==1)
       {
        if(OptimalX[ActionNumber]>0)
        {
          #Record it as an optimal decision
          StatesOptimalDecisions=c(StatesOptimalDecisions,NodesCanMoveTo[ActionNumber])
        }
       }
       else
       {
         if(OptimalX[sum(NumberOfActionsFromState[1:(StateNumber-1)])+ActionNumber]>0)
         {
           StatesOptimalDecisions=c(StatesOptimalDecisions,NodesCanMoveTo[ActionNumber])
         }
       }
      }
      OptimalDecision[StateNumber]=StatesOptimalDecisions 
    }
    else
    {
      #print("We have a transient State")
      #For transient states we look at y's
      
      #If the sum is postive then this state is recurrent and put into S_{x}
      StateType[StateNumber]=0
      
      StatesOptimalDecisions=vector(length=0)
      #We then run through and record all x>0
      for(ActionNumber in 1:NumberOfActionsFromState[StateNumber])
      {
        if(StateNumber==1)
        {
          if(OptimalY[ActionNumber]>0)
          {
            #Record it as an optimal decision
            StatesOptimalDecisions=c(StatesOptimalDecisions,NodesCanMoveTo[ActionNumber])
          }
        }
        else
        {
          if(OptimalY[sum(NumberOfActionsFromState[1:(StateNumber-1)])+ActionNumber]>0)
          {
            StatesOptimalDecisions=c(StatesOptimalDecisions,NodesCanMoveTo[ActionNumber])
          }
        }
      }
      OptimalDecision[StateNumber]=StatesOptimalDecisions 
    }
    
  }
  
  return(list(OptimalValue=Value,OptimalDecision=OptimalDecision,StateSpace=SVStateSpace,AlphaVec=CreatedDual$AlphaVec))
  
}

FindOptimalEquilibriumValuesByDual<-function(AdjMatrix,n,xVec,bVec,CostVec,LambdaVec,AlphaVec=NULL)
{
  BVec=ceiling(xVec)
  #We first solve the dual problem to find optimal actions
  print("Starting to solve Dual LP")
  Solved=OptimalDualDesicionPolicy(AdjMatrix,n,xVec,bVec,CostVec,LambdaVec)
  print("Dual LP solved")
  OptimalValue=Solved$OptimalValue
  OptimalDecision=Solved$OptimalDecision
  SVStateSpace=Solved$StateSpace
  if(is.null(AlphaVec))
  {
    AlphaVec=Solved$AlphaVec
  }
  
  
  #Know we look at creating a matrix of simultaneous equations for those actions under primal problem
  #We have 2xNumStates (for g's and y's)
  NumStates=nrow(SVStateSpace)
  ConstraintsMatrix=matrix(nrow=0,ncol=(2*NumStates))
  VectorValue=vector(length=0)
  #Form Type 1 constraints
  for(StateNumber in 1:NumStates)
  {
    State=SVStateSpace[StateNumber,]
    #For each state we get two constraints
    ConstraintRow=vector(length=(2*NumStates))
    
    #For type one constraint we get g(s)-\sum p(j|s,a)g(j)=0
    
    #Place a one in this states position
    ConstraintRow[StateNumber]=ConstraintRow[StateNumber]+1
    
    #Do the subtraction of proabilities
    NodeMovedTo=OptimalDecision[[StateNumber]][1] #Note if more than one optimal decision is stored we will use the first one
    New=NewSVState(State,NodeMovedTo,BVec,bVec,LambdaVec)
    NewStates=New$State
    NewStateProbs=New$Prob
    #print(NewStates)
    
    #For each of the new states we subtract the probability of ending up there
    for(i in 1:nrow(NewStates))
    {
      NewStateID=IdenityRow(NewStates[i,],SVStateSpace)
      ConstraintRow[NewStateID]=ConstraintRow[NewStateID]-NewStateProbs[i]
    }
    ConstraintsMatrix=rbind(ConstraintsMatrix,ConstraintRow)
    VectorValue=c(VectorValue,0)
  }
  #Form Type 2 Constraints
  for(StateNumber in 1:NumStates)
  {
    State=SVStateSpace[StateNumber,]
    #For each state we get two constraints
    ConstraintRow=vector(length=(2*NumStates))
    
    #For type one constraint we get g(s)+h(s)-\sum p(j|s,a)h(j)=c(s,a)
    
    #Place a one in this states position
    ConstraintRow[StateNumber]=ConstraintRow[StateNumber]+1
    ConstraintRow[NumStates+StateNumber]=ConstraintRow[NumStates+StateNumber]+1
    
    #Do the subtraction of proabilities
    NodeMovedTo=OptimalDecision[[StateNumber]][1] #Note if more than one optimal decision is stored we will use the first one
    New=NewSVState(State,NodeMovedTo,BVec,bVec,LambdaVec)
    NewStates=New$State
    NewStateProbs=New$Prob
    
    #For each of the new states we subtract the probability of ending up there
    for(i in 1:nrow(NewStates))
    {
      NewStateID=IdenityRow(NewStates[i,],SVStateSpace)
      ConstraintRow[NumStates+NewStateID]=ConstraintRow[NumStates+NewStateID]-NewStateProbs[i]
    }
    ConstraintsMatrix=rbind(ConstraintsMatrix,ConstraintRow)
    VectorValue=c(VectorValue,CostOfAction(State,NodeMovedTo,n,CostVec,xVec,LambdaVec))
  }
  
  print("Constraints are")
  print(ConstraintsMatrix)
  print("vector values are")
  print(VectorValue)
  
  #We now attempt to solve the problem (Still using LP solved with equal constraints)
  
  Objdir="max"
  Objective=c(AlphaVec,rep(0,NumStates))
  
  
  Constdir=rep("=",nrow(ConstraintsMatrix))
  
  
  print("Starting to solve")
  Solved=lp(Objdir,Objective,ConstraintsMatrix,Constdir,VectorValue)
  OptimaSol=Solved$solution
  
  return(list(OptimalEquilibrium=OptimaSol))
  
}



#Function to compare two policies (note. Assumed in same order as SVStateSpace)
ComparePolicies<-function(Policy1,Policy2,StateSpace=NULL)
{
  #We assume the Policies have the format of a vector of lists
  AgreeAt=vector(length=length(Policy1))
  for(i in 1:length(Policy1))
  {
    #print(Policy1[[i]])
    #print(Policy2[[i]])
    if(all(Policy1[[i]]==Policy2[[i]]))
    {
      AgreeAt[i]=1
    }
    else
    {
      print(paste("Policies disagree at state ",toString(StateSpace[i,])))
      print(paste("Policy 1 suggests: ",toString(Policy1[[i]]),"Policy 2 suggest: ",toString(Policy2[[i]])))
    }
  }

  return(AgreeAt)
}