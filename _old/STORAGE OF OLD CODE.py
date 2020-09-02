# -*- coding: utf-8 -*-
"""
Created on Tue Sep  1 15:01:54 2020

@author: eugen
"""

'''
def myloglik_grad(beta, ydata, nsamp, A, sens, spec):
    betaI = beta[0:A.shape[1]]
    betaJ = beta[A.shape[1]:]
    probs = ((1-invlogit(betaJ)) * np.squeeze(np.array(A @ invlogit(betaI)))  + invlogit(betaJ))
    probs_dirJ = -invlogit_grad(betaJ) * np.squeeze(np.array(A @ invlogit(betaI))) + invlogit_grad(betaJ)
    ilJ = np.array([invlogit(betaJ).T,]*(betaI.shape[0])).T
    ilK = np.array([invlogit_grad(betaI),]*(betaJ.shape[0]))
    probs_dirI =  np.multiply(np.multiply((1- ilJ ) , A) , np.squeeze(ilK)).T 
    probz_dirI = probs_dirI*sens - (probs_dirI) * (1-spec)
    probz_dirJ = probs_dirJ*sens - (probs_dirJ) * (1-spec)
    probsz = probs*sens + (1-probs) * (1-spec)
    negloglikI =  np.array(probz_dirI @ (-(np.asarray(ydata) / probsz  - (np.asarray(nsamp)-np.asarray(ydata)) / (1-probsz))).T)
    negloglikJ =  np.squeeze(np.array((probz_dirJ * (-(np.asarray(ydata) / probsz  - (np.asarray(nsamp)-np.asarray(ydata)) / (1-probsz)))).T))
    return -np.squeeze(np.concatenate([negloglikI,negloglikJ]))

def Est_LinearProjection(A,X): # Linear Projection
    # Uses the (estimated) transition matrix, A, and the (estimated) percentage SF
    # at each end node, X
    intProj = np.dot(np.linalg.inv(np.dot(A.T,A)),np.dot(A.T,X))
    endProj = np.subtract(X,np.dot(A,intProj))
    return np.ndarray.tolist(intProj.T)[0], np.ndarray.tolist(endProj.T)[0]

def Est_BernMLEProjection(A,X): #MLE OF BERNOULLI VARIABLE
    # USING ITERATIVELY REWEIGHTED LEAST SQUARES, SEE WIKIPEDIA FOR NOTATION
    A = np.array(A)
    X = np.array(X)
    np.seterr(all='print')
    currGap = 10
    tol = 1e-2
    n = A.shape[0] # Number of end nodes
    m = A.shape[1] # Number of intermediate nodes
    X = np.reshape(X,(n,1))
    w_k = np.zeros([m,1])
    while currGap > tol:
        mu_k = []
        for i in range(n):
            mu_k.append(float(1/(1+np.exp(-1*((np.dot(w_k.T,A[i])))))))
        Sdiag = []
        for i in range(n):
            Sdiag.append(mu_k[i]*(1-mu_k[i]))            
        mu_k = np.reshape(mu_k,(n,1))
        S_k = np.diag(Sdiag)
        w_k1 = np.dot(np.linalg.inv(np.dot(A.T,np.dot(S_k,A))),np.dot(A.T, np.subtract(np.add(np.dot(np.dot(S_k,A),w_k),X),mu_k)))
        currGap = np.linalg.norm(w_k-w_k1)
        w_k = np.copy(w_k1)
    # Now our importer SF rates are calculated; figure out variance + Wald statistics
    covarMat_Bern = np.linalg.inv(np.dot(A.T,np.dot(S_k,A)))
    w_Var = np.diag(covarMat_Bern)
    wald_stats = []
    for j in range(m):
        wald_stats.append(float((w_k[j]**2)/w_Var[j]))
    
    # Convert to intermediate and end node estimates
    intProj = np.ndarray.tolist(invlogit(w_k.T.tolist()[0]))
    errs_Bern = np.subtract(X,mu_k)
    endProj = errs_Bern.T.tolist()[0]    
    return intProj, endProj, covarMat_Bern, wald_stats

def PlumleeEstimates(ydata, numsamples, A, sens, spec, rglrWt = 0.1):
    ydata = np.array(ydata)
    numsamples = np.array(numsamples)
    beta0 = -6 * np.ones(A.shape[1]+A.shape[0])
    beta0 = beta0 + np.random.uniform(-1,1,np.size(beta0))
    def invlogit_INTERIOR(beta):
        return np.exp(beta)/(np.exp(beta)+1)
    def mynegloglik_INTERIOR(beta, ydata, numsamples, A, sens, spec):
        betaI = beta[0:A.shape[1]]
        betaJ = beta[A.shape[1]:]
        probs = (1-invlogit_INTERIOR(betaJ)) * np.matmul(A,invlogit_INTERIOR(betaI)) + invlogit_INTERIOR(betaJ)
        probsz = probs*sens + (1-probs) * (1-spec)
        return -np.sum(ydata * np.log(probsz) + (numsamples-ydata) * np.log(1-probsz)) \
            + rglrWt*np.sum(np.abs((betaJ - beta0[A.shape[1]:]))) #have to regularize to prevent problems
    
    bds = spo.Bounds(beta0-8, beta0+8)
    opval = spo.minimize(mynegloglik_INTERIOR, beta0+1,
                         args=(ydata, numsamples, A, sens, spec),
                         method='L-BFGS-B',
                         options={'disp': False},
                         bounds=bds)
    return invlogit_INTERIOR(opval.x)[0:A.shape[1]].tolist(), invlogit_INTERIOR(opval.x)[A.shape[1]:].tolist()
    
def testingScheduleGenerator(nodes = [], int_numDays=1000,
                             int_sampleBudget = 1000,
                             int_PolicyType = 0, 
                             arr_PolicyParameter = [0]
                             ):
    """
    Generates a testing schedule list for a desired number of time units, 
    with the following parameters entered:
            0) nodes: The node list generated from the generateNodeListsFromFile function
            1) int_numDays: Number of schedule days desired
            2) int_sampleBudget: Budget amount, in number of samples
            3) int_PolicyType: Desired policy type, one of the following
                0 = Deterministic-End Node; samples 1 end node each day
                ...
            4) arr_PolicyParameter: Array of the parameters required for generatng
                the test schedules
                ...
    
    Outputs a Python list with the following elements within each entry:
            1) Day: Simulation day of the test
            2) Node: Which node to test that day
    """
    #Initialize our output, a list with the above mentioned outputs
    sampleSchedule = []
    
    if int_PolicyType == 0: # Deterministic End-Node policy
        # Identify the list of end nodes
        endNodes = []
        for nodeInd in range(len(nodes)):
            if nodes[nodeInd][1] == 'TRUE': #We have an end node
                endNodes.append(nodeInd)
        # Generate a sampling schedule iterating through each end node
        nodeCount = 0
        currNode = endNodes[nodeCount]
        lastEndNode = endNodes[-1]            
        for i in range(int_numDays): 
            sampleSchedule.append([i,currNode])
            if currNode == lastEndNode:
                nodeCount = 0
                currNode = endNodes[nodeCount]
            else:
                nodeCount += 1
                currNode = endNodes[nodeCount]
        ### END DETERMINISTIC END-NODE POLICY
    
    elif int_PolicyType == 1: # Deterministic End-Node policy with multiple tests per day, set via int_sampleBudget
        # Identify the list of end nodes
        endNodes = []
        for nodeInd in range(len(nodes)):
            if nodes[nodeInd][1] == 'TRUE': #We have an end node
                endNodes.append(nodeInd)
        # Generate a sampling schedule iterating through each end node
        nodeCount = 0
        currNode = endNodes[nodeCount]
        lastEndNode = endNodes[-1]        
        for i in range(int_sampleBudget):
            day = np.mod(i,int_numDays)
            sampleSchedule.append([day,currNode])
            if currNode == lastEndNode:
                nodeCount = 0
                currNode = endNodes[nodeCount]
            else:
                nodeCount += 1
                currNode = endNodes[nodeCount]
        ### END DETERMINISTIC END-NODE POLICY WITH MULTIPLE TESTS
    
    elif int_PolicyType == 2: # Randomized End-Node policy with multiple tests per day, set via int_sampleBudget
        # Identify the list of end nodes
        endNodes = []
        for nodeInd in range(len(nodes)):
            if nodes[nodeInd][1] == 'TRUE': #We have an end node
                endNodes.append(nodeInd)
        numEndNodes = len(endNodes)
        
        # Generate a sampling schedule randomly sampling the list of end nodes
        for i in range(int_sampleBudget):
            day = np.mod(i,int_numDays)
            currEndInd = int(np.floor(np.random.uniform(low=0,high=numEndNodes,size=1)))
            currNode = endNodes[currEndInd]
            sampleSchedule.append([day,currNode])
        ### END RANODMIZED END-NODE POLICY WITH MULTIPLE TESTS
    
    else:
        print('Error generating the sampling schedule.')
    
    # Need to sort this list before passing it through
    sampleSchedule.sort(key=lambda x: x[0])
    return sampleSchedule

 ### END "testingScheduleGenerator" ###

def dynamicTestingGenerator(resultsList,
                             int_totalDays = 1000,
                             int_numDaysRemain = 1000,
                             int_totalBudget = 1000,
                             int_sampleBudgetRemain = 1000,
                             int_PolicyType = 0, 
                             arr_PolicyParameter = [0]
                             ):
    """
    Generates a dynamic testing schedule list for the following day, using the following:
            0) resultsList: A list of testing results so far, with each entry 
                formatted [Node ID, Num Samples, Num Positive, Positive Rate, [IntNodeSourceCounts]]
            1) int_totalDays: Number of simulation days started with
            2) int_numDaysRemain: Number of simulation days remaining
            3) int_totalBudget: Budget amount, in number of samples
            4) int_sampleBudgetRemain: Remaining budget, in number of samples
            5) int_PolicyType: Desired policy type, one of the following:
                0 = Epsilon-Greedy; with probability epsilon, exploit one of the
                nodes with the highest positive rate
                1 = Epsilon-Decreasing; (epsilon cubed) begin mostly exploring nodes, then overtime have a greater chance of exploitation node with highest positive rate (J curve of epsilon values) 
                3 = Epsilon-First; fully explore for eps*int_totalDays. Then exploit for the remaining (1-eps)*int_totalDays
                5 = Every other; explore when creating even sampling schedules, exploit when creating odd sampling schedules. 
                6 = Epsilon-Decreasing; (square root epsilon) begin mostly exploring nodes, then overtime have a greater chance of exploitation node with highest positive rate (rcurve of epsilon values) 
            6) arr_PolicyParameter: Array of the parameters required for generatng
                the test schedules
                ...
    
    Outputs a Python list with the following elements within each entry:
            1) Day: Simulation day of the test
            2) Node: Which node to test that day
    """
    #Initialize our output, a list with the above mentioned outputs
    sampleSchedule = []
    
    if int_PolicyType == 0: # Basic epsilon-greedy algorithm
        nextTestDay = int_totalDays - int_numDaysRemain # The day we are generating a schedule for
        eps = arr_PolicyParameter[0] # Our exploit parameter
        numToTest = int(np.floor(int_sampleBudgetRemain / int_numDaysRemain)) + min(int_sampleBudgetRemain % int_numDaysRemain,1) # How many samples to conduct in the next day
        # Generate a sampling schedule using the current list of results
        # First grab the pool of highest SF rate nodes
        maxSFRate = 0
        maxIndsList = []
        for rw in resultsList:
            if rw[3] > maxSFRate:
                maxSFRate = rw[3]
        for currInd in range(len(resultsList)):
            if resultsList[currInd][3] == maxSFRate:
                maxIndsList.append(currInd)
        for testNum in range(numToTest):
            # Explore or exploit?
            if np.random.uniform() < 1-eps: # Exploit
                exploitBool = True
            else:
                exploitBool = False
            # Based on the previous dice roll, generate a sampling point
            if exploitBool:
                testInd = np.random.choice(maxIndsList)
                NodeToTest = resultsList[testInd][0]
            else:
                testInd = np.random.choice(len(resultsList))
                NodeToTest = resultsList[testInd][0]

            sampleSchedule.append([nextTestDay,NodeToTest])
        ### END EPSILON-GREEDY STRATEGY    
    elif int_PolicyType==1:
        nextTestDay = int_totalDays - int_numDaysRemain # The day we are generating a schedule for
        eps = arr_PolicyParameter[nextTestDay] *arr_PolicyParameter[nextTestDay] *arr_PolicyParameter[nextTestDay]  # Our exploit parameter
        numToTest = int(np.floor(int_sampleBudgetRemain / int_numDaysRemain)) + min(int_sampleBudgetRemain % int_numDaysRemain,1) # How many samples to conduct in the next day
        # Generate a sampling schedule using the current list of results
        # First grab the pool of highest SF rate nodes
        maxSFRate = 0
        maxIndsList = []
        for rw in resultsList:
            if rw[3] > maxSFRate:
                maxSFRate = rw[3]
        for currInd in range(len(resultsList)):
            if resultsList[currInd][3] == maxSFRate:
                maxIndsList.append(currInd)
        for testNum in range(numToTest):
            # Explore or exploit?
            if np.random.uniform() < 1-eps: # Exploit
                exploitBool = True
            else:
                exploitBool = False
            # Based on the previous dice roll, generate a sampling point
            if exploitBool:
                testInd = np.random.choice(maxIndsList)
                NodeToTest = resultsList[testInd][0]
            else:
                testInd = np.random.choice(len(resultsList))
                NodeToTest = resultsList[testInd][0]

            sampleSchedule.append([nextTestDay,NodeToTest])
    
    elif int_PolicyType==2: #Exponential-decay epsilon
        nextTestDay = int_totalDays - int_numDaysRemain # The day we are generating a schedule for
        eps = np.exp(-1*(nextTestDay/int_totalDays)/arr_PolicyParameter[0])
        numToTest = int(np.floor(int_sampleBudgetRemain / int_numDaysRemain)) + min(int_sampleBudgetRemain % int_numDaysRemain,1) # How many samples to conduct in the next day
        # Generate a sampling schedule using the current list of results
        # First grab the pool of highest SF rate nodes
        maxSFRate = 0
        maxIndsList = []
        for rw in resultsList:
            if rw[3] > maxSFRate:
                maxSFRate = rw[3]
        for currInd in range(len(resultsList)):
            if resultsList[currInd][3] == maxSFRate:
                maxIndsList.append(currInd)
        for testNum in range(numToTest):
            # Explore or exploit?
            if np.random.uniform() < 1-eps: # Exploit
                exploitBool = True
            else:
                exploitBool = False
            # Based on the previous dice roll, generate a sampling point
            if exploitBool:
                testInd = np.random.choice(maxIndsList)
                NodeToTest = resultsList[testInd][0]
            else:
                testInd = np.random.choice(len(resultsList))
                NodeToTest = resultsList[testInd][0]

            sampleSchedule.append([nextTestDay,NodeToTest])
        ### END Exponential-decay epsilon policy ###
        
    elif int_PolicyType==3: #EPSILON FIRST POLICY
        nextTestDay = int_totalDays - int_numDaysRemain # The day we are generating a schedule for
        eps = arr_PolicyParameter[0] # Our exploit parameter
        numToTest = int(np.floor(int_sampleBudgetRemain / int_numDaysRemain)) + min(int_sampleBudgetRemain % int_numDaysRemain,1) # How many samples to conduct in the next day
        # Generate a sampling schedule using the current list of results
        # First grab the pool of highest SF rate nodes
        maxSFRate = 0
        maxIndsList = []
        for rw in resultsList:
            if rw[3] > maxSFRate:
                maxSFRate = rw[3]
        for currInd in range(len(resultsList)):
            if resultsList[currInd][3] == maxSFRate:
                maxIndsList.append(currInd)
        for testNum in range(numToTest):
            # Explore or exploit?
            if nextTestDay > (1-eps)*int_totalBudget: # Exploit
                exploitBool = True
            else:
                exploitBool = False
            # Based on the previous dice roll, generate a sampling point
            if exploitBool:
                testInd = np.random.choice(maxIndsList)
                NodeToTest = resultsList[testInd][0]
            else:
                testInd = np.random.choice(len(resultsList))
                NodeToTest = resultsList[testInd][0]

            sampleSchedule.append([nextTestDay,NodeToTest])
        ### EPSILON FIRST POLICY
    
    elif int_PolicyType == 4: # THOMPSON-SAMPLING
        nextTestDay = int_totalDays - int_numDaysRemain # The day we are generating a schedule for
        numToTest = int(np.floor(int_sampleBudgetRemain / int_numDaysRemain)) + min(int_sampleBudgetRemain % int_numDaysRemain,1) # How many samples to conduct in the next day
        # Generate a sampling schedule using the current list of results
        for testNum in range(numToTest):
            # Iterate through each end node, generating an RV according to the beta distribution of samples + positives
            betaSamples = []
            for rw in resultsList:
                alphaCurr = 1 + rw[2]
                betaCurr = 1 + (rw[1]-rw[2])
                sampleCurr = np.random.beta(alphaCurr,betaCurr)
                betaSamples.append(sampleCurr)
            # Select the highest variable
            maxSampleInd = betaSamples.index(max(betaSamples))
            NodeToTest = resultsList[maxSampleInd][0]
            sampleSchedule.append([nextTestDay,NodeToTest])
        
        ### END THOMPSON SAMPLING ###
    elif int_PolicyType == 5: #Every-other sampling. Explore on even days, exploit on odd days
        nextTestDay = int_totalDays - int_numDaysRemain # The day we are generating a schedule for
        eps = arr_PolicyParameter[0] # Our exploit parameter
        numToTest = int(np.floor(int_sampleBudgetRemain / int_numDaysRemain)) + min(int_sampleBudgetRemain % int_numDaysRemain,1) # How many samples to conduct in the next day
        # Generate a sampling schedule using the current list of results
        # First grab the pool of highest SF rate nodes
        maxSFRate = 0
        maxIndsList = []
        for rw in resultsList:
            if rw[3] > maxSFRate:
                maxSFRate = rw[3]
        for currInd in range(len(resultsList)):
            if resultsList[currInd][3] == maxSFRate:
                maxIndsList.append(currInd)
        for testNum in range(numToTest):
            # Explore or exploit?
            if nextTestDay%2  == 1: # Exploit if we are on an odd sampling schedule day
                exploitBool = True
            else:
                exploitBool = False
            # Based on the previous dice roll, generate a sampling point
            if exploitBool:
                testInd = np.random.choice(maxIndsList)
                NodeToTest = resultsList[testInd][0]
            else:
                testInd = np.random.choice(len(resultsList))
                NodeToTest = resultsList[testInd][0]

            sampleSchedule.append([nextTestDay,NodeToTest])
        
    elif int_PolicyType==6: #epsilon decreasing by r-curve (take square root of epsilon values)
        nextTestDay = int_totalDays - int_numDaysRemain # The day we are generating a schedule for
        eps = np.sqrt(arr_PolicyParameter[nextTestDay]) # Our exploit parameter
        numToTest = int(np.floor(int_sampleBudgetRemain / int_numDaysRemain)) + min(int_sampleBudgetRemain % int_numDaysRemain,1) # How many samples to conduct in the next day
        # Generate a sampling schedule using the current list of results
        # First grab the pool of highest SF rate nodes
        maxSFRate = 0
        maxIndsList = []
        for rw in resultsList:
            if rw[3] > maxSFRate:
                maxSFRate = rw[3]
        for currInd in range(len(resultsList)):
            if resultsList[currInd][3] == maxSFRate:
                maxIndsList.append(currInd)
        for testNum in range(numToTest):
            # Explore or exploit?
            if np.random.uniform() < 1-eps: # Exploit
                exploitBool = True
            else:
                exploitBool = False
            # Based on the previous dice roll, generate a sampling point
            if exploitBool:
                testInd = np.random.choice(maxIndsList)
                NodeToTest = resultsList[testInd][0]
            else:
                testInd = np.random.choice(len(resultsList))
                NodeToTest = resultsList[testInd][0]

            sampleSchedule.append([nextTestDay,NodeToTest])
    elif int_PolicyType==7: #epsilon changing by sin
        nextTestDay = int_totalDays - int_numDaysRemain # The day we are generating a schedule for
        eps = (np.sin(12.4*nextTestDay)) # Our exploit parameter
        numToTest = int(np.floor(int_sampleBudgetRemain / int_numDaysRemain)) + min(int_sampleBudgetRemain % int_numDaysRemain,1) # How many samples to conduct in the next day
        # Generate a sampling schedule using the current list of results
        # First grab the pool of highest SF rate nodes
        maxSFRate = 0
        maxIndsList = []
        for rw in resultsList:
            if rw[3] > maxSFRate:
                maxSFRate = rw[3]
        for currInd in range(len(resultsList)):
            if resultsList[currInd][3] == maxSFRate:
                maxIndsList.append(currInd)
        for testNum in range(numToTest):
            # Explore or exploit?
            if 0 < eps: # Exploit
                exploitBool = True
            else:
                exploitBool = False
            # Based on the previous dice roll, generate a sampling point
            if exploitBool:
                testInd = np.random.choice(maxIndsList)
                NodeToTest = resultsList[testInd][0]
            else:
                testInd = np.random.choice(len(resultsList))
                NodeToTest = resultsList[testInd][0]

            sampleSchedule.append([nextTestDay,NodeToTest])
    
    elif int_PolicyType == 8: # TS w NUTS, version 1
        # Grab intermediate and end node distros, then project onto end nodes
        # for different samples from the distro. Pick the largest projected 
        # SF estimate
        # arr_PolicyParameter = [# days to plan for, sensitivity, specificity, M,
        #                        Madapt, delta]
        # How many days to plan for?
        numDaysToSched = min(arr_PolicyParameter[0],int_numDaysRemain)
        usedBudgetSoFar = 0
        firstTestDay = int_totalDays - int_numDaysRemain
        
        if int_numDaysRemain == int_totalDays: # Our initial schedule should just be a distrubed exploration
            currNode = resultsList[0][0]
            for currDay in range(numDaysToSched):
                numToTest = int(np.floor((int_sampleBudgetRemain-usedBudgetSoFar) / (int_numDaysRemain-currDay))) +\
                            min((int_sampleBudgetRemain-usedBudgetSoFar) % (int_numDaysRemain-currDay),1) # How many samples to conduct in the next day
                for testInd in range(numToTest): # Iterate through our end nodes
                    if currNode > resultsList[len(resultsList)-1][0]:
                        currNode = resultsList[0][0]
                        sampleSchedule.append([firstTestDay+currDay,currNode])
                        currNode += 1
                    else:                        
                        sampleSchedule.append([firstTestDay+currDay,currNode])
                        currNode += 1
                    usedBudgetSoFar += 1
            
        else: # Generate NUTS sample using current results and use it to generate a new schedule
            ydata = []
            nSamp = []
            for rw in resultsList:
                ydata.append(rw[2])
                nSamp.append(rw[1])
            A = GenerateTransitionMatrix(resultsList)      
            sens, spec, M, Madapt, delta = arr_PolicyParameter[1:]
            NUTSsamples = GenerateNUTSsamples(ydata,nSamp,A,sens,spec,M,Madapt,delta)
            # Now pick from these samples to generate projections
            for currDay in range(numDaysToSched):
                numToTest = int(np.floor((int_sampleBudgetRemain-usedBudgetSoFar) / (int_numDaysRemain-currDay))) +\
                            min((int_sampleBudgetRemain-usedBudgetSoFar) % (int_numDaysRemain-currDay),1) # How many samples to conduct in the next day
                for testInd in range(numToTest):    
                    currSample = invlogit(NUTSsamples[random.randrange(len(NUTSsamples))])
                    probs = currSample[A.shape[1]:] + np.matmul(A,currSample[:A.shape[1]])
                    # Normalize? Or just pick largest value
                    highInd = [i for i,j in enumerate(probs) if j == max(probs)]
                    currNode = resultsList[highInd[0]][0]
                    sampleSchedule.append([firstTestDay+currDay,currNode])
                    usedBudgetSoFar += 1
  
        
    elif int_PolicyType == 9: # Enhancing exploration w NUTS
        # Grab intermediate and end node distros. Identify intermediate node
        # sample variances. Pick an intermediate node according to a weighting.
        # Pick an outlet from this int node's column in the transition matrix
        # A, again by a weighting (where 0% nodes have a non-zero probability
        # of being selected). log((p/1-p) + eps)?
        # How many days to plan for?
        numDaysToSched = min(arr_PolicyParameter[0],int_numDaysRemain)
        usedBudgetSoFar = 0
        firstTestDay = int_totalDays - int_numDaysRemain
        
        if int_numDaysRemain == int_totalDays: # Our initial schedule should just be a distrubed exploration
            currNode = resultsList[0][0]
            for currDay in range(numDaysToSched):
                numToTest = int(np.floor((int_sampleBudgetRemain-usedBudgetSoFar) / (int_numDaysRemain-currDay))) +\
                            min((int_sampleBudgetRemain-usedBudgetSoFar) % (int_numDaysRemain-currDay),1) # How many samples to conduct in the next day
                for testInd in range(numToTest): # Iterate through our end nodes
                    if currNode > resultsList[len(resultsList)-1][0]:
                        currNode = resultsList[0][0]
                        sampleSchedule.append([firstTestDay+currDay,currNode])
                        currNode += 1
                    else:                        
                        sampleSchedule.append([firstTestDay+currDay,currNode])
                        currNode += 1
                    usedBudgetSoFar += 1
            
        else: # Generate NUTS sample using current results and use it to generate a new schedule
            ydata = []
            nSamp = []
            for rw in resultsList:
                ydata.append(rw[2])
                nSamp.append(rw[1])
            A = GenerateTransitionMatrix(resultsList)      
            sens, spec, M, Madapt, delta = arr_PolicyParameter[1:]
            NUTSsamples = GenerateNUTSsamples(ydata,nSamp,A,sens,spec,M,Madapt,delta)
            # Store sample variances for intermediate nodes
            NUTSintVars = []
            for intNode in range(A.shape[1]):
                currVar = np.var(invlogit(NUTSsamples[:,intNode]))
                NUTSintVars.append(currVar)
            # Normalize sum of all variances to 1
            NUTSintVars = NUTSintVars/np.sum(NUTSintVars)
            
              
            # Now pick from these samples to generate projections
            for currDay in range(numDaysToSched):
                numToTest = int(np.floor((int_sampleBudgetRemain-usedBudgetSoFar) / (int_numDaysRemain-currDay))) +\
                            min((int_sampleBudgetRemain-usedBudgetSoFar) % (int_numDaysRemain-currDay),1) # How many samples to conduct in the next day
                for testInd in range(numToTest):    
                    # Pick an intermediate node to "target", with more emphasis on higher sample variances
                    rUnif = random.uniform(0,1)
                    for intInd in range(A.shape[1]):
                        if rUnif < np.sum(NUTSintVars[0:(intInd+1)]):
                            targIntInd = intInd
                            break
                    # Go through the same process with the column of A
                    # pertaining to this target intermediate node
                    AtargCol = [row[targIntInd] for row in A]
                    # Add a small epsilon, for 0 values, and normalize
                    AtargCol = np.add(AtargCol,1e-3)
                    AtargCol = AtargCol/np.sum(AtargCol)
                    rUnif = random.uniform(0,1)
                    for intEnd in range(A.shape[0]):
                        if rUnif < np.sum(AtargCol[0:(intEnd+1)]):
                            currInd = intEnd
                            break
                    currNode = resultsList[currInd][0]
                    sampleSchedule.append([firstTestDay+currDay,currNode])
                    usedBudgetSoFar += 1                   
    
    else:
        print('Error generating the sampling schedule.')
    
    # Need to sort this list before passing it through
    sampleSchedule.sort(key=lambda x: x[0])
    return sampleSchedule


 ### END "dynamicTestingGenerator" ###
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    '''