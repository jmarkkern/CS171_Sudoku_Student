import SudokuBoard
import Variable
import Domain as domain
import Trail
import Constraint
import ConstraintNetwork
import time
import random

class BTSolver:

    # ==================================================================
    # Constructors
    # ==================================================================

    def __init__ ( self, gb, trail, val_sh, var_sh, cc ):
        self.network = ConstraintNetwork.ConstraintNetwork(gb)
        self.hassolution = False
        self.gameboard = gb
        self.trail = trail

        self.varHeuristics = var_sh
        self.valHeuristics = val_sh
        self.cChecks = cc

    # ==================================================================
    # Consistency Checks
    # ==================================================================

    # Basic consistency check, no propagation done
    def assignmentsCheck ( self ):
        for c in self.network.getConstraints():
            if not c.isConsistent():
                return False
        return True

    """
        Part 1 TODO: Implement the Forward Checking Heuristic

        This function will do both Constraint Propagation and check
        the consistency of the network

        (1) If a variable is assigned then eliminate that value from
            the square's neighbors.

        Note: remember to trail.push variables before you assign them
        Return: a tuple of a dictionary and a bool. The dictionary contains all MODIFIED variables, mapped to their MODIFIED domain.
                The bool is true if assignment is consistent, false otherwise.
    """
    def forwardChecking ( self ):
        modified_vars = {}

        if len(self.trail.trailStack) == 0:
            #preprocessing
            for v in self.network.getVariables():
                if v.isAssigned():
                    for n in self.network.getNeighborsOfVariable(v):
                        n.removeValueFromDomain(v.getAssignment())
                        self.trail.push(n)
                    self.forwardChecking()
                    
        else:            
            trail_last = self.trail.trailStack[-1][0]
            neighbors = self.network.getNeighborsOfVariable(trail_last)
            for n in neighbors:
                if not n.isAssigned() and n.getDomain().isEmpty():
                    return (modified_vars,False)
                elif n.isChangeable and not n.isAssigned():
                    if n.getDomain().contains(trail_last.getAssignment()):
                        if n.getDomain().size() == 1:
                            return (modified_vars,False)
                        else:
                            self.trail.push(n)
                            n.removeValueFromDomain(trail_last.getAssignment())
                            modified_vars[n] = n.getDomain()
                    
                    if n.getDomain().size() == 1:
                        n_allowed = True
                        for nn in self.network.getNeighborsOfVariable(n):
                            if(nn.getDomain().contains(n.domain.values[0])):
                                n_allowed = False
                        if n_allowed:
                            self.trail.push(n)
                            n.assignValue(n.domain.values[0])
                    
        return (modified_vars,True)


    # =================================================================
	# Arc Consistency
	# =================================================================
    def arcConsistency( self ):
        assignedVars = []
        for c in self.network.constraints:
            for v in c.vars:
                if v.isAssigned():
                    assignedVars.append(v)
        while len(assignedVars) != 0:
            av = assignedVars.pop(0)
            for neighbor in self.network.getNeighborsOfVariable(av):
                if neighbor.isChangeable and not neighbor.isAssigned() and neighbor.getDomain().contains(av.getAssignment()):
                    neighbor.removeValueFromDomain(av.getAssignment())
                    if neighbor.domain.size() == 1:
                        neighbor.assignValue(neighbor.domain.values[0])
                        assignedVars.append(neighbor)

    
    """
        Part 2 TODO: Implement both of Norvig's Heuristics

        This function will do both Constraint Propagation and check
        the consistency of the network

        (1) If a variable is assigned then eliminate that value from
            the square's neighbors.

        (2) If a constraint has only one possible place for a value
            then put the value there.

        Note: remember to trail.push variables before you assign them
        Return: a pair of a dictionary and a bool. The dictionary contains all variables 
		        that were ASSIGNED during the whole NorvigCheck propagation, and mapped to the values that they were assigned.
                The bool is true if assignment is consistent, false otherwise.
    """
    def norvigCheck ( self ):
        modified_vars,firstHeuristic = self.forwardChecking()
        max_domain = self.gameboard.p * self.gameboard.q
        assigned_vars = {}
        
        for var,dom in modified_vars.items():
            neighbor_assignments = [n.getAssignment() for n in self.network.getNeighborsOfVariable(var) if n.isAssigned() or n.getDomain().size() == 1]
            if(len(set(neighbor_assignments))) >= max_domain:
                return ({},False)
            if var.getDomain().size() == 1 and dom.size() == 1 and not var.getDomain().values[0] in neighbor_assignments:
                self.trail.push(var)
                var.assignValue(var.getDomain().values[0])
                assigned_vars[var] = dom
                if not self.norvigCheck():
                    return ({},False)


        return (assigned_vars,firstHeuristic)

    """
         Optional TODO: Implement your own advanced Constraint Propagation

         Completing the three tourn heuristic will automatically enter
         your program into a tournament.
     """
    def getTournCC ( self ):
        return False

    # ==================================================================
    # Variable Selectors
    # ==================================================================

    # Basic variable selector, returns first unassigned variable
    def getfirstUnassignedVariable ( self ):
        for v in self.network.variables:
            if not v.isAssigned():
                return v

        # Everything is assigned
        return None

    """
        Part 1 TODO: Implement the Minimum Remaining Value Heuristic

        Return: The unassigned variable with the smallest domain
    """
    def getMRV ( self ):
        min_domain = self.gameboard.p * self.gameboard.q # largest domain = p * q
        mrv_var = None
        for v in self.network.getVariables():
            if not v.isAssigned():
                if v.getDomain().size() <= min_domain: # if v has smaller domain, becomes new minimum :)
                    min_domain = v.getDomain().size()
                    mrv_var = v
        return mrv_var

    """
        Part 2 TODO: Implement the Minimum Remaining Value Heuristic
                       with Degree Heuristic as a Tie Breaker

        Return: The unassigned variable with the smallest domain and affecting the  most unassigned neighbors.
                If there are multiple variables that have the same smallest domain with the same number of unassigned neighbors, add them to the list of Variables.
                If there is only one variable, return the list of size 1 containing that variable.
    """
    def MRVwithTieBreaker ( self ):
        min_domain = self.gameboard.p * self.gameboard.q # largest domain = p * q
        mrv_vars = []
        for v in self.network.getVariables():
            if not v.isAssigned():
                if v.getDomain().size() < min_domain: # if v has smaller domain, becomes new minimum :)
                    min_domain = v.getDomain().size()
                    mrv_vars = []
                    mrv_vars.append(v)
                elif v.getDomain().size() == min_domain:
                    mrv_vars.append(v)
        
        if len(mrv_vars) == 1:
            return mrv_vars
        elif len(mrv_vars) == 0:
            return [None]
        else:
            return_vars = []
            degree_max = 0
            for v in mrv_vars:
                degree = len(self.network.getConstraintsContainingVariable(v))
                if degree > degree_max:
                    return_vars = []
                    return_vars.append(v)
                elif degree == degree_max:
                    return_vars.append(v)
            return return_vars
        
            

    """
         Optional TODO: Implement your own advanced Variable Heuristic

         Completing the three tourn heuristic will automatically enter
         your program into a tournament.
     """
    def getTournVar ( self ):
        return None

    # ==================================================================
    # Value Selectors
    # ==================================================================

    # Default Value Ordering
    def getValuesInOrder ( self, v ):
        values = v.domain.values
        return sorted( values )

    """
        Part 1 TODO: Implement the Least Constraining Value Heuristic

        The Least constraining value is the one that will knock the least
        values out of it's neighbors domain.

        Return: A list of v's domain sorted by the LCV heuristic
                The LCV is first and the MCV is last
    """
    def getValuesLCVOrder ( self, v ):
        d_count = [] # empty list
        for i in range(0, len(v.domain.values)): # iterate thru domain of v by index
            count = 0
            for n in self.network.getNeighborsOfVariable(v):
                if v.domain.values[i] in n.domain.values: # count number of times value shows up in domain of neighbors
                    count += 1
            d_count.append((count, v.domain.values[i])) # append tuple of (# times on neighbor domain, value)
        output = [i[1] for i in sorted(d_count)] # output list of values sorted by # times on neighbor domain
        return output

    """
         Optional TODO: Implement your own advanced Value Heuristic

         Completing the three tourn heuristic will automatically enter
         your program into a tournament.
     """
    def getTournVal ( self, v ):
        return None

    # ==================================================================
    # Engine Functions
    # ==================================================================

    def solve ( self, time_left=600):
        if time_left <= 60:
            return -1

        start_time = time.time()
        if self.hassolution:
            return 0

        # Variable Selection
        v = self.selectNextVariable()

        # check if the assigment is complete
        if ( v == None ):
            # Success
            self.hassolution = True
            return 0

        # Attempt to assign a value
        for i in self.getNextValues( v ):

            # Store place in trail and push variable's state on trail
            self.trail.placeTrailMarker()
            self.trail.push( v )

            # Assign the value
            v.assignValue( i )

            # Propagate constraints, check consistency, recur
            if self.checkConsistency():
                elapsed_time = time.time() - start_time 
                new_start_time = time_left - elapsed_time
                if self.solve(time_left=new_start_time) == -1:
                    return -1
                
            # If this assignment succeeded, return
            if self.hassolution:
                return 0

            # Otherwise backtrack
            self.trail.undo()
        
        return 0

    def checkConsistency ( self ):
        if self.cChecks == "forwardChecking":
            return self.forwardChecking()[1]

        if self.cChecks == "norvigCheck":
            return self.norvigCheck()[1]

        if self.cChecks == "tournCC":
            return self.getTournCC()

        else:
            return self.assignmentsCheck()

    def selectNextVariable ( self ):
        if self.varHeuristics == "MinimumRemainingValue":
            return self.getMRV()

        if self.varHeuristics == "MRVwithTieBreaker":
            return self.MRVwithTieBreaker()[0]

        if self.varHeuristics == "tournVar":
            return self.getTournVar()

        else:
            return self.getfirstUnassignedVariable()

    def getNextValues ( self, v ):
        if self.valHeuristics == "LeastConstrainingValue":
            return self.getValuesLCVOrder( v )

        if self.valHeuristics == "tournVal":
            return self.getTournVal( v )

        else:
            return self.getValuesInOrder( v )

    def getSolution ( self ):
        return self.network.toSudokuBoard(self.gameboard.p, self.gameboard.q)
