import sys
import numpy as np
sys.path.append("../")
from opyrators.operators.operator import *

def make_term_from_indices(L, cdagsites, csites, nsites):
    # Sanity checks
    overlap = [x for x in cdagsites if x in csites]
    if len(overlap) > 0:
          print("Invalid cdagsites and csites have overlap (should be in nsites!)")
          return

    opstring = ["0"] * L

    for c in nsites:
      opstring[c] = "3"; #ns
    for c in cdagsites:
      opstring[c] = "1"; #cdags
    for c in csites:
      opstring[c] = "2"; #cs
    opstring = "".join(opstring)

    return opterm(1,opstring)

def generate_random_term(L):
    # Available sites
    sites = np.arange(L)

    # Randomly pick how many cdags, but at most L/2
    numCdags = np.random.randint(0, L/2)
    cdagsites = np.random.choice(sites, numCdags, replace=False)

    # Compute leftover sites
    leftovers = [s for s in sites if s not in cdagsites]
    # Pick c's (as many as cdags)
    csites = np.random.choice(leftovers, numCdags, replace = False)

    # leftovers for n
    leftovers = [s for s in sites if (s not in cdagsites and s not in csites)]
    numNs = np.random.randint(0,len(leftovers))
    nsites = np.random.choice(leftovers, numNs)

    # Construct the diagonal part
    coeff = np.random.uniform()
    return coeff, make_term_from_indices(L, cdagsites, csites, nsites)

class hoppingOperator:
    def __init__(self, diag, offdiag):
        self.diagonal    = diag
        self.offdiagonal = offdiag

    def addToDiagonal(self, diagonal):
        self.diagonal = self.diagonal + diagonal

    def getFullOperator(self):
        return self.diagonal * self.offdiagonal

    def getOffdiagonalRepresentation(self):
        rep = sorted([self.offdiagonal.opterms[0].string, self.offdiagonal.opterms[1].string])
        return "+".join(rep)

    def __str__(self):
        return self.getFullOperator().__str__()

class randomRangeH:
    def __init__(self, L, numTerms):

        # Store system size
        self.L = L
        self.verbose = False

        # All of the operators
        self.H = operator([])

        # Generate numTerms fully random terms (but conserving particle number)
        for n in range(numTerms):
            # Generate a new term
            coeff, newTerm = generate_random_term(self.L)
            self.H = self.H + coeff*0.5*(operator([newTerm]) + operator([newTerm]).conj())

        #-----------------
        # We want to group similar offdiagonals, so that we can base our choice
        # of which to rotate off of that
        #-----------------
        self.group_terms()

        #-----------------
        # Construct the transformation matrix that takes the operator basis
        # (which has a 1 whenever there is a density op on a site) to the
        # state-basis (computational basis, with a 1 if there is a particle
        # on a site)
        #-----------------
        mSingle = np.array([[1,0],[1,1]]) # Transformation for a single site
        self.m = mSingle
        for i in range(self.L-1):
            self.m = np.kron(self.m,mSingle)
        # Also store the inverse
        self.minv = np.linalg.inv(self.m)

    def group_terms(self):
        offdiagonalindices = {}
        diagonal = operator([])
        offdiagonal_list = []

        # Now we extract a list of all operators from the Hamiltonian
        already_seen_operators = []
        for term in self.H.opterms:

            # Check if we already considered the conjugate of this term
            if( term.conj().string in already_seen_operators ):
                continue
            already_seen_operators.append(term.string)

            if( term.isDiagonal() ):
                diagonal = diagonal + operator([term])
                continue

            # must be off-diagonal
            # Build hoppingOperator
            newop_diag = term.coeff*operator([term.getDiagonal()])
            newop_offdiag = operator([term.getOffDiagonal()]) + operator([term.getOffDiagonal()]).conj()
            ho = hoppingOperator(newop_diag, newop_offdiag)

            # See if we already have this
            strrep = ho.getOffdiagonalRepresentation()
            index = offdiagonalindices.get( strrep, -1 )

            if( index == -1 ):
                offdiagonalindices[ strrep ] = len(offdiagonal_list)
                offdiagonal_list.append( ho )
            else:
                offdiagonal_list[ index ].addToDiagonal( newop_diag )

        # Store diagonals and offdiagonals
        self.diagonals = diagonal
        self.offdiagonals = offdiagonal_list

    def getCoefficientDistributions(self, mean = False):

        # Zero out the dictionaries
        h = {}
        for r in range(self.L+1):
          h[r] = []

        J = {}
        for r in range(self.L+1):
          J[r] = []

        # This list only contains each conjugate once
        for term in self.H.opterms:

          # Extract the "range"
          r = term.range

          # Track diagonals
          if term.isDiagonal():
            h[r].append(np.abs(term.coeff))
            continue

          else:
            # Otherwise, add to the offdiags
            J[r].append(np.abs(term.coeff))

        if not mean:
            return h, J

       # Convert to means
        ha = []
        for r in range(self.L+1):
            if( len(h[r]) == 0 ):
              ha.append(0)
            else:
              ha.append(np.mean(h[r]))

        Ja = []
        for r in range(self.L+1):
            if( len(J[r]) == 0 ):
              Ja.append(0)
            else:
              Ja.append(np.mean(J[r]))

        return ha, Ja

    def invert(self, op, index):

        if( len(op.opterms) == 0):
            return None

        #print("  "*index + "Inverting operator:")
        #for term in op.opterms:
        #    print("  "*index + term.__str__())

        #print("%d"%index + "  "*index + "Considering index %d"%index)
        # If we get to the identity, we're done
        if len(op.opterms) == 1:
            if( op.opterms[0].diagonal_str == "0"*self.L ):
        #        print("  "*index + "Is identity, so we will return 1/coeff: %.3f"%(1/op.opterms[0].coeff))
                return operator([opterm(1/op.opterms[0].coeff,"0"*self.L)])

        # We loop over all possible local density operators
        #for site in range(self.L):
        # Construct density on this site
        opstring = ["0"] * self.L
        opstring[index] = "3";
        opstring = "".join(opstring)
        densop = opterm(1,opstring)
        #print("  "*index + "Checking for density term: ")
        #print("  "*index + densop.__str__())
        densop = operator([densop])

        # We haven't yet expaned this term
        expanded = False

        # now for each term in the operator, we see if the site bit is set.

        # Setting that density to 1 is the same as:
        # If it is, we replace it by a 0, and we just keep it
        # Setting that density to 0 is the same as:
        # If it is, we remove the term, and otherwise we just keep it
        newop1 = operator([])
        newop2 = operator([])
        encountered = False
        for term in op.opterms:

            if term.diagonal_str[index] == '3':
            #    print("  "*index + "This term has that density")
                # So we'll split off two operators
                # One in which we replace this density by an identity
                newopstring = term.diagonal_str[:index] + "0" + term.diagonal_str[index+1:]
                newop1 = newop1 + operator([opterm(term.coeff,newopstring)])
                # And one in which we set it to zero, so it disappears
                encountered = True
            else:
            #    print("  "*index + "This term does not have that density")
                # So we'll keep the term as-is
                newop1 = newop1 + operator([term])
                newop2 = newop2 + operator([term])


        if( not encountered ):
            # Skip this and go to the next index
            print("%d"%index + "  "*index + op.__str__())
            return self.invert(op, index+1)
        else:
            #print("%d"%index + "  "*index + "Current operators")
            #print("%d"%index + "  "*index + newop1.__str__())
            #print("%d"%index + "  "*index + newop2.__str__())
            newterm1 = self.invert(newop1, index+1)
            newterm2 = self.invert(newop2, index+1)

            result = newterm1*densop

            #if( len(newterm2.opterms) != 0):
            if( newterm2 != None ):
                result = result + newterm2*(operator([opterm(1,"0"*self.L)])-densop)

            print("%d"%index + "  "*index + result.__str__())
            return result

    def convertToStateBasis(self, op):
        operator_basis = np.zeros( 2**(self.L) )

        for term in op.opterms:
            binary_op_str = term.diagonal_str.replace('3', '1')
            integer_op = int(binary_op_str, 2)
            operator_basis[integer_op] = term.coeff

        # Convert to state basis
        state_basis = np.dot(self.m,operator_basis)
        return state_basis

    def convertToOperator(self, state_basis):
         # Convert to operator basis
        operator_basis = np.dot(self.minv,state_basis)

        if self.verbose:
            print("Convert to op: ")
            print(operator_basis)

        # Extract operators
        inverse_operator = operator([])
        for i in range(len(operator_basis)):
            if operator_basis[i] != 0:
                i_to_string_w_3s = format(i,'0%db'%self.L)
                i_to_string_w_3s = i_to_string_w_3s.replace('1','3')

                if self.verbose:
                    print("Nonzero: " + i_to_string_w_3s)
                    print("With val: ", operator_basis[i])
                newopterm = opterm(operator_basis[i], i_to_string_w_3s)
                inverse_operator = inverse_operator + operator([newopterm])

        return inverse_operator

    def computeInverse(self, op):

        state_basis = self.convertToStateBasis(op)
        # Invert
        state_basis = 1/state_basis

        return self.convertToOperator(state_basis)

    def rotateOut(self):
        # See if there is anything to rotate out; return True if we're diagonal
        if( len(self.offdiagonals) == 0 ):
            return True

        # Compute a list of all the coeffs of the offdiagonals
        amplitudes = []
        for o in self.offdiagonals:
            tmp = 0
            for i in o.diagonal.opterms:
                tmp += np.abs(i.coeff)**2
            amplitudes.append(tmp)

        # Pick the max
        index = np.argmax(amplitudes)
        op = self.offdiagonals[index]

        print("Picking %d"%index)

        # Get the diagonal part
        diagt = op.diagonal
        # Get the offdiagonal part, but the - version
        A = operator([op.offdiagonal.opterms[0]]) - operator([op.offdiagonal.opterms[1]])

        if self.verbose:
            print("Rotating out with:")
            print(A)
            print("Which has diagonal prefactor: ", diagt)

        # Construct unitary
        deltaEterms = (A*self.diagonals - self.diagonals*A).cleanup()

        if self.verbose:
            print("Commutator with diagonals: ")
            print(deltaEterms)

        #deltaEterms = [deltaEterms.opterms[i] for i in range(0,len(deltaEterms.opterms),2)]

        # Extract only the diagonal parts of all the operators
        deltaV = operator([])
        for d in deltaEterms.opterms:
            if( d.isDiagonal() ):
                deltaV = deltaV + operator([opterm(d.coeff, d.diagonal_str)])
            else:
                deltaV = deltaV + operator([opterm(d.coeff/2, d.diagonal_str)])

        if self.verbose:
            print("Delta V is: ")
            print(deltaV)

        delta_state_basis = self.convertToStateBasis(deltaV)
        t_state_basis = self.convertToStateBasis(diagt)

        if( delta_state_basis.all() == 0 ):
            r = np.ones_like(delta_state_basis)*np.pi/4
        else:
            r = np.arctan(2*t_state_basis/delta_state_basis)/2

        if self.verbose:
            print("vd: ", delta_state_basis)
            print("vt: ", t_state_basis)
            print("vq: ", r)

        sin_rotator = self.convertToOperator( np.sin(r) )
        cos_rotator = self.convertToOperator( 1 - np.cos(r) )

        if self.verbose:
            print("Sin rot")
            print(np.sin(r))
            print(self.convertToOperator(np.sin(r)))
            print(sin_rotator)
            print("Cos rot")
            print(cos_rotator)
            print(A)
            print(A*A)

        Sm = operator([opterm(1,"0"*A.length)]) - sin_rotator*A + cos_rotator*A*A
        Sp = operator([opterm(1,"0"*A.length)]) + sin_rotator*A + cos_rotator*A*A

        if self.verbose:
            print("transormation")
            print(Sm)

        # Rotate out
        newH = Sm * self.H * Sp
        # Update the Hamiltonian
        self.H = newH.cleanup()

        # Regroup
        self.group_terms()

        return self.H.isDiagonal()

np.random.seed(0)
L = 6
H = randomRangeH(L,4)

H.H = H.H + operator([opterm(0.1,"003000")])
H.H = H.H + operator([opterm(0.2,"000030")])
H.H = H.H + operator([opterm(0.3,"330000")])
H.H = H.H + operator([opterm(0.4,"000303")])


print(H.H)

h = []
J = []
hstep, Jstep = H.getCoefficientDistributions(mean=True)
h.append(hstep)
J.append(Jstep)

finished = False
#while not finished:
numSteps = 50
for step in range(numSteps):
    finished = H.rotateOut()
    hstep, Jstep = H.getCoefficientDistributions(mean=True)
    h.append(hstep)
    J.append(Jstep)

finished = H.rotateOut()

print(H.H)
