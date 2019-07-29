import sys
import numpy as np
sys.path.append("../../")
from operators.operator import *

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
            h[r].append(term.coeff)
            continue

          else:
            # Otherwise, add to the offdiags
            J[r].append(term.coeff)

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

        deltaEterms = [deltaEterms.opterms[i] for i in range(0,len(deltaEterms.opterms),2)]

        if self.verbose:
            print("IF THESE ARE NOT SORTED, CHECK")
            print("Resulting deltaEterms: ")

        # Invert the delta
        invDelta = operator([]) # Identity
        for i in deltaEterms:
            invDelta = invDelta + 1/i.coeff * operator([opterm(1,i.diagonal_str)])

        # Make identity if no operators
        if invDelta.length == 0:
            invDelta = 0

        if self.verbose:
            print("Inverse delta operator: ", invDelta)
            print("Rotator: ", rotator)

        rotator = (2 * diagt * invDelta).cleanup()

        sin_rotator = operator([])
        cos_rotator = operator([opterm(1,"0"*self.L)])

        for i in range(len(rotator.opterms)):

            if( rotator.opterms[0].coeff == 0 ):
                l = np.pi/4
            else:
                l = np.arctan(rotator.opterms[i].coeff)/2

            sin_rotator = sin_rotator + operator([opterm(np.sin(l), rotator.opterms[i].string)])
            cos_rotator = cos_rotator - operator([opterm(np.cos(l), rotator.opterms[i].string)])

        Sm = operator([opterm(1,"0"*A.length)]) - sin_rotator*A + cos_rotator*A*A
        Sp = operator([opterm(1,"0"*A.length)]) + sin_rotator*A + cos_rotator*A*A

        # Rotate out
        newH = Sm * self.H * Sp
        # Update the Hamiltonian
        self.H = newH

        # Regroup
        self.group_terms()
        return False


np.random.seed(0)
H = randomRangeH(4,10)

print(H.H)

#for o in H.offdiagonals:
#    print(o.getFullOperator())

#print("diagonals")
#for diags in H.diagonals.opterms:
#    print(diags)

h = []
J = []
hstep, Jstep = H.getCoefficientDistributions(mean=True)
h.append(hstep)
J.append(Jstep)

finished = False
while not finished:
    finished = H.rotateOut()
    hstep, Jstep = H.getCoefficientDistributions(mean=True)
    h.append(hstep)
    J.append(Jstep)

print(H.H)
