# Opyrators
Opyrators (pronounced the same way as the word 'operators') is a lightweight python package that represents many-body fermionic and spin operators as strings.

For example, the fermionic string "012003" stands for
a 6-site operator formed by a creation operator on site 2, an annihilation operator on site 3, a density operator on site 6 and identity operators on the rest. [See an example use below](#examples).

In this representation, operator manipulations such as addition and multiplication are easily implemented. That means commutation relations are easy too! For example:

<img src="https://latex.codecogs.com/svg.latex?[c^\dagger_2 n_3 c_4, c_1 n_2 c^\dagger_4] = \underbrace{[\texttt{0132},\texttt{2301}] = \texttt{2133}}_{\textrm{use opyrators!}} = c_1c^\dagger_2n_3n_4"/>

Opyrators has a few extra features, such as automatically keeping track of the extent/range of operator terms (the max distance between non-trivial operators), which is very useful for implementing flow equations based on operators.

## Spins
For spin-1/2 operators, the encoding works as follows:
* 0 = Identity operator
* 1 = Pauli-X operator
* 2 = Pauli-Y operator
* 3 = Pauli-Z operator

An operator with X on site 2, Y on site 4 and Z on site 5, in a system of 8 sites, hence would be "01023000".

## Fermions
For particles, the encoding works as follows:
* 0 = Identity operator
* 1 = Creation operator
* 2 = Annihilation operator
* 3 = Density operator

## Examples
```python
# Import the fermion operators.
from opyrators.fermions import *

# Create single operators, i.e. terms in what is possibly
# an operator with many terms.
A = opterm(1.3, "112233")  # Scalar prefactor 1.3
B = opterm(0.34, "112233")  # Scalar prefactor 0.34
# Construct two separate operators, both with only one
# term each
opA = operator([A])
opB = operator([B])
opC = opA * opB - opB * opA
print(opC)

# The output of the print operation shows that this complex
# operator consists of 6 terms, but is fully diagonal.
# Term 0: 0.442 303003
# Term 1: -0.442 303033
# Term 2: -0.442 303303
# Term 3: -0.442 000333
# Term 4: 0.442 003333
# Term 5: 0.442 300333
```
