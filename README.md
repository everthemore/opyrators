# Opyrators
Opyrators (pronounced the same way as the word 'operators') is a lightweight python package that represents many-body fermionic and spin operators as strings.

For example, the fermionic string "012003" stands for
a 6-site operator formed by a creation operator on site 2, an annihilation operator on site 3, a density operator on site 6 and identity operators on the rest.

In this representation, operator manipulations such as addition and multiplication are easily implemented. That means commutation relations are easy too! [See the example below](#examples).

## Example
Here is a quick example showing some of the basic features of fermionic opyrators.

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

## Encodings
### Spins
For spin-1/2 operators, the encoding works as follows:
* 0 = Identity operator
* 1 = Pauli-X operator
* 2 = Pauli-Y operator
* 3 = Pauli-Z operator

An operator with X on site 2, Y on site 4 and Z on site 5, in a system of 8 sites, hence would be "01023000".

### Fermions
For particles, the encoding works as follows:
* 0 = Identity operator
* 1 = Creation operator
* 2 = Annihilation operator
* 3 = Density operator
