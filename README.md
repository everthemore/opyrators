A python library that represents many-body fermionic operators as strings, e.g. "012003" for
an operator that has c-dagger on site 2, c on site 3 and a density operator on site 6.
Also has spin-1/2 operators implemented.

## Spins
For spin-1/2 operators, the encoding works as follows:
* 0 = Identity operator
* 1 = Pauli-X operator
* 2 = Pauli-Y operator
* 3 = Pauli-Z operator

An operator with X on site 2, Y on site 4 and Z on site 5, in a system of 6 sites, hence would be "010230".

## Fermions
For particles, the encoding works as follows:
* 0 = Identity operator
* 1 = Creation operator
* 2 = Annihilation operator
* 3 = Density operator


