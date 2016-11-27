# OllipticID

Olliptic is a multigrid elliptic sover for numerical relativity 

## Synopsis

Olliptic  computes initial data for black hole evolution.
See:

PHYSICAL REVIEW D 82, 024005 (2010)
Numerical evolution of multiple black holes with accurate initial data
Pablo Galaviz Bernd Bruegmann and Zhoujian Cao

for more details. 

## Installation

To compile olliptic:

`make`

run 

`make help` to see more options 

## Examples

### binary black hole

./Olliptic.x -eq Punctures -B Robin -BA 0.0 -Bn 1 -NumBH 2 -bhmp1 0.487209 -bhy1 5.5 -bhpx1 -0.0901099 -bhpy1 -0.000703975 -bhmp2 0.487209 -bhy2 -5.5 -bhpx2 0.0901099 -bhpy2 0.000703975 -L 1000 -N 65 -O 2 -l 9 -s 1e-11 -Dp 0.05 -nu2 20 -o TestBH -pc yes




