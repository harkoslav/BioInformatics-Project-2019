# BioInformatics-Project-2019
Bioinformatics project for course Bioinformatika FER
https://www.fer.unizg.hr/predmet/bio

## Compile & run:
`g++ -std=c++11 -Wall -I . ./*.cpp -o main.o `

`./main.o`


## Compile Parallel execution
`g++ -std=c++11 -fopenmp -Wall -I . ./*.cpp -o main.o`

`./main.o`

## Default values while executing
if no parameters are provided:
k=15
w=5
c_tres = 12
min_gap 15
ecoli = 0 (0 for lambda, 1 for ecoli)

example tu run from cmd for same default input:
`./main.o 15 5 12 15 0`


## Score checking
To check the score run:
`make lambda`
or
`make ecoli`
