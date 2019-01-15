# BioInformatics-Project-2019
University of Zagreb, Faculty of Electrical Engineering and Computing

Bioinformatics project for course [Bioinformatika](https://www.fer.unizg.hr/predmet/bio) 


## Compile & run:
`g++ -std=c++11 -Wall -I . ./*.cpp -o main.o `

`./main.o`


## Compile Parallel execution
`g++ -std=c++11 -fopenmp -Wall -I . ./*.cpp -o main.o`

`./main.o`

## Default values while executing
If no parameters are provided e.g. `./main.o`, the following default ones are used: 
```
k=15
w=5
c_tres = 12
min_gap 15
ecoli = 0 (0 for lambda or 1 for ecoli)
```

Example tu run from terminal for the same default input:
`./main.o 15 5 12 15 0`


## Score checking
To check the score run:
`make lambda`
or
`make ecoli`
