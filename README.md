# TodaOscillatorSim
Code to simulate the RLD-circuit used in the Toda oscillator experiment.

![bifurkation_sim](https://github.com/user-attachments/assets/137b6a30-d4f5-41eb-bc8b-39afde2f4655)


# Usage
Compile and run:
```
g++ main.cpp -fopenmp -o main && mkdir -p Bifurkation && ./main
```
The intervall and step width for the controll parameter A can be set via A0, A1 and dA respectively. In the bifurkation mode (activated by defining "BIFURKATIONMODE" in line 9 of main.cpp) only the current is saved ones per Period and only for the last 3/4-periods. To save the full simulation results, i.e. time, driving amplitude, current and charge, "BIFUKRATIONMODE" needs to be undefined, either by deleting or commenting line 9 (**This will lead to a large amount of data stored if many values for A are simulated!**)
