# drs4_read_events

A C++ program that reads the binary files produce by the DRS4 Oscilloscope program (based on an original program written by Stefan Ritt at PSI). 

** This code is very much in development. 

## Dependencies

So you'll Boost, and to link against Boost's Program Options library. 

## Compilation 

So to compile this program, please use the following: 

g++ read_events.cpp -o read_events -lm -lstdc++ -lboost_program_options
