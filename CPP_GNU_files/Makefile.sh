#!/bin/bash

# run solar system simulation
g++ -std=c++14  ./main.cpp ./add.cpp ./planet.cpp

./a.out

# remove output file
rm ./a.out




# generates a .png that shows the path of all planets (inclucding the Sun) within a certain timespan
gnuplot ./solarsystem.gnu


# generates a .png that shows the energy conservation
gnuplot ./energy_cons.gnu


# creates a GIF animation of the above mentioned path (needs 'animateSun.gnuplot')
gnuplot ./sunGIF.gnu


#creates a GIF animation of the planets' movement in a cartesian grid (needs 'animatePlanets.gnuplot')
gnuplot ./planetsGIF.gnu


# for temporary exclusion
: <<EOF

# insert here


EOF
