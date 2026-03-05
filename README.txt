RSolve suite, verion 1.0

1. General Overview

The code in this project is intended for creating possibly hierarchical,
filamentous networks and determining their sheet resistance. Utilities are
also provided for subtractive and additive manufacturing approaches to allow
such networks to be fabricated and tested. The project currently encompasses
the following executables:

-hnet_maker: A terminal based program that prompts the user for information
about a possibly hierarchical, filamentous network

-rsolve: A program for computing the resistance of a filamentous resistor 
network

-stlwrite: A program for creating STL files to allow for 3D simulations and
3D printing of networks

-svgwrite: A program for creating SVG files from filamentous networks, to
allow them to be cut with a laser cutter or cutting plotter.

In addition, files ending in ".lat" are included. Such files specify the
structure of a filamentous lattice. Upon learning the system by which these
files specify the structure of a lattice, a user can extend the functionality
of this software. At present, the following lattices are implemented by default:

1. square.lat - A square lattice

2. tri.lat - The triangular lattice

3. kag.lat - The Kagome lattice, which has the underlying rotational and
translational symmetries of the triangular lattice and a three site basis.
Here, nearest neighbor connections are assumed.

4. hcmb.lat - The honeycomb lattice, which has the translational and rotational
symmetries of the triangular lattice with a two site basis. Once more, nearest
neighbor connections are used, and the lattice is oriented such that the
"armchair" edge runs along the horizontal direction.

2. Installation

I have tried to write as much of this code as possible in standard C++. The
code does, however, have some external dependencies. Two of these, Triangle
and the Geometric Predicates libraries, both by Jonathan Richard Shewchuck,
are included with the code. Please note that, while professor Shewchuck has
released his geometric predicates library into the public domain, Triangle is
still under copyright. This code must not be used for any commercial purpose
without obtaining permission from Professor Shewchuck. In addition, this
code relies on SuiteSparse, written by Professor Timothy Alden Davis and
collaborators. While not necessary, rsolve can be sped up if you build
SuiteSparse with the optional GPU support.

In addition to the above software, you should also have a C++ compiler that
supports at least the C++ 17 standard. Once you have all dependencies installed,
make sure paths specified in the included Makefile are appropriate for your
system. Once you have done this, you should be able to build the entire
project by simply typing "make".

3. Acknowledgements

I began work on this project as a graduate student at the Georgia Institute of
Technology, where I completed my doctoral studies under the supervision of
Professor Peter Yunker. Professor Yunker proposed the project and provided
ample support. I have also incorporated some improvements based on techniques I
developed while working as a postdoctoral research fellow in the group of
Professor Moumita Das at RIT. This opportunity gave me the chance to hone
my skills and learn valuable techniques which I believe have made the software
more performant and robust. I was also alerted to numerous opportunities for
improvement during valuable discussions with Henry Murphy.

While this code contains a number of innovations of my own, I am also deeply
indebted to a community of experts in computational geometry and numerical
linear algebra. In addition to the libraries mentioned above, I have also
adapted code from Physically Based Rendering, Third Edition, by Matt Pharr,
Wenzel Jakob, and Greg Humphreys. These authors have generously made the
third edition of their remarkable text, and the accompanying code, available
free of charge on the Web. I have also benefitted considerably from the
excellent text Real-Time Collision Testing, by Christer Ericson, which is a gold
mine of information about intersection and proximity tests, among many other
topics.
