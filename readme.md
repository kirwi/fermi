### Mixed-dimensional fermions, oh my!

This is a collection of code I wrote in graduate school to solve a 
fermionic lattice problem. Two species of fermions interact on a lattice -
one species can hop around on a 1D chain, the other can hop around on a 2D
square lattice. Whenever the two species hop onto the same lattice site
they interact, repelling or attracting each other. These interactions drive
correlations in the 1D system, and manifest themselves in the form of an
"effective potential".

These codes calculate the effective potential experienced by the 1D system
mediated by interactions with the 2D system. Evolution of the correlation
functions betweem fermions on the 1D lattice are calculated by solving 
differential "flow equations" using the renormalization group (RG) field 
theoretic technique. As the correlation functions grow exponentially under
the renormalization process, the code checks for divergence and identifies
the phase of the system. 

The phase information is stored in a list and
plotted on the fly, allowing the user to input a list of initial conditions
with which to begin renormalizing (a typical part of the RG workflow).
Matplotlib is used to plot phase diagrams and visualizations of the 
effective interaction.
