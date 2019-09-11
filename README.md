# SpikeSynchrony.jl

Measuring distances, synchrony and correlation between spike trains. Most methods implemented here are summarized in the relevant [Scholarpedia entry](http://www.scholarpedia.org/article/Measures_of_spike_train_synchrony), but the documentation strings also cite the papers.

SpikeSynchrony.jl has a clear, small and easy-to-understand source code, which could be helpful in education purposes. Notice that the present code has not been optimized in any way for performance (as the functions are already fast enough that there is no reason to).

To install, press `]` in the Julia console to enter package mode and then
do
```
add https://github.com/Datseris/SpikeSynchrony.jl
```

To use it in your Julia sessions, do `using SpikeSynchrony`.

## Functionality
These are the exported functions from `SpikeSynchrony`. See their documentation strings for how to use them:
```
vanRossum # the van Rossum distance
SPIKE_distance # self-explanatory
SPIKE_distance_profile # the function S(t)
```

## Related software
The following software, implemented in other languages, have much more functionality than we do:

1. [SPIKY](http://wwwold.fi.isc.cnr.it/users/thomas.kreuz/Source-Code/SPIKY.html) : MATLAB with a Graphical User Interface
2. [PySpike](http://mariomulansky.github.io/PySpike/) : Python
3. [cSpike](http://wwwold.fi.isc.cnr.it/users/thomas.kreuz/Source-Code/cSPIKE.html) : MATLAB with MEX C++ backends

*(If you know more, add them here via a PR)*

## Acknowledgements
I would like to thank [Thomas Kreuz](http://wwwold.fi.isc.cnr.it/users/thomas.kreuz/) for support regarding the implementations here, as well as useful discussions.
