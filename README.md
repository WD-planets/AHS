'''
This code has different functions in terms of how many dips you see per transit. 
It's rudimentary. I tried to generalise it with lines 38 to 54 but failed. If you want
to explore that option what fails is the tuple to feed to curve fit

Procedure: I let the out of transit level, the depth factor scale, the
time at minimum flux and the ingress and egress time-scales as free parameters for the first time we fit 
a transit. If I want to track down a transit in different orbits: I fix the depth factor, ingress and 
egress and just let the out-of transit and time at minimum flux as free parameters 
(in order to "preserve") the shape of the transit.

The numbers in the first brackets are the values, and the second bracket correspond 
to uncertainties of the fit.
'''
