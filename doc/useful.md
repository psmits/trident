Tutz and Schmid book 
--------------------


One also speaks of grouped survival data or interval censoring. In grouped survival data there are typically some observations that have the same survival time. This phenomenon is usually referred to as “ties”. In continuous time, ties ideally should not occur. In fact, some models and estimation methods for continuous time even assume that there are no ties in the data. Nevertheless, in practical applications with continuous event times ties occur, which might be taken as a hint for underlying grouping. In some areas, for example in demography, discrete data are quite natural.


- intrinsically discrete measurements, where the measurements represent natural numbers
- grouped data, which represent events in underlying time intervals, and the response refers to an interval

While one has a constant hazard, h(t) = h, on the individual level, the hazard
on the population level is a decreasing function of t determined by the parameters  ̨
and ˇ. What one observes are in fact realizations of T on the population level, not on
the individual level. Therefore the (estimated) hazard function refers to the marginal
hazard rate, not to the individual hazard rate which is not observed. It should be
noted that even if the model for the marginal hazard fits the data well, this does not
mean that the model on the individual level holds. This is because similar population
models can be derived from quite different individual level models.


