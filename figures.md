# Figures

## 2016-11-01: Comparison of different methods

Here we show the addition and removal energies of Fei's IMSRG+QDPT (`qdpt`), Nathan's IMSRG+EOM (`eom`), and Sam's CC (`cc`) as a function of the number of shells in the single-particle basis.  Note that `eom_f` is still IMSRG+EOM but using Fei's IMSRG matrix elements with Nathan's EOM code.

We vary the following parameters:

  - frequency: $\omega = 0.28, 1.0$
  - number of particles: $N = 6, 12, 20$

In each case we compare the $N$-particle state with the $(N \pm 1)$-particle state constructed by adding a particle to the lowest shell above Fermi level / by removing a particle from the highest shell below Fermi level, while also satisfying $0 \le M_L \le 1$.  The $(N \pm 1)$-particle state is uniquely determined by these conditions.

For each figure there is a zoomed-in version that reveals the details on the convergence tail of each curve.

### Observations

**Effect due to difference in IM-SRG matrix elements.**  It is noteworthy that `eom` and `eom_f` differ by a very small amount ($\sim 0.01\%$), far than the discrepancy between the methods themselves.

**Overall trends.**  In the majority of cases, `eom` tends to be sandwiches between `qdpt` and `cc`.

**Outlier methods.**  Sometimes, one of the three methods becomes an "outlier", meaning that it is noticeably different from the other two.  In `\omega = 1.0` cases, the addition energy from `qdpt` is a bit of an outlier.  For all the removal energies, `cc` tends to be an outlier.

**Convergence.**  Convergence of all 3 methods is similar.  Convergence rate decreases with the number of particles as well as the frequency, as expected.  Convergence of removal energy is generally faster than that of addition energy.

**Drifting.**  While the converge initially appears to be rapid, in reality there is a residual drift even after the initial rapid drop.  The drift is small, and tends to be downward for removal energies, and upward for addition energies.  In some situations, the drifts may not even be in the same direction for different methods (for example, addition energy for 12 particles, $ML = 0$, $\omega = 0.28$).

![](FigureFiles/fig-compare-6-0.28-add-1.pdf){width=100%}

![](FigureFiles/fig-compare-12-0.28-add-0.pdf){width=100%}

![](FigureFiles/fig-compare-20-0.28-add-1.pdf){width=100%}

![](FigureFiles/fig-compare-6-1.0-add-1.pdf){width=100%}

![](FigureFiles/fig-compare-12-1.0-add-0.pdf){width=100%}

![](FigureFiles/fig-compare-20-1.0-add-1.pdf){width=100%}

![](FigureFiles/fig-compare-6-0.28-rm-1.pdf){width=100%}

![](FigureFiles/fig-compare-12-0.28-rm-0.pdf){width=100%}

![](FigureFiles/fig-compare-20-0.28-rm-1.pdf){width=100%}

![](FigureFiles/fig-compare-6-1.0-rm-1.pdf){width=100%}

![](FigureFiles/fig-compare-12-1.0-rm-0.pdf){width=100%}

![](FigureFiles/fig-compare-20-1.0-rm-1.pdf){width=100%}
