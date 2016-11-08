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

## 2016-11-08: Altered Coulomb interaction

(Author: Fei)

To investigate the drifts, Scott suggested a modified interaction (short-range, no cusp):

$$V_{\mathrm I}(r; \sigma, \mu) \equiv \left(1 - \mathrm e^{-r^2 / (2 \sigma)}\right)\frac{\mathrm e^{-\mu r}}{r}$$

where $\sigma$ and $\mu$ are chosen to be near $1$ in natural units.

Simen's OpenFCI code already supports a fairly general Hamiltonian consisting of a reciprocal power, polynomial factor, and a Gaussian factor.  It was therefore easier to use something more like this:

$$V_{\mathrm{II}}(r; \sigma_{\mathrm A}, \sigma_{\mathrm B}) \equiv C \left(1 - \mathrm e^{-r^2 / (2 \sigma_{\mathrm A}^2)}\right) \frac{\mathrm e^{-r^2 / (2 \sigma_{\mathrm B}^2)}}{r}$$

I choose $\sigma_{\mathrm A} = 0.5$ and $\sigma_{\mathrm B} = 4.0$ (in units of $a_0 / \sqrt{\omega}$ where $a_0$ is the Bohr radius) so that there would be some range remaining where the Coulomb interaction is still present.  I also re-scaled the envelope so that its peak is $1$.  The precise scaling factor is $C \equiv (1 + c)^{1 - 1/c} / c$ where $c \equiv \sqrt{\sigma_{\mathrm B} / \sigma_{\mathrm A}}$.

Judging from the fractional change of the last two data points (`rel_slope`)

$$\frac{E_{\mathtt{num\_shells}{=}15} - E_{\mathtt{num\_shells}{=}14}}{E_{\mathtt{num\_shells}{=}15}}$$

it does appear to converge faster: the fractional change of the modified interaction is about 2 orders of magnitude smaller than that of the standard Coulomb interaction (usually, at least.  It varies from case to case).  Keep in mind that the precision of the ODE solver in my calculations is set to $10^{-6}$ anyway, so you can't really expect any better.

![](FigureFiles/fig-compare-12-0.28-add-0_sigmaA=0.5_sigmaB=4.0.pdf){width=100%}
![](FigureFiles/fig-compare-12-0.28-rm-0_sigmaA=0.5_sigmaB=4.0.pdf){width=100%}
![](FigureFiles/fig-compare-12-1.0-add-0_sigmaA=0.5_sigmaB=4.0.pdf){width=100%}
![](FigureFiles/fig-compare-12-1.0-rm-0_sigmaA=0.5_sigmaB=4.0.pdf){width=100%}
![](FigureFiles/fig-compare-20-0.28-add-1_sigmaA=0.5_sigmaB=4.0.pdf){width=100%}
![](FigureFiles/fig-compare-20-0.28-rm-1_sigmaA=0.5_sigmaB=4.0.pdf){width=100%}
![](FigureFiles/fig-compare-20-1.0-add-1_sigmaA=0.5_sigmaB=4.0.pdf){width=100%}
![](FigureFiles/fig-compare-20-1.0-rm-1_sigmaA=0.5_sigmaB=4.0.pdf){width=100%}
![](FigureFiles/fig-compare-6-0.28-add-1_sigmaA=0.5_sigmaB=4.0.pdf){width=100%}
![](FigureFiles/fig-compare-6-0.28-rm-1_sigmaA=0.5_sigmaB=4.0.pdf){width=100%}
![](FigureFiles/fig-compare-6-1.0-add-1_sigmaA=0.5_sigmaB=4.0.pdf){width=100%}
![](FigureFiles/fig-compare-6-1.0-rm-1_sigmaA=0.5_sigmaB=4.0.pdf){width=100%}
