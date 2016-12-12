---
title: Figures
toc: true
---

## 2016-11-01: Comparison of different methods

Here we show the addition and removal energies of Fei's IMSRG+QDPT (`qdpt`), Nathan's IMSRG+EOM (`eom`), and Sam's CC (`cc`) as a function of the number of shells in the single-particle basis.  Note that `eom_f` is still IMSRG+EOM but using Fei's IMSRG matrix elements with Nathan's EOM code.

We vary the following parameters:

  - frequency: $\omega = 0.28, 1.0$
  - number of particles: $N = 6, 12, 20$

In each case we compare the $N$-particle state with the $(N \pm 1)$-particle state constructed by adding a particle to the lowest shell above Fermi level / by removing a particle from the highest shell below Fermi level, while also satisfying $0 \le M_L \le 1$.  The $(N \pm 1)$-particle state is uniquely determined by these conditions.

For each figure there is a zoomed-in version that reveals the details on the convergence tail of each curve.

### Observations

**Effect due to difference in IM-SRG matrix elements.**  It is noteworthy that `eom` and `eom_f` differ by a very small amount ($\sim 0.01\%$), far less than the discrepancy between the post-IM-SRG methods themselves.

**Overall trends.**  In the majority of cases, `eom` tends to be sandwiched between `qdpt` and `cc`.

**Outlier methods.**  Sometimes, one of the three methods becomes an "outlier", meaning that it is noticeably different from the other two.  In $\omega = 1.0$ cases, the addition energy from `qdpt` is a bit of an outlier.  For all the removal energies, `cc` tends to be an outlier.

**Convergence.**  Convergence of all 3 methods is similar.  Convergence rate decreases with the number of particles as well as the frequency, as expected.  Convergence of removal energy is generally faster than that of addition energy.

**Drifting.**  While the converge initially appears to be rapid, in reality there is a residual drift even after the initial rapid drop.  The drift is small.  It is downward for removal energies, but for addition energies, it can be either down or up, and can vary by the method too.

![](FigureFiles/fig-compare-6-0.28-add-0.pdf)
![](FigureFiles/fig-compare-6-1.0-add-0.pdf)
![](FigureFiles/fig-compare-6-0.28-rm-1.pdf)
![](FigureFiles/fig-compare-6-1.0-rm-1.pdf)
![](FigureFiles/fig-compare-12-0.28-add-1.pdf)
![](FigureFiles/fig-compare-12-1.0-add-1.pdf)
![](FigureFiles/fig-compare-12-0.28-rm-0.pdf)
![](FigureFiles/fig-compare-12-1.0-rm-0.pdf)
![](FigureFiles/fig-compare-20-0.28-add-0.pdf)
![](FigureFiles/fig-compare-20-1.0-add-0.pdf)
![](FigureFiles/fig-compare-20-0.28-rm-1.pdf)
![](FigureFiles/fig-compare-20-1.0-rm-1.pdf)

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

![](FigureFiles/fig-compare-6-0.28-add-0_sigmaA=0.5_sigmaB=4.0.pdf)
![](FigureFiles/fig-compare-6-0.28-rm-1_sigmaA=0.5_sigmaB=4.0.pdf)
![](FigureFiles/fig-compare-6-1.0-add-0_sigmaA=0.5_sigmaB=4.0.pdf)
![](FigureFiles/fig-compare-6-1.0-rm-1_sigmaA=0.5_sigmaB=4.0.pdf)
![](FigureFiles/fig-compare-12-0.28-add-1_sigmaA=0.5_sigmaB=4.0.pdf)
![](FigureFiles/fig-compare-12-0.28-rm-0_sigmaA=0.5_sigmaB=4.0.pdf)
![](FigureFiles/fig-compare-12-1.0-add-1_sigmaA=0.5_sigmaB=4.0.pdf)
![](FigureFiles/fig-compare-12-1.0-rm-0_sigmaA=0.5_sigmaB=4.0.pdf)
![](FigureFiles/fig-compare-20-0.28-add-0_sigmaA=0.5_sigmaB=4.0.pdf)
![](FigureFiles/fig-compare-20-0.28-rm-1_sigmaA=0.5_sigmaB=4.0.pdf)
![](FigureFiles/fig-compare-20-1.0-add-0_sigmaA=0.5_sigmaB=4.0.pdf)
![](FigureFiles/fig-compare-20-1.0-rm-1_sigmaA=0.5_sigmaB=4.0.pdf)

## 2016-11-11: Changes with respect to frequency

(Author: Fei)

Shown here are the plots against frequency from 0.1 to 1.0.  The changes are quite smooth and very roughly linear.  There is not a lot of interesting trends I see though.  Perhaps we should look at more frequencies beyond 1.0 or 0.1?

Erratum: All of the previous plots displayed the addition energy with an incorrect $M_{\mathrm L}$.  This has been corrected.  The trends are a little different now: there are now more cases where the addition energy drift in different directions for different methods.

I have also added a semilog-plot to show the absolute value of `rel_slope` (as defined previously) for different systems.  The x-axis is labeled by the system parameters: `(interaction, freq, num_particles)`.  The normal interaction is labeled `''` and the modified interaction is labeled `'_sigmaA=0.5_sigmaB=4.0'`.  From this plot we see that:

  - The convergence rate between removal and addition energy is not as dramatically different as we thought.  Addition energy has more unusual behavior (bending and then drifting upward), but the absolute rate seems to be not so different from that of removal energy.
  - For the normal interaction, the ordering of `rel_slope` for removal and addition energies depends on the frequency.
  - The softened interaction clearly improves convergences: the plot makes it very obvious now.  Removal energy is almost always has a higher `rel_slope`.  (There is still that strange spike for 20 particles with 0.28 frequency, but that is also the most difficult case we have looked at so far with only 15 shells available, so it's likely that we have not reached the asymptotics yet.)

![](FigureFiles/fig-by-freq-10.0-6-add-0.pdf)
![](FigureFiles/fig-by-freq-10.0-6-rm-1.pdf)
![](FigureFiles/fig-by-freq-10.0-6-add-0_sigmaA=0.5_sigmaB=4.0.pdf)
![](FigureFiles/fig-by-freq-10.0-6-rm-1_sigmaA=0.5_sigmaB=4.0.pdf)
![](FigureFiles/fig-rel-slopes.pdf)

## 2016-11-13: Model attempt #1

(Author: Fei)

A rough pre-Hartree-Fock approximation of the addition energy is given by

$$\varepsilon \approx f_{u, u} = \hat H^{\text o}_{u, u} + \sum_i V_{u, i, u, i}$$

where $u$ denotes the *harmonic oscillator* orbital that is being added.  If we assume the orbitals are harmonic oscillator ones, then we can dedimensionalize the Hamiltonian via the substitution $r' = \sqrt{\omega} r$ and $E' = E / \omega$.  Taking into account the form of the interaction $V \sim r^{-1}$, this results in

$$\varepsilon \approx \omega \hat H^{\text o\prime}_{u, u} + \sqrt{\omega} \sum_i V'_{u, i, u, i}$$

where the primed quantities are dedimensionalized and thus independent of $\omega$.  The linear coefficient $\hat H^{\text o\prime}_{u, u}$ is just the harmonic oscillator energy of $u$.  So if we have 6 particles and want to add another, this is just 3.

To see how well this matches up with the results, we re-write the equation and plot $\varepsilon / \sqrt\omega - 3 \sqrt\omega$ against $\sqrt{\omega}$.

The curve appears to approach a constant value as $\sqrt{\omega} \to \infty$, which is expected.  This value is $\sum_i V'_{u, i, u, i}$.  We can confirm this value by calculating it directly from the matrix elements:

    0.86165 * 2 - 0.23500 + (0.72457 * 2 - 0.13702) * 2 = 4.1125

This agrees with what is shown on the graph.

![](FigureFiles/fig-by-freq-model1-10.0-6-add-0.pdf)
![](FigureFiles/fig-by-freq-model1-zoomed-10.0-6-add-0.pdf)

## 2016-12-08: More altered interactions

Now we consider the cases where $(\sigma_{\mathrm{A}}, \sigma_{\mathrm{B}}) = (0.5, \infty)$ and $(\sigma_{\mathrm{A}}, \sigma_{\mathrm{B}}) = (0, 4.0)$.  The results are shown the figures below.

The results are also included in `fig-rel-slopes` (last figure of 2016-11-11).  Looking at this figure, one can see that making $\sigma_{\mathrm{A}}$ nonzero had a huge effect on improving the convergence, whereas making $\sigma_{\mathrm{B}}$ finite had little effect.

![](FigureFiles/fig-compare-6-0.28-add-0_sigmaA=0.5.pdf)
![](FigureFiles/fig-compare-6-0.28-rm-1_sigmaA=0.5.pdf)
![](FigureFiles/fig-compare-6-1.0-add-0_sigmaA=0.5.pdf)
![](FigureFiles/fig-compare-6-1.0-rm-1_sigmaA=0.5.pdf)
![](FigureFiles/fig-compare-12-0.28-add-1_sigmaA=0.5.pdf)
![](FigureFiles/fig-compare-12-0.28-rm-0_sigmaA=0.5.pdf)
![](FigureFiles/fig-compare-12-1.0-add-1_sigmaA=0.5.pdf)
![](FigureFiles/fig-compare-12-1.0-rm-0_sigmaA=0.5.pdf)
![](FigureFiles/fig-compare-20-0.28-add-0_sigmaA=0.5.pdf)
![](FigureFiles/fig-compare-20-0.28-rm-1_sigmaA=0.5.pdf)
![](FigureFiles/fig-compare-20-1.0-add-0_sigmaA=0.5.pdf)
![](FigureFiles/fig-compare-20-1.0-rm-1_sigmaA=0.5.pdf)

![](FigureFiles/fig-compare-6-0.28-add-0_sigmaB=4.0.pdf)
![](FigureFiles/fig-compare-6-0.28-rm-1_sigmaB=4.0.pdf)
![](FigureFiles/fig-compare-6-1.0-add-0_sigmaB=4.0.pdf)
![](FigureFiles/fig-compare-6-1.0-rm-1_sigmaB=4.0.pdf)
![](FigureFiles/fig-compare-12-0.28-add-1_sigmaB=4.0.pdf)
![](FigureFiles/fig-compare-12-0.28-rm-0_sigmaB=4.0.pdf)
![](FigureFiles/fig-compare-12-1.0-add-1_sigmaB=4.0.pdf)
![](FigureFiles/fig-compare-12-1.0-rm-0_sigmaB=4.0.pdf)
![](FigureFiles/fig-compare-20-0.28-add-0_sigmaB=4.0.pdf)
![](FigureFiles/fig-compare-20-0.28-rm-1_sigmaB=4.0.pdf)
![](FigureFiles/fig-compare-20-1.0-add-0_sigmaB=4.0.pdf)
![](FigureFiles/fig-compare-20-1.0-rm-1_sigmaB=4.0.pdf)
