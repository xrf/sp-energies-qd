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

## 2016-01-07: Extrapolations

Here, we are using a power law model:

$$E(K) = \alpha K^b + \gamma$$

where $K$ is the total number of shells.  Note: In the draft paper I used $-\beta \equiv b$ as the exponent.

If we plot the derivative of energy on a loglog plot, then it should show up as a straight line.  The draft paper goes into some more detail.

There are two horizontal dotted lines for each method in the plots.  These two represent the infinite-basis limit extrapolated using the fit.  The reason there are two is because I did the fits using 4 slightly different approaches:

  - `logderiv` is just a linear fit to the loglog energy derivative plot; this provides us with a good initial estimate for the exponent and coefficient;
  - `deriv` is a nonlinear power law fit to the energy derivative plot (not shown explicitly), using `logderiv` as a starting guess;
  - `fixedab` is a nonlinear power law fit to the energy plot, using `logderiv` as the parameters, but the coefficient and exponent are both *fixed* and only the constant term (i.e. the infinite basis limit) is allowed to vary;
  - `full` is like the above, but with all three parameters allowed to vary; it uses `fixedab` as the starting point.

Of these four methods, the last two give us the infinite basis limit, and the results are slightly different due to the difference in fitting strategy.  This actually gives us a *diagnostic*: discrepancy between the two strategies is a hint to the goodness of fit.

### Ground state energy

It seems that the extrapolated IMSRG energy is generally closer to DMC than MP2 for ground state energy.  However, for larger (`num_filled >= 6`) cases it is not sure if this holds.  More data is required for these cases we are not yet close to convergence.

HF has a hard time following the trend: the derivative data is very noisy, and moreover the derivative just seems to hit zero very quickly.  Its convergence *may* be exponential.  I removed them from the plots as a result.

### Addition/removal energy

The data here is similar, but there is a lot more stuff to look at.  Some interesting things:

  - Extrapolated QDPT results appear to disagree significantly in some cases: (add, 2, 0.28), (add, 6, 0.28), (add, 12, 0.28), and (add, 2, 1.0).  In these cases, QDPT energies sometimes *increase* rather than decrease wrt number of shells.  The turning of the curve can screw up the fit, since the loglog fit only sees the absolute value of the derivative, so the extrapolated result is probably not very trustworthy.  So what exactly causes these results to deviate so much?  Why does it seem to affect the *few particle* cases more?
  - In other cases, extrapolated QDPT is quite close to EOM, but further from CC.
  - In the few places where IMSRG-EOM with quadruples (`eom_quad`) data is available, it is closer to CC.
  - The agreement between DMC and all other methods is noticeably worse than for ground state energy.

![](FigureFiles/fig-fit-add-1-0.28.pdf)
![](FigureFiles/fig-fit-add-1-1.0.pdf)
![](FigureFiles/fig-fit-add-2-0.28.pdf)
![](FigureFiles/fig-fit-add-2-1.0.pdf)
![](FigureFiles/fig-fit-add-3-0.28.pdf)
![](FigureFiles/fig-fit-add-3-1.0.pdf)
![](FigureFiles/fig-fit-add-4-1.0.pdf)
![](FigureFiles/fig-fit-rm-1-0.28.pdf)
![](FigureFiles/fig-fit-rm-1-1.0.pdf)
![](FigureFiles/fig-fit-rm-2-0.28.pdf)
![](FigureFiles/fig-fit-rm-2-1.0.pdf)
![](FigureFiles/fig-fit-rm-3-0.28.pdf)
![](FigureFiles/fig-fit-rm-3-1.0.pdf)
![](FigureFiles/fig-fit-rm-4-1.0.pdf)
![](FigureFiles/fig-fit-ground-1-0.28.pdf)
![](FigureFiles/fig-fit-ground-1-1.0.pdf)
![](FigureFiles/fig-fit-ground-2-0.28.pdf)
![](FigureFiles/fig-fit-ground-2-1.0.pdf)
![](FigureFiles/fig-fit-ground-3-0.28.pdf)
![](FigureFiles/fig-fit-ground-3-1.0.pdf)
![](FigureFiles/fig-fit-ground-4-1.0.pdf)
![](FigureFiles/fig-fit-ground-5-0.28.pdf)
![](FigureFiles/fig-fit-ground-5-1.0.pdf)
![](FigureFiles/fig-fit-ground-6-1.0.pdf)

## 2016-02-14: Extrapolation uncertainty

To recap: we are performing a fit on the data using the power law:

$$E = a K^b + c$$

where $K$ is the number of shells and $E$ is ground, addition, or removal energy.  The constant term $c$ is the infinite-basis limit in this model, provided that $b < 0$.  For this section, we always use 5-point-consecutive fits.

The focus of this section is to answer these two questions:

  - How far are the extrapolated results from the actual results?
  - Can we predict the uncertainty of extrapolated results without knowledge of the infinite-basis result?

Unfortunately, due to the inability to perform calculations with an infinite basis, we must instead compromise by using an estimate.  Naturally, we would guess that the "best" result comes from extrapolated results near the maximum number of shells for which we have data (treat this as the *formal* definition of the word "best" in this section).  For systems with few particles, we have good confidence that the asymptotics have been reached, based on the quality of the fits.  But for larger cases, this is less certain.

To mitigate this deficiency, we would prefer to exclude the best fits from consideration if there is any suspicion that asymptotics have not been reached.  To investigate precisely what criterion to use, we performed a preliminary examination of each best fit by visual inspection and marked them as either "good" or "bad".  After some analysis, we found that this resulted in roughly three categories of fits:

  - A. Fits with good visual agreement on the log-log derivative plots.
  - B. Fits with apparently poor visual agreement in the log-log derivative
    plots, but relative sum of squared residuals (SSR) is unusually low.
    ("Relative" here means divided by $c$.)
  - C. Fits with poor visual agreement in both log-log derivative plots as
    well as the ordinary energy plots, and relative SSR is high.

Category B was somewhat unique in that it occurred predominantly for Hartree-Fock cases.  The convergence of Hartree-Fock tends to be very rapid, so the tails of the plots are often extremely flat, leading to visually poor power law fits but with good statistics (low SSR and uncertainty of $c$).

Even in category B fits are often marked "bad", we chose to accept them nonetheless as their uncertainties are often quite good.  Moreover, there are strong indications that HF does not have the power law tail that other methods do, in which judging based on the log-log derivative plot would be wrong.

After tallying up the results in a histogram of the relative SSR, we found that the boundary between category (A, B) and category C lies at roughly $10^{-6}$, which is coincidentally (?) the tolerance used in our IM-SRG calculations.  This is the final cutoff that we use to decide whether a fit was "good" or "bad".  Best fits with bad relative SSR are excluded from consideration in the remaining parts of this section.

Now that we have a selection of systems each of which has a best value for $c_{\mathrm{best}}$ that represents our infinite-basis limit, we proceed to perform fits on subsets of the data below the maximum number of shells (e.g. fit between 11-15, then 12-16, then 13-17, etc).  This allows us to simulate what happens had we stopped collecting data earlier than what we have now.

After doing this we gather all of the results and calculate the *discrepancy* $\varepsilon$:

$$\varepsilon = c - c_{\mathrm{best}}$$

This quantifies the "actual" error.  We divide the discrepancy by the fit uncertainty $\sigma_c$ as estimated by the Levenberg-Marquardt optimization algorithm:

$$\frac{\varepsilon}{\sigma_c}$$

This is plotted as the $y$-axis in all figures.  If the fit uncertainty is accurate, then we expect this quantity to form a standard normal distribution.  In reality, we find that the magnitude of this value can fluctuate anywhere between 1 and a few hundred.  However, the cases above 50 or so ("outliers") are actually quite rare, which suggests there's something that we aren't accounting for.  In the plots, outliers tend to occur when the fitting algorithm find a spurious minimum that seem far better than it really is.

After some trial and error, we found that the quantity
$$Q = \log_{10}(\sigma_c / \sqrt{\mathrm{RSS}})$$
is actually a useful statistic for excluding outliers.  The $Q$ number compares the uncertainty of the constant term with the residuals themselves: if the root-mean-squared residuals are larger than the fit uncertainty, that can be rather suspicious.

The $Q$ number is the $c$-component of $(\boldsymbol J^{\mathrm{T}} \boldsymbol J)^{-1/2}$ where $\boldsymbol J$ is the Jacobian of the model function.  This could also be interpreted as the $c$-component of an approximation of the inverse-square-root of the Hessian.

The plot against $Q$ is shown in first and third figures.  Note that we have split HF from non-HF methods.  HF has distinctively different behavior from the rest of the methods.  It was suggested that this is likely because Hartree-Fock fails to probe the Coulomb cusp, thus rendering the power law model highly inaccurate.  Nonetheless, despite the incorrect model, the discrepancy-to-uncertainty ratio is quite good.  This is probably due to the rapid HF convergence.

For non-HF, the results are not as great.  The outlier problem makes the plot difficult to read, but fortunately they occur only for $Q < 0$ or so.

Another parameter we could plot against is the relative uncertainty $\sigma_c / c$.  These are the second and fourth figures.  The outliers have been excluded by the $Q < 0$ criterion for the non-HF case.  For HF, we observe that HF uncertainties are generally very slightly overestimated.

For non-HF, there is a curious drift upward.  It seems that above around -3.5 or so the uncertainty quantifies the discrepancy well, but below -3.5 the uncertainty tends to underestimate the true discrepancy by about a factor of ten.  Moreover, the discrepancy tends to be skewed more toward the positive side, meaning that the extrapolated results tend to be above the true value.

In overall, the fit uncertainty is not a bad estimate of the error, albeit some corrections need to be applied for cases where $\log_{10}(\sigma_c / c) < -3.5$ for non-HF methods.

![](FigureFiles/fig-fit-predictiveness-contour-5-hessian-False.pdf)
![](FigureFiles/fig-fit-predictiveness-contour-5-err-True.pdf)
![](FigureFiles/fig-fit-predictiveness-contour-5-hessian-True.pdf)
![](FigureFiles/fig-fit-predictiveness-contour-5-err-False.pdf)
