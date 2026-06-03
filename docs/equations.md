## The physics behind

The temporal evolution of the phase space density $f$ can be described by the Fokker–Planck equation [(Subbotin and Shprits, 2012)](https://doi.org/10.1029/2011JA017467):

$$
\frac{\partial f}{\partial t} = \sum_{i, j} \frac{\partial}{\partial J_i} D_{J_i, J_j} \frac{\partial f}{\partial J_j}
$$

where $J_1 = \mu$, $J_2 = J$, and $J_3 = \Phi$ are the adiabatic invariants. $D_{J_i, J_j}$ are the diffusion coefficients due to wave–particle resonant interactions, which are typically calculated in terms of $(p, \alpha, L)$ . [Lyons et al. (1972)](https://doi.org/10.1029/JA077i019p03455) give the general form for bounce-averaged quasi-linear diffusion coefficients $\langle D_{\alpha \alpha} \rangle$ (pitch-angle diffusion), $\langle D_{\alpha p} \rangle$ (mixed diffusion), and $\langle D_{p p} \rangle$ (momentum diffusion) as

$$
\langle D_{\alpha \alpha} \rangle = \frac{1}{\tau_B}\int^{\tau_B}_{0} D_{\alpha \alpha}(\alpha) \left(\frac{\partial \alpha_{eq}}{\partial \alpha}\right)^2 \mathrm{d}t
$$

$$
\langle D_{\alpha p} \rangle = \frac{1}{\tau_B}\int^{\tau_B}_{0} D_{\alpha p}(\alpha) \left(\frac{\partial \alpha_{eq}}{\partial \alpha}\right) \mathrm{d}t
$$

$$
\langle D_{p p} \rangle = \frac{1}{\tau_B}\int^{\tau_B}_{0} D_{p p}(\alpha) \mathrm{d}t
$$

Here, $D_{\alpha \alpha}$, $D_{\alpha p}$, and $D_{p p}$ are local diffusion coefficients, $\alpha$ and $\alpha_{eq}$ are local and equatorial pitch angles, respectively, and $\tau_B$ is the electron bounce period.

The quasi-linear bounce-averaged diffusion coefficients are calculated using a dipole magnetic field in the form given by [Glauert and Horne (2005)](https://doi.org/10.1029/2004JA010851):

$$
\langle D_{\alpha \alpha} \rangle = \frac{1}{S(\alpha_{eq})}\int^{\lambda_m}_{0} D_{\alpha \alpha}(\alpha) \left(\frac{\cos \alpha}{\cos^2 \alpha_{eq}}\right)\cos^7 \lambda \,\mathrm{d}\lambda 
$$

$$
\langle D_{\alpha p} \rangle = \frac{1}{S(\alpha_{eq})}\int^{\lambda_m}_{0} D_{\alpha p}(\alpha) \left(\frac{\cos^4 \lambda (1 + 3 \sin^2 \lambda)^{\frac{1}{4}}}{\cos \alpha}\right) \mathrm{d}\lambda 
$$

$$
\langle D_{p p} \rangle = \frac{1}{S(\alpha_{eq})}\int^{\lambda_m}_{0} D_{p p}(\alpha) \left(\frac{\cos \lambda (1 + 3 \sin^3 \lambda)^{\frac{1}{2}}}{\cos \alpha}\right) \mathrm{d}\lambda
$$

where $\lambda_m$ is the mirror latitude and $S(\alpha_{eq})$ gives the variation of $\tau_B$ with $\alpha_{eq}$, approximated by $S(\alpha_{eq}) = 1.3 - 0.56 \sin \alpha_{eq}$ in a dipole magnetic field [(Hamlin et al., 1961)](https://doi.org/10.1029/JZ066i001p00001). $\alpha_{eq}$ is associated with the local pitch angle by $\sin^2 \alpha = \frac{\sqrt{4 - 3\cos^2\lambda}}{\cos^6\lambda}\sin^2\alpha_{eq}$ . $\lambda_m$ is the upper limit of the magnetic latitude determined either by the mirror latitude of the particles, which is found by solving the polynomial

$$
C_l^6 + 3 C_l \sin^4 \alpha_{eq} - 4 \sin^4 \alpha_{eq} = 0
$$

for $C_l = \cos^2 \lambda_m$, or by the maximum latitude of wave occurrence.

The local pitch-angle diffusion rate $D_{\alpha \alpha}$ is calculated following [Albert (2007)](https://doi.org/10.1029/2007JA012551):

$$
D_{\alpha \alpha} = \frac{\Omega_\sigma}{\gamma^2}\frac{B^2_{wave}}{B_0^2} \sum^\infty_{n=-\infty} \sum_\omega D^n_{\alpha \alpha} 
$$

$$
D^n_{\alpha \alpha} = \int^{\theta_{max}}_{\theta_{min}} \sin \theta \,\mathrm{d}\theta \,\Delta_N G_1 G_2 
$$

$$
\Delta_N(\omega, \theta) = \frac{\pi}{2} \frac{\sec \theta}{\left|\frac{v_{||}}{c}\right|^3} \Phi^2_n \frac{\left(-\sin^2 \alpha + \frac{n\Omega_\sigma}{\omega\gamma}\right)^2}{\left|1 - \frac{\left(\frac{\partial\omega}{\partial k_{||}}\right)_\theta}{v_{||}}\right|} 
$$

$$
G_1(\omega) = \frac{|\Omega_\sigma| B^2(\omega)}{\int^{\omega_{UC}}_{\omega_{LC}}B^2(\omega')\,\mathrm{d}\omega'} 
$$

$$
G_2(\omega, \theta) = \frac{g_w(\theta)}{N(\omega)} 
$$

$$
N(\omega) = \int^{\theta_{max}}_{\theta_{min}}\mathrm{d}\theta' \,\Gamma g_\omega (\theta') 
$$

$$
\Gamma = \mu^2\left| \mu + \omega \frac{\partial \mu}{\partial \omega} \right| 
$$

$$
\Phi^2_n = \left[ \left( \frac{D}{\mu^2 - S}\right)^2 \left( \frac{\mu^2 \sin^2 \theta - P}{\mu^2} \right)^2 + \left( \frac{P \cos^2 \theta}{\mu^2} \right) \right]^{-1} \left[ \frac{\mu^2 \sin^2\theta - P}{2\mu^2} \left( 1 + \frac{D}{\mu^2 - S} \right) J_{n+1} + \frac{\mu^2 \sin^2\theta - P}{2\mu^2} \left( 1 + \frac{D}{\mu^2 - S} \right) J_{n-1} + \cot\alpha\, \sin\theta \cos\theta\, J_n \right]
$$

Here, $\mu = \frac{kc}{\omega}$ is the refractive index given as a function of $(\omega, \theta)$ by cold plasma theory [(Stix, 1962)](https://ui.adsabs.harvard.edu/abs/1962tpw..book.....S). The term $\Phi^2_n$ accounts for the relationships between the components of the wave electric and magnetic fields as well as the gyroaveraging of the resonant wave–particle phase. It contains Bessel functions $J_n$ and $J_{n\pm1}$ with arguments $\frac{k_\perp p_\perp}{m_\sigma \Omega_\sigma}$ and the Stix coefficients $D$, $S$, and $P$ .

The function $B^2(\omega)$ describes the frequency distribution of wave power and is taken to be zero unless $\omega$ lies between the lower and upper cutoffs $\omega_{LC}$ and $\omega_{UC}$ . It can be specified in the input file and in the subroutine `Bw_wave.f90`.

Similarly, the distribution of wave power with normal angle $\theta$ is described by $g_\omega(\theta)$, which is zero unless $\theta$ lies between $\theta_{min}$ and $\theta_{max}$ . It can be specified in the input file and in the subroutine `Theta_wave.f90`.

$G_2(\omega, \theta)$ and $\Delta_N(\omega, \theta)$ are computed at the resonant frequency $\omega$ corresponding to $\theta$ and $n$, determined from the cyclotron resonance condition. There may be several such values of $\omega$ as simultaneous solutions to the cold plasma wave dispersion relation.

Local cross diffusion rate $D_{\alpha p}$ and momentum diffusion rate $D_{p p}$ can be subsequently obtained by

$$
D_{\alpha p} = D_{\alpha \alpha} \left[ \frac{\sin \alpha \cos \alpha}{\frac{n\Omega_\sigma}{\gamma \omega_k} - \sin^2\alpha} \right]
$$

$$
D_{p p} = D_{\alpha \alpha} \left[ \frac{\sin \alpha \cos \alpha}{\frac{n\Omega_\sigma}{\gamma \omega_k} - \sin^2\alpha} \right]^2
$$
