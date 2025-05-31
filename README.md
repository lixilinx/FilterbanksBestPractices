# Best practices with filterbanks

Filterbank is versatile but also tricky. Here, I’d like to share what I have learned from my many years of practice (mainly for audio processing). The math and design methods are from [my paper](https://ieeexplore.ieee.org/document/8304771). It is a straightforward time domain design framework that can cover most cases (DFT and DCT filterbanks, DWT and complex DWT, constraints like latency, symmetry and phase linearity, control on sidelobe tails, ...). 

Topics:

[What is a filterbank](https://github.com/lixilinx/FilterbanksBestPractices/tree/main#what-is-a-filterbank)

[MIRROR symmetry if latency is unconcerned](https://github.com/lixilinx/FilterbanksBestPractices/tree/main#mirror-symmetry-if-latency-is-unconcerned)

[SAME symmetry if low latency is crucial](https://github.com/lixilinx/FilterbanksBestPractices/tree/main#same-symmetry-if-low-latency-is-crucial)

[Symmetries and least squares reconstruction](https://github.com/lixilinx/FilterbanksBestPractices/tree/main#symmetries-and-least-squares-reconstruction)

[Designs with lower aliasing do not necessarily performs better for tasks like AEC](https://github.com/lixilinx/FilterbanksBestPractices/tree/main#designs-with-lower-aliasing-do-not-necessarily-performs-better-for-tasks-like-aec)

[Squeezing synthesis filter does not help low latency filterbank design](https://github.com/lixilinx/FilterbanksBestPractices/tree/main#squeezing-synthesis-filter-does-not-help-low-latency-filterbank-design)

[Replacing STFT with filterbank](https://github.com/lixilinx/FilterbanksBestPractices/tree/main#replacing-stft-with-filterbank)

[Spare that SRC](https://github.com/lixilinx/FilterbanksBestPractices/tree/main#spare-that-src)

[Frequency shift with analytic signal](https://github.com/lixilinx/FilterbanksBestPractices/tree/main#frequency-shift-on-analytic-signal-only)

[Remove the systemic phase biases](https://github.com/lixilinx/FilterbanksBestPractices/tree/main#remove-the-systemic-phase-biases)

[DCT filterbank for tasks like AEC?](https://github.com/lixilinx/FilterbanksBestPractices/tree/main#dct-filterbank-for-tasks-like-aec)

[Special designs: wavelet, complex wavelet, quadrature mirror filter (QMF), customized, nonsubsampled and nonuniform filterbanks](https://github.com/lixilinx/FilterbanksBestPractices/tree/main#special-designs-wavelet-complex-wavelet-quadrature-mirror-filter-qmf-customized-nonsubsampled-and-nonuniform-filterbanks)

[Optimization of discrete design freedoms](https://github.com/lixilinx/FilterbanksBestPractices/tree/main#optimization-of-discrete-design-freedoms)

[Exact phase linearity, a luxury or a must-have?](https://github.com/lixilinx/FilterbanksBestPractices/tree/main#exact-phase-linearity-a-luxury-or-a-must-have)

[Fine balance between resolution and leakage](https://github.com/lixilinx/FilterbanksBestPractices/tree/main#fine-balance-between-resolution-and-leakage)

Note: there are some issues rendering bold math symbols, say a bold 𝐚 maybe rendered as a or aa in the equation.   

### What is a filterbank

[This script](https://github.com/lixilinx/Best-practices-for-filterbanks/blob/main/what_is_a_filterbank.m) generates the following plot illustrating what is a DCT-IV modulated filterbank.

<img src="https://github.com/lixilinx/Best-practices-for-filterbanks/blob/main/what_is_a_filterbank.svg" width="500" />

From top to bottom, we have 1) the prototype filter; 2) the four cosine modulation sequences; 3) the four modulated filters, i.e., element-wise products between the prototype filter and modulation sequences; 4) frequency responses of the four modulated filters. We see that these four filters cover the whole frequency range nicely.

[These examples](https://github.com/lixilinx/FilterbanksBestPractices/tree/main#special-designs-wavelet-complex-wavelet-quadrature-mirror-filter-qmf-customized-nonsubsampled-and-nonuniform-filterbanks) illustrate a few special designs. The wavelet and complex wavelet examples also are good for explaining the concept of filterbank due to their simplicity.

Here, I mainly focus on the DFT modulated filterbanks. Still, all these modulated filterbanks share the same math, and only difference in modulation sequences. Notably, 
* The periodicity of modulation sequences induces the polyphase structure.
* Trigonometric modulation series make FFT like fast algorithms possible.

### MIRROR symmetry if latency is unconcerned

MIRROR symmetry imposes constraint $\pmb h = {\rm flip}(\pmb g)$, where $\pmb h$ and $\pmb g$ are the prototype analysis and synthesis filters, respectively. [This script](https://github.com/lixilinx/PracticalFilterbanks/blob/main/mirror_design_for_latency_insensitive_applications.m) generates the following design for applications insensitive to latency. With latency+1=filter length and symmetry setting either [1;0;0] or [0;0;0], the code finds this optimal design where analysis and synthesis filters mirror each other. Forcing symmetry=[-1;0;0] (analysis filter equals synthesis filter) leads to noticeable less sidelobe suppression as the resultant filter itself is symmetric, thus halves the continuous design freedoms.

MIRROR symmetry leads to the least aliasing and also the [least squares solution for signal reconstruction](https://github.com/lixilinx/FilterbanksBestPractices/tree/main#symmetries-and-least-squares-reconstruction), but brings two possible drawbacks: 
* Generally it is not good for low latency design as here, latency+1 must equal filter length.
* Typically neither the analysis nor the synthesis filter has good linear phase property (surely, their combined phase is linear as the analysis-synthesis filter is symmetric). 

<img src="https://github.com/lixilinx/Best-practices-for-filterbanks/blob/main/mirror_design_time.svg" width="400" />
<img src="https://github.com/lixilinx/Best-practices-for-filterbanks/blob/main/mirror_design_frequency.svg" width="400" />

### SAME symmetry if low latency is crucial

SAME symmetry imposes constraint $\pmb h=\pmb g$. [This script](https://github.com/lixilinx/PracticalFilterbanks/blob/main/low_latency_design_for_aec.m) generates the following low latency (latency+1 < filter length) design with symmetry=[-1;0;0], which suits realtime processings like acoustic echo cancellation (AEC), beamforming (BF), noise suppression (NS), etc. Clearly, MIRROR symmetry is incompatible with the low latency requirement. 

I gradually increase filter length until overshoot, i.e., bumpy mainlobe, arises. Another strategy is to start from a large filter length, and then gradually increase $\lambda$ to damp the overshoot if there is.

Setting symmetry=[0;0;0] releases all the design freedoms, and typically finds solutions with much lower aliasing. Lowering aliasing benefits certain applications like low, high or band pass filtering (LPF/HPF/BPF), sample-rate conversion (SRC), etc. However, symmetric design has nicer numerical properties and better fits other tasks like AEC. Another merit of the SAME symmetry is that both the analysis and synthesis filters must have approximately linear phase in the pass band as otherwise, nearly perfect reconstruction (NPR) is not possible. In the extreme case of latency+1 = filter length, the SAME symmetry leads to exactly linear phase prototype filters.  

<img src="https://github.com/lixilinx/Best-practices-for-filterbanks/blob/main/low_latency_design_for_aec.svg" width="400" /> 

### Symmetries and least squares reconstruction 
Symmetry is important due to the following numerical property:
* With the MIRROR symmetry, filterbank synthesis gives the least squares solution for fitting a given set of analysis results, whether consistent or not. 

Synthesis with the SAME symmetry generally only approximates the least squares solution for fitting a given set of analysis results due to the phase linearity in the pass band, not necessarily the stop band. If the phase also is linear in the stop band, then the DFT filterbank has both the MIRROR and SAME symmetries.
  
To show this property, let's consider a STFT with analysis and synthesis prototype filters $\pmb h=[h(0), h(1), \ldots]^T$ and $\pmb g=[g(0), g(1), \ldots]^T$, respectively. With $\pmb F$ being the normalized DFT matrix ($\pmb F^{H}\pmb F=\pmb I$), the STFT of a time signal $\pmb x=[x(0), x(1), \ldots]^T$ is given by 

$$
 \pmb X = \begin{bmatrix}
 \pmb F {\rm diag}({\rm flip}(\pmb h)) &  & \\
 & \pmb F {\rm diag}({\rm flip}(\pmb h)) & \\
 & & \ddots \\
\end{bmatrix} \pmb x  = \pmb H \pmb x, \qquad   ( {\rm what \\, shown \\, is \\, not\\, a \\, block \\,matrix, \\, just \\, the \\, layout \\, of \\, rows })
$$

The inverse STFT (ISTFT) on STFT results $\pmb Y$ is

$$
 \pmb y = \begin{bmatrix}
  {\rm diag}(\pmb g) \pmb F^{-1} &  & \\
 & {\rm diag}(\pmb g) \pmb F^{-1} & \\
 & & \ddots \\
\end{bmatrix} \pmb Y  = \pmb G \pmb Y, \qquad   ( {\rm what \\, shown \\, is \\, not\\, a \\, block \\,matrix, \\, just \\, the \\, layout \\, of \\, columns })
$$

By design, $\pmb G\pmb H = \pmb I$ (ignoring the delay and reconstruction errors). Now, let's consider the linear system  $\pmb H x = \pmb X$ with $\pmb x$ unknown. It is over-determined when $B < T $, and the least squares solution is given by $\hat{x}=(\pmb H^H \pmb H)^{-1} \pmb H^H \pmb X$. Is this $\hat{x}$ the same as $\pmb G \pmb X$? Not always. But, when $\pmb g = {\rm flip}(\pmb h)$, then indeed synthesis gives the least squares solution as shown by

$$ (\pmb H^H \pmb H)^{-1} H^H = (\pmb G\pmb H  )^{-1} \pmb G = \pmb G $$

In general, we can show that with the MIRROR symmetry, $\pmb G$ is the Moore–Penrose inverse of $\pmb H$, and synthesis $\hat{x} = \pmb G \pmb X$ gives the least squares solution for $\pmb x$ fitting equations $\pmb H \pmb x = \pmb X$. 

One application of this property is signal reconstruction from magnitude spectra $\pmb A \ge 0$. We are to solve the following optimization problem 

$$
\min_{\pmb s, \\, \pmb \phi} || \pmb H \pmb s - e^{\sqrt{-1} \pmb \phi} \pmb A ||^2
$$

with the Griffin-Lim algorithm, i.e., alternatively optimizing $\pmb s$ given $e^{\sqrt{-1} \pmb \phi} \pmb A$ (filterbank synthesis) and phase $\pmb \phi$ given $\pmb s$ (simply taking the phase of $\pmb H\pmb s$). Clearly, the Griffin-Lim algorithm requires MIRROR symmetry filterbanks to ensure its convergence. See this [audio stretching demo](https://github.com/lixilinx/FilterbanksBestPractices/blob/main/mag2sig.m) for illustration.     

### Designs with lower aliasing do not necessarily performs better for tasks like AEC

Lower aliasing is a necessary, but not sufficient, condition for better AEC performance. [This script](https://github.com/lixilinx/PracticalFilterbanks/blob/main/impulse_in_the_frequency_domain.m) shows what a naive time domain impulse looks like in the frequency domain. Intuitively, subband domain filters using designs with low aliasing but long prototype filters will need many taps to fit even simple time domain filters, causing slow convergence and high misadjustment. Symmetric and compact designs have better numerical condition number and are thus preferred. 

<img src="https://github.com/lixilinx/Best-practices-for-filterbanks/blob/main/impulse_in_the_frequency_domain.svg" width="400" />

### Squeezing synthesis filter does not help low latency filterbank design

Some manually designed low latency filterbanks have long and refined filters for analysis, and short and coarse filters for synthesis. Such a practice may not yield balanced and efficient designs as the analysis filters mainly focus on old samples, not the recent samples that are to be decomposed and re-synthesized. Aliasing caused by poor synthesis filters is another concern. Anyway, our time-frequency analysis cannot break the [uncertainty principle](https://en.wikipedia.org/wiki/Uncertainty_principle) The math is still the same. [This script](https://github.com/lixilinx/PracticalFilterbanks/blob/main/very_low_latency_design_for_hearing_aid.m) generates the following extremely low latency design suitable for applications like hearing aid. Note that since the oversampling ratio is high, I have halved the default mainlobe width by setting the cutoff frequency to $\omega_c = \pi/B/2$ to sharpen the frequency resolution.

Also see the [phase linearity discussions](https://github.com/lixilinx/FilterbanksBestPractices/tree/main#exact-phase-linearity-a-luxury-or-a-must-have) for more examples. 

<img src="https://github.com/lixilinx/Best-practices-for-filterbanks/blob/main/very_low_latency_design_for_hearing_aid.svg" width="400" />

### Replacing STFT with filterbank

STFT is a special DFT filterbank with prototype filter length equalling FFT size. Except for certain time-frequency analysis requiring nonnegative windows, generally we can replace STFT with filterbanks without much concern. [This code](https://github.com/lixilinx/PracticalFilterbanks/blob/main/replace_stft_with_fb.m) generates the following comparison results where both the filterbank and STFT are subject to the same latency, hop size and FFT size. Aliasings of the analysis (same for synthesis) filters of the STFT and filterbank are -14.5 dB and -35.7 dB, respectively. A virtually free design gain of 20+ dB!

The same argument holds for the [modified DCT (MDCT)](https://en.wikipedia.org/wiki/Modified_discrete_cosine_transform) as well. MDCT is a special DCT-IV filterbank with prototype filter length equalling half of the period of modulation sequences. We can get significantly better designs by increasing the prototype filter length while keeping everything else the same.    

<img src="https://github.com/lixilinx/PracticalFilterbanks/blob/main/replace_stft_with_fb.svg" width="400" />

### Spare that SRC

SRC can be expensive and invokes extra latency. But, a time domain SRC may be unnecessary in many cases. [This example](https://github.com/lixilinx/PracticalFilterbanks/blob/main/spare_that_SRC.m) shows how we can save a 16KHz $\rightarrow$ 48KHz SRC by analyzing with filter for 16 KHz design and synthesizing with the one for 48 KHz design. The trick is that [this code](https://github.com/lixilinx/PracticalFilterbanks/blob/main/low_latency_design_for_aec.m) designs the two filters from initial guess of the same shape, and thus they are compatible.

Actually, resampling in the frequency domain could be easier (with the help of a good FFT lib) than in the time domain when the conversion ratio is complicated, e.g., 44.1KHz $\leftrightarrow$ 16KHz. We just need to switch to the matched synthesis filters that are prepared in advance. The same argument holds for the analysis side too. Most likely the SRC before filterbank analysis can be saved as well, just by switching to the analysis filters matching the original sample rate.

### Frequency shift with analytic signal

Frequency shift is a common operation to cut off the acoustic feedback in systems like hearing aid and public address (PA) systems. One possible choice is to shift the signals in the frequency domain directly. This approximate way may introduce noticeable harmonic distortions when either the amount of shift or the hop size is large. [This code](https://github.com/lixilinx/PracticalFilterbanks/blob/main/frequency_shift_with_analytic_signal.m) shows an alternative practice: first split the signal into two analytic parts; then shift the low frequency part less, and high frequency part more. For this example, we see that the absolute normalized cross correlation between input and output reduces to about zero by tens Hz of shifting while without any audible artifacts.

### Remove the systemic phase biases

One possible neglect when dealing with the phases, say phase unwrapping, is to ignore the systemic or structural bias of the phases for frequency domain signals obtained by a DFT filterbank. I summarize these biases in the card below. Note that $E$ here takes ensemble average, or time average for a stationary ergodic process. Generally, these systemic biases cannot be put into a single smooth function of frame and bin indices as if assuming so, the resultant partial differential equation (PDE) is inconsistent with arbitrarily fine grids of time and frequency. This inconsistency is not a surprise: the phase representation is already succinct without assuming any further specific structures of the time domain signals, e.g., a pure tone or a certain type of chirp.

<img src="https://github.com/lixilinx/PracticalFilterbanks/blob/main/nature_of_phase.svg" width="400" />

[This script](https://github.com/lixilinx/PracticalFilterbanks/blob/main/structural_bias_of_phase.m) compares measured against predicted phase differences to generate the results as below. They do match well.

<img src="https://github.com/lixilinx/PracticalFilterbanks/blob/main/measured_and_predicted_phases.svg" width="400" />

[The same script](https://github.com/lixilinx/PracticalFilterbanks/blob/main/structural_bias_of_phase.m) also demonstrates how to exploit these structural biases for accurate phase unwrapping. Removal of these biases lets us deal with the inherent transient properties of phases, and leads to consistently and significantly cleaner phase unwrapping across all the frequencies.

<img src="https://github.com/lixilinx/PracticalFilterbanks/blob/main/phase_unwrapping_std.svg" width="400" />

### DCT filterbank for tasks like AEC?

DCT filterbanks like [the MDCT](https://en.wikipedia.org/wiki/Modified_discrete_cosine_transform) are widely used in subband codecs. Occasionally these critically decimated DCT filterbanks are used for other tasks like AEC. This may not work well as the aliasing due to downsampling is too high. Still, an oversampled DCT filterbank works for AEC as well as any good oversampled DFT filterbanks. [This script](https://github.com/lixilinx/FilterbanksBestPractices/blob/main/dct_filterbank_for_aec.m) compares a critically decimated and an oversampled DCT filterbanks to produce the following plot. The critically decimated one clearly suffers a lot from aliasing, while the oversampled one does not.

<img src="https://github.com/lixilinx/FilterbanksBestPractices/blob/main/dct_filterbank_for_aec.svg" width="400" />

### Special designs: wavelet, complex wavelet, quadrature mirror filter (QMF), customized, nonsubsampled and nonuniform filterbanks

Discrete wavelet (DWT) and QMF are DFT filterbanks with $T=2$, thus having the two periodic modulation sequences: $[1,\\,1]$ for low pass (LP) filters and $[1, \\,-1]$ for high pass (HF) filters. This observation lets us design all the possible flavors of QMFs: critically decimated ($B=2$) or oversampled ($B=1$), different symmetries, low latency and phase linearity. See the [QMF examples](https://github.com/lixilinx/FilterbanksBestPractices/blob/main/qmf_examples.m) for details.   

DWT is sensitive to signal shift. By replacing the two modulation sequences in DWT with $[1,\\,\sqrt{-1},\\,-1,\\,-\sqrt{-1}]$ and $[1,\\,-\sqrt{-1},\\,-1,\\,\sqrt{-1}]$, we obtain the complex wavelet which decomposes the signal into positive and negative frequency parts, where $\sqrt{-1}$ is the imaginary unit. The envelope of complex wavelet coefficients are less sensitive to signal shift as the phase compensates for it. Unlike the dual-tree complex DWT that can only be redundant, our complex DWT designs support critical sampling, see the [complex wavelet](https://github.com/lixilinx/FilterbanksBestPractices/blob/main/critically_sampled_complex_wavelet.m) example.  

Nonsubsampled filterbanks are the ones with $B=1$, i.e., no downsampling after analysis. See the [QMF examples](https://github.com/lixilinx/FilterbanksBestPractices/blob/main/qmf_examples.m) and the [nonuniform design example](https://github.com/lixilinx/FilterbanksBestPractices/blob/main/a_nonuniform_design.m).

Nonuniform filterbanks can be built from uniform ones by either splitting certain bands, e.g., the wavelet packet decomposition (WPD), or merging bands under certain conditions, e.g., shown in this [nonuniform design example](https://github.com/lixilinx/FilterbanksBestPractices/blob/main/a_nonuniform_design.m). 

The same method also supports customized designs. See the [nonuniform design example](https://github.com/lixilinx/FilterbanksBestPractices/blob/main/a_nonuniform_design.m), where the modulation sequences are that of DFT's shifted by $0.5$ bin to get rid of the real-valued DC and Nyquist bins in an ordinary DFT filterbanks so that all bins are complex. Following this demo code, we identify the period of the modulation sequences, the $\pmb\Gamma$ matrix, and then proceed as usual to design the filterbanks.    

The [QMF examples](https://github.com/lixilinx/FilterbanksBestPractices/blob/main/qmf_examples.m) generates the following sample designs:  

<img src="https://github.com/lixilinx/FilterbanksBestPractices/blob/main/qmf_examples.svg" width="400" />

The [complex wavelet](https://github.com/lixilinx/FilterbanksBestPractices/blob/main/critically_sampled_complex_wavelet.m) generate this sample design:

<img src="https://github.com/lixilinx/FilterbanksBestPractices/blob/main/critically_sampled_complex_wavelet.svg" width="400" />

The [nonuniform example](https://github.com/lixilinx/FilterbanksBestPractices/blob/main/a_nonuniform_design.m) generates this sample design:

<img src="https://github.com/lixilinx/FilterbanksBestPractices/blob/main/a_nonuniform_design.svg" width="600" />

### Optimization of discrete design freedoms

Filterbank design is an optimization problem loaded with local minima. We are almost sure to find the global minima for the continuous design freedoms by starting the filter coefficients from diverse random initial guesses. But, those discrete design freedoms are more tricky to optimize.    

Certain discrete design freedoms, say the phase shift pair $(i, j)$ in a DFT filterbank, have closed-form solutions. Other discrete design freedoms have clear boundaries, say latency+1 cannot be smaller than block size in any causal designs, latency+1 must equal filter length when MIRROR symmetry is enforced, $B$ cannot be smaller than $T$, ... The problem is how to find the optimal designs subject to certain constraints, say a maximum allowable latency. [This script](https://github.com/lixilinx/FilterbanksBestPractices/blob/main/optimization_of_discrete_design_freedoms.m) studies DFT filterbank designs with a fixed $T$ and filter length, and various $B$ and latency, to generate the following results. Nearly perfect reconstruction (NPR) always is feasible for each design point, but the quality of designs could jitter a lot with respect to latency. For the most common case of $T=2B$, local minima for latency are at points with mod(latency+1, $T$)=0. Generally, there are no simple recipes for the optimization of these discrete design freedoms. A good strategy is to sweep these discrete design freedoms on a scaled down design problem to find out the best recipe, and then apply the found optimal recipe to the full scale target design. 

<img src="https://github.com/lixilinx/FilterbanksBestPractices/blob/main/optimization_of_discrete_design_freedoms.svg" width="400" />

Let's consider a practical design example where the block size is $B=10$ ms so that our processing can connect with most audio codecs without an extra ping pong buffer. The FFT size is $T=16$ ms, resulting in an oversample ratio of $8/5$. With SAME symmetry and filter length $L=4T$, there are three optimal latencies: 16 ms, 40 ms, and 54 ms. Optimal latency 40 ms gives a good balance between performance and delay. [This script](https://github.com/lixilinx/FilterbanksBestPractices/blob/main/a_design_with_optimized_latency.m) generates the following comparison results. Design with the SAME symmetry clearly outperforms that with MIRROR symmetry when subjecting to the same latency.    

<img src="https://github.com/lixilinx/FilterbanksBestPractices/blob/main/a_design_with_optimized_latency.svg" width="400" />

### Exact phase linearity, a luxury or a must-have?

A linear phase filterbank suggests that all the modulated analysis and synthesis filters have linear phases, i.e., either symmetric or anti-symmetric. Hence, a linear phase prototype filter does not always imply a linear phase filterbank. For example, the [MDCT](https://en.wikipedia.org/wiki/Modified_discrete_cosine_transform), with either the rectangular or sine windows, is not a linear phase DCT filterbank, although the prototype filters, i.e.,  the rectangular or sine windows, are symmetric. Nevertheless, a DFT filterbank with a linear phase prototype filter is always a linear phase filterbank. 

For audio processing, our human ears are deaf to absolute phases. Thus, phase linearity might not be a must-have. For image processing, noticeable phase distortion is not acceptable. But, image processing is acausal as all the neighbors of any pixel are available for filtering. Thus, filters with large latency but linear phase are affordable and good to have there.  

[This script](https://github.com/lixilinx/FilterbanksBestPractices/blob/main/linear_phase_can_be_a_luxury.m) generates the following comparison results to show how it looks like to impose linear phase and low latency constraints at the same time. Here, I set symmetry=[0;1;1] to force both the analysis and synthesis filters to have linear phases, and different filter lengths and weights to get high resolution analysis filters and low resolution synthesis filters. The performance gap between linear phase and SAME symmetry designs is huge. Not a surprise. The SAME symmetry only forces phase linearity in the pass band, while symmetry=[0;1;1] forces the phase to be linear in both the pass and stop bands, thus wasteful.       

<img src="https://github.com/lixilinx/FilterbanksBestPractices/blob/main/linear_phase_can_be_a_luxury.svg" width="400" />

Exact phase linearity is such a strong constraint and sometimes cannot coexist with the perfect reconstruction (PR) condition. Two examples:
* By revisiting the frequency domain PR condition, e.g., [here](https://ccrma.stanford.edu/~jos/sasp/Linear_Phase_Quadrature_Mirror.html), we see that a critically decimated linear phase QMF filter with odd filter lengths (or even order) will have gain $0$ at $\omega=\pi/2$, which contradicts the PR condition. Still, odd or even order linear phase oversampled QMF always exists. 
* For a DCT-IV filterbank with filter length $T$, the first cosine modulation sequences can be symmetric or anti-symmetric with a proper shift $i$, but always has zero mean. Then with a symmetric prototype filter, the first modulated filter must have zero DC gain, and thus fails the PR condition at $\omega=0$.

Nevertheless, it is trivial to have approximately linear phase in the pass band for either the DFT or DCT filterbank. Luckily, that could be sufficient in most cases.  

### Fine balance between resolution and leakage

In the same way as balancing the spectral resolution and leakage in window designs, say Hann vs Hamming, we also provide such choices for filterbank design. The following card summarizes the common stopband losses for filterbank design. One also can come up with his own customized loss. The default choice is momentum $=0$, and two more choices I have implemented are $1$ and $2$.       

<img src="https://github.com/lixilinx/FilterbanksBestPractices/blob/main/smallest_momentum_designs.svg" width="400" />

[This script](https://github.com/lixilinx/FilterbanksBestPractices/blob/main/different_min_momentum_designs.m) generates the following comparison results showing the tradeoff between resolution and leakage. The default setting, i.e., minimizing the 0th momentum, has the narrowest mainlobe, but largest sidelobe tails. The design minimizing the 2nd momentum, i.e., the momentum of inertia, has the widest mainlobe, but significantly lower sidelobe tails.     

<img src="https://github.com/lixilinx/FilterbanksBestPractices/blob/main/different_min_momentum_designs.svg" width="500" />

### Resources 

1, [Periodic sequences modulated filter banks](https://ieeexplore.ieee.org/document/8304771), IEEE SPL, 2018. (Notation and math are from this paper. A typo in Table I: $\pmb J$ in the $\pmb \Gamma$ matrix for GDFT modulation should be $\pmb I$.)

2, The matlab filterbank design code is slightly adapted from [the original version here](https://sites.google.com/site/lixilinx/home/psmfb). I also have [Tensorflow and Pytorch implementations](https://github.com/lixilinx/Filterbank) of the DFT modulated filterbank, not polished but works, useful for end2end differentiable DSP designs.
