# Best practices with filterbanks

Filterbank is versatile but also tricky. I repeatedly see improper practices of it even from the most experienced DSP engineers. Here, I’d like to share what I have learned and am continuing to learn from my many years of practice. The math and design method are from [my paper](https://ieeexplore.ieee.org/document/8304771). Unlike the textbook styles, what makes life a little easier is that I treat all the uniform filterbanks, include special cases like QMF and nonsubsampled ones, with a single straightforward math framework. Topics:

[What is a filterbank](https://github.com/lixilinx/FilterbanksBestPractices/tree/main#what-is-a-filterbank)

[MIRROR symmetry if latency is unconcerned](https://github.com/lixilinx/FilterbanksBestPractices/tree/main#mirror-symmetry-if-latency-is-unconcerned)

[SAME symmetry if low latency is crucial](https://github.com/lixilinx/FilterbanksBestPractices/tree/main#same-symmetry-if-low-latency-is-crucial)

[Designs with lower aliasing do not necessarily performs better for tasks like AEC](https://github.com/lixilinx/FilterbanksBestPractices/tree/main#designs-with-lower-aliasing-do-not-necessarily-performs-better-for-tasks-like-aec)

[Squeezing synthesis filter does not help low latency filterbank design](https://github.com/lixilinx/FilterbanksBestPractices/tree/main#squeezing-synthesis-filter-does-not-help-low-latency-filterbank-design)

[Replacing STFT with filterbank](https://github.com/lixilinx/FilterbanksBestPractices/tree/main#replacing-stft-with-filterbank)

[Spare that SRC](https://github.com/lixilinx/FilterbanksBestPractices/tree/main#spare-that-src)

[Frequency shift on analytic signal only](https://github.com/lixilinx/FilterbanksBestPractices/tree/main#frequency-shift-on-analytic-signal-only)

[Remove the systemic phase biases](https://github.com/lixilinx/FilterbanksBestPractices/tree/main#remove-the-systemic-phase-biases)

[DCT filterbank for tasks like AEC?](https://github.com/lixilinx/FilterbanksBestPractices/tree/main#dct-filterbank-for-tasks-like-aec)

[Special designs: wavelet, quadrature mirror filter (QMF), nonsubsampled filterbank](https://github.com/lixilinx/FilterbanksBestPractices/blob/main/README.md#special-designs-wavelet-quadrature-mirror-filter-qmf-nonsubsampled-filterbank)

### What is a filterbank

Aside from many great resources like textbooks, [this script](https://github.com/lixilinx/Best-practices-for-filterbanks/blob/main/what_is_a_filterbank.m) generates the following plot illustrating what is a DCT-IV modulated filterbank.

<img src="https://github.com/lixilinx/Best-practices-for-filterbanks/blob/main/what_is_a_filterbank.svg" width="500" />

From top to bottom, we have 1) the prototype filter; 2) the four cosine modulation sequences; 3) the four modulated filters, i.e., element-wise products between the prototype filter and modulation sequences; 4) frequency responses of the four modulated filters. We see that these four filters cover the whole frequency range nicely.

[These examples](https://github.com/lixilinx/FilterbanksBestPractices/blob/main/README.md#special-designs-wavelet-quadrature-mirror-filter-qmf-nonsubsampled-filterbank) consider the simplest filterbanks, i.e., DFT ones with period T=2. They also are good for illustrating the concept of filterbank.

Here, I mainly focus on the DFT modulated filterbanks. Still, all these modulated filterbanks share the same math, and only difference in modulation sequences. Notably, 1) the periodicity of modulation sequences induces the polyphase structure; 2) trigonometric modulation series make FFT like fast algorithms possible.

### MIRROR symmetry if latency is unconcerned

[This script](https://github.com/lixilinx/PracticalFilterbanks/blob/main/mirror_design_for_latency_insensitive_applications.m) generates the following design for applications insensitive to latency. With symmetry setting either [1;0;0] or [0;0;0], the code finds this optimal design where analysis and synthesis filters mirror each other. Forcing symmetry=[-1;0;0] (analysis filter equals synthesis filter) leads to noticeable worse design as the resultant filter is symmetric, thus halves the design freedoms.

<img src="https://github.com/lixilinx/Best-practices-for-filterbanks/blob/main/mirror_design_time.svg" width="400" />
<img src="https://github.com/lixilinx/Best-practices-for-filterbanks/blob/main/mirror_design_frequency.svg" width="400" />

### SAME symmetry if low latency is crucial

[This script](https://github.com/lixilinx/PracticalFilterbanks/blob/main/low_latency_design_for_aec.m) generates the following low latency design with symmetry=[-1;0;0], which suits realtime processings like acoustic echo cancellation (AEC), beamforming (BF), noise suppression (NS), etc.

I gradually increase filter length until overshoot, i.e., bumpy mainlobe, arises. Another strategy is to start from a large filter length, and then gradually increase lambda to damp the overshoot if there is.

Setting symmetry=[0;0;0] releases all the design freedoms, and typically finds solutions with much lower aliasing (around 20 dB lower sidelobes in this example!). Lowering aliasing benefits certain applications like low, high or band pass filtering (LPF/HPF/BPF), sample-rate conversion (SRC), etc. However, symmetric design has nicer numerical properties and better fits other tasks like AEC.

<img src="https://github.com/lixilinx/Best-practices-for-filterbanks/blob/main/low_latency_design_for_aec.svg" width="400" /> 

### Designs with lower aliasing do not necessarily performs better for tasks like AEC

Lower aliasing is a necessary, but not sufficient, condition for better AEC performance. [This script](https://github.com/lixilinx/PracticalFilterbanks/blob/main/impulse_in_the_frequency_domain.m) shows what a naive time domain impulse looks like in the frequency domain. Intuitively, subband domain filters using designs with low aliasing but long prototype filters will need many taps to fit even simple time domain filters, causing slow convergence and high misadjustment. Symmetric and compact designs have better numerical condition number and are thus preferred. Nevertheless, aliasing alone is a good performance index for applications like LPF, HPF, BPF, SRC, etc.

<img src="https://github.com/lixilinx/Best-practices-for-filterbanks/blob/main/impulse_in_the_frequency_domain.svg" width="400" />

### Squeezing synthesis filter does not help low latency filterbank design

Many manually designed low latency filterbanks have long and refined filters for analysis, and short and coarse filters for synthesis. Such a practice may not yield balanced and efficient designs as the analysis filters mainly focus on old samples, not the recent samples that are to be decomposed and re-synthesized. Anyway, our time-frequency analysis cannot break the [uncertainty principle](https://en.wikipedia.org/wiki/Uncertainty_principle) The math is still the same. [This script](https://github.com/lixilinx/PracticalFilterbanks/blob/main/very_low_latency_design_for_hearing_aid.m) generates the following extremely low latency design suitable for applications like hearing aid. Note that since the oversampling ratio is high, I have halved the default mainlobe width by setting the cutoff frequency to pi/B/2 to sharpen the frequency resolution.

<img src="https://github.com/lixilinx/Best-practices-for-filterbanks/blob/main/very_low_latency_design_for_hearing_aid.svg" width="400" />

### Replacing STFT with filterbank

STFT is a special filterbank with prototype filter length equalling FFT size. Except for certain time-frequency analysis requiring nonnegative windows, generally we can replace STFT with filterbanks without much concern. [This code](https://github.com/lixilinx/PracticalFilterbanks/blob/main/replace_stft_with_fb.m) generates the following comparison results where both the filterbank and STFT are subject to the same latency, hop size and FFT size. Aliasings of the analysis (same for synthesis) filters of the STFT and filterbank are -14.5 dB and -35.7 dB, respectively. A virtually free design gain of 20+ dB!

<img src="https://github.com/lixilinx/PracticalFilterbanks/blob/main/replace_stft_with_fb.svg" width="400" />

### Spare that SRC

SRC can be expensive and invokes extra latency. But, a time domain SRC may be unnecessary in many cases. [This example](https://github.com/lixilinx/PracticalFilterbanks/blob/main/spare_that_SRC.m) shows how we can save a 16KHz-to-48KHz SRC by analyzing with filter for 16 KHz design and synthesizing with the one for 48 KHz design. The trick is that [this code](https://github.com/lixilinx/PracticalFilterbanks/blob/main/low_latency_design_for_aec.m) designs the two filters from initial guess of the same shape, and thus they are compatible.

Actually, resampling in the frequency domain could be easier (with a good FFT lib) than in the time domain when the conversion ratio is complicated, e.g., 44.1KHz-to-16KHz. We just need to switch to the matched synthesis filters that are prepared in advance. The same argument applies to the analysis side. Most likely the SRC before filterbank analysis can be saved as well, just by switching to the analysis filters matching the original sample rate.

### Frequency shift on analytic signal only

Frequency shift is a basic operation to cut off the acoustic feedback in systems like hearing aid and public address (PA) systems. I saw that many implementations shift the signals in the frequency domain. This practice is improper as the synthesis filters are not shifted accordingly, and this misalignment may cause artifacts like beats even with shift of a few Hz. [This code](https://github.com/lixilinx/PracticalFilterbanks/blob/main/frequency_shift_with_analytic_signal.m) shows a baseline practice: first split the signal into two analytic parts; then shift the low frequency part less, and high frequency part more. For this example, we see that the absolute normalized cross correlation between input and output reduces to about 0.1 for (2, 20) Hz shift and 0.0 for (3, 30) Hz shift, while without any audible artifacts.

### Remove the systemic phase biases

One common neglect when dealing with the phases, say phase unwrapping, is to ignore the systemic or structural bias of phases for most natural signals. I summarize these biases in the figure below. Note that E here takes ensemble average, or time average for a stationary ergodic process. Generally, these systemic biases cannot be put into a single smooth function of frame and bin indices as the resultant partial differential equation (PDE) is inconsistent with arbitrarily fine grids of time and frequency. This inconsistency should not be a surprise: the phase representation is already succinct without assuming any further structures of the signals.

<img src="https://github.com/lixilinx/PracticalFilterbanks/blob/main/nature_of_phase.svg" width="400" />

[This script](https://github.com/lixilinx/PracticalFilterbanks/blob/main/structural_bias_of_phase.m) compares measured against predicted phase differences to generate the results as below. They do match well.

<img src="https://github.com/lixilinx/PracticalFilterbanks/blob/main/measured_and_predicted_phases.svg" width="400" />

[The same script](https://github.com/lixilinx/PracticalFilterbanks/blob/main/structural_bias_of_phase.m) also demonstrates how to exploit this knowledge for accurate phase unwrapping. Removal of these biases lets us deal with the inherent transient properties of phases, and leads to consistently and significantly cleaner unwrapping across all the frequencies.

<img src="https://github.com/lixilinx/PracticalFilterbanks/blob/main/phase_unwrapping_std.svg" width="400" />

### DCT filterbank for tasks like AEC?

DCT filterbanks like lapped transform are widely used in subband codec. One common mistake for beginners is to use those critically decimated DCT filterbanks, designed for subband codec, for other tasks like AEC. This cannot work well as the aliasing due to downsampling is too high. Then there is another myth stating that DCT filterbanks cannot be used for AEC. Neither is this true. Actually, an oversampled DCT filterbank works for AEC as well as any good oversampled DFT filterbanks. [This script](https://github.com/lixilinx/FilterbanksBestPractices/blob/main/dct_filterbank_for_aec.m) compares critically decimated and oversampled DCT filterbanks to produce the following plot. The critically decimated one clearly suffers a lot from aliasing, while the oversampled one does not.

<img src="https://github.com/lixilinx/FilterbanksBestPractices/blob/main/dct_filterbank_for_aec.svg" width="400" />

### Special designs: wavelet, quadrature mirror filter (QMF), nonsubsampled filterbank

Discrete wavelet and QMF are filterbanks with the two periodic modulation sequences: [...,1,1,1,1,...] and [...,1,-1,1,-1,...]. Thus, they are DFT modulated filterbanks with period T=2. Nonsubsampled filterbanks are the ones with B=1, i.e., no downsampling. All these special filterbanks can be designed with the same code, providing the same knobs for tweaking. No extra math is needed.

[This script](https://github.com/lixilinx/FilterbanksBestPractices/blob/main/some_special_designs.m) demonstrates such designs to have the following figure. When T=2, I only show the low pass (LP) analysis filter. The high pass (HP) filter is obtained just by modulating the LP one with sequence [...,1,-1,1,-1;...], i.e., alternative sign change. Hence, their frequency responses mirror each other around pi/2. The nonsubsampled one is free of aliasing. These filterbanks are the building blocks for wavelet packet decomposition (WPD), with or without decimation at any resolution level. 

<img src="https://github.com/lixilinx/FilterbanksBestPractices/blob/main/some_special_designs.svg" width="400" />

### Refs

1, [Periodic sequences modulated filter banks](https://ieeexplore.ieee.org/document/8304771), IEEE SPL, 2018.

2, The matlab filterbank design code is the same as [here](https://sites.google.com/site/lixilinx/home/psmfb). I also have [Tensorflow and Pytorch implementations](https://github.com/lixilinx/Filterbank) of the DFT modulated filterbank, not polished but works, useful for end2end differential designs.
