# Best practices for filterbanks

Filterbank can be tricky. I repeatedly saw improper practices of it even from the most experienced DSP engineers. Hence, I’d like to share what I have learned and am continuing to learn from my many years of  experiences. The math and design tool are from [my paper](https://ieeexplore.ieee.org/document/8304771).  

### What is a filterbank

[This script](https://github.com/lixilinx/Best-practices-for-filterbanks/blob/main/what_is_a_filterbank.m) generates the following plot showing what is a DCT-IV modulated filterbank.

<img src="https://github.com/lixilinx/Best-practices-for-filterbanks/blob/main/what_is_a_filterbank.svg" width="500" /> 

From top to bottom, we have 1) the prototype filter; 2) the four cosine modulation sequences; 3) the four modulated filters, i.e., element-wise products between the prototype filter and modulation sequences; 4) frequency responses of the four modulated filters. We see that these four filters cover the whole frequency range evenly, i.e., a filterbank. 

Here, I mainly focus on the DFT modulated filterbanks. Still, all these modulated filterbanks share the same math, and only diff by modulated sequences. 

### MIRROR symmetry if latency is unconcerned  

[This script](https://github.com/lixilinx/PracticalFilterbanks/blob/main/mirror_design_for_latency_insensitive_applications.m) generates the following designed filterbank for applications where latency is unconcerned, e.g., vocoder. With symmetry setting either [1;0;0] or [0;0;0], the code finds this optimal design where analysis and synthesis filters mirror each other. Forcing symmetry=[-1;0;0] (analysis filter equals synthesis filter) leads to noticeable worse design as the resultant filter is symmetric and half of the design freedoms are gone.         

<img src="https://github.com/lixilinx/Best-practices-for-filterbanks/blob/main/mirror_design_time.svg" width="400" /> 
<img src="https://github.com/lixilinx/Best-practices-for-filterbanks/blob/main/mirror_design_frequency.svg" width="400" /> 

### SAME symmetry if low latency is crucial 

[This script](https://github.com/lixilinx/PracticalFilterbanks/blob/main/mirror_design_for_latency_insensitive_applications.m) generates the following low latency design with symmetry=[-1;0;0], which suits realtime processings like acoustic echo cancellation (AEC), beamforming (BF), noise suppression (NS), etc. 

I gradually increase filter length until overshoot, i.e., bumpy mainlobe, arises. Another strategy is to start from a large filter length and gradually increase lambda to damp the overshoot.

Setting symmetry=[0;0;0] releases all the design freedoms, and typically finds designs with much lower aliasing (around 20 dB lower sidelobes for this example!). Lower aliasing benefits simple applications like low or high pass filtering (LPF/HPF), sample-rate conversion (SRC), etc. However, a symmetric design has better numerical properties and better fits more complicated applications like AEC.    

<img src="https://github.com/lixilinx/Best-practices-for-filterbanks/blob/main/low_latency_design_for_aec.svg" width="400" /> 

### Designs with lower aliasing do not necessarily performs better for tasks like AEC 

Lower aliasing is a necessary, but not sufficient, condition for better AEC performance. [This script](https://github.com/lixilinx/PracticalFilterbanks/blob/main/impulse_in_the_frequency_domain.m) shows what a naive time domain impulse looks like in the frequency domain. Frequency domain filters using designs with low aliasing but long prototype filters will need many taps to fit even simple time domain filters, causing slow convergence and high misadjustment. Symmetric and compact designs have better numerical condition number and are thus preferred. Nevertheless, applications like LPF, HPF, SRC, etc., benefit from low aliasing.   

<img src="https://github.com/lixilinx/Best-practices-for-filterbanks/blob/main/impulse_in_the_frequency_domain.svg" width="400" /> 

### Squeezing synthesis filter does not help low latency filterbank design

Many manual low latency designs use a long and refined filter for analysis, and a short and coarse filter for synthesis. This practice cannot give a balanced design as the analysis focuses on old samples, but ignores recent samples, thus wasteful. [This script](https://github.com/lixilinx/PracticalFilterbanks/blob/main/very_low_latency_design_for_hearing_aid.m) generates the following extremely low latency design suitable for applications like hearing aid. Note that since the oversampling ratio is high, I narrow the mainlobe by setting the cutoff frequency to pi/B/2.      

<img src="https://github.com/lixilinx/Best-practices-for-filterbanks/blob/main/very_low_latency_design_for_hearing_aid.svg" width="400" /> 

### Spare that SRC

SRC can be expensive and invokes extra latency. But, a time domain SRC is unnecessary in most cases. [This example](https://github.com/lixilinx/PracticalFilterbanks/blob/main/spare_that_SRC.m) shows how we can save a 16KHz-to-48KHz SRC by analyzing with filter for 16 KHz design and synthesizing with the one for 48 KHz design. The trick is that [this code](https://github.com/lixilinx/PracticalFilterbanks/blob/main/mirror_design_for_latency_insensitive_applications.m) design the two filters from initial guess with the same shape, and thus they are compatible. 

Actually, resampling in the frequency domain could be easier than in the time domain when the conversion ratio is complicated, e.g., 44.1KHz-to-16KHz. We just need to switch the synthesis filters that are prepared in advance.        

### Frequency shift on analytic signal only 

Frequency shift is a basic operation to cut off the acoustic feedback in systems like hearing aid and public address (PA) systems. I saw that many implementations shift the signals in the frequency domain. This practice is improper as the synthesis filters are not shifted accordingly, and thus causing artifacts like beats even with shifting of a few Hz. [This code](https://github.com/lixilinx/PracticalFilterbanks/blob/main/frequency_shift_with_analytic_signal.m) shows baseline practice: first split the signal into two analytic parts; then shift the low frequency part less, and high frequency part more. For this example, we see that tens of Hz of shifting is enough to reduce the absolute normalized correlation between input and output to about 0.1 without audible artifacts.        

