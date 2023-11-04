# Best practices for filterbanks

Filterbank can be tricky. I repeatedly saw improper practices of it even from the most experienced DSP engineers. Hence, I’d to share what I have learned and am keeping learning from my years of work experiences. The math and design tool are from [my paper](https://ieeexplore.ieee.org/document/8304771).  

### What is a filterbank

[This script](https://github.com/lixilinx/Best-practices-for-filterbanks/blob/main/what_is_a_filterbank.m) generates the following plot showing what is a DCT-IV modulated filterbank.

<img src="https://github.com/lixilinx/Best-practices-for-filterbanks/blob/main/what_is_a_filterbank.svg" width="500" /> 

From top to bottom, we have 1) the prototype filter; 2) the four cosine modulation sequences; 3) the four modulated filters, i.e., element-wise products between the prototype filter and modulation sequences; 4) frequency responses of the four modulated filters. We see that these four filters cover the whole frequency range evenly, i.e., a filterbank. 

Here, I mainly focus on the DFT modulated filterbanks. Still, all these modulated filterbanks share the same math, and only diff by modulated sequences. 

### MIRROR symmetry if latency is unconcerned  

<img src="https://github.com/lixilinx/Best-practices-for-filterbanks/blob/main/mirror_design_time.svg" width="400" /> 
<img src="https://github.com/lixilinx/Best-practices-for-filterbanks/blob/main/mirror_design_frequency.svg" width="400" /> 

### SAME symmetry if low latency is crucial 

<img src="https://github.com/lixilinx/Best-practices-for-filterbanks/blob/main/low_latency_design_for_aec.svg" width="400" /> 

### Designs with lower aliasing do not necessarily performs better for tasks like AEC 

<img src="https://github.com/lixilinx/Best-practices-for-filterbanks/blob/main/impulse_in_the_frequency_domain.svg" width="400" /> 

### Squeezing synthesis filter does not help low latency filterbank design

<img src="https://github.com/lixilinx/Best-practices-for-filterbanks/blob/main/very_low_latency_design_for_hearing_aid.svg" width="400" /> 

### Spare that SRC
### Frequency shift on analytic signal only 

