function fb = FilterBankStruct( )
% Written by Xi-Lin Li, lixilinx@gmail.com
% Returns an empty filter bank structure with the following fields
%
% Gamma: the Gamma matrix \tilde{W}*W
% T: minimum period of modulation sequences 
% B: decimation ratio, or block size, or frame size, or hop size in STFT/MDCT
% tau0: system delay
% i: circular shift amount in S(i-1)
% j: circular shift amount in S(j)
% h: analysis prototype filter
% g: synthesis prototype filter
% w_cut: cut-off frequency for prototype filters
% zeta: relative design weight on synthesis filter's stop band energy
% symmetry: symmetry code with meaning:
%           symmetry(1)=1 means g(t)=h(L+1-t) when g and h have the same length L
%           symmetry(1)=-1 means g(t)=h(t) when g and h have the same length L
%           symmetry(2)=1 means h(t)=h(length(h)+1-t)
%           symmetry(3)=1 means g(t)=g(length(g)+1-t)
%           set it to 0 if no symmetry constraint is required
%
fb.Gamma=[]; 
fb.T=[]; 
fb.B=[]; 
fb.tau0=[]; 
fb.i=[]; 
fb.j=[]; 
fb.h=[]; 
fb.g=[]; 
fb.w_cut=[]; 
fb.zeta=[]; 
fb.symmetry=[]; 