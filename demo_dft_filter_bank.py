"""Random irregular DFT modulated filterbank example. 
"""
import numpy as np
from filterbank_design import design_filter_bank, DesignResult, FilterBank
from matplotlib import pyplot as plt 

rng = np.random.default_rng()
T = int(rng.integers(10, 21))
B = int(rng.integers(T//2, 2*T//3)) 
Lh = int(rng.integers(2*T, 3*T))
Lg = int(rng.integers(2*T, 3*T))
tau0 = int(rng.integers(T - 1, min(Lh, Lg)))
momentum = int(rng.integers(0, 3))
eta = 1e5
lambda_ = 1.0

best_result: DesignResult | None = None
for trial in range(1, 11):
    h, g = rng.uniform(size=Lh), rng.uniform(size=Lg) # rand initialization is better than fbd_random_initial_guess.m
    result = design_filter_bank(FilterBank(T=T, B=B, tau0=tau0, h=h, g=g, momentum=momentum), 
                                eta=eta, lambda_=lambda_, max_iterations=100)
    print(f"Trial: {trial}; design cost: {result.cost}")
    if best_result is None or result.cost < best_result.cost:
        best_result = result 

refined = design_filter_bank(best_result.filter_bank, eta=eta, lambda_=lambda_, max_iterations=1000)
print(f"Refined cost: {refined.cost}; reconstruction error: {refined.reconstruction_error}")

fb = refined.filter_bank
h, g = fb.h, fb.g
shift_i, shift_j = fb.i, fb.j
x = rng.normal(size=10*B)
x[np.abs(x) < 1] = 0
x = np.sign(x)
padded_h = np.pad(h, (0, (-Lh) % T))
padded_g = np.pad(g, (0, (-Lg) % T))
ana_bfr = np.zeros_like(padded_h)
syn_bfr = np.zeros_like(padded_g)
y = np.zeros_like(x)
for block_start in range(0, x.size - B + 1, B):
    ana_bfr[:-B] = ana_bfr[B:]
    ana_bfr[-B:] = x[block_start : block_start + B]
    windowed = padded_h[::-1] * ana_bfr
    bar_x = windowed.reshape(-1, T).sum(axis=0)
    subbands = np.fft.rfft(np.roll(bar_x, 1 - shift_i))
    subbands = 1.0 * subbands + 0.0 # subband processing here
    syn_vector = np.roll(np.fft.irfft(subbands, T), -shift_j)
    syn_bfr += padded_g * np.tile(syn_vector, len(padded_g) // T)
    y[block_start : block_start + B] = syn_bfr[:B]
    syn_bfr[:-B] = syn_bfr[B:]
    syn_bfr[-B:] = 0.0

fft_size = 32768
omega = np.linspace(0.0, np.pi, fft_size // 2, endpoint=False)
h_db = 20 * np.log10(np.abs(np.fft.rfft(h, fft_size)[:fft_size // 2]) + 1e-9)
g_db = 20 * np.log10(np.abs(np.fft.rfft(g, fft_size)[:fft_size // 2]) + 1e-9)

fig, axes = plt.subplots(3, 1, figsize=(10, 9), constrained_layout=True)
axes[0].plot(h, label="ana filter")
axes[0].plot(g, label="syn filter")
axes[0].set(xlabel="time", ylabel="impulse response")
axes[0].legend()

axes[1].plot(omega, h_db, label="ana filter")
axes[1].plot(omega, g_db, label="syn filter")
axes[1].set(xlabel="rad/sample", ylabel="magnitude (dB)")
axes[1].set_xlim([0, np.pi])
axes[1].legend()

axes[2].stem(x, linefmt="C0", label="original")
axes[2].stem(y, linefmt="C1", label=f"reconstructed, {tau0 - B + 1}-samples lag")
axes[2].set(xlabel="time", ylabel="signal")
axes[2].legend()
plt.show() 
