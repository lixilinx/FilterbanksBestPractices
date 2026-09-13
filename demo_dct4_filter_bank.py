"""Random irregular DCT4 modulated filterbank example.  
"""
import numpy as np
from filterbank_design import design_filter_bank, DesignResult, FilterBank
from matplotlib import pyplot as plt 
from scipy.fft import dct, idct 

rng = np.random.default_rng()
T = 24 # must have T % 4 = 0
B = T//4 # typically T//4 % B = 0 
Lh = int(rng.integers(T, 2*T))
Lg = int(rng.integers(T, 2*T))
tau0 = int(rng.integers(T - 1, min(Lh, Lg)))
momentum = int(rng.integers(0, 3))
eta = 1e5
lambda_ = 1.0
band_count = T // 4

gamma_quadrant = np.eye(T // 2) - np.fliplr(np.eye(T // 2))
Gamma = np.block([[gamma_quadrant, -gamma_quadrant], [-gamma_quadrant, gamma_quadrant]])

best_result: DesignResult | None = None
for trial in range(1, 101):
    h, g = rng.uniform(size=Lh), rng.uniform(size=Lg) # rand initialization is better than fbd_random_initial_guess.m
    shift_i = int(rng.integers(0, T)) # also search for (i, j)-shift pair
    shift_j = (-tau0 - shift_i) % T
    result = design_filter_bank(FilterBank(Gamma=Gamma, T=T, B=B, tau0=tau0, i=shift_i, j=shift_j, h=h, g=g, momentum=momentum), 
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
    shifted_bar_x = np.roll(bar_x, 1 - shift_i)
    folded = (shifted_bar_x[:band_count] - shifted_bar_x[band_count:2*band_count][::-1]
              -shifted_bar_x[2*band_count:3*band_count] + shifted_bar_x[3*band_count:][::-1])
    subbands = dct(folded, type=4)
    subbands = 1.0 * subbands + 0.0 # subband processing here
    inverse_dct = idct(subbands, type=4)
    syn_vector = np.roll(np.concatenate([inverse_dct, -inverse_dct[::-1], -inverse_dct, inverse_dct[::-1]]), -shift_j)
    syn_bfr += padded_g * np.tile(syn_vector, len(padded_g) // T)
    y[block_start : block_start + B] = syn_bfr[:B]
    syn_bfr[:-B] = syn_bfr[B:]
    syn_bfr[-B:] = 0.0

fft_size = 32768
omega = np.linspace(0.0, np.pi, fft_size // 2, endpoint=False)
fig, axes = plt.subplots(2, 1, figsize=(10, 9), constrained_layout=True)
time_indices = np.arange(Lh)
for band in range(band_count):
    modulated_h = h * np.cos(np.pi / band_count * ((-time_indices - shift_i) + 0.5) * (band + 0.5))
    h_db = 20 * np.log10(np.abs(np.fft.rfft(modulated_h, fft_size)[:fft_size // 2]) + 1e-9)
    axes[0].plot(omega, h_db, color="black", linestyle="-" if band%2==0 else "--")
axes[0].set(xlabel="rad/sample", ylabel="Ana-filters magnitude (dB)")

axes[1].stem(x, linefmt="C0", label="original")
axes[1].stem(y, linefmt="C1", label=f"reconstructed, {tau0 - B + 1}-samples lag")
axes[1].set(xlabel="time", ylabel="signal")
axes[1].legend()
plt.show() 
