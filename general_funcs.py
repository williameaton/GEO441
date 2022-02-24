import matplotlib.pyplot as plt
import numpy as np
import matplotlib.animation as animation
import matplotlib.pyplot as plt


def calc_fft_grad(f, dx):
    # Take gradient and calculate freq/wavenumbers (fk):
    fk = np.fft.fftfreq(len(f), dx)
    F = np.fft.fft(f)

    grad = 2*np.pi*np.fft.ifft(1j * fk * F)

    return grad

def normalise(x):
    return x/np.amax(np.abs(x))


def homo_ofsize(val, ofsize):
    return ofsize*0 + val

def hetero_ofsize(bounds, vals, x):
    assert(len(bounds)==len(vals)+1)

    bounds = np.array(bounds)
    vals  = np.array(vals)

    out = np.zeros(len(x))
    for i in range(len(vals)):
        mask = [x>=bounds[i]]
        out[mask] = vals[i]
    return out



def save_anim_mp4(anim, fname, fps):
    # Animate and save:
    writervideo = animation.FFMpegWriter(fps=fps)
    print("Saving video:")
    anim.save(fname, writer=writervideo)
    print(f"Written to {fname}")