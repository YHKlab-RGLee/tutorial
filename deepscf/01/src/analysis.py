"""Small plotting helpers for the periodic-polymer tutorial."""
from pathlib import Path
import re
import numpy as np
import matplotlib.pyplot as plt
import sisl


def density(path):
    # rhoatom and *.RHO share SIESTA's binary grid format.
    reader = type(sisl.get_sile('dummy.RHO'))
    grid = reader(str(path)).read_grid()
    assert np.isfinite(grid.grid).all()
    return grid  # sisl returns electron density in e / Ang^3.


def plot_density(reference, atomic=None, prediction=None):
    ref = density(reference)
    grids = [('SCF', ref)]
    if atomic is not None:
        grids.insert(0, ('Atomic', density(atomic)))
    if prediction is not None:
        grids.insert(-1, ('DeepSCF', density(prediction)))
    for _, grid in grids:
        assert grid.shape == ref.shape
        assert np.allclose(grid.cell, ref.cell)
    cell = np.diag(ref.cell)
    assert np.allclose(ref.cell, np.diag(cell))
    dy = cell[1] / ref.shape[1]
    projected = [(label, g.grid.sum(axis=1) * dy) for label, g in grids]
    if prediction is not None:
        diff = density(prediction).grid - ref.grid
        projected.append(('DeepSCF - SCF', diff.sum(axis=1) * dy))
    fig, axes = plt.subplots(1, len(projected), figsize=(4*len(projected), 4), constrained_layout=True, squeeze=False)
    vmax = max(values.max() for label, values in projected if label != 'DeepSCF - SCF')
    for ax, (label, values) in zip(axes[0], projected):
        delta = label == 'DeepSCF - SCF'
        limit = max(abs(values.min()), abs(values.max()), 1e-12) if delta else vmax
        im = ax.imshow(values.T, origin='lower', extent=(0, cell[0], 0, cell[2]),
                       cmap='RdBu_r' if delta else 'viridis', vmin=-limit if delta else 0, vmax=limit,
                       aspect='equal')
        ax.set(xlabel='x (Angstrom)', ylabel='z (Angstrom)', title=label, xlim=(7,13))
        fig.colorbar(im, ax=ax, label=r'$\int\rho\,dy$ (e / $\mathrm{\AA}^2$)', shrink=0.8)
    electrons = ref.grid.sum() * ref.dvolume
    print('Integrated SCF charge: %.6f electrons' % electrons)
    for label, grid in grids:
        if label != 'SCF':
            error = 100 * np.abs(grid.grid-ref.grid).sum() / ref.grid.sum()
            print('%s density error (full 3D grid): %.6f %%' % (label,error))
    plt.show()
    return fig


def read_bands(directory):
    reader = sisl.get_sile(str(Path(directory) / 'Molecular.bands'))
    ticks, k, values = reader.read_data()
    assert values.shape[1] == 1, 'This example is non-spin-polarized.'
    return np.asarray(k), values[:,0,:] + reader.read_fermi_level()


def plot_bands(scf, non_scf=None):
    k, reference = read_bands(scf)
    # C8H16: 48 valence electrons, 24 doubly occupied bands.
    nocc = 24
    zero = reference[:, nocc-1].max()
    fig, ax = plt.subplots(figsize=(7,6), constrained_layout=True)
    cases = [('SCF', k, reference, 'black', '-')]
    if non_scf is not None:
        km, ml = read_bands(non_scf)
        assert ml.shape == reference.shape and np.allclose(km,k)
        cases.append(('DeepSCF / non-SCF', km, ml, 'tab:red', '--'))
    for label, axis, energy, color, style in cases:
        # Display 12 occupied and 12 unoccupied bands around the gap.
        ax.plot(axis, energy[:, nocc-12:nocc+12] - zero, color=color, ls=style, lw=1.2)
        ax.plot([], [], color=color, ls=style, label=label)
        gap = energy[:,nocc].min() - energy[:,nocc-1].max()
        print('%s gap along sampled Gamma-X path: %.6f eV' % (label,gap))
    ax.set(xticks=[k[0],k[-1]], xticklabels=[r'$\Gamma$', 'X'],
           xlim=(k[0],k[-1]), ylabel=r'$E-E_{\mathrm{VBM}}^{\mathrm{SCF}}$ (eV)')
    ax.axhline(0,color='0.6',lw=0.5)
    ax.legend()
    plt.show()
    return fig


def compare_energy(scf, non_scf):
    # Use the same final Kohn-Sham total-energy field in both logs.
    def energy(directory):
        text = (Path(directory)/'stdout.txt').read_text(errors='replace')
        values = re.findall(r'^siesta:\s+Total\s*=\s*([-+0-9.EeDd]+)',text,re.M)
        if not values:
            raise ValueError('Final siesta: Total energy not found in '+str(directory))
        return float(values[-1].replace('D','E').replace('d','e'))
    a,b=energy(scf),energy(non_scf)
    print('SCF total energy:      %.8f eV' % a)
    print('non-SCF total energy:  %.8f eV' % b)
    print('|Delta E| / atom:      %.6f meV/atom' % (abs(b-a)/24*1000))
