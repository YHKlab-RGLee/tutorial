"""Predict with the repository model and preserve the reference RHO header."""
import argparse
from pathlib import Path
import sys
import time
import h5py
import numpy as np
import torch
from omegaconf import OmegaConf

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('project', type=Path)
p.add_argument('work', type=Path)
a = p.parse_args()
project, work = a.project.resolve(), a.work.resolve()
sys.path.insert(0, str(project))
from model import ResUNet
from utils.siestaio import readGrid, writeGrid

config = OmegaConf.load(str(project / 'configs/predict/examples/r_plus_s_x3.yaml'))
scale = dict(line.split(' = ') for line in (project / 'datasets/scale/scale_r_plus_s_x3.txt').read_text().splitlines() if line.strip())
assert float(scale['mean3']) == config.model.mean
assert float(scale['std3']) == config.model.std
# CPU is the default: no CUDA dependency for this tutorial.
torch.set_num_threads(2)
model = ResUNet(config.model).cpu()
checkpoint = torch.load(str(project / config.load.model), map_location='cpu')
model.load_state_dict(checkpoint['model'] if 'model' in checkpoint else checkpoint)
model.eval()
with h5py.File(str(work / 'dataset/polymer.h5'), 'r') as f:
    feature, target = f['feature'][:], f['target'][:]
cell, mesh, reference = readGrid(str(work / 'scf/polymer/Molecular.RHO'))
assert feature.shape == (3, *mesh) and target.shape == tuple(mesh)
assert tuple(mesh) == (80, 80, 40)
assert np.all(mesh % 8 == 0)
assert np.isfinite(feature).all() and np.isfinite(target).all()
assert np.allclose(target, reference[0])
start = time.monotonic()
with torch.no_grad():
    rho = model(torch.from_numpy(feature).float().unsqueeze(0)).numpy()
assert np.isfinite(rho).all()
output = work / 'predict'
output.mkdir(exist_ok=True)
writeGrid(str(output / 'polymer.RHO'), cell, mesh, rho)
# Check the saved density, not a second inference.
c2, m2, r2 = readGrid(str(output / 'polymer.RHO'))
assert np.allclose(c2, cell) and np.array_equal(m2, mesh)
assert np.allclose(r2, rho)
initial = feature[1] * config.model.std + config.model.mean
for label, values in [('Atomic', initial), ('ML', r2[0])]:
    error = 100 * np.abs(values - target).sum() / target.sum()
    print('%s density error: %.6f %%' % (label, error))
print('Prediction time: %.1f s' % (time.monotonic()-start))
print('Saved:', output / 'polymer.RHO')
