"""L2 计算和五张图：python L2.py
同目录放 L2_data.csv、L2_meta.json、已有的 results.json。
依赖：numpy matplotlib。不使用 Git，不生成其他辅助文件。
"""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

folder = Path(__file__).resolve().parent
meta = json.loads((folder/'L2_meta.json').read_text(encoding='utf-8-sig'))
raw = np.genfromtxt(folder/'L2_data.csv', delimiter=',', names=True)
T = float(meta['T_K'])
kT = 8.617333e-5*T
indices = np.unique(raw['measurement_index']).astype(int)
repeats = []
for i in indices:
    rows = raw[raw['measurement_index'] == i]
    rows = rows[np.argsort(rows['repeat'])]
    assert len(rows) == 3 and len(np.unique(rows['repeat'])) == 3
    repeats.append(np.column_stack([rows['R_H_raw'], rows['R_xx_raw']]))
repeats = np.array(repeats)
assert np.all(np.isfinite(repeats)) and np.all(repeats[:, :, 1] > 0)
means = repeats.mean(axis=1)

# 1. 每个通道分别归一化，保留 Hall 正负号。
def normalise(pairs):
    return pairs / np.array([np.max(np.abs(pairs[:, 0])), np.max(pairs[:, 1])])

# 2. 按指南的固定网格建立 cone 和 hard-gap 模型库。
# e、mu、NA 和 c 的公共正比例常数在各自归一化后消去，故此处设为 1。
energy = np.linspace(-2.5, 2.5, 4000)
ef_grid = np.linspace(-8*kT, 9*kT, 400)
gaps = np.r_[0.0, np.arange(0.5, 8.001, 0.25)]
f = 1/(np.exp(np.clip((energy[None, :]-ef_grid[:, None])/kT, -60, 60))+1)
# 梯形权重与 np.trapezoid(np.where(...), energy) 相同。
weights = np.full(4000, energy[1]-energy[0])
weights[[0, -1]] *= 0.5
n_kernel = np.where(energy[None, :] >= ef_grid[:, None], f, 0)*weights
p_kernel = np.where(energy[None, :] <= ef_grid[:, None], 1-f, 0)*weights
library = []
for gap in gaps:
    dos = np.maximum(np.abs(energy)-gap*kT/2, 0)
    n, p = n_kernel@dos, p_kernel@dos
    rh = (p-n)/(p+n)**2
    rs = 1/(p+n)
    library.append(normalise(np.column_stack([rh, rs])))
library = np.array(library)

# 3. 在二维归一化空间中找最近点，再用 15% 规则选择模型。
def fit(pairs):
    xy = normalise(pairs)
    d2 = ((library[:, None, :, :]-xy[None, :, None, :])**2).sum(axis=-1)
    nearest = d2.argmin(axis=2)
    rms = np.sqrt(d2.min(axis=2).mean(axis=1))
    best_gap = 1+np.argmin(rms[1:])
    chosen = best_gap if rms[best_gap] < 0.85*rms[0] and gaps[best_gap] >= 1 else 0
    return chosen, ef_grid[nearest[chosen]], rms, best_gap, nearest[chosen]

chosen, trajectory, rms, best_gap, nearest = fit(means)
model = 'hard_gap' if chosen else 'cone'

# 4. 200 次重抽样：每个测量位置随机取一整对重复数据，重新归一化及拟合。
rng = np.random.default_rng(20260923)
boot_gaps, boot_ef = [], []
for _ in range(200):
    pairs = repeats[np.arange(len(indices)), rng.integers(0, 3, len(indices))]
    b, ef, _, _, _ = fit(pairs)
    boot_gaps.append(gaps[b])
    boot_ef.append(ef)
lo, median, hi = np.percentile(boot_gaps, [16, 50, 84])
ef_lo, ef_hi = np.percentile(boot_ef, [16, 84], axis=0)
fraction = float(np.mean(np.array(boot_gaps) > 0))
assert len(trajectory) == len(indices) and np.all(np.isfinite(trajectory))

# 5. 更新三个结果字段和完整 EF 轨迹；保留 L1、L3 及其他已有字段。
path = folder/'results.json'
results = json.loads(path.read_text(encoding='utf-8-sig'))
results['L2'] = dict(model=model, Eg_kT=float(gaps[chosen]),
                    Eg_uncertainty_kT=float((hi-lo)/2),
                    EF_trajectory_eV=trajectory.tolist())
path.write_text(json.dumps(results, indent=2, ensure_ascii=False, allow_nan=False)+'\n', encoding='utf-8')

# 6. 画图：每张直接保存在代码所在文件夹。
plt.rcParams.update({'font.size': 11, 'axes.spines.top': False, 'axes.spines.right': False})
def figure(title, xlabel, ylabel):
    fig, ax = plt.subplots(figsize=(9, 6), layout='constrained')
    ax.set(title=title, xlabel=xlabel, ylabel=ylabel)
    ax.grid(alpha=.2)
    return fig, ax

def save(fig, ax, name):
    ax.legend(fontsize=9)
    fig.savefig(folder/name, dpi=250)
    plt.close(fig)

# 图 1：最优带隙拟合，和无带隙基准比较。
fig, ax = figure(f'G11 | Selected {model}: Eg/kBT = {gaps[chosen]:g}',
                 'Normalised Hall coefficient', 'Normalised sheet resistance')
ax.plot(*library[0].T, '--', color='gray', label=f'Cone: RMS={rms[0]:.6f}')
ax.plot(*library[best_gap].T, color='#167b8b', label=f'Best gap: RMS={rms[best_gap]:.6f}')
xy = normalise(means)
points = ax.scatter(*xy.T, c=indices, cmap='plasma', s=25, label='G11 repeat means')
fig.colorbar(points, ax=ax, label='Measurement index')
save(fig, ax, '最优带隙.png')

# 图 2：带隙搜索；2 kBT 只作示例候选，不是指南规定的初始值。
fig, ax = figure('G11 | Band-gap search', 'Eg / kBT', 'Nearest-curve RMS')
ax.plot(gaps[1:], rms[1:], 'o-', ms=4, color='#167b8b', label='Hard-gap candidates')
ax.axhline(rms[0], color='gray', ls='--', label='Cone baseline')
ax.axhline(.85*rms[0], color='gray', ls=':', label='0.85 x cone RMS')
example = int(np.argmin(abs(gaps-2)))
for j, label in [(example, 'Example'), (best_gap, 'Best')]:
    ax.plot(gaps[j], rms[j], '*', ms=12)
    ax.annotate(f'{label}: {gaps[j]:g} kBT\nRMS={rms[j]:.6f}', (gaps[j], rms[j]),
                xytext=(10, 30), textcoords='offset points', fontsize=9,
                arrowprops=dict(arrowstyle='-', color='gray'))
save(fig, ax, '带隙优化准确性图.png')

# 图 3：EF 轨迹及逐点 bootstrap 区间，单位 eV。
fig, ax = figure('G11 | Fermi-level trajectory', 'Measurement index', 'EF (eV)')
ax.fill_between(indices, ef_lo, ef_hi, color='#167b8b', alpha=.25, label='16th-84th percentile')
ax.plot(indices, trajectory, 'o-', ms=4, label='Fit to repeat means')
ax.axhline(0, color='gray', ls=':')
save(fig, ax, '费米能级变化图.png')

# 图 4、5：RH 和 Rs 对 EF；每条曲线按自身最大值归一化。
for column, name, ylabel in [(0, 'RH和EF图.png', 'RH / max|RH|'),
                              (1, 'Rs和EF图.png', 'Rs / max(Rs)')]:
    fig, ax = figure('G11 | '+('RH' if column==0 else 'Rs')+' versus EF', 'EF / kBT', ylabel)
    for j, label, color in [(0, 'Cone: Eg=0', 'gray'),
                             (example, 'Example: Eg=2 kBT', '#d18c35'),
                             (chosen, f'Selected: Eg={gaps[chosen]:g} kBT', '#167b8b')]:
        ax.plot(ef_grid/kT, library[j, :, column], color=color, label=label)
    ax.axvline(0, color='gray', ls=':')
    if column==0: ax.axhline(0, color='gray', ls=':')
    save(fig, ax, name)

print(f'Model: {model}; Eg/kBT: {gaps[chosen]:g}; Eg: {gaps[chosen]*kT:.6f} eV')
print(f'Cone RMS: {rms[0]:.6f}; best-gap RMS: {rms[best_gap]:.6f}')
print(f'Bootstrap gap interval: [{lo:g}, {hi:g}] kBT; hard-gap fraction: {fraction:.1%}')
print('结果已更新，五张图已保存。区间半宽为零不代表物理不确定度为零。')
