from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# 1. 读取与代码放在同一个文件夹中的数据和参数
folder = Path(r"C:\Users\Administrator\Qxu6028_Lab")
files = [folder / name for name in
         ['L1_data.csv', 'L1_data.xls', 'L1_data.xlsx', 'L1_data']
         if (folder / name).exists()]
if len(files) != 1:
    raise ValueError("未找到数据文件，或者找到了多个数据文件，请检查文件夹。")
data_file = files[0]
# 按实际内容识别格式，也能读取扩展名为 .xls 的 CSV 文本。
header = data_file.read_bytes()[:8]
if header.startswith(b'PK'):
    data = pd.read_excel(data_file, engine='openpyxl')
elif header.startswith(bytes.fromhex('D0CF11E0A1B11AE1')):
    data = pd.read_excel(data_file, engine='xlrd')
else:
    data = pd.read_csv(data_file, encoding='utf-8-sig')
data.columns = data.columns.str.strip()
with open(folder / 'L1_params.json', encoding='utf-8-sig') as f:
    params = json.load(f)

EF = data['EF_eV'].to_numpy(dtype=float)
n_band = data['n_band_cm2'].to_numpy(dtype=float)
rho = data['rho_ohm_sq'].to_numpy(dtype=float)
T, c = params['T'], params['c']
EF_op, rho_op = params['EF_op_eV'], params['rho_ohm_sq']
kB, NA, e = 8.617333e-5, 3.816e15, 1.602176634e-19

# 2. 积分计算载流子浓度、修正系数和迁移率
# lower=0 为带边计数；lower=EF 为费米参考计数。
def count(ef, lower):
    energy = np.linspace(lower, 2.5, 20001)
    fermi = 1 / (np.exp(np.clip((energy - ef) / (kB * T), -60, 60)) + 1)
    y = c * energy * fermi
    integral = np.trapezoid(y, energy) if hasattr(np, 'trapezoid') else np.trapz(y, energy)
    return float(NA * integral)

n_fermi = np.array([count(ef, ef) for ef in EF])
C = n_band / n_fermi
mobility = 1 / (e * n_fermi * rho)
nb_op = count(EF_op, 0)
nf_op = count(EF_op, EF_op)
C_op = nb_op / nf_op
mu_op = 1 / (e * nf_op * rho_op)

# 在扫描范围内寻找 C=2 和 C=4 的交点。
def crossing(target):
    low, high = float(EF.min()), float(EF.max())
    def ratio(ef):
        return count(ef, 0) / count(ef, ef)
    if not ratio(low) <= target <= ratio(high):
        return None
    for _ in range(50):
        mid = (low + high) / 2
        if ratio(mid) < target:
            low = mid
        else:
            high = mid
    return (low + high) / 2

# 3. 画四联图，星号表示精确工作点
fig, axes = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)
curves = [n_band / 1e12, n_fermi / 1e12, C, mobility / 1e4]
points = [nb_op / 1e12, nf_op / 1e12, C_op, mu_op / 1e4]
titles = ['Band-edge carrier density', 'Fermi-referenced carrier density',
          'Correction factor', 'Corrected mobility']
labels = [r'$n_{band}$ ($10^{12}$ cm$^{-2}$)',
          r'$n_F$ ($10^{12}$ cm$^{-2}$)', r'$C = n_{band}/n_F$',
          r'$\mu$ ($10^4$ cm$^2$ V$^{-1}$ s$^{-1}$)']
colors = ['#236c9e', '#238579', '#7557ab', '#bb7740']
for ax, y, point, title, label, color in zip(axes.flat, curves, points, titles, labels, colors):
    ax.plot(EF, y, 'o-', color=color, markerfacecolor='white')
    ax.axvline(EF_op, color='gray', linestyle=':', linewidth=1)
    ax.plot(EF_op, point, '*', color='black', markersize=12)
    ax.set(title=title, xlabel=r'$E_F$ (eV)', ylabel=label)
    ax.text(0.04, 0.92, f'Operating point: {point:.6f}', transform=ax.transAxes,
            bbox=dict(facecolor='white', edgecolor='none', alpha=0.85))
    ax.grid(alpha=0.25)

for target in [2, 4]:
    x = crossing(target)
    if x is not None:
        ax = axes[1, 0]
        ax.axhline(target, color='gray', linestyle='--', linewidth=0.8)
        ax.plot(x, target, 'o', color=colors[2])
        ax.annotate(f'C={target}: {x:.5f} eV', (x, target),
                    xytext=(8, -24 if target == 2 else 24), textcoords='offset points',
                    fontsize=9, bbox=dict(facecolor='white', edgecolor='none', alpha=0.85),
                    arrowprops=dict(arrowstyle='-', color='gray'))
fig.suptitle(f"{params.get('group', 'L1')} | T={T:g} K | "
             f'EF={EF_op:.8f} eV | C={C_op:.6f}', fontsize=14)
fig.savefig(folder / 'L1 FIGURE.png', dpi=300)
plt.close(fig)

# 4. 只更新已有 results.json 中的三个 L1 数值，保留其他字段
with open(folder / 'results.json', encoding='utf-8-sig') as f:
    results = json.load(f)
results.setdefault('L1', {}).update({
    'correction_factor': C_op,
    'n_fermi_op_cm2': nf_op,
    'mobility_corrected_cm2_Vs': mu_op
})
with open(folder / 'results.json', 'w', encoding='utf-8') as f:
    json.dump(results, f, indent=2, ensure_ascii=False, allow_nan=False)
    f.write('\n')
print(f'读取：{data_file.name}、L1_params.json')
print('已保存：L1 FIGURE.png、results.json')
print(json.dumps(results['L1'], indent=2))
