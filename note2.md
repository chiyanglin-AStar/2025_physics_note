##  ICTP note

- [Official repo for the ICTP/School on Parallel Programming and Parallel Architecture for High Performance Computing](https://github.com/Sera91/SMR3935-2024/tree/main)

- [Materials  for MaX  QE-2021 online school](https://gitlab.com/QEF/materials-for-max-qe2021-online-school.git)

- [Materials-for-MaX-CoE-ENCCS-workshop-2022](https://gitlab.com/QEF/materials-for-max-coe-enccs-workshop-2022.git)

- [Material for Ljubljana QE summer school 2019](https://gitlab.com/QEF/material_ncc_czechia_qe_school)

- [Material for the MaX NCC Czechia QE-2024 online school](https://gitlab.com/QEF/material_ncc_czechia_qe_school.git)

- []()

- []()

- []()

## ASE & Quantum Espresso installation
- [Quantum Espresso installation](https://pranabdas.github.io/espresso/setup/install)

- [ase install](https://wiki.fysik.dtu.dk/ase/install.html)

- [ase install 2](https://pypi.org/project/ase/)

- [ase Tutorials](https://wiki.fysik.dtu.dk/ase/tutorials/tutorials.html)

- [Hands-on • Quantum Espresso Tutorial](https://pranabdas.github.io/espresso/category/hands-on/)

- [Quantum Espresso Tutorial ref.](https://pranabdas.github.io/espresso/category/hands-on/)

- [ASE sources gitlab](https://gitlab.com/ase/ase.git)

- [ASE Calculator](https://wiki.fysik.dtu.dk/ase/ase/calculators/calculators.html)

### ASE examples:
```shell
pip install --upgrade ase

pip install --upgrade ase[test]

pip install ase quantum-espresso
```
```shell
sudo apt update && sudo apt upgrade

sudo apt install --no-install-recommends   libfftw3-dev quantum-espresso
```
### ASE test

```shell
ase test  # takes 1 min.
```
```python
from ase import Atoms
from ase.calculators.emt import EMT

atom = Atoms('N')
atom.calc = EMT()
e_atom = atom.get_potential_energy()

d = 1.1
molecule = Atoms('2N', [(0.0, 0.0, 0.0), (0.0, 0.0, d)])
molecule.calc = EMT()
e_molecule = molecule.get_potential_energy()

e_atomization = e_molecule - 2 * e_atom

print('Nitrogen atom energy: %5.2f eV' % e_atom)
print('Nitrogen molecule energy: %5.2f eV' % e_molecule)
print('Atomization energy: %5.2f eV' % -e_atomization)
```

```python
from ase import Atoms
from ase.optimize import BFGS
from ase.calculators.nwchem import NWChem
from ase.io import write
h2 = Atoms('H2',positions=[[0, 0, 0],[0, 0, 0.7]])
h2.calc = NWChem(xc='PBE')
opt = BFGS(h2, trajectory='h2.traj')
opt.run(fmax=0.02)
write('H2.xyz', h2)
h2.get_potential_energy()  # ASE's units are eV and Ang
```

```python
from ase.build import bulk
from ase.visualize import view

# 建立冰 Ih 的六方晶系結構 (簡化模型)
# 實際研究中會使用更大的超胞
a = 4.5  # Å
c = 7.3  # Å
ice = bulk('H2O', 'ice', a=a, c=c)
