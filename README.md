# Chi_Package
This Python script uses "**sympy**" to calculate the linear and higher-order nonlinear susceptibility tensor ($\chi$) for a given material symmetry. It constructs symmetry operations (rotation, inversion, mirror) as Kronecker‑product matrices acting on the input and output fields, sets up the linear system of equations that enforce invariance under these operations, and solves for the tensor components. 

## Features
Build symmetry transformation matrices for:

1. Rotation about an arbitrary axis (with optional inversion).

2. Inversion (parity).

3. Mirror reflection across a plane with arbitrary normal.

Handle arbitrary input order (e.g., second‑order, third‑order) and output order (first‑order, second‑order, etc.). Automatically construct the susceptibility matrix with symbolic components. Impose symmetry constraints by requiring $\chi = H_{out} · \chi · H_{in}^T$ for each symmetry operation. Solve the resulting linear system and display the simplified tensor. Optionally enforce Kleinman symmetry, which makes the tensor fully symmetric under permutation of all indices.

## Dependencies
1. Python 3.x

2. sympy (tested with version 1.10+)

3. itertools (standard library)

Install sympy via pip if needed: **pip install sympy**
