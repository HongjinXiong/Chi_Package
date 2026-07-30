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

2. sympy

3. itertools

Install sympy via pip if needed: **pip install sympy**

## Usage

### Example

For a 3m crystal, the symmetry elements are a 3‑fold rotation about z and three mirror planes. The code below computes the reduced $\chi$ tensor for second‑order input (order_input=2) and first‑order output (order_output=1).

```
order_input=2
order_output=1
# C3 rotation axis z, angle 2π/3
axis1 = sp.Matrix([0, 0, 1])
theta1 = 2 * sp.pi / 3
r1 = H_Rotation(order_input, order_output, axis1, theta1, Rotation_Inversion=False)

# Mirror planes (normals in the xy‑plane)
axis2 = sp.Matrix([1, sp.sqrt(3), 0])
axis3 = sp.Matrix([1, -sp.sqrt(3), 0])
axis4 = sp.Matrix([1, 0, 0])

m1 = H_Mirror(order_input, order_output, axis2)
m2 = H_Mirror(order_input, order_output, axis3)
m3 = H_Mirror(order_input, order_output, axis4)

# Compute the reduced tensor (without Kleinman)
Calculation_Chi(order_input, order_output, r1, m1, m2, m3, Kleinman_Symmetry=False)
```

## Output
All output is \Latex format.

1. The reduced susceptibility matrix with free parameters.

2. The solution for each independent component, expressed in terms of free variables.

3. Operation matrix for each symmetric.
