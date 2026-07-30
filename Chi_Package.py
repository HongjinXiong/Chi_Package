import sympy as sp
from itertools import product

def H_Inversion(order_input,order_output):

    Inversion = sp.Matrix([
        [-1, 0, 0],
        [0, -1, 0],
        [0, 0, -1]
    ])

    Inversion_input = sp.Matrix([1])
    Inversion_output = sp.Matrix([1])

    for i in range(order_input):
        Inversion_input = sp.kronecker_product(Inversion_input, Inversion)
    for i in range(order_output):
        Inversion_output = sp.kronecker_product(Inversion_output, Inversion)

    print(f'Inversion_input_matrix (latex format):\n {sp.latex(Inversion_input)}\n')
    print(f'Inversion_output_matrix (latex format):\n {sp.latex(Inversion_output)}\n')
    print()
    return Inversion_input,Inversion_output

def H_Rotation(order_input,order_output,rotate_axis,theta,Rotation_Inversion):

    s = sp.sin(theta)
    c = sp.cos(theta)

    R_z = sp.Matrix([
        [c, -s, 0],
        [s, c, 0],
        [0, 0, 1]
    ])
    if rotate_axis [0,0] !=0 or rotate_axis [1,0] !=0:
        axis_xy = sp.Matrix([rotate_axis[0, 0], rotate_axis[1, 0]])
        length_axis = (sp.transpose(rotate_axis) * rotate_axis) ** 0.5
        length_axisxy = (sp.transpose(axis_xy) * axis_xy) ** 0.5
        cos_theta_axis_to_z =sp.nsimplify( rotate_axis[2, 0] / length_axis[0, 0])
        sin_theta_axis_to_z =sp.nsimplify( sp.sqrt(rotate_axis[0, 0] ** 2 + rotate_axis[1, 0] ** 2) / length_axis[0, 0])
        cos_axisxy_to_x = sp.nsimplify(rotate_axis[0, 0] / length_axisxy[0, 0])
        sin_axisxy_to_x = sp.nsimplify(rotate_axis[1, 0] / length_axisxy[0, 0])


        z_rotate = sp.Matrix([
            [cos_axisxy_to_x, sin_axisxy_to_x, 0],
            [-sin_axisxy_to_x, cos_axisxy_to_x, 0],
            [0, 0, 1]
        ])
        y_rotate=sp.Matrix([
            [cos_theta_axis_to_z, 0, -sin_theta_axis_to_z],
            [0, 1, 0],
            [sin_theta_axis_to_z, 0, cos_theta_axis_to_z]
        ])

        O=sp.transpose(y_rotate*z_rotate)*R_z*y_rotate*z_rotate

        H_input = sp.Matrix([1])
        H_output = sp.Matrix([1])
        for i in range(order_input):
            H_input=sp.kronecker_product(H_input,O)
        for i in range(order_output):
            H_output = sp.kronecker_product(H_output, O)

        Rotation_input=H_input
        Rotation_output=H_output
        if Rotation_Inversion==True:
            Inversion=H_Inversion(order_input,order_output)
            Rotation_input = Inversion[0]*Rotation_input
            Rotation_output =Inversion[1]*Rotation_output
    else:
        H_Rz_input = sp.Matrix([1])
        H_Rz_output = sp.Matrix([1])
        for i in range(order_input):
            H_Rz_input = sp.kronecker_product(H_Rz_input, R_z)
        for i in range(order_output):
            H_Rz_output = sp.kronecker_product(H_Rz_output, R_z)
        Rotation_input = H_Rz_input
        Rotation_output =  H_Rz_output

        if Rotation_Inversion==True:
            Inversion=H_Inversion(order_input,order_output)
            Rotation_input = Inversion[0]*Rotation_input
            Rotation_output =Inversion[1]*Rotation_output

    print(f'Rotation axis (latex format): {sp.latex(rotate_axis)}\n')
    print(f'Rotation_input_matrix (latex format):\n {sp.latex(Rotation_input)}\n')
    print(f'Rotation_output_matrix (latex format):\n {sp.latex(Rotation_output)}\n')
    print()
    return Rotation_input, Rotation_output

def H_Mirror(order_input,order_output,Mirror_axis):

    Mirror_xy = sp.Matrix([
        [1, 0, 0],
        [0, 1, 0],
        [0, 0, -1]
    ])


    if Mirror_axis[0, 0] != 0 or Mirror_axis[1, 0] != 0:
        axis_xy = sp.Matrix([Mirror_axis[0, 0], Mirror_axis[1, 0]])
        length_axis = (sp.transpose(Mirror_axis) * Mirror_axis) ** 0.5
        length_axisxy = (sp.transpose(axis_xy) * axis_xy) ** 0.5
        cos_theta_axis_to_z = sp.nsimplify(Mirror_axis[2, 0] / length_axis[0, 0])
        sin_theta_axis_to_z = sp.nsimplify(
            sp.sqrt(Mirror_axis[0, 0] ** 2 + Mirror_axis[1, 0] ** 2) / length_axis[0, 0])
        cos_axisxy_to_x = sp.nsimplify(Mirror_axis[0, 0] / length_axisxy[0, 0])
        sin_axisxy_to_x = sp.nsimplify(Mirror_axis[1, 0] / length_axisxy[0, 0])

        z_rotate = sp.Matrix([
            [cos_axisxy_to_x, sin_axisxy_to_x, 0],
            [-sin_axisxy_to_x, cos_axisxy_to_x, 0],
            [0, 0, 1]
        ])
        y_rotate = sp.Matrix([
            [cos_theta_axis_to_z, 0, -sin_theta_axis_to_z],
            [0, 1, 0],
            [sin_theta_axis_to_z, 0, cos_theta_axis_to_z]
        ])

        O = sp.transpose(y_rotate * z_rotate) * Mirror_xy * y_rotate * z_rotate

        H_input = sp.Matrix([1])
        H_output = sp.Matrix([1])
        for i in range(order_input):
            H_input = sp.kronecker_product(H_input, O)
        for i in range(order_output):
            H_output = sp.kronecker_product(H_output, O)


        Mirror_input = H_input
        Mirror_output = H_output
    else:
        Mirror_input = sp.Matrix([1])
        Mirror_output = sp.Matrix([1])
        for i in range(order_input):
            Mirror_input = sp.kronecker_product(Mirror_input, Mirror_xy)
        for i in range(order_output):
            Mirror_output = sp.kronecker_product(Mirror_output, Mirror_xy)

        Mirror_input = Mirror_input
        Mirror_output =  Mirror_output

    print(f'Normal axis of mirror plane (latex format): {sp.latex(Mirror_axis)}\n')
    print(f'Mirror_input_matrix (latex format):\n {sp.latex(Mirror_input)}\n')
    print(f'Mirror_output_matrix (latex format):\n {sp.latex(Mirror_output)}\n')
    print()
    return Mirror_input,Mirror_output




def Calculation_Chi(order_input, order_output, *args,Kleinman_Symmetry):
    indices = ['x', 'y', 'z']
    input_combos = [''.join(p) for p in product(indices, repeat=order_input)]
    output_combos = [''.join(p) for p in product(indices, repeat=order_output)]

    if Kleinman_Symmetry==True:
        def elem_name(i, j):
            combined = ''.join(sorted(output_combos[i] + input_combos[j]))
            return sp.symbols(f'\chi_{{{combined}}}')
    else:
        def elem_name(i, j):
            return sp.symbols(f'\chi_{{{output_combos[i]}{input_combos[j]}}}')

    chi_matrix = sp.Matrix(len(output_combos), len(input_combos), elem_name)


    all_vars = []
    seen = set()
    for i in range(chi_matrix.rows):
        for j in range(chi_matrix.cols):
            var = chi_matrix[i, j]
            if var not in seen:
                seen.add(var)
                all_vars.append(var)
    equations=[]
    for inp, outp in args:
        eq_matrix = chi_matrix - outp * chi_matrix * sp.Transpose(inp)

        equations.append(eq_matrix)




    solution = sp.solve(equations, all_vars)
    subs_dict = {}
    for var in all_vars:
        if var in solution:
            subs_dict[var] = solution[var]
        else:
            subs_dict[var] = var

    substituted_matrix = chi_matrix.subs(subs_dict)


    print("Reduced Chi Matrix After Calculation (latex format):")
    print(sp.latex(substituted_matrix),'\n')

    print("\nSolution (latex format):")

    #all solutions

    for var in all_vars:
        if var in solution:
            expr = sp.simplify(solution[var])
            print(f"{var} = {expr}")
        else:
            print(f"{var} = {var} (free variable)")

    return solution, substituted_matrix

order_input=2
order_output=1

#Exsample: for a 3m symmetric crystal,

# its C3 rotation axis is z axis [0,0,1], thus:
axis1=sp.Matrix([0,0,1])
theta1=2*sp.pi/3
r1=H_Rotation(order_input,order_output,axis1,theta1,Rotation_Inversion=False)
#If the material has rotation inversion symmetry, the "Rotation_Inversion" is Ture

# its normal axes of mirror planes are:
axis2=sp.Matrix([1,sp.sqrt(3),0])
axis3=sp.Matrix([1,-sp.sqrt(3),0])
axis4=sp.Matrix([1,0,0])

m1=H_Mirror(order_input,order_output,axis2)
m2=H_Mirror(order_input,order_output,axis3)
m3=H_Mirror(order_input,order_output,axis4)

# To calculate Chi
Calculation_Chi(order_input,order_output,r1,m1,m2,m3,Kleinman_Symmetry=False)
#If consider the Kleinman symmetry, the "Kleinman_Symmetry" is Ture
