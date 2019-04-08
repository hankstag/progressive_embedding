from sympy import *
from sympy.matrices import *
import os
import re
import argparse

# local
import pretty_print

if __name__ == "__main__":

    dims = [2, 3]

    cpp = ""
    hpp = "#pragma once\n\n#include <Eigen/Dense>\n\n"
    cpp = cpp + "namespace autogen{\n"
    hpp = hpp + "namespace autogen{\n"

    lambdaa = Symbol('lambda', real=True)

    V = zeros(3, 2)
    for i in range(0, 3):
        for j in range(0, 2):
            V[i, j] = Symbol('V[' + str(i) + ',' + str(j) + ']', real=True)

    G = zeros(2, 3)
    for i in range(0, 2):
        for j in range(0, 3):
            G[i, j] = Symbol('G[' + str(i) + ',' + str(j) + ']', real=True)

    J = G*V
    Ji = J**-1
    J_F = J[0, 0]*J[0, 0] + J[0, 1]*J[0, 1] +J[1, 0]*J[1, 0] +J[1, 1]*J[1, 1]
    Ji_F = Ji[0, 0]*Ji[0, 0] + Ji[0, 1]*Ji[0, 1] +Ji[1, 0]*Ji[1, 0] +Ji[1, 1]*Ji[1, 1];

    # Energy
    E = J_F + Ji_F

    # Gradient
    grad = zeros(1, 2)
    for j in range(0, 2):
        grad[0,j] = diff(E,V[0,j])

    # Hessian
    hess = zeros(2, 2)
    for i in range(0, 2):
        for j in range(0, 2):
            hess[i,j] = diff(grad[0,i],V[0,j])

    # lambdas = simplify(lambdas)
    # E = simplify(E)

    # Energy
    c99 = pretty_print.C99_print(E)
    c99 = re.sub(r"V\[(\d{1}),(\d{1})\]", r'V(\1,\2)', c99)
    c99 = re.sub(r"G\[(\d{1}),(\d{1})\]", r'G(\1,\2)', c99)
    signature = "template<typename Scalar>\nScalar sd_energy(const Eigen::Matrix<Scalar,3,2>& V, const Eigen::Matrix<Scalar,2,3>& G) "
    hpp = hpp + signature + " {\n" + "Scalar result_0;\n" + c99  + "\nreturn result_0;\n}\n"

    # Gradient
    c99 = pretty_print.C99_print(grad)
    c99 = re.sub(r"V\[(\d{1}),(\d{1})\]", r'V(\1,\2)', c99)
    c99 = re.sub(r"G\[(\d{1}),(\d{1})\]", r'G(\1,\2)', c99)
    c99 = re.sub(r"result_0\[(\d{1})\]", r'res(\1)', c99)
    signature = "template<typename Scalar>\nvoid sd_grad(const Eigen::Matrix<Scalar,3,2>& V, const Eigen::Matrix<Scalar,2,3>& G, Eigen::Matrix<Scalar,1,2>& res) "
    hpp = hpp + signature + " {\n" + c99 + "\n}\n"

    # Hessian
    c99 = pretty_print.C99_print(hess)
    c99 = re.sub(r"V\[(\d{1}),(\d{1})\]", r'V(\1,\2)', c99)
    c99 = re.sub(r"G\[(\d{1}),(\d{1})\]", r'G(\1,\2)', c99)
    c99 = re.sub(r"result_0\[(\d{1})\]", r'res(\1)', c99)
    signature = "template<typename Scalar>\nvoid sd_hess(const Eigen::Matrix<Scalar,3,2>& V, const Eigen::Matrix<Scalar,2,3>& G, Eigen::Matrix<Scalar,2,2>& res) "
    hpp = hpp + signature + " {\n" + c99 + "\n}\n"

    # Hessian

    cpp = cpp + "\n}\n"
    hpp = hpp + "\n}\n"

    path = os.getcwd()

    print("saving...")
    with open(os.path.join(path, "auto_grad.cpp"), "w") as file:
        file.write(cpp)

    with open(os.path.join(path, "auto_grad.hpp"), "w") as file:
        file.write(hpp)

    print("done!")
