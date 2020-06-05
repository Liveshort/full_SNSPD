from network_elements import Nanowire, Capacitor, Resistor, Inductor, CircuitPiece, Bias, KVL, KCL, CapEq
from network_functions import string_lhs_equations, string_rhs_equations, print_lhs_equations, print_rhs_equations, print_equations, write_equations

numberOfBranches = 2
numberOfCaps = 1
numberOfBias = 1
numberOfNanowires = 1

# Bias (name, num = None)
# Nanowire (indName = "XW", resName = "R_w", indNum = None, resNum = None)
# Inductor (name, num = None)
# Resistor (name = "R", group = None, num = None, timeDependent = False)
# Capacitor (name, num, replacementResistance, equation, currents, directions)

p = []
p.append(CircuitPiece([0], ["pro"], [Nanowire(indNum=0, resNum=0), Resistor(group=0, num=0)]))
p.append(Capacitor(num=0, equation=numberOfBranches+0, currents=[1], directions=["pro"]))
p.append(CircuitPiece([1], ["pro"], [Resistor(group=1, num=0)]))

# CircuitPiece (currents, directions, elements)
# KVL (CircuitPieces, directions, numberOfBranches)
# KCL (currents, directions, biasCurrents)
# CapEq (capacitor, numberOfBranches)

equations = []
equations.append(KVL(p[0:3], ["contra", "pro", "pro"], numberOfBranches))
equations.append(KCL([0, 1], ["pro", "pro"], [Bias("Iv_b", 0)]))
equations.append(CapEq(p[1], numberOfBranches))

write_equations(equations, fileName="../C/electrical.c")
