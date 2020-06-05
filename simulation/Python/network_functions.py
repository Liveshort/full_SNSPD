import os
import shutil

def print_equations(equations, lhsName="Am", rhsName="b"):
    print_lhs_equations(equations, lhsName)
    print_rhs_equations(equations, rhsName)

def string_lhs_equations(equations, name="Am"):
    str = ""
    for i, eq in enumerate(equations):
        str += eq.string_lhs(name, i)
        str += "\n"

    return str

def string_rhs_equations(equations, name="b"):
    str = ""
    for i, eq in enumerate(equations):
        str += eq.string_rhs(name, i)

    return str

def print_lhs_equations(equations, name="Am"):
    print(string_lhs_equations(equations, name), end="")

def print_rhs_equations(equations, name="b"):
    print(string_lhs_equations(equations, name), end="")

def write_equations(equations, lhsName="Am", rhsName="b", fileName="../C/electrical.c"):
    # make a backup of the current electrical.c file to electrical.c.BU
    try:
        os.remove("{}.BU".format(fileName))
    except FileNotFoundError as e:
        print("error {}, nothing to see here, moving on...".format(e))
    shutil.copyfile(fileName, "{}.BU".format(fileName))

    # copy the template electrical.c.TEMPLATE file to electrical.c
    try:
        os.remove(fileName)
    except FileNotFoundError as e:
        print("error {}, nothing to see here, moving on...".format(e))
    shutil.copyfile("{}.TEMPLATE".format(fileName), fileName)

    # search for the correct lines to replace, and add the equations to the correct lines
    lines = []
    with open(fileName, "r") as f:
        lines = f.readlines()
        # print(lines)

        sizeTemplate = "#### TEMPLATE LINE #### PARAMETER m ####"
        sizeIdx = [i for i, l in enumerate(lines) if sizeTemplate in l]
        lines.pop(sizeIdx[0])
        lines.insert(sizeIdx[0], "    m = {};\n".format(len(equations)))

        matrixTemplate = "#### TEMPLATE LINE #### PARAMETER Am ####"
        matrixIdx = [i for i, l in enumerate(lines) if matrixTemplate in l]
        lines.pop(matrixIdx[0])
        for lhsLine in reversed(string_lhs_equations(equations, lhsName).splitlines()):
            lines.insert(matrixIdx[0], "    {}\n".format(lhsLine))

        vectorTemplate = "#### TEMPLATE LINE #### PARAMETER b ####"
        vectorIdx = [i for i, l in enumerate(lines) if vectorTemplate in l]
        lines.pop(vectorIdx[0])
        for rhsLine in reversed(string_rhs_equations(equations, rhsName).splitlines()):
            lines.insert(vectorIdx[0], "    {}\n".format(rhsLine))

        # print(lines)

    # close the file in reading mode, reopen in writing mode, write the lines
    with open(fileName, "w") as f:
        f.writelines(lines)

    print("Write successful!")
