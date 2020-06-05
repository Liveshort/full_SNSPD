class Nanowire:
    def __init__(self, indName = "XW", resName = "R_w", indNum = None, resNum = None):
        self.indName = indName
        self.resName = resName
        self.indNum = indNum
        self.resNum = resNum

    def n(self, direction):
        id = self.indName
        if self.indNum is not None:
            id += "[{}]".format(self.indNum)
        rn = self.resName
        if self.resNum is not None:
            rn += "[{}]".format(self.resNum)

        if direction == "pro":
            return "- {} - {}[n]".format(id, rn)
        else:
            return "{} + {}[n]".format(id, rn)

    def nm1(self, direction):
        id = self.indName
        if self.indNum is not None:
            id += "[{}]".format(self.indNum)
        rn = self.resName
        if self.resNum is not None:
            rn += "[{}]".format(self.resNum)

        if direction == "pro":
            return "- {} + {}[n-1]".format(id, rn)
        else:
            return "{} - {}[n-1]".format(id, rn)


class Inductor:
    def __init__(self, name, num = None):
        self.name = name
        self.num = num

    def n(self, direction):
        name = self.name
        if self.num is not None:
            name += "[{}]".format(self.num)

        if direction == "pro":
            return "- {}[n]".format(name)
        else:
            return "{}[n]".format(name)

    def nm1(self, direction):
        name = self.name
        if self.num is not None:
            name += "[{}]".format(self.num)

        if direction == "pro":
            return "- {}[n-1]".format(name)
        else:
            return "{}[n-1]".format(name)


class Resistor:
    def __init__(self, name = "R", group = None, num = None, timeDependent = False):
        self.name = name
        self.group = group
        self.num = num
        self.timeDependent = timeDependent

    def n(self, direction):
        name = self.name
        if self.group is not None:
            name += "[{}]".format(self.group)
        if self.num is not None:
            name += "[{}]".format(self.num)

        if direction == "pro":
            res = "- {}".format(name)
        else:
            res = "{}".format(name)

        if self.timeDependent:
            res += "[n]"

        return res

    def nm1(self, direction):
        name = self.name
        if self.group is not None:
            name += "[{}]".format(self.group)
        if self.num is not None:
            name += "[{}]".format(self.num)

        if direction == "pro":
            res = "{}".format(name)
        else:
            res = "- {}".format(name)

        if self.timeDependent:
            res += "[n-1]"

        return res


class Capacitor:
    def __init__(self, name = "V_c", num=None, replacementResistance = "Y", equation=None, currents=None, directions=None):
        self.name = name
        self.num = num
        self.replacementResistance = replacementResistance
        self.equation = equation
        self.currents = currents
        self.directions = directions

    def n(self, direction):
        name = self.name
        if self.num is not None:
            name += "[{}]".format(self.num)

        if direction == "pro":
            return "- {}[n]".format(name)
        else:
            return "{}[n]".format(name)

    def nm1(self, direction):
        name = self.name
        if self.num is not None:
            name += "[{}]".format(self.num)

        if direction == "pro":
            return "{}[n-1]".format(name)
        else:
            return "- {}[n-1]".format(name)

    def lhs_kvl(self):
        res = dict()
        res[self.equation] = -1

        return res

    def rhs_kvl(self):
        res = dict()
        res[self.equation] = self.nm1("pro")

        return res

    def lhs_cap(self):
        replacementResistance = self.replacementResistance
        if self.num is not None:
            replacementResistance += "[{}]".format(self.num)

        res = dict()
        for current, direction in zip(self.currents, self.directions):
            if direction == "pro":
                res[current] = "-{}".format(replacementResistance)
            else:
                res[current] = "{}".format(replacementResistance)
        res[self.equation] = 1

        return res

    def rhs_cap(self):
        replacementResistance = self.replacementResistance
        if self.num is not None:
            replacementResistance += "[{}]".format(self.num)

        res = dict()
        for current, direction in zip(self.currents, self.directions):
            if direction == "pro":
                res[current] = "{}".format(replacementResistance)
            else:
                res[current] = "-{}".format(replacementResistance)
        res[self.equation] = self.nm1("pro")

        return res


class Bias:
    def __init__(self, name, num = None):
        self.name = name
        self.num = num

    def n(self, direction):
        name = self.name
        if self.num is not None:
            name += "[{}]".format(self.num)

        if direction == "pro":
            return "- {}[n]".format(name)
        else:
            return "{}[n]".format(name)

    def nm1(self, direction):
        name = self.name
        if self.num is not None:
            name += "[{}]".format(self.num)

        if direction == "pro":
            return "{}[n-1]".format(name)
        else:
            return "- {}[n-1]".format(name)

    def static(self, direction):
        name = self.name
        if self.num is not None:
            name += "[{}]".format(self.num)

        if direction == "pro":
            return "{}".format(name)
        else:
            return "- {}".format(name)


class CircuitPiece:
    def __init__(self, currents, directions, elements):
        self.currents = currents
        self.directions = directions
        self.elements = elements

    def lhs_kvl(self):
        res = dict()
        for current, dir in zip(self.currents, self.directions):
            res[current] = ""
            for i, element in enumerate(self.elements):
                if i == 0:
                    res[current] += "{}".format(element.n(dir))
                else:
                    res[current] += " + {}".format(element.n(dir))

        return res

    def rhs_kvl(self):
        res = dict()
        for current, dir in zip(self.currents, self.directions):
            res[current] = ""
            for i, element in enumerate(self.elements):
                if i == 0:
                    res[current] += "{}".format(element.nm1(dir))
                else:
                    res[current] += " + {}".format(element.nm1(dir))

        return res


class KVL:
    def __init__(self, CircuitPieces, directions, numberOfBranches):
        self.CircuitPieces = CircuitPieces
        self.directions = directions
        self.numberOfBranches = numberOfBranches

    def string_lhs(self, lhsName, rowNumber):
        str = ""
        res = dict()

        for vd, dir in zip(self.CircuitPieces, self.directions):
            subres = vd.lhs_kvl()

            for cur, val in subres.items():
                if cur not in res:
                    if dir == "pro":
                        res[cur] = val
                    else:
                        res[cur] = "-({})".format(val)
                else:
                    if dir == "pro":
                        res[cur] += " + {}".format(val)
                    else:
                        res[cur] += " - ({})".format(val)

        for cur, val in res.items():
            str += "{}[{}][{}] = {};\n".format(lhsName, rowNumber, cur, val)

        return str

    def string_rhs(self, rhsName, rowNumber):
        str = ""
        res = dict()

        for vd, dir in zip(self.CircuitPieces, self.directions):
            subres = vd.rhs_kvl()

            for cur, val in subres.items():
                if cur not in res:
                    if dir == "pro":
                        res[cur] = val
                    else:
                        res[cur] = "-({})".format(val)
                else:
                    if dir == "pro":
                        res[cur] += " + {}".format(val)
                    else:
                        res[cur] += " - ({})".format(val)

        str += "{}[{}] = ".format(rhsName, rowNumber)
        for i, (cur, val) in enumerate(res.items()):
            if cur < self.numberOfBranches:
                if i == 0:
                    str += "Iv[{}][n-1]*({})".format(cur, val)
                else:
                    str += " + Iv[{}][n-1]*({})".format(cur, val)
            else:
                if i == 0:
                    str += "{}".format(val)
                else:
                    str += " + {}".format(val)
        str += ";\n"

        return str

    def print_lhs(self, lhsName, rowNumber):
        print(self.string_lhs(lhsName, rowNumber), end="")

    def print_rhs(self, rhsName, rowNumber):
        print(self.string_rhs(rhsName, rowNumber), end="")


class KCL:
    def __init__(self, currents, directions, biasCurrents):
        self.currents = currents
        self.directions = directions
        self.biasCurrents = biasCurrents

    def string_lhs(self, lhsName, rowNumber):
        str = ""
        for cur, dir in zip(self.currents, self.directions):
            if dir == "pro":
                str += "{}[{}][{}] = {};\n".format(lhsName, rowNumber, cur, 1)
            else:
                str += "{}[{}][{}] = {};\n".format(lhsName, rowNumber, cur, -1)

        return str

    def string_rhs(self, rhsName, rowNumber):
        str = ""
        str += "{}[{}] = ".format(rhsName, rowNumber)
        for i, bc in enumerate(self.biasCurrents):
            if i == 0:
                str += "{}".format(bc.static("pro"))
            else:
                str += " + {}".format(bc.static("pro"))
        str += ";\n"

        return str

    def print_lhs(self, lhsName, rowNumber):
        print(self.string_lhs(lhsName, rowNumber), end="")

    def print_rhs(self, rhsName, rowNumber):
        print(self.string_rhs(rhsName, rowNumber), end="")


class CapEq:
    def __init__(self, capacitor, numberOfBranches):
        self.capacitor = capacitor
        self.numberOfBranches = numberOfBranches

    def string_lhs(self, lhsName, rowNumber):
        str = ""
        for cur, val in self.capacitor.lhs_cap().items():
            str += "{}[{}][{}] = {};\n".format(lhsName, rowNumber, cur, val)

        return str

    def string_rhs(self, rhsName, rowNumber):
        str = ""
        str += "{}[{}] = ".format(rhsName, rowNumber)
        for i, (cur, val) in enumerate(self.capacitor.rhs_cap().items()):
            if cur < self.numberOfBranches:
                if i == 0:
                    str += "{}*Iv[{}][n-1]".format(val, cur)
                else:
                    str += " + {}*Iv[{}][n-1]".format(val, cur)
            else:
                if i == 0:
                    str += "{}".format(val)
                else:
                    str += " + {}".format(val)
        str += ";\n"

        return str

    def print_lhs(self, lhsName, rowNumber):
        print(self.string_lhs(lhsName, rowNumber), end="")

    def print_rhs(self, rhsName, rowNumber):
        print(self.string_rhs(rhsName, rowNumber), end="")
