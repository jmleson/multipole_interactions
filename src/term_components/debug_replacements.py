from src.term_components.Multipole import MultipoleMoment



def get_replacements_according_to_traceless_conditions(molecule:str):
    if molecule not in ["A", "B"]:
        raise Exception("unknown molecule")
    replacements = []

    mx = MultipoleMoment(["x", "x"])
    mx.molecule = molecule
    my = MultipoleMoment(["y", "y"])
    my.molecule = molecule
    mz = MultipoleMoment(["z", "z"])
    mz.molecule = molecule
    replacements.append({"to_be_replaced": mz, "replacements": [mx, my] })

    # m = MultipoleMoment(["z", "z", "y"])
    # m.molecule = molecule
    # m2 = MultipoleMoment(["y", "y", "y"])
    # m2.molecule = molecule
    # m3 = MultipoleMoment(["x", "x", "y"])
    # m3.molecule = molecule
    # replacements.append({"to_be_replaced": m, "replacements": [m2, m3]})

    m = MultipoleMoment(["z", "z", "y"])
    m.molecule = molecule
    m2 = MultipoleMoment(["y", "y", "y"])
    m2.molecule = molecule
    m3 = MultipoleMoment(["x", "x", "y"])
    m3.molecule = molecule
    replacements.append({"to_be_replaced": m3, "replacements": [m, m2]})

    # x double:
    mx = MultipoleMoment(["x", "x", "x", "x"])
    mx.molecule = molecule
    my = MultipoleMoment(["x", "x", "y", "y"])
    my.molecule = molecule
    mz = MultipoleMoment(["x", "x", "z", "z"])
    mz.molecule = molecule
    replacements.append({"to_be_replaced": mz, "replacements": [mx, my] })

    # y double:
    mx = MultipoleMoment(["y", "y", "x", "x"])
    mx.molecule = molecule
    my = MultipoleMoment(["y", "y", "y", "y"])
    my.molecule = molecule
    mz = MultipoleMoment(["y", "y", "z", "z"])
    mz.molecule = molecule
    replacements.append({"to_be_replaced": mz, "replacements": [mx, my] })

    # z double:
    mx = MultipoleMoment(["z", "z", "x", "x"])
    mx.molecule = molecule
    my = MultipoleMoment(["z", "z", "y", "y"])
    my.molecule = molecule
    mz = MultipoleMoment(["z", "z", "z", "z"])
    mz.molecule = molecule
    replacements.append({"to_be_replaced": mz, "replacements": [mx, my] })

    return replacements