from term_components.Multipole import MultipoleMoment



def get_replacements_according_to_traceless_conditions(molecule:str):
    if molecule not in ["A", "B"]:
        raise Exception("unknown molecule")
    replacements = []

    m = MultipoleMoment(["x", "x"])
    m.molecule = molecule
    m2 = MultipoleMoment(["y", "y"])
    m2.molecule = molecule
    m3 = MultipoleMoment(["z", "z"])
    m3.molecule = molecule
    replacements.append({"to_be_replaced": m, "replacements": [m2, m2] })


    return replacements