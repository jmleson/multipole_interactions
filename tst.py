
orders = {}
for a in ["x","y", "z"]:
    for b in ["x","y", "z"]:
        for c in ["x","y", "z"]:
            string = "".join(sorted([a,b,c]))
            if string in orders:
                orders[string] += 1
            else:
                orders[string] = 1
print(orders)
