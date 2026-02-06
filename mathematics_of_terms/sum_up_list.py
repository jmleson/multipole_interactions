
####### dummy functions for testing : #####

def addable_fictitious(i:int, j:int):
    if i < 0 or j < 0:
        return False
    return True

def add_fictitious(i:int, j:int) -> list[int]:
    return [i+j]




####### actual relevant functions: ######

def combine_general(acc, x, adding_function:callable, addable_function:callable):
    for i, a in enumerate(acc):
        if addable_function(a, x):
            # since addable -> result of adding_function has to have length 1
            acc[i] = adding_function(x, acc[i])[0]
            return acc
    acc.append(x)
    return acc

def sum_up_list(to_be_summed_up_list: list, adding_function:callable, addable_function:callable):
    result = [to_be_summed_up_list[0]]
    for x in to_be_summed_up_list[1:]:
        result = combine_general(result, x, adding_function=adding_function, addable_function=addable_function)
    return result


