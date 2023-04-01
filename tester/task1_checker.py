import sys
from typing import List, Set

def read_trusses(truss_file):
    truss_file = open(truss_file, 'r')
    all_trusses = []
    line = truss_file.readline().strip()
    while line:
        truss_exists = int(line)
        if truss_exists == 0:
            line = truss_file.readline().strip()
            continue
        n = int(truss_file.readline().strip())
        trusses = set()
        for i in range(n):
            trusses.add(truss_file.readline().strip())

        all_trusses.append(trusses)
        line = truss_file.readline().strip()

    return all_trusses

def compare_trusses(gold: list[set[str]], pred: list[set[str]]):
    for (gold_trusses, pred_trusses) in zip(gold, pred):
        if gold_trusses != pred_trusses:
            print("Unequal trusses")
            print(f"{pred_trusses.difference(gold_trusses)} present in pred but not in gold")
            print(f"{gold_trusses.difference(pred_trusses)} present in gold but not in pred")
            return False
    return True

if __name__ == "__main__":

    gold_trusses = read_trusses(sys.argv[1])
    pred_trusses = read_trusses(sys.argv[2])

    if compare_trusses(gold_trusses, pred_trusses):
        print("All trusses equal")
