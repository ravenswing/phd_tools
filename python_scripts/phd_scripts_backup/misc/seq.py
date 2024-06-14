import pickle
import matplotlib.pyplot as plt
import os
import sys

if len(sys.argv) > 4:
    print("\nERROR: Too many inputs\n")
    sys.exit()
# Print all sequences
if sys.argv[1] in ["p", "print", "PRINT"]:
    pickle_in = open("sequences.dict", "rb")
    seq = pickle.load(pickle_in)
    for s in seq:
        print(s)
        print("  " + "".join(seq[s]))
# Add a sequence
elif sys.argv[1] in ["a", "add", "ADD"]:
    pickle_in = open("sequences.dict", "rb")
    seq = pickle.load(pickle_in)
    seq[sys.argv[2]] = list(sys.argv[3])
    pickle_out = open("sequences.dict", "wb")
    pickle.dump(seq, pickle_out)
    pickle_out.close()
# Remove a sequence
elif sys.argv[1] in ["r", "remove", "REMOVE"]:
    pickle_in = open("sequences.dict", "rb")
    seq = pickle.load(pickle_in)
    if sys.argv[2] in list(seq.keys()):
        del seq[sys.argv[2]]
    else:
        print("Unable to remove: Molecule not in sequences.dict")
        sys.exit()
    pickle_out = open("sequences.dict", "wb")
    pickle.dump(seq, pickle_out)
    pickle_out.close()
