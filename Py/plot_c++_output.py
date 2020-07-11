import numpy as np
import matplotlib.pyplot as plt
import sys

#try:
#    infile = sys.argv[1]
#except:
#    print("Need input file from C++ calculation of RDF.")
#    sys.exit(1)

rdf_data = np.loadtxt("Py/rdf_c++_output.txt")
ax1 = plt.figure().gca()
ax1.plot(rdf_data[:,0], rdf_data[:,1])
plt.xlabel("Pair Distance ($\AA$)")
plt.ylabel("g(r)")
plt.show()