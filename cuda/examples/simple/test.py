import cuda_nblist_py as nblist
import numpy as np
import time

shape = (50000, 3)
inputs = np.random.random(shape) * 10

box = nblist.Box([50.0, 50.0, 50.0])
nb = nblist.createNeighborList("celllist",box, 2.0, 0.0)
t1 = time.time()
nb.build(inputs)
for i in range(3):
    nb.update(inputs)
t2 = time.time()
print("build time: ", t2-t1)
nb.out()
