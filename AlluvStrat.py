import base
reload(base) # For interactiveness

delta = base.structure()
delta.initialize()
for t in range(a.nt):
  delta.update()

from matplotlib import pyplot as plt
from numpy import flipud
plt.imshow(flipud(a.space),interpolation='nearest')
plt.colorbar()
plt.show()
   
#a.finalize()

