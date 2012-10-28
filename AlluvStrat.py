import base
reload(base) # For interactiveness

delta = base.structure()
delta.initialize()
for t in range(delta.nt):
  delta.update()

from matplotlib import pyplot as plt
from numpy import flipud
plt.imshow(flipud(delta.space),interpolation='nearest')
plt.colorbar()
plt.show()
   
#delta.finalize()

