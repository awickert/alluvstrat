import base
reload(base) # For interactiveness

a = base.structure()
a.initialize()
for t in range(a.nt):
  a.update()

from matplotlib import pyplot as plt
from numpy import flipud
plt.imshow(flipud(a.space),interpolation='nearest')
plt.colorbar()
plt.show()
   
#a.finalize()

