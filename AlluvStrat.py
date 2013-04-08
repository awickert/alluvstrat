#! /usr/bin/python

import base
#reload(base) # For interactiveness

delta = base.structure()
delta.initialize()
for t in range(delta.nt):
  delta.update()

from matplotlib import pyplot as plt
from numpy import flipud
fig_width = 14.
fig_height = (fig_width / delta.space.shape[-1]) * delta.space.shape[0] + fig_width / 10. # last bit for labels
plt.figure(figsize = (fig_width, fig_height) )
plt.imshow(flipud(delta.space),interpolation='nearest')
plt.xlabel(r'Cross-valley (along-strike) distance, $x$ [m]', fontsize=16)
plt.ylabel(r'Vertical distance, $z$ [m]', fontsize=16)
plt.title('Green channel deposit, blue overbank deposit, red above ground surface\n$b$ = '+str(delta.b)+' m; $h$ = '+str(delta.h)+' m; $\dot{\eta}_{ch}$ = '+str(delta.etadot_ch)+' m/yr; $\dot{\zeta}$ = '+str(delta.zetadot)+' m/yr; $x^*_{\dot{\eta}_{ob}}$ = '+str(delta.etadot_ob_xstar)+' m\n', fontsize=16)
#plt.colorbar()
plt.show()
   
#delta.finalize()

