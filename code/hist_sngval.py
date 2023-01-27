import matplotlib.pyplot as plt 
import numpy as np
import pandas as pd

df=pd.read_csv("sngval.txt")
print(df)

#ax = df.plot.hist(bins=1000)
#ax.set_ylim(-0.01,100)

plt.subplot(211)
hist, bins, _ = plt.hist(df, bins=1000)
plt.ylim(-0.01,200)


logbins = np.logspace(np.log10(bins[0]),np.log10(bins[-1]),len(bins))
plt.subplot(212)
plt.hist(df, bins=logbins)
plt.xscale('log')

plt.show()
