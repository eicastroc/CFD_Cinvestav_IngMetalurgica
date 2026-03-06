#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 14 12:12:36 2024

@author: ecastro
"""


import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse


centroids = True
saveFig = True

X = [[1, 3, 5], [1, 3, 5]]
Y = [[0.75, 0.75, 0.75], [2.25, 2.25, 2.25]]

#%%
fig, ax = plt.subplots(figsize=(4, 4))


# centroids
ax.scatter(X, Y, marker='o', s=150, 
           c='lightblue', alpha=0.5, edgecolors='k')

# ellipse
ell = Ellipse(
    xy=(3, 1.5),
    width=4,
    height=0.5
)
ax.add_artist(ell)
ell.set_facecolor('gray')

if centroids == True:
    
    # centroids
    ax.annotate("1", xy=[1, 0.75], fontweight='bold', ha='center', va='center')
    ax.annotate("2", xy=[3, 0.75], fontweight='bold', ha='center', va='center')
    ax.annotate("3", xy=[5, 0.75], fontweight='bold', ha='center', va='center')
    ax.annotate("4", xy=[1, 2.25], fontweight='bold', ha='center', va='center')
    ax.annotate("5", xy=[3, 2.25], fontweight='bold', ha='center', va='center')
    ax.annotate("6", xy=[5, 2.25], fontweight='bold', ha='center', va='center')

    
# grid  
ax.hlines([0, 1.5, 3], xmin=[0, 0, 0], xmax=[6, 6, 6],
          linestyle='-', color='gray')
ax.vlines([0, 2, 4, 6], ymin=[0, 0, 0, 0], ymax=[3, 3, 3, 3],
          linestyle='-', color='gray')


ax.spines[['left', 'right', 'top', 'bottom']].set_visible(False) 
# grid
ax.set_xticklabels([])
ax.set_xlim(-1, 7)
ax.set_yticklabels([])
ax.set_ylim(-1.5, 4.5)


fig.tight_layout()

plt.show() 

if saveFig == True:
    name = "mesh_with_6_elements"
    fig.savefig(".".join([name,"pdf"]), dpi=300)
    fig.savefig(".".join([name,"png"]), dpi=300)