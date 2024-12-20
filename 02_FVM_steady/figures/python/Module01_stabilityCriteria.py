# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


fig, ax = plt.subplots(figsize=(6, 4))

# y-coord for plotitng
yPoint = 0

# plot horizontal line (x-axis)
ax.axhline(0, ls='-', color='k')

# plot nodes on the FV grid
cNodes = [-1, 0, 1]
cLabels = ['W', 'P', 'E']
for xPoint, label in zip(cNodes, cLabels):
    ax.plot(xPoint, yPoint, ls='', marker='o', color='k', markerfacecolor='w')
    ax.annotate(label, xy=(xPoint, yPoint-0.01), ha='center', va='bottom')

# Plot distances between nodes 
for xPoint in cNodes:
    ax.plot([xPoint, xPoint], [0.01, 0.04], ls='-', color='k')

arrWP = mpatches.FancyArrowPatch((-1, 0.025), (0, 0.025),
                            arrowstyle='<->,head_width=.15', mutation_scale=20)
arrPE = mpatches.FancyArrowPatch((0, 0.025), (1, 0.025),
                            arrowstyle='<->,head_width=.15', mutation_scale=20)
ax.add_patch(arrWP)
ax.add_patch(arrPE)
ax.annotate(r"$(\delta x)_w$", xy=(-0.5, 0.03), ha='center', va='bottom')
ax.annotate(r"$(\delta x)_e$", xy=(0.5, 0.03), ha='center', va='bottom')


# plot positions of surfaces between FV cells
aNodes = [-0.5, 0.5]
aLabels = ['w', 'e']
for xPoint, label in zip(aNodes, aLabels):
    ax.plot(xPoint, yPoint, ls='', marker='.', color='k')
    ax.annotate(label, xy=(xPoint, yPoint-0.006), ha='left', va='bottom')

# plot surfaces between FV cells
aLabels=['Aw', 'Ae']
for xPoint, label in zip(aNodes, aLabels):
    ax.plot([xPoint, xPoint], [-0.02, 0.02], ls='--', color='red')
    ax.annotate(label, xy=(xPoint, yPoint+0.005), 
                ha='right', va='top', color='red')

# plot the delta-x distance
ax.plot([-0.5, 0.5], [-0.02, -0.02], ls='--', color='k')
ax.plot([-0.5, 0.5], [0.02, 0.02], ls='--', color='k')
arrdx = mpatches.FancyArrowPatch((-0.5, -0.025), (0.5, -0.025),
                            arrowstyle='<->,head_width=.15', mutation_scale=20)
ax.add_patch(arrdx)
ax.annotate(r"\Delta x")

# finishing touches
ax.set_xlim(-1.5, 1.5)
ax.set_ylim(-0.05, 0.05)


