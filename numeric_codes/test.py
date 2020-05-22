# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from numpy import exp, pi, sin

import os

##########CONFIG###########
# 動画の保存形式を選択
GIF = 0
MP4 = 1
PLT = 0

TEST = 1
kappa = 1.0
R = 10.0

def analytic(t, r):
  return R/pi*exp(-kappa*pi*pi/R/R*t)*sin(pi*r/R)/r

# 初期画像を設定
def set_ax(ax):
  ax.set_ylabel("u", fontsize=21)
  ax.set_xlabel("x", fontsize=21)
  ax.set_title("1D Diffusion", fontsize=24)
  ax.grid(linestyle="dotted")
  ax.xaxis.set_tick_params(direction='in')
  ax.yaxis.set_tick_params(direction='in')
  ax.tick_params(labelsize=21)

r = np.linspace(0,R,1000+1)
r[0] = 1e-20

# 基本的な部品を宣言  
fig = plt.figure(figsize=(8, 8))
fig.subplots_adjust(left=0.15)
ax = fig.add_subplot(111)
set_ax(ax)
ax.plot(r,analytic(0,r))

plt.show()

