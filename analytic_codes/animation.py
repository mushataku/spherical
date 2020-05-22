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
###########################

#--------------def-------------
def readcsv(datafilename):
  df=pd.read_csv(datafilename,header=None)
  dft = np.array(df).T
  dftr = dft[:,1:]
  dftru = np.array(dftr)
  return dftru

################### PARAMETER ##################
FILE_PATH = '../figs/animation'
df_time = pd.read_csv("../data/output_time.csv")
T = df_time["time"]
df_F = pd.read_csv("../data/conservation.csv")
F = df_F["F"]
co = readcsv("../data/condition.csv")
Lr = float(co[0])
Lmax = float(co[1])
kappa = float(co[4])
################### PARAMETER ##################

# 解析解
def analytic(t, r):
  return Lmax/pi*exp(-kappa*pi*pi/Lmax/Lmax*t)*sin(pi*r/Lmax)/r

#########################描画のための関数#########################

# frame 枚目の写真を作るデータを取得
def get_f(frame):
  fp = "../data/f/%.3f.csv" % T[frame]
  df_f = pd.read_csv(fp)
  r = df_f["r"]
  f = df_f["f"]
  return r,f

# 初期画像を設定
def init_f(ax):
  ax.set_ylabel("f", fontsize=21)
  ax.set_xlabel("r", fontsize=21)
  ax.set_title("Spherical Diffusion", fontsize=24)
  ax.grid(linestyle="dotted")
  ax.xaxis.set_tick_params(direction='in')
  ax.yaxis.set_tick_params(direction='in')
  ax.tick_params(labelsize=21)
  # ax.set_yscale("log")
  r,f = get_f(0)
  if TEST == 0:
    im_f, = ax.plot(r,f, "r", label="numeric")
    ax.legend(fontsize=20)
    return im_f
  else:
    im_f_analytic, = ax.plot(r,analytic(T[0],r), "b", label="analytic")
    im_f, = ax.plot(r,f, "r", label="numeric")
    ax.legend(fontsize=20)
    return im_f, im_f_analytic

# データ更新
def reset_f(im,frame):
  r,f = get_f(frame)
  im.set_data(r,f)

def reset_f_analytic(im,frame):
  r,_ = get_f(frame)
  im.set_data(r,analytic(T[frame],r))

#########################描画のための関数#########################



################################################################
########################### main ###############################
################################################################

# 基本的な部品を宣言  
fig = plt.figure(figsize=(8, 8))
fig.subplots_adjust(left=0.15)
ax_f = fig.add_subplot(111)
time_text = fig.text(0.01, 0.99, '', size=20, color="white", horizontalalignment='left',
            verticalalignment='top', backgroundcolor='black')
# kappa 出力 後でやる
fig.text(0, 0.01, r"$\kappa$="+str(kappa),
          backgroundcolor="black",color="white", size=20)

#### アニメの初期画像生成
if TEST == 0:
  im_f = init_f(ax_f)
else:
  im_f, im_f_analytic = init_f(ax_f)
time_text.set_text("time = 0.000\n"+r"$\int$fdr = %.3f"%F[0])

#### 画像更新用関数
def animate(frame):
  if(frame == 0):
    return
  print("frame:",frame)
  reset_f(im_f, frame)
  if(TEST):
    reset_f_analytic(im_f_analytic, frame)
  
  time_text.set_text("time = %.3f\n"%T[frame]+r"$\int$fdr = %.3f"%F[frame])

ani = FuncAnimation(fig, animate, frames=int(len(T))
              , interval=200, repeat=True, blit=False)


if(GIF == 1):
    ani.save(FILE_PATH+".gif", writer='pillow')
if(MP4 == 1):
    ani.save(FILE_PATH+".mp4", writer="ffmpeg", fps=5)
if(PLT == 1):
    plt.show()

