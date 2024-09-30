import torch
import torch.nn as nn
import numpy as np
import pandas as pd


class GravNN(nn.Module):
      def __init__(self,N):
          super(GravNN,self).__init__()
          self.N = N
          fc1 = nn.Linear(2,self.N)
          self.fc1 = fc1
          fc2 = nn.Linear(self.N,1)
      def forward(self,fi,al):
          x = self.fc1(fi,al)
          x = torch.sigmoid(x)
          x = self.fc2(x)
          return x

crd = np.loadtxt('testgrid.grid')
N_fi = np.unique(crd[:, 0]).shape[0]
N_al = np.unique(crd[:, 1]).shape[0]
df = df_user_key_word_org = pd.read_csv("new",
                                        sep="\s+|;|:",
                                        engine="python")
# net = GeoNet(10)
v = df['P_mod'].values
fi = df['fi'].values
lb = df['lb'].values
# t = np.concatenate((fi.reshape(fi.shape[0], 1), lb.reshape(fi.shape[0], 1)), axis=1)
# t = torch.from_numpy(t)
# vt = net.forward(t.float())
# optimizer = torch.optim.Adam(net.parameters(),lr = 0.01)
# lf = torch.ones(1)*1e9
v_torch = torch.from_numpy(v)

from read_harmonics import read_gfc

koef = read_gfc()

g = GravNN(10)

fi = 45.0
lb = 45.0

fifi = np.arange(-90,90,0.25)
lblb = np.arange(0,360,0.75)
i = np.where(fifi == fi)
k = np.where(lblb == lb)

from spher_harm import get_spherical_harm
V  = get_spherical_harm(0,0,1.0,fi,lb,1.0,koef)


loss = (g(fi,al) - Pmn(sin(fi)))

