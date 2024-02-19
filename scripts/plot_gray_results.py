import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def main():
   df_x = pd.read_csv("../build/tests/gray-test-x.csv", sep=',', header=None)
   df_phi = pd.read_csv("../build/tests/gray-test-phi.csv", delim_whitespace=True, header=None)
   df_F = pd.read_csv("../build/tests/gray-test-F.csv", delim_whitespace=True, header=None)

   num_G = df_phi.shape[0]
   num_N = df_phi.shape[1]

   print("num groups: ", num_G)
   print("num_cells: ", num_N)

   x_arr = df_x.to_numpy()[:,0] # delimeter for x.csv was newline char
   phi_arr = df_phi.to_numpy()[0]
   
   for g in range(num_G):
      F_arr_g = df_F.to_numpy()[g]
      _label="g=" + str(g)
      plt.scatter(x_arr, F_arr_g,label=_label)

   # print(df_phi)

   # plt.plot(x_arr, phi_arr)
   plt.legend()
   plt.title("Radiative Flux")
   plt.show()

   return 

# then we put main at the bottom to run everything
main()