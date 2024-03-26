import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def main():
   ### Code from C++
   # // Plot each contribution using matplotlib
   # plt::figure_size(1200, 780);

   # // Create vectors for python to use
   # vector<double> y_py(ph.get_N());
   # vector<double> F_py(ph.get_N());
   # vector<double> x_py(x.data(), x.data() + x.rows() * x.cols());

   # for (int i = 0; i < ph.get_M(); i++)
   # {
   #    double mu = ctv::G_x[i];
   #    for (int j = 0; j < ph.get_N(); j++)
   #    {
   #       y_py[j] = psi_mat(i,j);
   #    }
   #    string plot_tag = "psi_mat for mu = " + to_string(mu);
   #    plt::named_plot(plot_tag, x_py, y_py);
   #    // plt::scatter(x, y_py);
   # }

   # for (int i = 0; i < ph.get_N(); i++)
   # {
   #    y_py[i] = phi[i];
   #    F_py[i] = F[i];
   # }
   # plt::named_plot("phi", x_py, y_py);
   # plt::named_plot("Radiative flux", x_py, F_py);

   # // plt::axis("on");
   # plt::xlabel("x");
   # plt::ylabel("y");

   # plt::legend();
   # plt::title("Testing");
   # plt::show();
   df_x = pd.read_csv("../build/x.csv", sep=',', header=None)
   df_phi = pd.read_csv("../build/phi.csv", sep='\\s+', header=None)
   df_phi_plus = pd.read_csv("../build/phi_plus.csv", sep='\\s+', header=None)
   df_F = pd.read_csv("../build/F.csv", sep='\\s+', header=None)
   df_psi = pd.read_csv("../build/psi.csv", sep='\\s+', header=None)

   df_eave  = pd.read_csv("../build/e_ave.csv", sep=",", header=None)
   df_rends = pd.read_csv("../build/right_ends.csv", sep=",", header=None).astype(float)
   df_lends = pd.read_csv("../build/left_ends.csv", sep=",", header=None).astype(float)

   num_G = df_phi.shape[0]
   num_N = df_phi.shape[1]
   num_M = df_psi.shape[0]

   print("num groups: ", num_G)
   print("num_cells: ", num_N)
   print("num_dir: ", num_M)

   x_arr = df_x.to_numpy()[:,0] # delimeter for x.csv was newline char
   phi_arr = df_phi.to_numpy()[0]
   
   # F
   for g in range(num_G):
      F_arr_g = df_F.to_numpy()[g]
      # print("F for g: ", F_arr_g)
      _label="F, g=" + str(g)
      plt.plot(x_arr, F_arr_g,label=_label)

   # plt.plot(x_arr, phi_arr)
   plt.legend()
   plt.savefig("F.png")
   plt.clf()

   # phi
   for g in range(num_G):
      phi_arr_g = df_phi.to_numpy()[g]
      _label = "phi, g=" + str(g)
      plt.plot(x_arr, phi_arr_g, label=_label)
   
   plt.legend()
   plt.savefig("phi.png")
   plt.clf()

   # phi plus
   for g in range(num_G):
      phi_plus_arr_g = df_phi_plus.to_numpy()[g]
      _label = "phi_plus, g=" + str(g)
      plt.plot(x_arr, phi_plus_arr_g, label=_label)
   
   plt.legend()
   plt.savefig("phi_plus.png")
   plt.clf()

   # psi 
   for m in range(num_M):
      psi_arr_m = df_psi.to_numpy()[m]
      for g in range(num_G):
         psi_arr_m_g = psi_arr_m[g::num_G]
         _label = "psi, m=" + str(m) + ", g = " + str(g)
         plt.scatter(x_arr, psi_arr_m_g, label=_label, s=15)
   
   plt.legend()
   plt.savefig("psi.png")
   # plt.show()
   plt.clf()

   # plot left and right ends
   print(df_eave.to_numpy()[:,0])
   print(df_lends.to_numpy()[:,0])
   # plt.plot(df_eave.to_numpy()[:,0], df_lends.to_numpy()[:,0], label="left ends")
   plt.plot(df_eave.to_numpy()[:,0], df_rends.to_numpy()[:,0], label="right ends")
   plt.legend()
   plt.savefig("ends.png")
   plt.clf()

   return 

# then we put main at the bottom to run everything
main()
