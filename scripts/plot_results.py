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
   df_phi = pd.read_csv("../build/phi.csv", delim_whitespace=True, header=None)
   df_F = pd.read_csv("../build/F.csv", delim_whitespace=True, header=None)

   num_G = df_phi.shape[0]
   num_N = df_phi.shape[1]

   print("num groups: ", num_G)
   print("num_cells: ", num_N)

   x_arr = df_x.to_numpy()[:,0] # delimeter for x.csv was newline char
   phi_arr = df_phi.to_numpy()[0]
   
   for g in range(num_G):
      F_arr_g = df_F.to_numpy()[g]
      print("F for g: ", F_arr_g)
      _label="F, g=" + str(g)
      plt.plot(x_arr, F_arr_g,label=_label)

   # print(df_phi)

   # plt.plot(x_arr, phi_arr)
   plt.legend()
   plt.show()

   return 

# then we put main at the bottom to run everything
main()