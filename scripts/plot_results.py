import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def main():
   ### Code from C++
   # // Plot each contribution using matplotlib
   # plt::figure_size(1200, 780);

   # // Create vectors for python to use
   # vector<double> y_py(ctv::N);
   # vector<double> F_py(ctv::N);
   # vector<double> x_py(x.data(), x.data() + x.rows() * x.cols());

   # for (int i = 0; i < ctv::M; i++)
   # {
   #    double mu = ctv::G_x[i];
   #    for (int j = 0; j < ctv::N; j++)
   #    {
   #       y_py[j] = psi_mat(i,j);
   #    }
   #    string plot_tag = "psi_mat for mu = " + to_string(mu);
   #    plt::named_plot(plot_tag, x_py, y_py);
   #    // plt::scatter(x, y_py);
   # }

   # for (int i = 0; i < ctv::N; i++)
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
   x_arr = pd.read_csv("../build/x.csv", sep=',', header=None).to_numpy()
   phi_arr = pd.read_csv("../build/phi.csv", sep=',', header=None).to_numpy()
   psi_mat = pd.read_csv("../build/psi.csv", sep=',', header=None).to_numpy()
   F_arr = pd.read_csv("../build/F.csv", sep=',', header=None).to_numpy()
   plt.plot(x_arr, phi_arr)
   plt.show()

   return 

# then we put main at the bottom to run everything
main()