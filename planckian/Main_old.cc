#include "NDArray.h"
#include "Energy.h"
#include "Planck.h"
#include "GLQuad.h"

#include <vector>
#include <iostream>
#include <iomanip>
#include <cmath>

using std::cout;
using std::endl;
using std::setw;
using std::left;
using std::showpos;
using std::noshowpos;
using std::cin;
using std::vector;
using ndarray::array;
using std::pow;

int main()
{
	int num_dir;                        // Number of directions
	double sum_wt=Constants::FOUR_PI;   
		
	cout <<  "Enter number of directions" << endl;
	cin >> num_dir;
	cout << "Number of directions = " << num_dir << "\n" <<endl;
	
  GLQuad quad(num_dir, sum_wt);
  const auto &mu = quad.mu();
  const auto &wt = quad.wt();
  cout << setw(16) << left << "Mu" << setw(16) << left << "Wt" << endl
       << setw(16) << left << "--" << setw(16) << left << "--" << endl;
       
  for(int i = 0;i < num_dir;++i)
  {
    cout << showpos << setw(16) << left << mu(i) << setw(16) << left << wt(i) << endl;
  }
    cout << noshowpos << endl;
 
//------------------------------
// ENTER AND PRINT SPATIAL INPUT
//------------------------------
  double thik;
  int num_cells;
  
  cout << "Enter Thickness (cm)" << endl;  // Spatial domain thickness (cm))
  cin >> thik;      
  cout << "Thickness = " << thik << "\n" << endl;            
  cout << "Enter Number of Cells" << endl; // Number of spatial cells
  cin >> num_cells; 
  cout << "Number of Cells = " << num_cells << "\n" << endl;
  		
//----------------------------------
// GENERATE AND PRINT SPATIAL ARRAYS
//----------------------------------
  
  double dx;  //cell thickness
  array<double, 1> x_edge(num_cells+1), x_ave(num_cells); //cell edge and average coordinates
  
  dx = thik/num_cells;
  
  cout << left << setw(12) << "Cell Index"
       << left << setw(11) << "Average x" 
       << left << setw(14) << "Right Edge x" << endl;
  cout << left << setw(12) << "--------"
       << left << setw(11) << "---------" 
       << left << setw(14) << "------------" << endl;
  		
  x_edge(0) = 0.0;
  for(int i = 0; i < num_cells; ++i)
  {
    x_edge(i+1) = x_edge(i)+dx;
    x_ave(i) = 0.5*(x_edge(i)+x_edge(i+1));
    
    cout << left << setw(12) << i
         << left << setw(11) << x_ave(i) 
         << left << setw(14) << x_edge(i+1) << endl;
  }
  cout << " " << endl;
//--------------------------------------
// READ AND PRINT ENERGY GRID PARAMETERS 
//--------------------------------------
  
  double efirst;   // Right edge energy for lowest energy group in keV
  double elast;    // Right edge energy for highest energy group in keV
  int num_groups;  // Number of groups
  
  cout << "Enter number of groups" << endl;
  cin >> num_groups;  
  cout << "Number of groups = " << num_groups << "\n" << endl;
  cout << "Enter right edge energy for first group (keV)" << endl;
  cin >> efirst; 
  cout << "Right edge energy for first group = " << efirst << " (keV)" << "\n" << endl;
  cout << "Enter right edge energy for last group (keV)" << endl;
  cin >> elast;
  cout << "Right edge energy for last group = " << elast << " (keV)" << "\n" << endl;   

//------------------------------------------------
// GENERATE ENERGY GROUP EDGES AND AVERAGES IN KEV
//------------------------------------------------
  array<double, 1> ekev_edge(num_groups+1), ekev_ave(num_groups);  //Group edge and average energies in keV
  double logfac;                                                   //logarithmic factor for group generation
  
  logfac = (log(elast)-log(efirst))/(num_groups-1.0);
  logfac = exp(logfac);

  ekev_edge(0) = 0.0;
  ekev_edge(1) = efirst;
  ekev_ave(0) = 0.5*(ekev_edge(0)+ekev_edge(1));
  
  for(int g = 1;g < num_groups;++g)
  {
    ekev_edge(g+1) = ekev_edge(g)*logfac;
    ekev_ave(g) = 0.5*(ekev_edge(g)+ekev_edge(g+1));
  }
  
//---------------------------------------------
// PRINT ENERGY GROUP EDGES AND AVERAGES IN KEV
//---------------------------------------------

  cout << left << setw(13) << "Group Index" 
    << left << setw(22) << "Average energy (keV)" 
  	<< left << setw(20) << "Upper energy (keV)" << endl;
  cout << left << setw(13) << "-----------"
  	<< left << setw(22) << "--------------------" 
    << left << setw(20) << "------------------" << endl;
    	
  for(int g = 0;g < num_groups;++g)
  {
  cout << left << setw(13) << g
  	<< left << setw(22) << ekev_ave(g) 
    << left << setw(20) << ekev_edge(g+1) << endl;
  }
  cout << " " << endl;	
//------------------------------------
// READ AND PRINT MATERIAL TEMPERATURE 
//------------------------------------
  double T; // Temperature in keV
  
  cout << "Enter Temperature in keV" << endl;
  cin >> T;
  cout << "Temperature (keV) = " << T << " (keV)" << "\n" << endl;
  	
//---------------------------------------------------------------------------------
// FILL ENERGY BOUND ARRAYS FOR PLANCK INTEGRATION ROUTINES (MUST CONVERT TO JERKS)
//---------------------------------------------------------------------------------
  double k = Constants::BOLTZMANN_CONSTANT;
  vector<Group> energy_discretization;  // Creates a vector of pairs of doubles for the Planck Routines
  energy_discretization.reserve(num_groups);
  	
  for(int g = 0;g < num_groups;++g)
  {
    energy_discretization.emplace_back(Group(ekev_edge(g)*k,ekev_edge(g+1)*k));
  }
  
//--------------------------
// GENERATE PLANCK INTEGRALS
//--------------------------

  Planck planck;
  ndarray::array<double, 1> B(num_groups), dBdT(num_groups);
  planck.get_Planck(T, energy_discretization, B, dBdT);

  cout << left << setw(16) << "Group"
       << left << setw(16) << "EG Min"
       << left << setw(16) << "EG Max"
       << left << setw(16) << "B"
       << left << setw(16) << "dBdT" << endl;

  cout << left << setw(16) << "-----"
       << left << setw(16) << "------"
       << left << setw(16) << "------"
       << left << setw(16) << "-"
       << left << setw(16) << "----" << endl;

  for(int g = 0;g < num_groups;++g)
  {
    cout << left << setw(16) << g 
	 << left << setw(16) << energy_discretization[g].eg_min 
	 << left << setw(16) << energy_discretization[g].eg_max
	 << left << setw(16) << B(g)
	 << left << setw(16) << dBdT(g) << endl;  
  }
  cout << " " << endl;
  
  double bsum, dbsum;
  
  bsum=0.0;
  dbsum=0.0;
  for ( int g=0;g < num_groups;g++)
  {
  bsum += B(g);
  dbsum += dBdT(g);
  }
  double ac = Constants::RADIATION_CONSTANT_A*Constants::SPEED_OF_LIGHT;
  double acT4=ac*pow(T,4);
  double dacT4=4.0*ac*pow(T,3);
  
  cout << "acT^4 = " << acT4 << " B sum = " << bsum << endl;
  cout << "4acT^3 = " << dacT4 << " dBdT sum = " << dbsum << "\n" << endl;

  return 0;
}
