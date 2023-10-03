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
	
//--------------------------
// PLANCK FUNCTION PROTOTYPE
//--------------------------
  double pf (double, double);
//----------------
// PLANCK FUNCTION
//----------------
//------------------------------------------------------
// PHOTON ENERGY E AND TEMPERATURE T ARE IN UNITS OF keV
// THE FUNCTION HAS UNITS OF jk/(cm^2-sh-keV-steradian)
//------------------------------------------------------
  double pf (double E, double T)
  {
  double h = Constants::PLANCK_CONSTANT;
  double c = Constants::SPEED_OF_LIGHT;
  double k = Constants::BOLTZMANN_CONSTANT_JPK;
  double B;

  B = k*pow(E,3)/(pow(h,3)*pow(c,2)*(exp(E/T)-1.0));
  return B;
  }

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
  cout << "\n" << endl;
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
  array<double, 1> e_edge(num_groups+1), e_ave(num_groups), de_ave(num_groups);  //Group edge and average energies and average group widths in kev
  double logfac;                                                                 //Logarithmic factor for group generation
  
  logfac = (log(elast)-log(efirst))/(num_groups-1.0);
  logfac = exp(logfac);

  e_edge(0) = 0.0;
  e_edge(1) = efirst;
  e_ave(0) = 0.5*(e_edge(0)+e_edge(1));
  de_ave(0) = e_edge(1) - e_edge(0);
  
  for(int g = 1;g < num_groups;++g)
  {
    e_edge(g+1) = e_edge(g)*logfac;
    e_ave(g) = 0.5*(e_edge(g)+e_edge(g+1));
    de_ave(g) = e_edge(g+1) - e_edge(g);
  }
  
//---------------------------------------------
// PRINT ENERGY GROUP EDGES AND AVERAGES IN KEV
//---------------------------------------------

  cout << left << setw(13) << "Group Index" 
       << left << setw(16) << "Average Energy" 
       << left << setw(14) << "Upper Energy" 
       << left << setw(13) << "Group Width" << endl;
  cout << left << setw(13) << "-----------"
       << left << setw(16) << "(keV)---------" 
       << left << setw(14) << "(keV)-------" 
    	 << left << setw(13) << "(keV)------" << endl;
  for(int g = 0;g < num_groups;++g)
  {
    cout << left << setw(13) << g
    	   << left << setw(16) << e_ave(g) 
         << left << setw(14) << e_edge(g+1) 
         << left << setw(13) << de_ave(g) << endl;
  }
  cout << "\n" << endl;	
//------------------------------------
// READ AND PRINT MATERIAL TEMPERATURE 
//------------------------------------
  double T; // Temperature in keV
  
  cout << "Enter Temperature in keV" << endl;
  cin >> T;
  cout << "Temperature (keV) = " << T << " (keV)" << "\n" << endl;
  	
//---------------------------------------------------------
// FILL ENERGY BOUND ARRAYS FOR PLANCK INTEGRATION ROUTINES
//---------------------------------------------------------

  vector<Group> energy_discretization;  // Creates a vector of pairs of group energies edges for the Planck routines
  energy_discretization.reserve(num_groups);
  	
  for(int g = 0;g < num_groups;++g)
  {
    energy_discretization.emplace_back(Group(e_edge(g),e_edge(g+1)));
  }
  
//----------------------------------------------------
// GENERATE PLANCK INTEGRALS IN UNITS OF KEV/(CM^2-SH)
//----------------------------------------------------

  Planck planck;
  ndarray::array<double, 1> B(num_groups), dBdT(num_groups);
  planck.get_Planck(T, energy_discretization, B, dBdT);
  
//---------------------------------------------
// CHANGE PLANCK INTEGRAL UNITS TO JK/(CM^2-SH)
//---------------------------------------------

  double kcon = Constants::BOLTZMANN_CONSTANT_JPK;
  	
  for(int g = 0;g < num_groups;++g)
  {
    B(g)    = kcon*B(g);
    dBdT(g) = kcon*dBdT(g);
  }
  
  

  cout << left << setw(7)  << "Group"
       << left << setw(14) << "EG Min"
       << left << setw(14) << "EG Max"
       << left << setw(14) << "B"
       << left << setw(22) << "dBdT" << endl;

  cout << left << setw(7)  << "-----"
       << left << setw(14) << "(keV)-"
       << left << setw(14) << "(keV)-"
       << left << setw(14) << "(jk/cm^2-sh)"
       << left << setw(18) << "(jk/cm^2-sh-keV)" << endl;

  for(int g = 0;g < num_groups;++g)
  {
   cout << left << setw(7)  << g 
	      << left << setw(14) << energy_discretization[g].eg_min 
	      << left << setw(14) << energy_discretization[g].eg_max
	      << left << setw(14) << B(g)
	      << left << setw(18) << dBdT(g) << endl;  
  }
  cout << "\n" << endl;
  
//------------------------------------------
// CHECK FOR PROPER SUMS OF PLANCK INTEGRALS
//------------------------------------------
  
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
  
//------------------------------
// GENERATE MULTIGROUP OPACITIES
//------------------------------
  
  ndarray::array<double, 1> ukappa(num_groups), kappa(num_groups);  // Unnormalized and final opacities respectively
  ndarray::array<double, 1> ckappa(num_groups);                     // Opacities evaluated at group center energies
  double kappa_nfac;                                                // Opacity normalization factor
  double kappa_grey;                                                // Grey opacity used to define total emission at a given temperature
  double fourpi = Constants::FOUR_PI;    
  double hcon = Constants::PLANCK_CONSTANT;
  double ccon = Constants::SPEED_OF_LIGHT;
  double tmp1;                                                      // Temporary variable
  double tmp2;                                                      // Temporary variable
  	                               
//----------------------------
// READ AND PRINT GREY OPACITY                                   
//----------------------------                                               
  
  cout << "Enter grey opacity in cm^2/g" << endl;
  cin >> kappa_grey;
  cout << "grey opacity  = " << kappa_grey << "(cm^2/g)" << "\n" << endl; 
  
//-----------------------------------------------------------
// COMPUTE UNNOMALIZED OPACITIES AND GROUP-CENTERED OPACITIES
//-----------------------------------------------------------
  
  for(int g = 0;g < num_groups;++g)
  {
    tmp1 = 2.0*fourpi*kcon*T/(pow(hcon,3)*pow(ccon,2));
    tmp2 = exp(-e_edge(g)/T) - exp(-e_edge(g+1)/T); 
    ukappa(g) = tmp1*tmp2/B(g);
    ckappa(g) = (1.0-exp(-e_ave(g)/T))/pow(e_ave(g),3);
  }
  
//----------------------------------------------------------
// PRINT UNNORMALIZED OPACITIES AND GROUP-CENTERED OPACITIES
//----------------------------------------------------------

  cout << left << setw(7)  << "Group"
       << left << setw(14) << "Energy"
       << left << setw(25) << "Planck-Averaged Opacity"
       << left << setw(24) << "Group-Centered Opacity" << endl; 

  cout << left << setw(7)  << "-----"
       << left << setw(14) << "(keV)-"
       << left << setw(25) << "(cm^2/g)---------------"
       << left << setw(24) << "(cm^2/g)--------------" << endl;
       
  for(int g = 0;g < num_groups;++g)
  {
  cout << left << setw(7)  << g 
       << left << setw(14) << e_ave(g)
       << left << setw(25) << ukappa(g)
       << left << setw(24) << ckappa(g) << endl; 
  }
  cout << "\n" << endl;
//-------------------------------
// COMPUTE NORMALIZATION CONSTANT
//-------------------------------
  double emis_tot = 0.0;

  for(int g = 0;g < num_groups;g++)
  {
  emis_tot += ukappa(g)*B(g);
  }
  kappa_nfac = acT4/emis_tot;
  
//--------------------------------------------------------
// COMPUTE AND PRINT FINAL OPACITIES AND EMISSION SPECTRUM
//--------------------------------------------------------

  cout << left << setw(7)  << "Group"
       << left << setw(14) << "Energy"
       << left << setw(15) << "Final Opacity"
       << left << setw(15) << "Emission Spec" << endl; 

  cout << left << setw(7)  << "-----"
       << left << setw(14) << "(keV)-"
       << left << setw(15) << "(cm^2/g)-----"
       << left << setw(15) << "(jk-g-sh)----" << endl;
       

  ndarray::array<double, 1> emis_spec(num_groups);

  for(int g = 0;g < num_groups;g++)
  {
  kappa(g)= kappa_grey*kappa_nfac*ukappa(g);
  emis_spec(g) = kappa(g)*B(g);
  
  cout << left << setw(7)  << g 
       << left << setw(14) << e_ave(g)
       << left << setw(15) << kappa(g)
       << left << setw(15) << emis_spec(g) << endl; 
  }
  cout << "\n" << endl;
//-----------------------------
//CHECK FOR PROPER EMISSION SUM
//-----------------------------

  double sigacT4 = kappa_grey*acT4;
  
  emis_tot = 0.0;
  for(int g = 0;g < num_groups;g++)
  {
    emis_tot += kappa(g)*B(g);
  }

  cout << "Total Emission = " << emis_tot 
       << " kappa*acT^4 = " << sigacT4 << "\n"  << endl;
       
//----------------------------------------
// COMPUTE COMPONENTS OF CORRECTION SOURCE
//----------------------------------------

  ndarray::array<double, 1> dEB(num_groups), dkapEdE(num_groups), dkapEB(num_groups); 

//---------------------------------
// COMPUTE ENERGY DIFFERENCES OF EB
//---------------------------------

  dEB(0) = e_edge(1)*pf(e_edge(1),T);

  for(int g = 1;g < num_groups-1;g++)
  {
    dEB(g) = e_edge(g+1)*pf(e_edge(g+1),T) - e_edge(g)*pf(e_edge(g),T);
  }
  dEB(num_groups-1) = - e_edge(num_groups-1)*pf(e_edge(num_groups-1),T);
  
//-----------------------------------
//  CHECK FOR ZERO SUM OF DIFFERENCES
//-----------------------------------
    double sum = 0.0;
    
    for(int g = 0;g < num_groups;g++)
    {
    	sum += dEB(g);
    }
//-----------
//  PRINT SUM
//-----------
  cout << "Sum EB = " << sum << "\n" << endl;
    
//-------------------------------
// PRINT ENERGY DIFFERENCES OF EB
//-------------------------------
  cout << left << setw(7)  << "Group"
       << left << setw(14) << "Energy" 
       << left << setw(14) << "dEB" << endl;
       
  cout << left << setw(7)  << "-----"
       << left << setw(14) << "------" 
       << left << setw(14) << "-----" << endl;
       
   
  for(int g = 0;g < num_groups;g++)
  {
    cout << left << setw(7)  << g 
         << left << setw(14) << e_ave(g) 
         << left << setw(14) << dEB(g) << endl;
  }
    cout << "\n" << endl;
    
//-----------------------------
// COMPUTE GROUP-EDGE OPACITIES
//-----------------------------
  ndarray::array<double, 1> kappa_edge(num_groups+1);  // Edge opacities
  double wgt_L, wgt_R;                                 // Left and right weights for edge opacity generation
                                     
  kappa_edge(0) = kappa(0);
  
  for(int g = 1;g < num_groups;g++)
  {
  	wgt_L = (e_ave(g)-e_edge(g))/(e_ave(g) - e_ave(g-1));
  	wgt_R = (e_edge(g)-e_ave(g-1))/(e_ave(g) - e_ave(g-1));
    kappa_edge(g) = kappa(g-1)*wgt_L + kappa(g)*wgt_R;
  }  
  kappa_edge(num_groups) = kappa(num_groups-1);
  
//---------------------------
// PRINT GROUP EDGE OPACITIES
//---------------------------
  cout << left << setw(12) << "Edge Index"
       << left << setw(14) << "Energy" 
       << left << setw(14) << "Opacity" << endl;
       
  cout << left << setw(12) << "---------"
       << left << setw(14) << "------" 
       << left << setw(14) << "-------" << endl;
       
   
  for(int g = 0;g < num_groups+1;g++)
  {
    cout << left << setw(12)  << g 
         << left << setw(14) << e_edge(g) 
         << left << setw(14) << kappa_edge(g) << endl;
  }
    cout << "\n" << endl;
//--------------------------------------
// COMPUTE ENERGY DERIVATIVES OF KAPPA E
//--------------------------------------
  dkapEdE(0) = kappa_edge(1)*e_edge(1)/de_ave(0);
  for(int g = 1;g < num_groups-1;g++)
  {
    dkapEdE(g) = (kappa_edge(g+1)*e_edge(g+1) - kappa_edge(g)*e_edge(g))/de_ave(g);
  }
  dkapEdE(num_groups-1) = -kappa_edge(num_groups)*e_edge(num_groups)/de_ave(num_groups-1);
//------------------------------------------
//  CHECK INTEGRAL OF DERIVATIVES OF KAPPA E 
//------------------------------------------ 
    double sumabs=0.0;  
    for(int g = 0;g < num_groups-1;g++)
    {
    	sum += dkapEdE(g)*de_ave(g);
    	sumabs += fabs(dkapEdE(g))*de_ave(g);
    }
    sum=sum/sumabs;
//----------------
//  PRINT INTEGRAL
//----------------
  cout << "Integral dkapEdE/Integral |dkapEdE| = " << sum << "\n" << endl;
  
//-------------------------------------
//  PRINT ENERGY DERIVATIVES OF KAPPA E
//-------------------------------------
  cout << left << setw(7)  << "Group"
       << left << setw(14) << "Energy" 
       << left << setw(14) << "dkapEdE" << endl;
       
  cout << left << setw(7)  << "-----"
       << left << setw(14) << "------" 
       << left << setw(14) << "-------" << endl;
       
  for(int g = 0;g < num_groups;g++)
  {
    cout << left << setw(7)  << g 
         << left << setw(14) << e_ave(g) 
         << left << setw(14) << dkapEdE(g) << endl;
  }
  cout << "\n" << endl;
  
//---------------------------------------
//  COMPUTE ENERGY DIFERENCES OF KAPPA*EB 
//---------------------------------------
     
  dkapEB(0) = kappa_edge(1)*e_edge(1)*pf(e_edge(1),T);
  for(int g = 1;g < num_groups-1;g++)
  {  
    dkapEB(g) = kappa_edge(g+1)*e_edge(g+1)*pf(e_edge(g+1),T) - kappa_edge(g)*e_edge(g)*pf(e_edge(g),T);
  }
  dkapEB(num_groups-1) = - kappa_edge(num_groups-1)*e_edge(num_groups-1)*pf(e_edge(num_groups-1),T);
  
//-----------------------------------
//  CHECK FOR ZERO SUM OF DIFFERENCES
//-----------------------------------
    sum = 0.0;   
    for(int g = 0;g < num_groups;g++)
    {
    	sum += dkapEB(g);
    	sumabs += fabs(dkapEB(g));
    }
  sum=sum/sumabs;
//-----------
//  PRINT SUM
//-----------
  cout << "Sum dkapEB/Sum |dkapEB| = " << sum << "\n" << endl;
  
//------------------------------------
// PRINT ENERGY DIFERENCES OF KAPPA*EB
//------------------------------------  

  cout << left << setw(7)  << "Group"
       << left << setw(14) << "Energy" 
       << left << setw(14) << "dkapEB" << endl;
       
  cout << left << setw(7) << "-----"
       << left << setw(14) << "------" 
       << left << setw(14) << "------" << endl;
       
   
  for(int g = 0;g < num_groups;g++)
  {
    cout << left << setw(7)  << g 
         << left << setw(14) << e_edge(g) 
         << left << setw(14) << dkapEB(g) << endl;
  }
  cout << "\n" << endl;
  
//--------------------------------------------
//  GENERATE COMPONENTS OF THE CORRECTION TERM
//--------------------------------------------
  ndarray::array<double, 1> cor1(num_groups), cor2(num_groups), cor3(num_groups); // correction terms 1, 2, and 3

   for(int g = 0;g < num_groups;g++)
  {
     cor1(g) = dkapEdE(g);
     cor2(g) = 3.0*kappa(g)*B(g) - dkapEB(g);
     cor3(g) = cor1(g)*(4.0*B(g)-dEB(g));
  }
//
//  PRINT CORRECTION TERMS
// 
  cout << left << setw(7)  << "Group"
       << left << setw(14) << "Energy" 
       << left << setw(14) << "cor1" 
       << left << setw(14) << "cor2" 
       << left << setw(14) << "cor3" << endl;
 
       
  cout << left << setw(7)  << "-----"
       << left << setw(14) << "------" 
       << left << setw(14) << "----" 
       << left << setw(14) << "----" 
       << left << setw(14) << "----" << endl;
       
   
  for(int g = 0;g < num_groups;g++)
  {
    cout << left << setw(7)  << g 
         << left << setw(14) << e_ave(g) 
         << left << setw(14) << cor1(g) 
         << left << setw(14) << cor2(g) 
         << left << setw(14) << cor3(g) << endl;
  }
  cout << "\n" << endl;
  
  return 0;
}
