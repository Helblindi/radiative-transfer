#include "NDArray.h"

#include <cmath>
#include <algorithm>
#include <iostream>
#include <iomanip>

using ndarray::array;
using std::fill;
using std::cout;
using std::endl;

void oned(array<double, 1> &a, array<double, 1> &b, int N) 
{
  double val = 0.0;
  for(int i = 0;i < N;++i)
  {
    a(i) = val;
    b(i) = a(i)*2.0;
    val += 0.1;
  }
}

void twod(array<double, 2> &c, int M, int N)
{
  double val = 5.0;
  for(int i = 0;i < M;++i)
  {
    for(int j = 0;j < N;++j)
      c(i,j) = 5.0*val;
  }
}

void threed(array<double, 3> &d, int M, int N, int K)
{
  for(int i = 0;i < M;++i)
  {
    for(int j = 0;j < N;++j)
    { 
      for(int k = 0;k < K;++k)
      {	        
	if(i == M/2)
          d(i, j, k) = 2.0;
	else
          d(i, j, k) = 1.0;
      }
    }
  }
}

int main()
{
  array<double, 1> a(10), b(10);
  oned(a, b, 10);
  for(auto n : a)
    cout << n << " ";
  cout << endl;
  for(auto n : b)
    cout << n << " ";
  cout << endl << endl;

  array<double, 2> c(3, 5);
  twod(c, 3, 5);
  for(int i = 0;i < 3;++i)
  {
    for(int j = 0;j < 5;++j)
      cout << c(i,j) << " ";
    cout << endl;
  }
  cout << endl;

  array<double, 3> d(5, 2, 5);
  threed(d, 5, 2, 5);
  for(int i = 0;i < 5;++i)
  {
    cout << "Slice: " << i << std::endl;
    for(int j = 0;j < 2;++j)
    {
      for(int k =0;k < 5;++k)
        cout << d(i,j,k) << " ";
      cout << endl;
    }
  }
  cout << endl;

  // This is quick way to fill an array of any dimension if all values are the same
  fill(a.begin(), a.end(), 0.0);
  fill(d.begin(), d.end(), 1.0);

  return 0;
}
