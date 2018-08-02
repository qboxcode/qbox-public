#include<iostream>
#include "Function3d.h"
using namespace std;

int main()
{
  Function3d f;
  f.read("test.xml");

  cout << "function name: " << f.name_ << endl;
  cout << "function base64 size: " << f.str_.size() << endl;
  cout << "function array size: " << f.val_.size() << endl;

  double sum = 0.0;
  for ( int i = 0; i < f.val_.size(); i++ )
  {
    sum += f.val_[i]*f.val_[i];
  }
  cout << "function norm2: " << sum / f.val_.size() << endl;
}
