#include <iostream>
#include <vector>
#include <string>
#include <regex>
#include <stdexcept>
#include <cmath>
#include <Python.h>
#include "matplotlibcpp.h"

using namespace std;
namespace plt = matplotlibcpp;

// Parsing Values
struct ParsedResult {
    double alpha;
    double beta;
    int expx;
    int expy;
};



double fx (double a , double b)
{
    return (a*a) + (b*b) ;
}


int main() {
    string input;
    int ro = 4;

  PyObject * ax = plt::chart(111);
    plt::Clear3DChart(ax);

    int n = 60;
    double t0 = -4.0, t1 = 4.0;
    double dT = (t1 - t0)/(n - 1);

    std::vector<std::vector<double>> x, y, z;
    std::vector<double> tx, ty, tz;

    double rx, ry;

    for(int i = 0; i < n; ++i){
        tx.clear();
        ty.clear();
        tz.clear();
        rx = t0 + i*dT;
        for(int j = 0; j < n; ++j){
            ry = t0 + j*dT;
            tx.push_back(rx);
            ty.push_back(ry);
            tz.push_back(fx(rx, ry));
        }
        x.push_back(tx);
        y.push_back(ty);
        z.push_back(tz);
    }


    plt::surface3D(ax, x, y, z, "blue", 0.7);
    plt::show();
    
    return 0;
}
