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

/*
ParsedResult parseExpression(const string& input) {
    regex pattern(R"(\(([-+]?\d*\.?\d+)\(x\^(\d+)\)\+([-+]?\d*\.?\d+)\(y\^(\d+)\)\))");
    smatch matches;

    if (regex_match(input, matches, pattern)) {
        ParsedResult result;
        result.alpha = stod(matches[1].str());
        result.expx = stoi(matches[2].str());
        result.beta = stod(matches[3].str());
        result.expy = stoi(matches[4].str());
        return result;
    } else {
        throw invalid_argument("Input does not match the required format.");
    }
}
*/

// double fx(int ro, double theta, double  phi) {
   
// return ro*sin(phi)*cos(theta);
// }

// double fy(int ro, double theta, double  phi) {
   
// return ro*sin(theta)*sin(phi);
// }

// double fz(int ro, double  phi) {
   
// return ro*cos(phi);
// }

double fx (double a , double b)
{
    return (a*a) + (b*b) ;
}


int main() {
    string input;
    int ro = 4;
// Taking Inputs

/*
    cout << "Enter the expression in the format (alpha(x^expx)+beta(y^expy)): ";
    getline(cin, input);

    ParsedResult params;
    try {
        params = parseExpression(input);
        cout << "Alpha: " << params.alpha << "\n";
        cout << "Beta: " << params.beta << "\n";
        cout << "ExpX: " << params.expx << "\n";
        cout << "ExpY: " << params.expy << "\n";
    } catch (const invalid_argument& e) {
        cerr << "Error: " << e.what() << "\n";
        return 1;
    }
*/
// Plotting Graph

    // PyObject* ax = plt::chart(111);
    // plt::Clear3DChart(ax);

    // int n = 60;
    // double t0 = 0, t1 = 2*3.1428;
    // double dT1 = (t1 - t0) / (n - 1);

    // double t2 = 3.1428;
    // double dT2 = (t2 - t0) / (n - 1);

    // vector<vector<double>> x, y, z;
    // vector<double> tx, ty, tz;

    // double theta, phi;

    // for (int i = 0; i < n; ++i) {
    //     tx.clear();
    //     ty.clear();
    //     tz.clear();
    //     theta = t0 + i * dT1;

    //     for (int j = 0; j < n; ++j) {
    //         phi = t0 + j * dT2;
    //         tx.push_back(fx(ro,theta,phi));
    //         ty.push_back(fy(ro,theta,phi));
    //         tz.push_back(fz(ro,phi));
    //     }
    //     x.push_back(tx);
    //     y.push_back(ty);
    //     z.push_back(tz);
    // }

    // plt::surface3D(ax, x, y, z, "blue", 0.7);
    // plt::show();

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
