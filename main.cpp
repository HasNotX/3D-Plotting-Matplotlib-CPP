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


double func (double a , double b)
{
    return a*a + b*b;
}



class ThreeD {
    private:
        string equation;
    public:

        ThreeD(){  // Constructor
            //equation = s;
        }

        virtual void input_equation() {};
        virtual void plot() {};
};

class z1 : public ThreeD{

    private:
        vector<vector<double>> x,y,z; // Coordiante 2-D Matrices Storing coordinates
        vector<double> tx,ty,tz; // 1-D temporary coordinates matrices

        int n = 60; // No. of times points are generated
        double t0 = -4.0 , t1 = 4.0; // Range of points to be generated
        double dT = (t1 - t0) / (n - 1); // Small increment value 
        double rx, ry; // Temporary Coordinates

    public:
         z1 () { // Constructor

         }
        //Function to Calculate coordinates
        double func (double a , double b){
            return a*a + b*b;
        }

        void plot () {
        PyObject* ax = plt::chart(111);
        plt::Clear3DChart(ax);

        for (int i = 0; i < n; ++i) {
            tx.clear();
            ty.clear();
            tz.clear();
            rx = t0 + i * dT;

            for (int j = 0; j < n; ++j) {
                ry = t0 + j * dT;
                tx.push_back(rx);
                ty.push_back(rx);
                tz.push_back(func(rx,ry));
            }
            x.push_back(tx);
            y.push_back(ty);
            z.push_back(tz);
    }

    plt::surface3D(ax, x, y, z, "blue", 0.7);
    plt::show();
    }
};

class z2 : public ThreeD{
private:
    double t0 , t1 , t2 , t3;
    double dT1 , dT2;
    int n = 60;
    
public:
    int a , b , c ;
    vector<vector<double>> x,y,z,z2; // Coordiante 2-D Matrices Storing coordinates
    vector<double> tx,ty,tz,tz2; // 1-D temporary coordinates matrices
    double theta,phi;
    virtual double x_calc(int param, double theta, double phi) = 0;
    virtual double y_calc(int param, double theta, double phi) = 0;
    virtual double z_calc(int param, double phi) = 0;
    
//Function to plot
    void calc_bounds (double l , double m , double o , double p)
    {
     t0 = l, t1 = m;
     dT1 = (t1 - t0) / (n - 1);

     t2 = o , t3 = p;
     dT2 = (t2 - t3) / (n - 1);
    }

    void plot(double l , double m , double o , double p, int a , int b , int c){
        
        calc_bounds (l,m,o,p);

        PyObject* ax = plt::chart(111);
        plt::Clear3DChart(ax);

        

        for (int i = 0; i < n; ++i) {
            tx.clear();
            ty.clear();
            tz.clear();
            theta = t0 + i * dT1;

           for (int j = 0; j < n; ++j) {
                phi = t2 + j * dT2;
                tx.push_back(x_calc(a,theta,phi));
                ty.push_back(y_calc(b,theta,phi));
                tz.push_back(z_calc(c,phi));
            }
            x.push_back(tx);
            y.push_back(ty);
            z.push_back(tz);
        }

        plt::surface3D(ax, x, y, z, "blue", 0.7);
        plt::show();
    }
};


class sphere : public z2{

public:

// Functions to calculate coordinates
    double x_calc(int ro, double theta, double  phi) {
        return ro*sin(phi)*cos(theta);
    }

    double y_calc(int ro, double theta, double  phi) {
        return ro*sin(theta)*sin(phi);
    }

    double z_calc(int ro, double  phi) {
        return ro*cos(phi);
    }


};



class ellipsoid : public z2{
     // Denominator values
    private:

    public:

    // Functions to calculate coordinates
        double x_calc(int denom, double theta, double  phi) {
            return denom*sin(phi)*cos(theta);
        }

        double y_calc(int denom, double theta, double  phi) {
            return denom*sin(theta)*sin(phi);
        }

        double z_calc(int denom, double  phi) {
            return denom*cos(phi);
        }

};


class cone : public z2{
    private:

    public:

    double x_calc(int denom, double theta, double  phi) {
        return denom*cosh(phi)*cos(theta);
    }

    double y_calc(int denom, double theta, double  phi) {
        return denom*sin(theta)*cosh(phi);
    }

    double z_calc(int denom, double  phi) {
        return denom*sinh(phi);
    }

};

class twoSheetHyperboloid : public z2{
    private:

    int a = 2;
    int b = 4;
    int c = 3;

    int n = 60;
    double t0 = 0, t1 = 2*3.1428;
    double dT1 = (t1 - t0) / (n - 1);

    double t2 = 4;
    double dT2 = (t2 - t0) / (n - 1);

    public:

    double x_calc(int denom, double theta, double  phi) {
        return denom*sinh(phi)*cos(theta);
    }

    double y_calc(int denom, double theta, double  phi) {
        return denom*sin(theta)*sinh(phi);
    }

    double z_calc(int denom, double  phi) {
        return denom*cosh(phi);
    }
    
    void plot(){
        PyObject* ax = plt::chart(111);
        plt::Clear3DChart(ax);

        

        for (int i = 0; i < n; ++i) {
            tx.clear();
            ty.clear();
            tz.clear();
            tz2.clear();
            theta = t0 + i * dT1;

            for (int j = 0; j < n; ++j) {
                phi = t0 + j * dT2;
                tx.push_back(x_calc(a,theta,phi));
                ty.push_back(y_calc(b,theta,phi));
                tz.push_back(z_calc(c,phi));
                tz2.push_back(z_calc(-c,phi));
            }
            x.push_back(tx);
            y.push_back(ty);
            z.push_back(tz);
            z2.push_back(tz2);
        }
        
        plt::surface3D(ax, x, y, z, "blue", 0.7);
        plt::surface3D(ax, x, y, z2, "blue", 0.7);
        plt::show();
    }
};


int main()
{
    
    // sphere surface;

    // surface.plot(0, M_PI, 0, 2 * M_PI, 1, 1, 1);
    
    // ellipsoid ell;
    // ell.plot(0, M_PI, 0, 2 * M_PI, 2, 1, 1); // Plot an ellipsoid

     
cone c1;
c1.plot(0, 2 * M_PI, -4, 4, 2, 4, 3);

    return 0;
}
