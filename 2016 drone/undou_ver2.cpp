#include <stdio.h>
#include <math.h>
#include <iostream>
#include <fstream>
using namespace std;

double LM(double l,double T1,double T2);
double MM(double l,double T3,double T4);
double NM(double l,double T1,double T2,double T3,double T4);
double Zp(double T1,double T2,double T3,double T4);
double Udot(double q,double W0,double th,double th0,double fx,double g,double m)
double Vdot(double p,double r,double U0,double W0,double phi,double th0,double fy,double g,double m)
double Wdot(double q,double W0,double th,double th0,double fz,double g,double m)
double pdot(double Ixx,doubel Izz,double Ixz,double drtL,double drtN)
double qdot(double Iyy,double drtM)
double rdot(double Ixx,doubel Izz,double Ixz,double drtL,double drtN)
double phidot(double p,doubel r,double th0)
double thdot(double q)
double yawdot(double r,double th0) 



double LM(double l,double T1,double T2)
{
	return l * (T2 - T1);
}

double MM(double l,double T3,double T4)
{
	return l * (T3 - T4);
}

double NM(double k,double T1,double T2,double T3,double T4)
{
	return k * (T1 + T2 - (T3 + T4));
}

double Zp(double T1,double T2,double T3,double T4)
{
	return T1+T2+T3+T4;
}

double Udot(double q,double W0,double th,double th0,double fx,double g,double m)
{
	return -th * g * cos(th0) (fx / m) - q * W0;
}

double Vdot(double p,double r,double U0,double W0,double phi,double th0,double fy,double g,double m)
{
	return phi * g * cos(th0) + (fy / m) - r * U0 + p * W0;
}

double Wdot(double q,double W0,double th,double th0,double fz,double g,double m)
{ 
	return -th * g * cos(th0) + (fz / m) + q * W0;
}

double pdot(double Ixx,doubel Izz,double Ixz,double drtL,double drtN)
{
	return ((drtL * Izz) + (drtN * Ixz)) / ((Ixx * Izz) - (Ixz * Ixz));
}

double qdot(double Iyy,double drtM)
{
	return drtM / Iyy;
}

double rdot(double Ixx,doubel Izz,double Ixz,double drtL,double drtN)
{
	return ((drtL * Ixz) + (drtN * Ixx)) / ((Ixx * Izz) - (Ixz * Ixz));
}

double phidot(double p,doubel r,double th0)
{
	return p + r * tan(th0);
}

double thdot(double q)
{
	return q;
}

double yawdot(double r,double th0) 
{
	return r / cos(th0);
}

double Fu(double q,double W0,double th,double th0,double fx,double g,double m)
{
	


}

double Vdot(double p,double r,double U0,double W0,double phi,double th0,double fy,double g,double m)
double Wdot(double q,double W0,double th,double th0,double fz,double g,double m)
double pdot(double Ixx,doubel Izz,double Ixz,double drtL,double drtN)
double qdot(double Iyy,double drtM)
double rdot(double Ixx,doubel Izz,double Ixz,double drtL,double drtN)
double phidot(double p,doubel r,double th0)
double thdot(double q)
double yawdot(double r,double th0) 
