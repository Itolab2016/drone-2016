#include <stdio.h>
#include <math.h>
#include <iostream>
#include <fstream>
using namespace std;

//b = 33.3 * a -50.0

double delL(double l,double delT1,double delT2);
double delM(double l,double delT3,double delT4);
double delN(double kt1,double kt2,double kt3,double kt4,double pwm1,double pwm2,double pwm3,double pwm4);
double Zp(double delT1,double delT2,double delT3,double delT4,double m);
double Udot(double th);
double Vdot(double phi);
double Wdot(double Z,double m);
double Pdot(double a,double b,double k1,double k2,double kt1,double kt2,double kt3,double kt4,double pwm1,double pwm2,double pwm3,double pwm4,double l);
double Qdot(double M,double Iy);
double Rdot(double c,double d,double k1,double k2,double kt1,double kt2,double kt3,double kt4,double pwm1,double pwm2,double pwm3,double pwm4,double l);
double Phidot(double p);
double Thdot(double q);
double Yawdot(double r);
double Fu(double th,double U,double h);
double Fv(double phi,double V,double h);
double Fw(double Z,double m,double W,double h);
double Fp(double a,double b,double k1,double k2,double kt1,double kt2,double kt3,double kt4,double pwm1,double pwm2,double pwm3,double pwm4,double l,double P,double h);
double Fq(double M,double Iy,double Q,double h);
double Fr(double c,double d,double k1,double k2,double kt1,double kt2,double kt3,double kt4,double pwm1,double pwm2,double pwm3,double pwm4,double l,double R,double h);
double Fphi(double p,double Phi,double h);
double Fth(double q,double Th,double h);
double Fyaw(double r,double Yaw,double h);

main()
{
	double X,Y,Z,m,p,q,r,th,phi,yaw,W,V,U,P,Q,R,Phi,Th,Yaw,Ix,Iy,Iz,Ixz,Izx,L,M,N,h,l,k,af,ab,ar,al,bf,bb,br,bl;
	double delT1,delT2,delT3,delT4;
	double k1,k2,k3,k4,kt1,kt2,kt3,kt4,pwm1,pwm2,pwm3,pwm4;
	double a,b,c,d;
	double udot,vdot,wdot,pdot,qdot,rdot,phidot,thdot,yawdot;
	//double k1[9],k2[9],k3[9],k4[9];

	X   = 0.0;	//Xの初期値
	Y   = 0.0;	//Yの初期値
	Z   = 0.0;	//Zの初期値
	m   = 0.737;	//機体の重量
	p   = 0.0;	//x方向の回転角の初期位置
	q   = 0.0;	//y方向の回転角の初期位置
	r   = 0.0;	//z方向の回転角の初期位置
	W   = 0.0;	//x方向の速度の初期位置
	V   = 0.0;	//Y方向の速度の初期位置
	U   = 0.0;	//z方向の速度の初期位置
	th  = 0.0;	//thの初期角度
	phi = 0.0;	//phiの初期角度
	yaw = 0.0;	//yawの初期角度
	Ix  = 0.010694;	//Ixの初期値
	Iy  = 0.010895;	//Iyの初期値
	Iz  = 0.012944;	//Izの初期値
	Ixz = 0.0000361921;	//Ixzの初期値
	Izx = 0.0000361921;	//Izxの初期値
	L   = 0.0;	//Lの初期値
	M   = 0.0;	//Mの初期値
	N   = 0.0;	//Nの初期値
	h   = 0.01;	//刻み幅
	l   = 0.136;       //腕の長さ
	k1  = 4.0857;
	k2  = 3.5421;
	k3  = 3.5284;
	k4  = 3.7375;
	kt1 = 2;
	kt2 = 2;
	kt3 = 2;
	kt4 = 2;
	pwm1   = 1.5;
	pwm2   = 1.5;
	pwm3   = 1.5;
	pwm4   = 1.5;
	delT1  = -k1 * pwm1; 
	delT2  = -k2 * pwm2;
	delT3  = -k3 * pwm3;
	delT4  = -k4 * pwm4;
	a   = 77.259;
	b   = 0.26152;
	c   = 0.26152;
	d   = 93.51;
	k   = 20.0;

	char filename[] = "flightlog.txt";
	char outstr[255];
	std::ofstream fs(filename);

	for (short time=0;time<=1000;time++){
		printf("U=%2.3f V=%2.3f W=%2.3f p=%2.3f q=%2.3f r=%2.3f phi=%2.3f th=%2.3f yaw=%2.3f\n"
			,U,V,W,p,q,r,phi,th,yaw);

		sprintf(outstr,"%f %f %f %f %f %f %f %f %f %d",
			U,V,W,
			p,q,r,
			phi,th,yaw,
			time);
		fs << outstr << endl;		

		L = delL(l,delT1,delT2);
		M = delM(l,delT3,delT4);
		N = delN(kt1,kt2,kt3,kt4,pwm1,pwm2,pwm3,pwm4);
		Z = Zp(delT1,delT2,delT3,delT4,m);

		//udot   = Udot(&X,&m,&q,&r,&W,&V,&th);
		//cout<<"Udot="<<udot<<" ";
		//vdot   = Vdot(&Y,&m,&r,&p,&U,&W,&th,&phi);
		//cout<<"Vdot="<<pdot<<" ";
		//wdot   = Wdot(&Z,&m,&p,&q,&V,&U,&th,&phi);
		//cout<<"Wdot="<<wdot<<" ";
		//pdot   = Pdot(&Ix,&Iy,&Iz,&Ixz,&Izx,&p,&q,&r,&L,&N);
		//cout<<"pdot="<<pdot<<" ";
		//qdot   = Qdot(&Ix,&Iy,&Iz,&Ixz,&p,&r,&M);
		//cout<<"qdot="<<qdot<<" ";
		//rdot   = Rdot(&Ix,&Iy,&Iz,&Ixz,&Izx,&p,&q,&r,&L,&N);
		//cout<<"rdot="<<Rdot(&Ix,&Iy,&Iz,&Ixz,&Izx,&p,&q,&r,&L,&N)<<" ";
		//phidot = Phidot(&p,&q,&r,&phi,&th);
		//cout<<"Phidot="<<phidot<<" ";
		//thdot  = Thdot(&q,&r,&phi);
		//cout<<"Thdot="<<Thdot(&q,&r,&phi)<<" ";
		//yawdot = Yawdot(&q,&r,&phi,&th);
		//cout<<"yawdot="<<yawdot<<" ";
		U   = Fu(th,U,h);
		//cout<<"U="<<U<<" ";
		V   = Fv(phi,V,h);
		//cout<<"V="<<V<<" ";
		W   = Fw(Z,m,W,h);
		//cout<<"W="<<W<<" ";
		p   = Fp(a,b,k1,k2,kt1,kt2,kt3,kt4,pwm1,pwm2,pwm3,pwm4,l,P,h);
		//cout<<"p="<<p<<" ";
		q   = Fq(M,Iy,Q,h);
		//cout<<"q="<<q<<" ";
		r   = Fr(c,d,k1,k2,kt1,kt2,kt3,kt4,pwm1,pwm2,pwm3,pwm4,l,R,h);
		//cout<<"r="<<r<<" ";
		phi = Fphi(p,Phi,h);
		//cout<<"phi="<<phi<<" ";
		th  = Fth(q,Th,h);
		//cout<<"th="<<th<<" ";
		yaw = Fyaw(r,Yaw,h);
		//cout<<"yaw="<<yaw<<" ";



	}

	fs.close();

} 
	 
double delL(double l,double delT1,double delT2)
{
	return l * (delT2 - delT1);
}

double delM(double l,double delT3,double delT4)
{
	return l * (delT3 - delT4);
}

double delN(double kt1,double kt2,double kt3,double kt4,double pwm1,double pwm2,double pwm3,double pwm4)
{
	return (kt1 * pwm1 + kt2 * pwm2) - (kt3 * pwm1 + kt4 * pwm4);
}

double Zp(double delT1,double delT2,double delT3,double delT4,double m)
{
	return (4 * delT1 / m)+(4 * delT2 / m)+(4 * delT3 / m)+(4 * delT4 / m);
}

double Udot(double th)
{

	return -th * 9.80665;
}

double Vdot(double phi)
{

	return phi * 9.80665;
}

double Wdot(double Z,double m)
{
	return Z/m;
}

double Pdot(double a,double b,double k1,double k2,double kt1,double kt2,double kt3,double kt4,double pwm1,double pwm2,double pwm3,double pwm4,double l)
{
	return (-a*l*k1+b*k1)*pwm1 + (a*l*k2+b*kt2)*pwm2 - b*kt3*pwm3 - b*kt4*pwm4;
}

double Qdot(double M,double Iy)
{ 
	return M / Iy;
}

double Rdot(double c,double d,double k1,double k2,double kt1,double kt2,double kt3,double kt4,double pwm1,double pwm2,double pwm3,double pwm4,double l)
{
	return (-c*l*k1+d*kt1)*pwm1 + (c*l*k2+d*kt2)*pwm2 -d*kt3*pwm3 - d*kt4*pwm4;
}

double Phidot(double p)
{
	return p;
}

double Thdot(double q)
{

	return q;
}

double Yawdot(double r)
{

	return r;
}

double  Fu(double th,double U,double h)
{
	double k[4];
	double Futh;
	k[0] = Udot(th);
	for(int i=0;i<3;i++)
	{
		Futh=((th)+h*k[i]*0.5);		
		
		k[i+1] = Udot(Futh);
	}

		Futh=((th)+h*k[2]);		
		
		k[3] = Udot(Futh);

	U = U + (h * (k[0] + 2*k[1] + 2*k[2] + k[3])) / 6;
	return U;
}


double  Fv(double phi,double V,double h)
{
	double k[4];
	double Fvphi;
	k[0] = Vdot(phi);
	for(int i=0;i<3;i++)
	{
		Fvphi=(phi+h*k[i]*0.5);
		
		k[i+1] = Vdot(Fvphi);
	}

		Fvphi=(phi+h*k[2]);
		
		k[3] = Vdot(Fvphi);

	V = V + (h * (k[0] + 2*k[1] + 2*k[2] + k[3])) / 6;
	return V;
}


double  Fw(double Z,double m,double W,double h)
{
	double k[4];
	double FwZ;
	k[0] = Wdot(Z,m);
	for(int i=0;i<3;i++)
	{
		FwZ=(Z+h*k[i]*0.5);
		
		k[i+1] = Wdot(FwZ,m);
	}

		FwZ=(Z+h*k[2]);
		
		k[3] = Wdot(FwZ,m);

	W = W + (h * (k[0] + 2*k[1] + 2*k[2] + k[3])) / 6;
	return W;
}

double Fp(double a,double b,double k1,double k2,double kt1,double kt2,double kt3,double kt4,double pwm1,double pwm2,double pwm3,double pwm4,double l,double P,double h)
{
	double Fpa,Fpb,Fpk1,Fpk2,Fpkt1,Fpkt2,Fpkt3,Fpkt4,Fppwm1,Fppwm2,Fppwm3,Fppwm4,Fpl;
	double k[4];
	k[0] = Pdot(a,b,k1,k2,kt1,kt2,kt3,kt4,pwm1,pwm2,pwm3,pwm4,l);

	for(int i=0;i<3;i++)
	{
	 	Fpa = a+h*k[i]*0.5;
		Fpb = b+h*k[i]*0.5;
		Fpk1 = k1+h*k[i]*0.5;
		Fpk2 = k2+h*k[i]*0.5;
		Fpkt1 = kt1+h*k[i]*0.5;
		Fpkt2 = kt2+h*k[i]*0.5;
		Fpkt3 = kt3+h*k[i]*0.5;
		Fpkt4 = kt4+h*k[i]*0.5;
		Fppwm1 = pwm1+h*k[i]*0.5;
		Fppwm2 = pwm2+h*k[i]*0.5;
		Fppwm3 = pwm3+h*k[i]*0.5;
		Fppwm4 = pwm4+h*k[i]*0.5;
		Fpl = l+h*k[i]*0.5;	
		
		k[i+1] = Pdot(Fpa,Fpb,Fpk1,Fpk2,Fpkt1,Fpkt2,Fpkt3,Fpkt4,Fppwm1,Fppwm2,Fppwm3,Fppwm4,Fpl);
	}

		Fpa = a+h*k[2]*0.5;
		Fpb = b+h*k[2]*0.5;
		Fpk1 = k1+h*k[2]*0.5;
		Fpk2 = k2+h*k[2]*0.5;
		Fpkt1 = kt1+h*k[2]*0.5;
		Fpkt2 = kt2+h*k[2]*0.5;
		Fpkt3 = kt3+h*k[2]*0.5;
		Fpkt4 = kt4+h*k[2]*0.5;
		Fppwm1 = pwm1+h*k[2]*0.5;
		Fppwm2 = pwm2+h*k[2]*0.5;
		Fppwm3 = pwm3+h*k[2]*0.5;
		Fppwm4 = pwm4+h*k[2]*0.5;
		Fpl = l+h*k[2]*0.5;

		k[3] = Pdot(Fpa,Fpb,Fpk1,Fpk2,Fpkt1,Fpkt2,Fpkt3,Fpkt4,Fppwm1,Fppwm2,Fppwm3,Fppwm4,Fpl);

	P = P + (h * (k[0] + 2*k[1] + 2*k[2] +k[3])) / 6;
	return P;
}

double Fq(double M,double Iy,double Q,double h)
{
	double k[4];
	double FqM,FqIy;
	k[0] = Qdot(M,Iy);

	for(int i=0;i<3;i++)
	{
		FqM = M+h*k[i]*0.5;
		FqIy = Iy+h*k[i]*0.5;
		
		k[i+1] = Qdot(FqM,FqIy);
	}

		FqM = M+h*k[2];
		FqIy = Iy+h*k[2];
		
		k[3] = Qdot(FqM,FqIy);

	Q = Q + (h * (k[0] + 2*k[1] + 2*k[2] + k[3])) / 6;
	return Q;

}

double Fr(double c,double d,double k1,double k2,double kt1,double kt2,double kt3,double kt4,double pwm1,double pwm2,double pwm3,double pwm4,double l,double R,double h)
{
	double k[4];
	double Frc,Frd,Frk1,Frk2,Frk3,Frk4,Frkt1,Frkt2,Frkt3,Frkt4,Frpwm1,Frpwm2,Frpwm3,Frpwm4,Frl;
	k[0] = Rdot(c,d,k1,k2,kt1,kt2,kt3,kt4,pwm1,pwm2,pwm3,pwm4,l);
	for(int i=0;i<3;i++)
	{
		Frc = c+h*k[i]*0.5;
		Frd = d+h*k[i]*0.5;
		Frk1 = k1+h*k[i]*0.5;
		Frk2 = k2+h*k[i]*0.5;
		Frkt1 = kt1+h*k[i]*0.5;
		Frkt2 = kt2+h*k[i]*0.5;
		Frkt3 = kt3+h*k[i]*0.5;
		Frkt4 = kt4+h*k[i]*0.5;
		Frpwm1 = pwm1+h*k[i]*0.5;
		Frpwm2 = pwm2+h*k[i]*0.5;
		Frpwm3 = pwm3+h*k[i]*0.5;
		Frpwm4 = pwm4+h*k[i]*0.5;
		Frl = l+h*k[i]*0.5;
	
		k[i+1] = Rdot(Frc,Frd,Frk1,Frk2,Frkt1,Frkt2,Frkt3,Frkt4,Frpwm1,Frpwm2,Frpwm3,Frpwm4,Frl);
	}
	
		Frc = c+h*k[2]*0.5;
		Frd = d+h*k[2]*0.5;
		Frk1 = k1+h*k[2]*0.5;
		Frk2 = k2+h*k[2]*0.5;
		Frkt1 = kt1+h*k[2]*0.5;
		Frkt2 = kt2+h*k[2]*0.5;
		Frkt3 = kt3+h*k[2]*0.5;
		Frkt4 = kt4+h*k[2]*0.5;
		Frpwm1 = pwm1+h*k[2]*0.5;
		Frpwm2 = pwm2+h*k[2]*0.5;
		Frpwm3 = pwm3+h*k[2]*0.5;
		Frpwm4 = pwm4+h*k[2]*0.5;
		Frl = l+h*k[2]*0.5;

		k[3] = Rdot(Frc,Frd,Frk1,Frk2,Frkt1,Frkt2,Frkt3,Frkt4,Frpwm1,Frpwm2,Frpwm3,Frpwm4,Frl);

	R = R + (h * (k[0] + 2*k[1] + 2*k[2] +k[3])) / 6;
	return R;
}

double Fphi(double p,double Phi,double h)
{
	double k[4];
	double Fphi_p;
	k[0] = Phidot(p);
	for(int i=0;i<3;i++)
	{
		Fphi_p = p+h*k[i]*0.5;

		k[i+1] = Phidot(Fphi_p);	
	}

		Fphi_p = p+h*k[2];

		k[3] = Phidot(Fphi_p);

	Phi = Phi + (h * (k[0] + 2*k[1] + 2*k[2] + k[3])) / 6;
	return Phi;
}


double Fth(double q,double Th,double h)
{
	double k[4];
	double Fth_q;
	k[0] = Thdot(q);
	for(int i=0;i<3;i++)
	{
		Fth_q = q+h*k[i]*0.5;
		k[i+1] = Thdot(Fth_q);
	}
	
		Fth_q = q+h*k[2];
		k[3] = Thdot(Fth_q);
	Th = Th + (h * (k[0] + 2*k[1] + 2*k[2] + k[3])) / 6;
	return Th;
}

double Fyaw(double r,double Yaw,double h)
{
	double k[4];
	double Fyaw_r;
	k[0] = Yawdot(r);
	for(int i=0;i<3;i++)
	{
		Fyaw_r = r+h*k[i]*0.5;
		
		k[i+1] = Yawdot(Fyaw_r);
	}

		Fyaw_r = r+h*k[2];

		k[3] = Yawdot(Fyaw_r);

	Yaw = Yaw + (h * (k[0] + 2*k[1] + 2*k[2] + k[3])) / 6;
	return Yaw;
}


