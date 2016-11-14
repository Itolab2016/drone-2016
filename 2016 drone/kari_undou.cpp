#include <stdio.h>
#include <math.h>
using namespace std;

//double Fu(double X,double m,double q,double r,double W,double V,double th);
//double Fv(double Y,double m,double r,double p,double U,double W,double th,double phi);
//double Fw(double W,double m,double p,double q,double V,double U,double th,double phi);
//double Fp(double Ix,double Iy,double Iz,double Ixz,double Izx,double p,double q,double r,double L,double N);
//double Fq(double Ix,double Iy,double Iz,double Ixz,double p,double r,double M);
//double Fr(double Ix,double Iy,double Iz,double Ixz,double Izx,double p,double q,double r,double L,double N);
//double Fphi(double p,double q,double r,double phi,double th);
//double Fth(double q,double r,double phi);
//double Fyaw(double q,double r,double phi,double th);


	double X,Y,Z,m,p,q,r,W,V,U,th,phi,yaw,Ix,Iy,Iz,Ixz,Izx,L,M,N,h,LM,l,MM,NM;
	double T1,T2,T3,T4;
	double udot,vdot,wdot,pdot,qdot,rdot,phidot,thdot,yawdot;
	//double k1[9],k2[9],k3[9],k4[9];

	X   = 0.0;	//Xの初期値
	Y   = 0.0;	//Yの初期値
	Z   = 0.0;	//Zの初期値
	m   = 617.0;    //機体の重量
	p   = 0.0;	//x方向の回転角の初期位置
	q   = 0.0;	//y方向の回転角の初期位置
	r   = 0.0;	//z方向の回転角の初期位置
	W   = 0.0;	//x方向の速度の初期位置
	V   = 0.0;	//Y方向の速度の初期位置
	U   = 0.0;	//z方向の速度の初期位置
	th = 0.0;	//thの初期角度
	phi =0.0;	//phiの初期角度
	yaw = 0.0;	//yawの初期角度
	Ix  = 0.0;	//Ixの初期値
	Iy = 0.0;	//Iyの初期値
	Iz = 0.0;	//Izの初期値
	Ixz = 0.0;	//Ixzの初期値
	Izx = 0.0;	//Izxの初期値
	L   = 0.0;	//Lの初期値
	M   = 0.0;	//Mの初期値
	N   = 0.0;	//Nの初期値
	h   = 0.01;	//刻み幅


main(){
	
	while (1){
		printf("%5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f"
			,U,V,W,p,q,r,phi,th,yaw);

		LM = LM(l,T1,T2);
		MM = MM(l,T1,T2);
		NM = NM(l,T1,T2,T3,T4);
		
		udot   = Udot(X,m,q,r,W,V,th);
		vdot   = Vdot(Y,m,r,p,U,W,th,phi);
		wdot   = Wdot(Z,m,p,q,V,U,th,phi);
		pdot   = Pdot(Ix,Iy,Iz,Ixz,Izx,p,q,r,L,N);
		qdot   = Qdot(Ix,Iy,Iz,Ixz,p,r,M);
		rdot   = Rdot(Ix,Iy,Iz,Ixz,Izx,p,q,r,L,N);
		phidot = Phidot(p,q,r,phi,th);
		thdot  = Thdot(q,r,phi);
		yawdot = Yawdor(q,r,phi,th);

		U   = Fu(X,m,p,r,W,V,th);
		V   = Fv(Y,m,p,r,W,U,th,phi);
		W   = Fw(Z,m,p,q,V,U,th,phi);
		p   = Fp(Ix,Iy,Iz,Ixz,Izx,p,q,r,L,N);
		q   = Fq(Ix,Iy,Iz,Ixz,p,r,M);
		r   = Fr(Ix,Iy,Iz,Ixz,Izx,p,q,r,L,N);
		phi = Fphi(p,q,r,th,phi);
		th  = Fth(q,r,phi);
		yaw = Fyaw(q,r,th,phi);
	}
} 
	 
double LM(double *l,double *T1,double *T2)
{
	return (*l) * ((*T1) - (*T2));
}

double MM(double *l,double *T3,double *T4)
{
	return(*l) * ((*T3) - (*T4));
}

double NM(double *l,double *T1,double *T2,double *T3,double *T4)
{
	return (*l) * ((*T1) + (*T2) - ((*T3) + (+T4)));
}

double Udot(double *X,double *m,double *q,double *r,double *W,double *V,double *th)
{
	double udot;	
	udot = (*X)/(*m) - 9.80665*sin(*th) - (*q)*(*W) + (*r)*(*V);
	return udot;
}

double Vdot(double *Y,double *m,double *r,double *p,double *U,double *W,double *th,double *phi)
{
	double vdot;
	vdot = (*Y)/(*m) + 9.80665*cos(*th)*sin(*phi) - (*r)*(*U) + (*p)*(*W);
	return vdot;
}

double Wdot(double *Z,double *m,double *p,double *q,double *V,double *U,double *th,double *phi)
{
	double wdot;
	wdot = (*Z)/(*m )+ 9.80665*cos(*th)*sin(*phi) - (*p)*(*V) + (*q)*(*U);
	return wdot;
}

double Pdot(double *Ix,double *Iy,double *Iz,double *Ixz,double *Izx,double *p,double *q,double *r,double *L,double *N)
{
	double pdot,pbunbo,pbunshi;
	pbunbo = ((*Iz)*((*L) - (*q)*(*r)*((*Iz)-(*Iy)) + (*p)*(*q)*(*Ixz)) + (*Izx)*((*N) - (*p)*(*q)*((*Iy)-(*Ix)) - (*q)*(*r)*(*Ixz)));
	pbunshi = (((*Ix)*(*Iz)) - ((*Ixz)*(*Izx)));
	pdot = pbunbo / pbunshi;
	return pdot;
}

double Qdot(double *Ix,double *Iy,double *Iz,double *Ixz,double *r,double *p,double *M)
{ 
	double qdot;
	qdot = ((*M) - (*r)*(*p)*((*Ix)-(*Iz)) - (*Ixz)(((*r)*(*r))-((*p)*(*p)))) / Iy;
	return qdot;
}

double Rdot(double *Ix,double *Iy,double *Iz,double *Ixz,double *p,double *q,double *r,double *N)
{
	double rdot,rbunbo,rbunshi;
	rbunbo = ((*Ixz)*((*L)- (*q)*(*r)*((*Iz)-(*Iy)) + (*p)*(*q)*(*Ixz)) + (*Ix)*((*N) - (*p)*(*q)*((*Iz)-(*Ix)) - (*q)*(*r)*(*Ixz)));
	rbunshi = (((*Ix)*(*Iz)) - ((*Izx)*(*Ixz)));
	rdot = rbunbo / rbunshi;
	return rdot;
}

double Phidot(double *p,double *q,double *r,double *th,double *phi)
{
	double phidot;
	phidot = (*p) + (*q)*sin(*phi)*tan(*th) + (*r)*cos(*phi)*tan(*th);
	return phidot;
}

double Thdot(double *q,double *r,double *phi)
{
	double thdot;
	thdot = (*q)*cos(*phi) - (*r)*sin(*phi);
	return thdot;
}

double Yawdot(double *q,double *r,double *phi,double *th)
{
	double yawdot;
	yawdot = (*q)*sin(*phi)*(1/cos(*th)) + (*r)*cos(*phi)*(1/cos(*th));
	return yawdot;
}

double  Fu(double *X,double *m,double *q,double *r,double *W,double *V,double *th)
{
	double k1,k2,k3,k4;
	k1 = udot(X,m,q,r,W,V,th);
	k2 = udot((X)+(h)*k1*0.5,(*m)+(*h)*k1*0.5,(*q)+(*h)*k1*0.5,(*r)+(*h)*k1*0.5,(*W)+(*h)*k1*0.5,(*V)+(*h)*k1*0.5,(*th)+(*h)*k1*0.5);
	k3 = udot((*X)+(*h)*k2*0.5,(*m)+(*h)*k2*0.5,(*q)+(*h)*k2*0.5,(*r)+(*h)*k2*0.5,(*W)+(*h)*k2*0.5,(*V)+(*h)*k2*0.5,(*th)+(*h)*k2*0.5);
	k4 = udot((*X)+(*h)*k3*0.5,(*m)+(*h)*k3*0.5,(*q)+(*h)*k3*0.5,(*r)+(*h)*k3*0.5,(*W)+(*h)*k3*0.5,(*V)+(*h)*k3*0.5,(*th)+(*h)*k3*0.5);
	u = u + ((*h) * (k1 + 2*k2 + 2*k3 + k4)) / 6;
	return U;
}

double  Fv(double *Y,double *m,double *p,double *r,double *W,double *U,double *th,double *phi)
{
	double k1,k2,k3,k4;
	k1 = vdot(*Y,*m,*p,*r,*W,*U,*th,*phi);
	k2 = vdot((*Y)+(*h)*k1*0.5,(*m)+(*h)*k1*0.5,(*p)+(*h)*k1*0.5,(*r)+(*h)*k1*0.5,(*W)+(*h)*k1*0.5,(*U)+(*h)*k1*0.5,(*th)+(*h)*k1*0.5,(*phi)+(*h)*k1*0.5);
	k3 = udot((*Y)+(*h)*k2*0.5,(*m)+(*h)*k2*0.5,(*p)+(*h)*k2*0.5,(*r)+(*h)*k2*0.5,(*W)+(*h)*k2*0.5,(*U)+(*h)*k2*0.5,(*th)+(*h)*k2*0.5,(*phi)+(*h)*k2*0.5);
	k4 = udot((*Y)+(*h)*k3*0.5,(*m)+(*h)*k3*0.5,(*p)+(*h)*k3*0.5,(*r)+(*h)*k3*0.5,(*W)+(*h)*k3*0.5,(*U)+(*h)*k3*0.5,(*th)+(*h)*k3*0.5,(*phi)+(*h)*k3*0.5);
	v = v + ((*h) * (k1 + 2*k2 + 2*k3 + k4)) / 6;
	return V;
}

double  Fw(double *Z,double *m,double *p,double *q,double *V,double *U,double *th,double *phi)
{
	double k1,k2,k3,k4;
	k1 = vdot(*Z,*m,*p,*q,*V,*U,*th,*phi);
	k2 = vdot((*Z)+(*h)*k1*0.5,(*m)+(*h)*k1*0.5,(*p)+(*h)*k1*0.5,(*q)+(*h)*k1*0.5,(*V)+(*h)*k1*0.5,(*U)+(*h)*k1*0.5,(*th)+(*h)*k1*0.5,(*phi)+(*h)*k1*0.5);
	k3 = udot((*Z)+(*h)*k2*0.5,(*m)+(*h)*k2*0.5,(*p)+(*h)*k2*0.5,(*q)+(*h)*k2*0.5,(*V)+(*h)*k2*0.5,(*U)+(*h)*k2*0.5,(*th)+(*h)*k2*0.5,(*phi)+(*h)*k2*0.5);
	k4 = udot((*Z)+(*h)*k3*0.5,(*m)+(*h)*k3*0.5,(*p)+(*h)*k3*0.5,(*q)+(*h)*k3*0.5,(*V)+(*h)*k3*0.5,(*U)+(*h)*k3*0.5,(*th)+(*h)*k3*0.5,(*phi)+(*h)*k3*0.5);
	v = v + ((*h) * (k1 + 2*k2 + 2*k3 + k4)) / 6;
	return W;
}

double Fp(double *Ix,double *Iy,double *Iz,double *Ixz,double *Izx,double *p,double *q,double *r,double *L,double *N)
{
	double k1,k2,k3,k4;
	k1 = pdot(*Ix,*Iy,*Iz,*Ixz,*Izx,*p,*q,*r,*L,*N);
	k2 = pdot((*Ix)+(*h)*k1*0.5,(*Iy)+(*h)*k1*0.5,(*Iz)+(*h)*k1*0.5,(*Ixz)+(*h)*k1*0.5,(*Izx)+(*h)*k1*0.5,(*p)+(*h)*k1*0.5,(*q)+(*h)*k1*0.5,(*r)+(*h)*k1*0.5,(*L)+(*h)*k1*0.5,(*N)+(*h)*k1*0.5);	
	k3 = pdot((*Ix)+(*h)*k2*0.5,(*Iy)+(*h)*k2*0.5,(*Iz)+(*h)*k2*0.5,(*Ixz)+(*h)*k2*0.5,(*Izx)+(*h)*k2*0.5,(*p)+(*h)*k2*0.5,(*q)+(*h)*k2*0.5,(*r)+(*h)*k2*0.5,(*L)+(*h)*k2*0.5,(*N)+(*h)*k2*0.5);
	k4 = pdot((*Ix)+(*h)*k3*0.5,(*Iy)+(*h)*k3*0.5,(*Iz)+(*h)*k3*0.5,(*Ixz)+(*h)*k3*0.5,(*Izx)+(*h)*k3*0.5,(*p)+(*h)*k3*0.5,(*q)+(*h)*k3*0.5,(*r)+(*h)*k3*0.5,(*L)+(*h)*k3*0.5,(*N)+(*h)*k3*0.5);
	p = p + ((*h) * (k1 + 2*k2 + 2*k3 +k4)) / 6;
	return p;
}

double Fq(double *Ix,double *Iy,double *Iz,double *Ixz,double *p,double *r,double *M)
{
	double k1,k2,k3,k4;
	k1 = qdot(*Ix,*Iy,*Iz,*Ixz,*p,*r,*M);
	k2 = qdpt((*Ix)+(*h)*k1*0.5,(*Iy)+(*h)*k1*0.5,(*Iz)+(*h)*k1*0.5,(*Ixz)+(*h)*k1*0.5,(*p)+(*h)*k1*0.5,(*r)+(*h)*k1*0.5,(*M)+(*h)*k1*0.5);
	k3 = qdpt((*Ix)+(*h)*k2*0.5,(*Iy)+(*h)*k2*0.5,(*Iz)+(*h)*k2*0.5,(*Ixz)+(*h)*k2*0.5,(*p)+(*h)*k2*0.5,(*r)+(*h)*k2*0.5,(*M)+(*h)*k2*0.5);
	k4 = qdpt((*Ix)+(*h)*k3*0.5,(*Iy)+(*h)*k3*0.5,(*Iz)+(*h)*k3*0.5,(*Ixz)+(*h)*k3*0.5,(*p)+(*h)*k3*0.5,(*r)+(*h)*k3*0.5,(*M)+(*h)*k3*0.5);
	q = q + ((*h) * (k1 + 2*k2 + 2*k3 + k4)) / 6;
	return q;
}

double Fr(double *Ix,double *Iy,double *Iz,double *Ixz,double *Izx,double *p,double *q,double *r,double *L,double *N)
{
	double k1,k2,k3,k4;
	k1 = rdot(*Ix,*Iy,*Iz,*Ixz,*Izx,*p,*q,*r,*L,*N);
	k2 = rdot((*Ix)+(*h)*k1*0.5,(*Iy)+(*h)*k1*0.5,(*Iz)+(*h)*k1*0.5,(*Ixz)+(*h)*k1*0.5,(*Izx)+(*h)*k1*0.5,(*p)+(*h)*k1*0.5,(*q)+(*h)*k1*0.5,(*r)+(*h)*k1*0.5,(*L)+(*h)*k1*0.5,(*N)+(*h)*k1*0.5);	
	k3 = rdot((*Ix)+(*h)*k2*0.5,(*Iy)+(*h)*k2*0.5,(*Iz)+(*h)*k2*0.5,(*Ixz)+(*h)*k2*0.5,(*Izx)+(*h)*k2*0.5,(*p)+(*h)*k2*0.5,(*q)+(*h)*k2*0.5,(*r)+(*h)*k2*0.5,(*L)+(*h)*k2*0.5,(*N)+(*h)*k2*0.5);
	k4 = rdot((*Ix)+(*h)*k3*0.5,(*Iy)+(*h)*k3*0.5,(*Iz)+(*h)*k3*0.5,(*Ixz)+(*h)*k3*0.5,(*Izx)+(*h)*k3*0.5,(*p)+(*h)*k3*0.5,(*q)+(*h)*k3*0.5,(*r)+(*h)*k3*0.5,(*L)+(*h)*k3*0.5,(*N)+(*h)*k3*0.5);
	r = r + ((*h) * (k1 + 2*k2 + 2*k3 +k4)) / 6;
	return r;
}

double Fphi(double *p,double *q,double *r,double *th,double *phi)
{
	double k1,k2,k3,k4;
	k1 = phidot(*p,*q,*r,*th,*phi);
	k2 = phidot((*p)+(*h)*k1*0.5,(*q)+(*h)*k1*0.5,(*r)+(*h)*k1*0.5,(*th)+(*h)*k1*0.5,(*phi)+(*h)*k1*0.5);
	k3 = phidot((*p)+(*h)*k2*0.5,(*q)+(*h)*k2*0.5,(*r)+(*h)*k2*0.5,(*th)+(*h)*k2*0.5,(*phi)+(*h)*k2*0.5);
	k4 = phidot((*p)+(*h)*k3*0.5,(*q)+(*h)*k3*0.5,(*r)+(*h)*k3*0.5,(*th)+(*h)*k3*0.5,(*phi)+(*h)*k3*0.5);
	phi = phi + ((*h) * (k1 + 2*k2 + 2*k3 + k4)) / 6;
	return phi;
}

double Fth(double *q,double *r,double *phi)
{
	double k1,k2,k3,k4;
	k1 = thdot(*q,*r,*phi);
	k2 = thdot((*q)+(*h)*k1*0.5,(*r)+(*h)*k1*0.5,(*phi)+(*h)*k1*0.5);
	k3 = thdot((*q)+(*h)*k2*0.5,(*r)+(*h)*k2*0.5,(*phi)+(*h)*k2*0.5);
	k4 = thdot((*q)+(*h)*k3*0.5,(*r)+(*h)*k3*0.5,(*phi)+(*h)*k3*0.5);
	th = th + ((*h) * (k1 + 2*k2 + 2*k3 + k4)) / 6;
	return th;
}

double Fyaw(double *q,double *r,double *th,double *phi)
{
	double k1,k2,k3,k4;
	k1 = yawdot(*q,*r,*th,*phi);
	k2 = yawdot((*q)+(*h)*k1*0.5,(*r)+(*h)*k1*0.5,(*th)+(*h)*k1*0.5,(*phi)+(*h)*k1*0.5);
	k3 = yawdot((*q)+(*h)*k2*0.5,(*r)+(*h)*k2*0.5,(*th)+(*h)*k2*0.5,(*phi)+(*h)*k2*0.5);
	k4 = yawdot((*q)+(*h)*k3*0.5,(*r)+(*h)*k3*0.5,(*th)+(*h)*k3*0.5,(*phi)+(*h)*k3*0.5);
	yaw = yaw + ((*h) * (k1 + 2*k2 + 2*k3 + k4)) / 6;
	return yaw;
}




