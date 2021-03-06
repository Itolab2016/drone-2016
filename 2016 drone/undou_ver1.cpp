#include <stdio.h>
#include <math.h>
#include <iostream>
using namespace std;


double LM(double *l,double *T1,double *T2);
double MM(double *l,double *T3,double *T4);
double NM(double *l,double *T1,double *T2,double *T3,double *T4);
double Udot(double *X,double *m,double *q,double *r,double *W,double *V,double *th);
double Vdot(double *Y,double *m,double *r,double *p,double *U,double *W,double *th,double *phi);
double Wdot(double *Z,double *m,double *p,double *q,double *V,double *U,double *th,double *phi);
double Wdot(double *Z,double *m,double *p,double *q,double *V,double *U,double *th,double *phi);
double Pdot(double *Ix,double *Iy,double *Iz,double *Ixz,double *Izx,double *p,double *q,double *r,double *L,double *N);
double Qdot(double *Ix,double *Iy,double *Iz,double *Ixz,double *r,double *p,double *M);
double Rdot(double *Ix,double *Iy,double *Iz,double *Ixz,double *Izx,double *p,double *q,double *r,double *L,double *N);
double Phidot(double *p,double *q,double *r,double *th,double *phi);
double Thdot(double *q,double *r,double *phi);
double Yawdot(double *q,double *r,double *phi,double *th);
double Fu(double *X,double *m,double *q,double *r,double *W,double *V,double *th,double *h,double *U);
double Fv(double *Y,double *m,double *p,double *r,double *W,double *U,double *th,double *h,double *phi,double *V);
double Fw(double *Z,double *m,double *p,double *q,double *V,double *U,double *th,double *h,double *phi,double *W);
double Fp(double *Ix,double *Iy,double *Iz,double *Ixz,double *Izx,double *p,double *q,double *r,double *L,double *N,double *h,double *P);
double Fq(double *Ix,double *Iy,double *Iz,double *Ixz,double *p,double *r,double *M,double *h,double *Q);
double Fr(double *Ix,double *Iy,double *Iz,double *Ixz,double *Izx,double *p,double *q,double *r,double *L,double *N,double *h,double *R);
double Fphi(double *p,double *q,double *r,double *th,double *phi,double *h,double *Phi);
double Fth(double *q,double *r,double *phi,double *h,double *Th);
double Fyaw(double *q,double *r,double *th,double *phi,double *h,double *Yaw);
main(){
	double X,Y,Z,m,p,q,r,phi,W,V,U,P,Q,R,th,Phi,Th,Yaw,yaw,Ix,Iy,Iz,Ixz,Izx,L,M,N,h,l;
	double T1,T2,T3,T4;
	double udot,vdot,wdot,pdot,qdot,rdot,phidot,thdot,yawdot;
	//double k1[9],k2[9],k3[9],k4[9];

	X   = 0.0;	//Xの初期値
	Y   = 0.0;	//Yの初期値
	Z   = 0.0;	//Zの初期値
	m   = 737.0;	//機体の重量
	p   = 0.0;	//x方向の回転角の初期位置
	q   = 0.0;	//y方向の回転角の初期位置
	r   = 0.0;	//z方向の回転角の初期位置
	W   = 0.0;	//x方向の速度の初期位置
	V   = 0.0;	//Y方向の速度の初期位置
	U   = 0.0;	//z方向の速度の初期位置
	th  = 0.0;	//thの初期角度
	phi = 0.0;	//phiの初期角度
	yaw = 0.0;	//yawの初期角度
	Ix  = 10.389059;	//Ixの初期値
	Iy  = 10.570814;	//Iyの初期値
	Iz  = 12.683685;	//Izの初期値
	Ixz = 0.0361921;	//Ixzの初期値
	Izx = 9.0;	//Izxの初期値
	L   = 0.0;	//Lの初期値
	M   = 0.0;	//Mの初期値
	N   = 0.0;	//Nの初期値
	h   = 0.01;	//刻み幅
	l   = 136.0;       //腕の長さ
	T1  = 0.0;
	T2  = 0.0;
	T3  = 0.0;
	T4  = 0.0;


	for (short AS=0;AS<=100;AS++){
		printf("U=%2.3f V=%2.3f W=%2.3f p=%2.3f q=%2.3f r=%2.3f phi=%2.3f th=%2.3f yaw=%2.3f\n"
			,U,V,W,p,q,r,phi,th,yaw);
		cout<<endl;
		L = LM(&l,&T1,&T2);
		M = MM(&l,&T3,&T4);
		N = NM(&l,&T1,&T2,&T3,&T4);

/*
		udot   = Udot(&X,&m,&q,&r,&W,&V,&th);
		cout<<"Udot="<<udot<<" ";
		vdot   = Vdot(&Y,&m,&r,&p,&U,&W,&th,&phi);
		cout<<"Vdot="<<pdot<<" ";
		wdot   = Wdot(&Z,&m,&p,&q,&V,&U,&th,&phi);
		cout<<"Wdot="<<wdot<<" ";
		pdot   = Pdot(&Ix,&Iy,&Iz,&Ixz,&Izx,&p,&q,&r,&L,&N);
		cout<<"pdot="<<pdot<<" ";
		qdot   = Qdot(&Ix,&Iy,&Iz,&Ixz,&p,&r,&M);
		cout<<"qdot="<<qdot<<" ";
		rdot   = Rdot(&Ix,&Iy,&Iz,&Ixz,&Izx,&p,&q,&r,&L,&N);
		cout<<"rdot="<<Rdot(&Ix,&Iy,&Iz,&Ixz,&Izx,&p,&q,&r,&L,&N)<<" ";
		phidot = Phidot(&p,&q,&r,&phi,&th);
		cout<<"Phidot="<<phidot<<" ";
		thdot  = Thdot(&q,&r,&phi);
		cout<<"Thdot="<<Thdot(&q,&r,&phi)<<" ";
		yawdot = Yawdot(&q,&r,&phi,&th);
		cout<<"yawdot="<<yawdot<<" ";
*/
		U   = Fu(&X,&m,&p,&r,&W,&V,&th,&h,&U);
		//cout<<"U="<<U<<" ";
		V   = Fv(&Y,&m,&p,&r,&W,&U,&th,&h,&phi,&V);
		//cout<<"V="<<V<<" ";
		W   = Fw(&Z,&m,&p,&q,&V,&U,&th,&h,&phi,&W);
		//cout<<"W="<<W<<" ";
		p   = Fp(&Ix,&Iy,&Iz,&Ixz,&Izx,&p,&q,&r,&L,&N,&h,&P);
		//cout<<"p="<<p<<" ";
		q   = Fq(&Ix,&Iy,&Iz,&Ixz,&p,&r,&M,&h,&Q);
		//cout<<"q="<<q<<" ";
		r   = Fr(&Ix,&Iy,&Iz,&Ixz,&Izx,&p,&q,&r,&L,&N,&h,&R);
		//cout<<"r="<<r<<" ";
		phi = Fphi(&p,&q,&r,&th,&phi,&h,&Phi);
		//cout<<"phi="<<phi<<" ";
		th  = Fth(&q,&r,&phi,&h,&Th);
		//cout<<"th="<<th<<" ";
		yaw = Fyaw(&q,&r,&th,&phi,&h,&Yaw);
		//cout<<"yaw="<<yaw<<" ";

	}
} 
	 
double LM(double *l,double *T1,double *T2)
{
	return (*l) * ((*T2) - (*T1));
}

double MM(double *l,double *T3,double *T4)
{
	return(*l) * ((*T3) - (*T4));
}

double NM(double *l,double *T1,double *T2,double *T3,double *T4)
{
	return (*l) * ((*T1) + (*T2) - ((*T3) + (*T4)));
}

double Udot(double *X,double *m,double *q,double *r,double *W,double *V,double *th)
{

	return (*X)/(*m) - 9.80665*sin(*th) - (*q)*(*W) + (*r)*(*V);
}

double Vdot(double *Y,double *m,double *r,double *p,double *U,double *W,double *th,double *phi)
{

	return (*Y)/(*m) + 9.80665*cos(*th)*sin(*phi) - (*r)*(*U) + (*p)*(*W);
}

double Wdot(double *Z,double *m,double *p,double *q,double *V,double *U,double *th,double *phi)
{
	return (*Z)/(*m )+ 9.80665*cos(*th)*cos(*phi) - (*p)*(*V) + (*q)*(*U);
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

	return ((*M) - ((*r)*(*p)*((*Ix)-(*Iz))) - ((*Ixz)*((*r)*(*r) - (*p)*(*p)))) / (*Iy);
}

double Rdot(double *Ix,double *Iy,double *Iz,double *Ixz,double *Izx,double *p,double *q,double *r,double *L,double *N)
{
	double rbunbo,rbunshi;
	rbunbo = ((*Ixz)*((*L) - (*q)*(*r)*((*Iz)-(*Iy)) + (*p)*(*q)*(*Ixz)) + (*Ix)*((*N) - (*p)*(*q)*((*Iz)-(*Ix)) - (*q)*(*r)*(*Ixz)));
	rbunshi = (((*Ix)*(*Iz)) - ((*Izx)*(*Ixz)));
	return rbunbo / rbunshi;
}

double Phidot(double *p,double *q,double *r,double *th,double *phi)
{
	return (*p) + (*q)*sin(*phi)*tan(*th) + (*r)*cos(*phi)*tan(*th);
}

double Thdot(double *q,double *r,double *phi)
{

	return ((*q)*cos(*phi)) - ((*r)*sin(*phi));
}

double Yawdot(double *q,double *r,double *phi,double *th)
{

	return (*q)*sin(*phi)*(1/cos(*th)) + (*r)*cos(*phi)*(1/cos(*th));
}

double  Fu(double *X,double *m,double *q,double *r,double *W,double *V,double *th,double *h,double *U)
{
	double k[4];
	double Fux,Fum,Fuq,Fur,FuW,FuV,Futh;
	k[0] = Udot(X,m,q,r,W,V,th);
	for(int i=0;i<3;i++)
	{
		Fux=((*X)+(*h)*k[i]*0.5);
		Fum=((*m)+(*h)*k[i]*0.5);
		Fuq=((*q)+(*h)*k[i]*0.5);
		Fur=((*r)+(*h)*k[i]*0.5);
		FuW=((*W)+(*h)*k[i]*0.5);
		FuV=((*V)+(*h)*k[i]*0.5);
		Futh=((*th)+(*h)*k[i]*0.5);		
		
		k[i+1] = Udot(&Fux,&Fum,&Fuq,&Fur,&FuW,&FuV,&Futh);
	}

		Fux=((*X)+(*h)*k[2]);
		Fum=((*m)+(*h)*k[2]);
		Fuq=((*q)+(*h)*k[2]);
		Fur=((*r)+(*h)*k[2]);
		FuW=((*W)+(*h)*k[2]);
		FuV=((*V)+(*h)*k[2]);
		Futh=((*th)+(*h)*k[2]);		
		
		k[3] = Udot(&Fux,&Fum,&Fuq,&Fur,&FuW,&FuV,&Futh);

	*U = *U + ((*h) * (k[0] + 2*k[1] + 2*k[2] + k[3])) / 6;
	return *U;
}


double  Fv(double *Y,double *m,double *p,double *r,double *W,double *U,double *th,double *h,double *phi,double *V)
{
	double k[4];
	double Fvy,Fvm,Fvp,Fvr,FvW,FvU,Fvth,Fvphi;
	k[0] = Vdot(Y,m,p,r,W,U,th,phi);
	for(int i=0;i<3;i++)
	{
		Fvy=((*Y)+(*h)*k[i]*0.5);
		Fvm=((*m)+(*h)*k[i]*0.5);
		Fvp=((*p)+(*h)*k[i]*0.5);
		Fvr=((*r)+(*h)*k[i]*0.5);
		FvW=((*W)+(*h)*k[i]*0.5);
		FvU=((*U)+(*h)*k[i]*0.5);
		Fvth=((*th)+(*h)*k[i]*0.5);		
		Fvphi=(*phi)+(*h)*k[i]*0.5;

		k[i+1] = Vdot(&Fvy,&Fvm,&Fvp,&Fvr,&FvW,&FvU,&Fvth,&Fvphi);
	}

		Fvy=((*Y)+(*h)*k[2]);
		Fvm=((*m)+(*h)*k[2]);
		Fvp=((*p)+(*h)*k[2]);
		Fvr=((*r)+(*h)*k[2]);
		FvW=((*W)+(*h)*k[2]);
		FvU=((*U)+(*h)*k[2]);
		Fvth=((*th)+(*h)*k[2]);		
		Fvphi=(*phi)+(*h)*k[2];

		k[3] = Vdot(&Fvy,&Fvm,&Fvp,&Fvr,&FvW,&FvU,&Fvth,&Fvphi);

	*V = *V + ((*h) * (k[0] + 2*k[1] + 2*k[2] + k[3])) / 6;
	return *V;
}


double  Fw(double *Z,double *m,double *p,double *q,double *V,double *U,double *th,double *h,double *phi,double *W)
{
	double k[4];
	double Fwy,Fwm,Fwp,Fwq,FwV,FwU,Fwth,Fwphi;
	k[0] = Wdot(Z,m,p,q,V,U,th,phi);
	for(int i=0;i<3;i++)
	{
		Fwy=((*Z)+(*h)*k[i]*0.5);
		Fwm=((*m)+(*h)*k[i]*0.5);
		Fwp=((*p)+(*h)*k[i]*0.5);
		Fwq=((*q)+(*h)*k[i]*0.5);
		FwV=((*W)+(*h)*k[i]*0.5);
		FwU=((*U)+(*h)*k[i]*0.5);
		Fwth=((*th)+(*h)*k[i]*0.5);		
		Fwphi=(*phi)+(*h)*k[i]*0.5;

		k[i+1] = Wdot(&Fwy,&Fwm,&Fwp,&Fwq,&FwV,&FwU,&Fwth,&Fwphi);
	}

		Fwy=((*Z)+(*h)*k[2]);
		Fwm=((*m)+(*h)*k[2]);
		Fwp=((*p)+(*h)*k[2]);
		Fwq=((*q)+(*h)*k[2]);
		FwV=((*W)+(*h)*k[2]);
		FwU=((*U)+(*h)*k[2]);
		Fwth=((*th)+(*h)*k[2]);		
		Fwphi=(*phi)+(*h)*k[2];

		k[3] = Wdot(&Fwy,&Fwm,&Fwp,&Fwq,&FwV,&FwU,&Fwth,&Fwphi);

	*W = *W + ((*h) * (k[0] + 2*k[1] + 2*k[2] + k[3])) / 6;
	return *W;
}

double Fp(double *Ix,double *Iy,double *Iz,double *Ixz,double *Izx,double *p,double *q,double *r,double *L,double *N,double *h,double *P)
{
	double FpIx,FpIy,FpIz,FpIxz,FpIzx,Fpp,Fpq,Fpr,FpL,FpN;
	double k[4];
	k[0] = Pdot(Ix,Iy,Iz,Ixz,Izx,p,q,r,L,N);

	for(int i=0;i<3;i++)
	{
	 	FpIx = (*Ix)+(*h)*k[i]*0.5;
		FpIy = (*Iy)+(*h)*k[i]*0.5;
		FpIz = (*Iz)+(*h)*k[i]*0.5;
		FpIxz = (*Ixz)+(*h)*k[i]*0.5;
		FpIzx = (*Izx)+(*h)*k[i]*0.5;
		Fpp = (*p)+(*h)*k[i]*0.5;
		Fpq = (*q)+(*h)*k[i]*0.5;
		Fpr = (*r)+(*h)*k[i]*0.5;
		FpL = (*L)+(*h)*k[i]*0.5;
		FpN = (*N)+(*h)*k[i]*0.5;	
		
		k[i+1] = Pdot(&FpIx,&FpIy,&FpIz,&FpIxz,&FpIzx,&Fpp,&Fpq,&Fpr,&FpL,&FpN);
	}

		FpIx = (*Ix)+(*h)*k[2];
		FpIy = (*Iy)+(*h)*k[2];
		FpIz = (*Iz)+(*h)*k[2];
		FpIxz = (*Ixz)+(*h)*k[2];
		FpIzx = (*Izx)+(*h)*k[2];
		Fpp = (*p)+(*h)*k[2];
		Fpq = (*q)+(*h)*k[2];
		Fpr = (*r)+(*h)*k[2];
		FpL = (*L)+(*h)*k[2];
		FpN = (*N)+(*h)*k[2];

		k[3] = Pdot(&FpIx,&FpIy,&FpIz,&FpIxz,&FpIzx,&Fpp,&Fpq,&Fpr,&FpL,&FpN);

	*P = *P + ((*h) * (k[0] + 2*k[1] + 2*k[2] +k[3])) / 6;
	return *P;
}

double Fq(double *Ix,double *Iy,double *Iz,double *Ixz,double *p,double *r,double *M,double *h,double *Q)
{
	double k[4];
	double FqIx,FqIy,FqIz,FqIxz,Fqp,Fqr,FqM;
	k[0] = Qdot(Ix,Iy,Iz,Ixz,p,r,M);

	for(int i=0;i<3;i++)
	{
		FqIx = (*Ix)+(*h)*k[i]*0.5;
		FqIy = (*Iy)+(*h)*k[i]*0.5;
		FqIz = (*Iz)+(*h)*k[i]*0.5;
		FqIxz = (*Ixz)+(*h)*k[i]*0.5;
		Fqp = (*p)+(*h)*k[i]*0.5;
		Fqr = (*r)+(*h)*k[i]*0.5;
		FqM = (*M)+(*h)*k[i]*0.5;
		
		k[i+1] = Qdot(&FqIx,&FqIy,&FqIz,&FqIxz,&Fqp,&Fqr,&FqM);
	}

		FqIx = (*Ix)+(*h)*k[2];
		FqIy = (*Iy)+(*h)*k[2];
		FqIz = (*Iz)+(*h)*k[2];
		FqIxz = (*Ixz)+(*h)*k[2];
		Fqp = (*p)+(*h)*k[2];
		Fqr = (*r)+(*h)*k[2];
		FqM = (*M)+(*h)*k[2];

		k[3] = Qdot(&FqIx,&FqIy,&FqIz,&FqIxz,&Fqp,&Fqr,&FqM);

	*Q = *Q + ((*h) * (k[0] + 2*k[1] + 2*k[2] + k[3])) / 6;
	return *Q;

}

double Fr(double *Ix,double *Iy,double *Iz,double *Ixz,double *Izx,double *p,double *q,double *r,double *L,double *N,double *h,double *R)
{
	double k[4];
	double FrIx,FrIy,FrIz,FrIxz,FrIzx,Frp,Frq,Frr,FrL,FrN;
	k[0] = Rdot(Ix,Iy,Iz,Ixz,Izx,p,q,r,L,N);
	for(int i=0;i<3;i++)
	{
		FrIx = (*Ix)+(*h)*k[i]*0.5;
		FrIy = (*Iy)+(*h)*k[i]*0.5;
		FrIz = (*Iz)+(*h)*k[i]*0.5;
		FrIxz = (*Ixz)+(*h)*k[i]*0.5;
		FrIzx = (*Izx)+(*h)*k[i]*0.5;
		Frp = (*p)+(*h)*k[i]*0.5;
		Frq = (*q)+(*h)*k[i]*0.5;
		Frr = (*r)+(*h)*k[i]*0.5;
		FrL = (*L)+(*h)*k[i]*0.5;
		FrN = (*N)+(*h)*k[i]*0.5;
	
		k[i+1] = Rdot(&FrIx,&FrIy,&FrIz,&FrIxz,&FrIzx,&Frp,&Frq,&Frr,&FrL,&FrN);
	}
	
		FrIx = (*Ix)+(*h)*k[2];
		FrIy = (*Iy)+(*h)*k[2];
		FrIz = (*Iz)+(*h)*k[2];
		FrIxz = (*Ixz)+(*h)*k[2];
		FrIzx = (*Izx)+(*h)*k[2];
		Frp = (*p)+(*h)*k[2];
		Frq = (*q)+(*h)*k[2];
		Frr = (*r)+(*h)*k[2];
		FrL = (*L)+(*h)*k[2];
		FrN = (*N)+(*h)*k[2];

		k[3] = Rdot(&FrIx,&FrIy,&FrIz,&FrIxz,&FrIzx,&Frp,&Frq,&Frr,&FrL,&FrN);

	*R = *R + ((*h) * (k[0] + 2*k[1] + 2*k[2] +k[3])) / 6;
	return *R;
}

double Fphi(double *p,double *q,double *r,double *th,double *phi,double *h,double *Phi)
{
	double k[4];
	double Fphi_p,Fphi_q,Fphi_r,Fphi_th,Fphi_phi;
	k[0] = Phidot(p,q,r,th,phi);
	for(int i=0;i<3;i++)
	{
		Fphi_p = (*p)+(*h)*k[i]*0.5;
		Fphi_q = (*q)+(*h)*k[i]*0.5;
		Fphi_r = (*r)+(*h)*k[i]*0.5;
		Fphi_th = (*th)+(*h)*k[i]*0.5;
		Fphi_phi = (*phi)+(*h)*k[i]*0.5;

		k[i+1] = Phidot(&Fphi_p,&Fphi_q,&Fphi_r,&Fphi_th,&Fphi_phi);	
	}

		Fphi_p = (*p)+(*h)*k[2];
		Fphi_q = (*q)+(*h)*k[2];
		Fphi_r = (*r)+(*h)*k[2];
		Fphi_th = (*th)+(*h)*k[2];
		Fphi_phi = (*phi)+(*h)*k[2];

		k[3] = Phidot(&Fphi_p,&Fphi_q,&Fphi_r,&Fphi_th,&Fphi_phi);

	*Phi = *Phi + ((*h) * (k[0] + 2*k[1] + 2*k[2] + k[3])) / 6;
	return *Phi;
}


double Fth(double *q,double *r,double *phi,double *h,double *Th)
{
	double k[4];
	double Fth_q,Fth_r,Fth_phi;
	k[0] = Thdot(q,r,phi);
	for(int i=0;i<3;i++)
	{
		Fth_q = (*q)+(*h)*k[i]*0.5;
		Fth_r = (*r)+(*h)*k[i]*0.5;
		Fth_phi = (*phi)+(*h)*k[i]*0.5;
		k[i+1] = Thdot(&Fth_q,&Fth_r,&Fth_phi);
	}
	
		Fth_q = (*q)+(*h)*k[2];
		Fth_r = (*r)+(*h)*k[2];
		Fth_phi = (*phi)+(*h)*k[2];
		k[3] = Thdot(&Fth_q,&Fth_r,&Fth_phi);
	*Th = *Th + ((*h) * (k[0] + 2*k[1] + 2*k[2] + k[3])) / 6;
	return *Th;
}

double Fyaw(double *q,double *r,double *th,double *phi,double *h,double *Yaw)
{
	double k[4];
	double Fyaw_q,Fyaw_r,Fyaw_th,Fyaw_phi;
	k[0] = Yawdot(q,r,th,phi);
	for(int i=0;i<3;i++)
	{
		Fyaw_q = (*q)+(*h)*k[i]*0.5;
		Fyaw_r = (*r)+(*h)*k[i]*0.5;
		Fyaw_th = (*th)+(*h)*k[i]*0.5;
		Fyaw_phi = (*phi)+(*h)*k[i]*0.5;

		k[i+1] = Yawdot(&Fyaw_q,&Fyaw_r,&Fyaw_th,&Fyaw_phi);
	}

		Fyaw_q = (*q)+(*h)*k[2];
		Fyaw_r = (*r)+(*h)*k[2];
		Fyaw_th = (*th)+(*h)*k[2];
		Fyaw_phi = (*phi)+(*h)*k[2];

		k[3] = Yawdot(&Fyaw_q,&Fyaw_r,&Fyaw_th,&Fyaw_phi);

	*Yaw = *Yaw + ((*h) * (k[0] + 2*k[1] + 2*k[2] + k[3])) / 6;
	return *Yaw;
}




