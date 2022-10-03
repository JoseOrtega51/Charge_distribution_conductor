#include<iostream>
#include<stdio.h>
#include<math.h>

#define r 1.0
#define R 3.0

static const int Nu=200,Nv=200;
static const double pi=3.141592;
static double X[Nu][Nv], Y[Nu][Nv], Z[Nu][Nv];
static double X1[Nu][Nv], Y1[Nu][Nv], Z1[Nu][Nv];
static double VectorOrtogonal[Nu][Nv][3],Sactual[Nu][Nv],Santerior[Nu][Nv],NormadS[Nu][Nv];
static double VectorOrtogonal2[Nu][Nv][3],Sactual2[Nu][Nv],Santerior2[Nu][Nv],NormadS1[Nu][Nv];

///Devuelve el valor de un parámetro u o v dado su indice
double conversionIndividual (int n, int N, double inicio, double fin)
{
    double total=fin-inicio;
    double Delta=total/N;
    return inicio+Delta*n;
}

///Devuelve el valor de (x,y,z) dado el valor de (u,v)
void cCoord(double u, double v, double &x, double &y, double &z, int geometria)
{
    switch(geometria)
    {
    case 1:
        x=r*cos(u)*sin(v);
        y=r*sin(u)*sin(v);
        z=r*cos(v);
        break;
    case 2:
        x=(R+r*cos(u))*cos(v);
        y=(R+r*cos(u))*sin(v);
        z=r*sin(u);
        break;
    case 3:
        x=u*cos(v);
        y=u*sin(v);
        z=u;
        break;
    case 4:
        x=u*cos(v);
        y=u*sin(v);
        z=1;
    }
}

///Devuelve (x,y,z) dado el índice de u y v
void conversion (int nu, int nv,int geometria, double &x, double &y, double &z)
{

    switch(geometria)
    {
    case 1:
        cCoord(conversionIndividual(nu,Nu,0,2*pi),conversionIndividual(nv,Nv-1,0,pi),x,y,z,geometria);
        break;
    case 2:
        cCoord(conversionIndividual(nu,Nu,0,2*pi),conversionIndividual(nv,Nv,0,2*pi),x,y,z,geometria);
        break;
    case 3:
        cCoord(conversionIndividual(nu,Nu,0.001,1),conversionIndividual(nv,Nv,0,2*pi),x,y,z,geometria);
        break;
    case 4:
        cCoord(conversionIndividual(nu,Nu-1,0.001,1),conversionIndividual(nv,Nv,0,2*pi),x,y,z,geometria);


    }
}


double prodEscalar(double ux, double uy, double uz, double vx, double vy, double vz)
{
    return(ux*vx+uy*vy+uz*vz);
}

double norma(double x,double y, double z)
{
    return(sqrt(prodEscalar(x,y,z,x,y,z)));
}

double coseno(double ux,double uy, double uz, double vx, double vy, double vz)
{
    double cs=prodEscalar(ux,uy,uz,vx,vy,vz)/(norma(ux,uy,uz)*norma(vx,vy,vz));
    return cs;
}

void prodVectorial(double ux,double uy, double uz, double vx, double vy, double vz, double & wx, double & wy, double & wz)
{
    wx=uy*vz-uz*vy;
    wy=-(ux*vz-uz*vx);
    wz=ux*vy-uy*vx;
}

///Devuelve (ux,uy,uz)-(vx,vy,vz)
void vectDiferencia(double ux,double uy, double uz, double vx, double vy, double vz, double & wx, double & wy, double & wz)
{
    wx=ux-vx;
    wy=uy-vy;
    wz=uz-vz;
}

///Devuelve el vector r (R-R') dados los indices de (u,v)
void vectDiferenciaIndices(int nu, int nv, int nu2, int nv2, double &vx,double &vy, double &vz)
{
    vectDiferencia(X[nu][nv],Y[nu][nv],Z[nu][nv],X[nu2][nv2],Y[nu2][nv2],Z[nu2][nv2],vx,vy,vz);
}

void igualar(int geometria)
{
    int i,j;
    for (i=0; i<Nu; i++)
    {
        for(j=0; j<Nv; j++)
        {
            Santerior[i][j]=Sactual[i][j];
        }
    }
    if(geometria==3){
       for (i=0; i<Nu; i++)
    {
        for(j=0; j<Nv; j++)
        {
            Santerior2[i][j]=Sactual2[i][j];
        }
    }
    }
}

///Devuelve true si la desviación media respecto a la iteración anterior es menor que la precision fijada
bool compara(double precision)
{
    double suma=0;
    int i=0, j=0;
    bool condicion;
    for ( i=0; i<Nu; i++)
    {
        for( j=0; j<Nv; j++)
        {
            suma=suma + (Sactual[i][j]-Santerior[i][j])*(Sactual[i][j]-Santerior[i][j]);
        }
    }
    suma=suma/((Nu)*(Nv));
    printf("salir for %lf %lf el booleano %d\n",suma,precision,suma<=precision);
    if(suma<=precision)
    {
        condicion=true;
    }
    else
    {
        condicion=false;
    }
    return condicion;
}

///Crea los arrays X,Y,Z según la geometría
void definirXYZ(int geometria)
{
    double x,y,z;
    int i=0,j=0;

    for(i=0; i<Nu; i++)
    {
        for(j=0; j<Nv; j++)
        {
            conversion(i,j,geometria,x,y,z);
            X[i][j]=x;
            Y[i][j]=y;
            Z[i][j]=z;
        }
    }

    if(geometria==3){
        for(i=0; i<Nu; i++)
        {
            for(j=0; j<Nv; j++)
            {
                conversion(i,j,4,x,y,z);
                X1[i][j]=x;
                Y1[i][j]=y;
                Z1[i][j]=z;
            }
        }
    }

}

///Devuelve el valor de (u,v) dados sus indices
void cCoordDoble(int nu, int nv, double &u, double &v, int geometria)
{
    switch(geometria)
    {
    case 1:
        u=conversionIndividual(nu,Nu,0,2*pi);
        v=conversionIndividual(nv,Nv-1,0,pi);
        break;

    case 2:
        u=conversionIndividual(nu,Nu,0,2*pi);
        v=conversionIndividual(nv,Nv,0,2*pi);

        break;
    case 3:
        u=conversionIndividual(nu,Nu,0.001,1);
        v=conversionIndividual(nv,Nv,0,2*pi);
    }
}

///Crea un array que contiene los vectores ortogonales para cada (x,y,z)
void crearVectorOrtogonal(int tipo)
{
    double xu,xv,yu,yv,zu,zv,u,v;
    double du,dv;
    switch (tipo)
    {
    case 1:
        du=2*pi/(Nu);
        dv=(pi)/(Nv-1);
        for(int i=0; i<Nu; i++)
        {
            for(int j=0; j<Nv; j++)
            {
                cCoordDoble(i,j,u,v,tipo);
                xu=-r*sin(v)*sin(u);
                xv=r*cos(u)*cos(v);
                yu=r*cos(u)*sin(v);
                yv=r*sin(u)*cos(v);
                zu=0;
                zv=-r*sin(v);
                prodVectorial(xv,yv,zv,xu,yu,zu,VectorOrtogonal[i][j][0],VectorOrtogonal[i][j][1],VectorOrtogonal[i][j][2]);
                NormadS[i][j]=r*sin(v)*du*dv;
                ///NormadS[i][j]=norma(VectorOrtogonal[i][j][0],VectorOrtogonal[i][j][1],VectorOrtogonal[i][j][2])*du*dv;
                //printf("(%lf %lf) %lf %lf %lf, %lf %lf %lf\n",u,v, xu,yu,zu, xv, yv, zv);
            }
        }
        break;

    case 2:
        du=2*pi/(Nu);
        dv=2*pi/(Nv);
        for(int i=0; i<Nu; i++)
        {
            for(int j=0; j<Nv; j++)
            {
                cCoordDoble(i,j,u,v,tipo);
                xu=-cos(v)*r*sin(u);
                xv=-(R+r*cos(u))*sin(v);
                yu=-r*sin(v)*sin(u);
                yv=(R+r*cos(u))*cos(v);
                zu=r*cos(u);
                zv=0;
                prodVectorial(xv,yv,zv,xu,yu,zu,VectorOrtogonal[i][j][0],VectorOrtogonal[i][j][1],VectorOrtogonal[i][j][2]);
                NormadS[i][j]=norma(VectorOrtogonal[i][j][0],VectorOrtogonal[i][j][1],VectorOrtogonal[i][j][2])*du*dv;

            }
        }
        break;
    case 3:
        du=1./(Nu-1);
        dv=2*pi/(Nv);
        for(int i=0; i<Nu; i++)
        {
            for(int j=0; j<Nv; j++)
            {
                cCoordDoble(i,j,u,v,tipo);
                xu=cos(v);
                xv=-u*sin(v);
                yu=sin(v);
                yv=u*cos(v);
                zu=1;
                zv=0;
                prodVectorial(xv,yv,zv,xu,yu,zu,VectorOrtogonal[i][j][0],VectorOrtogonal[i][j][1],VectorOrtogonal[i][j][2]);
                NormadS[i][j]=norma(VectorOrtogonal[i][j][0],VectorOrtogonal[i][j][1],VectorOrtogonal[i][j][2])*du*dv;

            }
        }
        for(int i=0; i<Nu; i++)
        {
            for(int j=0; j<Nv; j++)
            {
                cCoordDoble(i,j,u,v,4);
                xu=cos(v);
                xv=-u*sin(v);
                yu=sin(v);
                yv=u*cos(v);
                zu=0;
                zv=0;
                prodVectorial(xu,yu,zu,xv,yv,zv,VectorOrtogonal2[i][j][0],VectorOrtogonal2[i][j][1],VectorOrtogonal2[i][j][2]);
                NormadS1[i][j]=norma(VectorOrtogonal2[i][j][0],VectorOrtogonal2[i][j][1],VectorOrtogonal2[i][j][2])*du*dv;

            }
        }
    }
}



void imprimirFichero(double matriz[][Nv],FILE*f)
{
    int i;
    for( i=0; i<Nu; i++)
    {
        for(int j=0; j<Nv; j++)
        {
            fprintf(f,"%lf\t",matriz[i][j]);
        }
        // fprintf(f,"%lf\t",matriz[i][0]);
        fprintf(f,"\n");
    }
    i=0;
    /**for(int j=0;j<Nv;j++){
                fprintf(f,"%lf\t",matriz[i][j]);
        }
             fprintf(f,"%lf\t",matriz[i][0]);
        fprintf(f,"\n");**/
}

void integral(int nu, int nv,double Ex, double Ey, double Ez, int geometria)
{
    int i,j;
    double aux,rx,ry,rz;
    Sactual[nu][nv]=0.0;
    //printf("%lf",Sactual[nu][nv]);
    switch(geometria)
    {
    case 1:
        if((nv==0||nv==(Nv-1))&&nu!=0)
        {
            Sactual[nu][nv]=Sactual[0][nv];

        }
        else if((nv==0||nv==(Nv-1))&&nu==0)
        {
            for(i=0; i<Nu; i++)
            {
                for(j=0; j<Nv; j++)
                {
                    if((i!=nu||j!=nv)&&(j!=0&&j!=(Nv-1)))
                    {
                        vectDiferenciaIndices(nu,nv,i,j,rx,ry,rz);
                        aux=Sactual[nu][nv];

                        Sactual[nu][nv]=aux+1/(36*pi)*pow(10,-9)*prodEscalar(Ex,Ey,Ez,VectorOrtogonal[nu][nv][0],VectorOrtogonal[nu][nv][1],VectorOrtogonal[nu][nv][2])+1/(2*pi*(norma(rx,ry,rz)*norma(rx,ry,rz)))*coseno(X[nu][nv],Y[nu][nv],Z[nu][nv],rx,ry,rz)*NormadS[i][j]*Santerior[i][j];

                        //printf("Sactual(%d)(%d)=%lf aux:%lf %lf x %lf y %lf z %lf %lf coseno: %lf  ; %lf\n",i,j,Sactual[nu][nv],aux,norma(rx,ry,rz),VectorOrtogonal[nu][nv][0],VectorOrtogonal[nu][nv][1],VectorOrtogonal[nu][nv][2],NormadS[i][j],coseno(rx,ry,rz,VectorOrtogonal[nu][nv][0],VectorOrtogonal[nu][nv][1],VectorOrtogonal[nu][nv][2]),1/(2*pi*(norma(rx,ry,rz)*norma(rx,ry,rz)))*coseno(rx,ry,rz,VectorOrtogonal[nu][nv][0],VectorOrtogonal[nu][nv][1],VectorOrtogonal[nu][nv][2])*NormadS[i][j]);
                    }
                    else if(((nv==0&&j==(Nv-1))||(nv==(Nv-1)&&j==0))&&i==0)
                    {
                        vectDiferenciaIndices(nu,nv,i,j,rx,ry,rz);
                        aux=Sactual[nu][nv];

                        Sactual[nu][nv]=aux+1/(36*pi)*pow(10,-9)*prodEscalar(Ex,Ey,Ez,VectorOrtogonal[nu][nv][0],VectorOrtogonal[nu][nv][1],VectorOrtogonal[nu][nv][2])+1/(2*pi*(norma(rx,ry,rz)*norma(rx,ry,rz)))*coseno(X[nu][nv],Y[nu][nv],Z[nu][nv],rx,ry,rz)*NormadS[i][j]*Santerior[i][j];

                    }

                }
            }
        }
        else
        {
            for(i=0; i<Nu; i++)
            {
                for(j=0; j<Nv; j++)
                {
                    if((i!=nu||j!=nv)&&(j!=0&&j!=(Nv-1)))
                    {
                        vectDiferenciaIndices(nu,nv,i,j,rx,ry,rz);
                        aux=Sactual[nu][nv];

                        Sactual[nu][nv]=aux+1/(36*pi)*pow(10,-9)*prodEscalar(Ex,Ey,Ez,VectorOrtogonal[nu][nv][0],VectorOrtogonal[nu][nv][1],VectorOrtogonal[nu][nv][2])+1/(2*pi*(norma(rx,ry,rz)*norma(rx,ry,rz)))*coseno(X[nu][nv],Y[nu][nv],Z[nu][nv],rx,ry,rz)*NormadS[i][j]*Santerior[i][j];

                        //printf("Sactual(%d)(%d)=%lf aux:%lf %lf x %lf y %lf z %lf %lf coseno: %lf  ; %lf\n",i,j,Sactual[nu][nv],aux,norma(rx,ry,rz),VectorOrtogonal[nu][nv][0],VectorOrtogonal[nu][nv][1],VectorOrtogonal[nu][nv][2],NormadS[i][j],coseno(rx,ry,rz,VectorOrtogonal[nu][nv][0],VectorOrtogonal[nu][nv][1],VectorOrtogonal[nu][nv][2]),1/(2*pi*(norma(rx,ry,rz)*norma(rx,ry,rz)))*coseno(rx,ry,rz,VectorOrtogonal[nu][nv][0],VectorOrtogonal[nu][nv][1],VectorOrtogonal[nu][nv][2])*NormadS[i][j]);
                    }
                    else if(i==0&&(j==0||j==(Nv-1)))
                    {
                        //printf("u=%d,v=%d, i=%d, j=%d, S=%lf \n",nu,nv,i,j,Sactual[nu][nv]);
                        vectDiferenciaIndices(nu,nv,i,j,rx,ry,rz);
                        aux=Sactual[nu][nv];

                        Sactual[nu][nv]=aux+1/(36*pi)*pow(10,-9)*prodEscalar(Ex,Ey,Ez,VectorOrtogonal[nu][nv][0],VectorOrtogonal[nu][nv][1],VectorOrtogonal[nu][nv][2])+1/(2*pi*(norma(rx,ry,rz)*norma(rx,ry,rz)))*coseno(X[nu][nv],Y[nu][nv],Z[nu][nv],rx,ry,rz)*NormadS[i][j]*Santerior[i][j];

                        //printf("u=%d,v=%d, i=%d, j=%d, S=%lf \n",nu,nv,i,j,Sactual[nu][nv]);
                    }

                }
            }
        }
        break;

    case 2:
        for(i=0; i<Nu; i++)
        {
            for(j=0; j<Nv; j++)
            {
                if(i!=nu||j!=nv)
                {
                    vectDiferenciaIndices(nu,nv,i,j,rx,ry,rz);
                    aux=Sactual[nu][nv];

                    Sactual[nu][nv]=aux+1/(36*pi)*pow(10,-9)*prodEscalar(Ex,Ey,Ez,VectorOrtogonal[nu][nv][0],VectorOrtogonal[nu][nv][1],VectorOrtogonal[nu][nv][2])+1/(2*pi*(norma(rx,ry,rz)*norma(rx,ry,rz)))*coseno(VectorOrtogonal[nu][nv][0],VectorOrtogonal[nu][nv][1],VectorOrtogonal[nu][nv][2],rx,ry,rz)*NormadS[i][j]*Santerior[i][j];

                    //printf("Sactual(%d)(%d)=%lf aux:%lf %lf x %lf y %lf z %lf %lf coseno: %lf  ; %lf\n",i,j,Sactual[nu][nv],aux,norma(rx,ry,rz),VectorOrtogonal[nu][nv][0],VectorOrtogonal[nu][nv][1],VectorOrtogonal[nu][nv][2],NormadS[i][j],coseno(rx,ry,rz,VectorOrtogonal[nu][nv][0],VectorOrtogonal[nu][nv][1],VectorOrtogonal[nu][nv][2]),1/(2*pi*(norma(rx,ry,rz)*norma(rx,ry,rz)))*coseno(rx,ry,rz,VectorOrtogonal[nu][nv][0],VectorOrtogonal[nu][nv][1],VectorOrtogonal[nu][nv][2])*NormadS[i][j]);
                }
            }
        }
        break;
    case 3:

        for(i=0; i<Nu; i++)
        {
            for(j=0; j<Nv; j++)
            {
                if(i!=nu||j!=nv)
                {
                    vectDiferenciaIndices(nu,nv,i,j,rx,ry,rz);
                    aux=Sactual[nu][nv];

                    Sactual[nu][nv]=aux+1/(36*pi)*pow(10,-9)*prodEscalar(Ex,Ey,Ez,VectorOrtogonal[nu][nv][0],VectorOrtogonal[nu][nv][1],VectorOrtogonal[nu][nv][2])+1/(2*pi*(norma(rx,ry,rz)*norma(rx,ry,rz)))*coseno(VectorOrtogonal[nu][nv][0],VectorOrtogonal[nu][nv][1],VectorOrtogonal[nu][nv][2],rx,ry,rz)*NormadS[i][j]*Santerior[i][j];

                    //printf("Sactual(%d)(%d)=%lf aux:%lf %lf x %lf y %lf z %lf %lf coseno: %lf  ; %lf\n",i,j,Sactual[nu][nv],aux,norma(rx,ry,rz),VectorOrtogonal[nu][nv][0],VectorOrtogonal[nu][nv][1],VectorOrtogonal[nu][nv][2],NormadS[i][j],coseno(rx,ry,rz,VectorOrtogonal[nu][nv][0],VectorOrtogonal[nu][nv][1],VectorOrtogonal[nu][nv][2]),1/(2*pi*(norma(rx,ry,rz)*norma(rx,ry,rz)))*coseno(rx,ry,rz,VectorOrtogonal[nu][nv][0],VectorOrtogonal[nu][nv][1],VectorOrtogonal[nu][nv][2])*NormadS[i][j]);
                }
            }
        }

        for(i=0; i<Nu; i++)
        {
            for(j=0; j<Nv; j++)
            {
                if(i!=nu||j!=nv)
                {
                    vectDiferencia(X1[i][j],Y1[i][j],Z1[i][j],X[nu][nv],Y[nu][nv],Z[nu][nv],rx,ry,rz);
                    aux=Sactual[nu][nv];

                    Sactual[nu][nv]=aux+1/(36*pi)*pow(10,-9)*prodEscalar(Ex,Ey,Ez,VectorOrtogonal2[nu][nv][0],VectorOrtogonal2[nu][nv][1],VectorOrtogonal2[nu][nv][2])+1/(2*pi*(norma(rx,ry,rz)*norma(rx,ry,rz)))*coseno(VectorOrtogonal2[nu][nv][0],VectorOrtogonal2[nu][nv][1],VectorOrtogonal2[nu][nv][2],rx,ry,rz)*NormadS1[i][j]*Santerior2[i][j];

                    //printf("Sactual(%d)(%d)=%lf aux:%lf %lf x %lf y %lf z %lf %lf coseno: %lf  ; %lf\n",i,j,Sactual[nu][nv],aux,norma(rx,ry,rz),VectorOrtogonal[nu][nv][0],VectorOrtogonal[nu][nv][1],VectorOrtogonal[nu][nv][2],NormadS[i][j],coseno(rx,ry,rz,VectorOrtogonal[nu][nv][0],VectorOrtogonal[nu][nv][1],VectorOrtogonal[nu][nv][2]),1/(2*pi*(norma(rx,ry,rz)*norma(rx,ry,rz)))*coseno(rx,ry,rz,VectorOrtogonal[nu][nv][0],VectorOrtogonal[nu][nv][1],VectorOrtogonal[nu][nv][2])*NormadS[i][j]);
                }
            }
        }
    }
    //printf("Sactual: %lf\n",Sactual[nu][nv]);

}

int main()
{

    int geometria=1;

    ///CAMPO EXTERNO:
    double Ex=0,Ey=0,Ez=0;
    FILE *fx=fopen("coordenadasX.txt","w"),*fy=fopen("coordenadasY.txt","w"),*fz=fopen("coordenadasZ.txt","w"),*fs=fopen("valoresS.txt","w"),*fDx=fopen("vectDx.txt","w"),*fDy=fopen("vectDy.txt","w"),*fDz=fopen("vectDz.txt","w"),
            *fxx=fopen("coordenadasX1.txt","w"),*fyy=fopen("coordenadasY1.txt","w"),*fzz=fopen("coordenadasZ1.txt","w");
    definirXYZ(geometria);
    crearVectorOrtogonal(geometria);
    imprimirFichero(X,fx);
    imprimirFichero(Y,fy);
    imprimirFichero(Z,fz);
    imprimirFichero(X1,fxx);
    imprimirFichero(Y1,fyy);
    imprimirFichero(Z1,fzz);
    int i,j;
    for(i=0; i<Nu; i++)
    {
        for(j=0; j<Nv; j++)
        {
            fprintf(fDx,"%lf\t",VectorOrtogonal2[i][j][0]);
        }
        fprintf(fDx,"\n");
    }
    for(i=0; i<Nu; i++)
    {
        for(j=0; j<Nv; j++)
        {
            fprintf(fDy,"%lf\t",VectorOrtogonal2[i][j][1]);
        }
        fprintf(fDy,"\n");
    }
    for(i=0; i<Nu; i++)
    {
        for(j=0; j<Nv; j++)
        {
            fprintf(fDz,"%lf\t",VectorOrtogonal2[i][j][2]);
        }
        fprintf(fDz,"\n");
    }
    int f=0;
    double Qtotal=10;
    printf("for");

    double densidadInicial=Qtotal/(4*pi*r*r);

    for(i=0; i<Nu; i++)
    {
        for(j=0; j<Nv; j++)
        {
            Santerior[i][j]=0;
            Sactual[i][j]=densidadInicial;

            if(geometria==3){
               Santerior2[i][j]=0;
                Sactual2[i][j]=densidadInicial;
            }
        }
    }
    while(!compara(0.001)||f==0)
        //for(int k=0;k<1;k++)
    {
        printf("vuelta %d\n",f);
        igualar(geometria);
        for(i=0; i<Nu; i++)
        {
            for(j=0; j<Nv; j++)
            {
                integral(i,j,Ex,Ey,Ez,geometria);
            }
        }
        f++;
    }

    imprimirFichero(Sactual,fs);


    fclose(fx);
    fclose(fy);
    fclose(fz);
    fclose(fs);

    fclose(fxx);
    fclose(fyy);
    fclose(fzz);
}
