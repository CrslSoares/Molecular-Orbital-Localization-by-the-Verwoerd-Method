#include <stdio.h>
#define max 20
#define TOL 1e-12
#include <math.h>
#include <stdlib.h>

void jacobi(double a[max][max], double b[max][max], int n)
{
  int i,j,k,l,iter=0;
  double M,alpha,t,h,c,s,r,u,v;

  for (i=1;i<=n;i++)
     { for (j=1;j<=n;j++) b[i][j]=0;
       b[i][i]=1;
     }

  do { M=0;
       for (i=1;i<n;i++)
         for (j=i+1;j<=n;j++)
            if (fabs(a[i][j])>M) { k=i; l=j; M=fabs(a[i][j]);}
       if (M<TOL) break;
       iter++;
       alpha=0.5*(a[l][l]-a[k][k])/a[k][l];
       t=sqrt(1+alpha*alpha);
       if (alpha<0) t=-t;
       t-=alpha;
       h=t*a[k][l];
       c=1/sqrt(1+t*t);
       s=t*c;
       r=s/(1+c);
       for (j=1;j<=n;j++)
          { u=b[j][k];
            v=b[j][l];
            b[j][k]-=s*(v+r*u);
            b[j][l]+=s*(u-r*v);
            if (j!=k && j!=l)
              { u=a[k][j];
                v=a[l][j];
                a[k][j]-=s*(v+r*u);
                a[l][j]+=s*(u-r*v);
                a[j][k]=a[k][j];
                a[j][l]=a[l][j];
              }
          }
       a[k][k]-=h;
       a[l][l]+=h;
       a[k][l]=a[l][k]=0;
     } while (M>TOL && iter<1000);
}

void transposer(double a[max][max], double b[max][max], int c, int d)
{
    int i,j;
    for(i=1;i<= c;i++)
    {
        for(j=1;j<=d;j++)
        {
            b[j][i]=a[i][j];
        }
    }
}

void multiplication(double a[max][max], double b[max][max], double c[max][max], int d)
{
   int i,j,k1;
   for(i=1;i<=d;i++)
    {
        for(j=1;j<=d;j++)
        {
            c[i][j]=0.0;
            for(k1=1;k1<=d+1;k1++)
            {
                c[i][j] += a[i][k1]*b[k1][j];
            }
        }
    }
}

int main()

{
/** #### Declaration Des Variables #### **/

    double E[max][max],Et[max][max];
    double R1[max][max],R2[max][max],R3[max][max],R4[max][max];
    double T1[max][max],T2[max][max],T3[max][max],T4[max][max];
    double P1[max][max],P2[max][max],P3[max][max],P4[max][max];
    double D1[max][max],D2[max][max],D3[max][max],D4[max][max];
    double B1[max][max],B2[max][max],B3[max][max],B4[max][max];
    double X1[max][max],X2[max][max],X3[max][max],X4[max][max];
    double L1[max][max],L2[max][max],L3[max][max],L4[max][max];
    double temp1,temp2,temp3,temp4;
    int p1,p2,p3,p4;
    int i, j, k1, k2;
    int N_C,N_CH,N_H,n_OM,n_OC;
    char line[200];
    char nom[1];

/** #### Lecture du Fichier #### **/

    FILE *pf;
    FILE *pf1;
    FILE *pf2;

    pf=fopen("CH4.mgf","r");
    pf1=fopen("CH4 sortir.mgf", "w");
    pf2=fopen("CH4.out","r");

     //Lecture du Fichier CH4.out
    for(i=1; i<32;i++)
    {
        fgets(line, sizeof(line), pf2);
    }

    do
    {
        j=0;
        fscanf(pf2, "%d       %s", &j, &nom);
        fgets(line, sizeof(line), pf2);
    }
    while(nom[0] == 'C');

    N_C = j -1; // Nombre de Carbons.

    printf("Nombre de  C = %d\n", N_C);
    printf("\n");
    printf("\n");

    N_H = 4*N_C; // Nombre de Hidrogens.
    N_CH = N_C + N_H; // Nombre de atomes.

    //Pour lire CH4.mgf est necessaire de savoir le nombre de lignes avant les donnees.
    // 1 ligne de title et 2 x N_CH lignes de donnees que ne sont pas important.

    for(i=1;i<=2*N_CH+1;i++)
    {
       fgets(line, sizeof(line), pf);
    }
    //Je n'ai pas autre fichier, par example C2H8, pour savoir comme des donnees seraient organizer.
    //Mais mon idée serait pour aprendre le modèle et faire le algorithme généraliste.
    //J'ai choisi pour lire seulement les orbitales que sont occupé.
    fgets(line, sizeof(line), pf);
    fscanf(pf, "%lf %lf %lf %lf %lf", &E[1][1], &E[1][2], &E[1][3], &E[1][4], &E[1][5]);
    fgets(line, sizeof(line), pf);
    fscanf(pf, "%lf %lf %lf", &E[1][6], &E[1][7], &E[1][8]);
    fgets(line, sizeof(line), pf);
    fgets(line, sizeof(line), pf);
    fscanf(pf, "%lf %lf %lf %lf %lf", &E[2][1], &E[2][2], &E[2][3], &E[2][4], &E[2][5]);
    fgets(line, sizeof(line), pf);
    fscanf(pf, "%lf %lf %lf", &E[2][6], &E[2][7], &E[2][8]);
    fgets(line, sizeof(line), pf);
    fgets(line, sizeof(line), pf);
    fscanf(pf, "%lf %lf %lf %lf %lf", &E[3][1], &E[3][2], &E[3][3], &E[3][4], &E[3][5]);
    fgets(line, sizeof(line), pf);
    fscanf(pf, "%lf %lf %lf", &E[3][6], &E[3][7], &E[3][8]);
    fgets(line, sizeof(line), pf);
    fgets(line, sizeof(line), pf);
    fscanf(pf, "%lf %lf %lf %lf %lf", &E[4][1], &E[4][2], &E[4][3], &E[4][4], &E[4][5]);
    fgets(line, sizeof(line), pf);
    fscanf(pf, "%lf %lf %lf", &E[4][6], &E[4][7], &E[4][8]);

    fclose(pf);

    n_OM = 4*N_C + N_H; // Nombre de Orbitales Moleculaires
    n_OC = n_OM/2; // Nombre de Orbitales Occupés

/** #### Matrice R #### **/

        for(i=1;i<=n_OC;i++)
        {
                for(j=1;j<=n_OC;j++)
                {
                    R1[i][j] = E[i][j];
                    R2[i][j] = E[i][j];
                    R3[i][j] = E[i][j];
                    R4[i][j] = E[i][j];
                }
        }
        //R1 = Liaison C-H1
        for(i=1;i <= n_OC; i++)
        {
            R1[i][5] = E[i][5];
        }
        //R2 = Liaison C-H2
        for(i=1; i <= n_OC; i++)
        {
            R2[i][5] = E[i][6];
        }
        //R3 = Liaison C-H3
        for(i=1; i <= n_OC; i++)
        {
            R3[i][5] = E[i][7];
        }
        //R4 = = Liaison C-H4
        for(i=1; i <= n_OC; i++)
        {
            R4[i][5] = E[i][8];
        }

/** #### Transposée Matrice R(matrice T) #### **/

            //T1
            transposer(R1,T1,n_OC,n_OC+1);
            //T2
            transposer(R2,T2,n_OC,n_OC+1);
            //T3
            transposer(R3,T3,n_OC,n_OC+1);
            //T4
            transposer(R4,T4,n_OC,n_OC+1);

/** #### Multiplication RxT(matrice P) #### **/

        //P1
        multiplication(R1,T1,P1,n_OC);
        //P2
        multiplication(R2,T2,P2,n_OC);
        //P3
        multiplication(R3,T3,P3,n_OC);
        //P4
        multiplication(R4,T4,P4,n_OC);

/** #### Diagnolize P(Matrice D et B) #### **/

        for(i=1; i <= n_OC; i++)
        {
            for(j=1; j <= n_OC; j++)
            {
                D1[i][j] = P1[i][j];
                D2[i][j] = P2[i][j];
                D3[i][j] = P3[i][j];
                D4[i][j] = P4[i][j];
            }
        }

        jacobi(D1,B1,n_OC);
        jacobi(D2,B2,n_OC);
        jacobi(D3,B3,n_OC);
        jacobi(D4,B4,n_OC);

/** #### Matrice X #### **/

    temp1=0;
    p1=1;
    for(i=1;i<=n_OC;i++)
    {
        if(fabs(1-temp1)>fabs(1-D1[i][i]))
        {
                temp1=D1[i][i];
                p1=i;
        }
    }
      for(j=1;j<=n_OC;j++)
    {
        X1[j][1]=B1[j][p1];
    }

    temp2=0;
    p2=1;
    for(j=1; j <=n_OC; j++)
    {
        if(fabs(1-temp2)>fabs(1-D2[j][j]))
        {
                temp2=D2[j][j];
                p2=j;
        }
    }
    for(j=1; j<=n_OC; j++)
    {
        X2[j][1]=B2[j][p2];
    }

    temp3=0;
    p3=1;
    for(k1=1; k1 <=n_OC; k1++)
    {
        if(fabs(1-temp3)>fabs(1-D3[k1][k1]))
        {
                temp3=D3[k1][k1];
                p3=k1;
        }
    }
    for(j=1; j<=n_OC; j++)
    {
        X3[j][1]=B3[j][p3];
    }

    temp4=0;
    p4=1;
    for(k2=1; k2 <=n_OC; k2++)
    {
        if(fabs(1-temp4)>fabs(1-D4[k2][k2]))
        {
                temp4=D4[k2][k2];
                p4=k2;
        }
    }
    for(j=1; j<=n_OC; j++)
    {
        X4[j][1]=B4[j][p4];
    }

/** #### Transposée Matrice E #### **/

    transposer(E,Et,n_OC,n_OM);

/** #### Matrice L #### **/

     //L1
     for(i=1; i <=n_OM; i++)
        {
            L1[i][1] = 0.0;
            for(k1=1; k1 <=n_OM; k1++)
            {
                L1[i][1] += Et[i][k1]*X1[k1][1];
            }
        }

    //L2
     for(i=1; i <=n_OM; i++)
        {
            L2[i][1] = 0.0;
            for(k1=1; k1 <=n_OM; k1++)
            {
                L2[i][1] += Et[i][k1]*X2[k1][1];
            }
        }

    //L3
     for(i=1; i <=n_OM; i++)
        {
            L3[i][1] = 0.0;
            for(k1=1; k1 <=n_OM; k1++)
            {
                L3[i][1] += Et[i][k1]*X3[k1][1];
            }
        }

     //L4
     for(i=1; i <=n_OM; i++)
        {
            L4[i][1] = 0.0;
            for(k1=1; k1 <=n_OM; k1++)
            {
                L4[i][1] += Et[i][k1]*X4[k1][1];
            }

        }

/** #### Print File CH4 sortir.mgf #### **/
/* OBS:
Le fichier "CH4 sortir,mgf" n’a pas été utilisé car la mise en forme des données n’est pas cohérente.
Le fichier "CH4end.mgf" est le quie a été utilisé sur Jmol.
Le fichier "CH4end.mgf" a été obéi en modifiant manuellement le format des données imprimées sur "CH4 sortir,mgf".
*/

    pf=fopen("CH4.mgf","r");
    for(i=1; i <= 2*N_CH+1;i++)
    {
       fgets(line, sizeof(line), pf);
       fprintf(pf1,"%s",line);
    }
    fgets(line, sizeof(line), pf);
    fprintf(pf1,"%s",line);
    fprintf(pf1, " %15.8E %15.8E %15.8E %15.8E %15.8E\n", L1[1][1], L1[2][1], L1[3][1], L1[4][1], L1[5][1]);
    fprintf(pf1, " %15.8E %15.8E %15.8E\n", L1[6][1], L1[7][1], L1[8][1]);
    fgets(line, sizeof(line), pf);
    fgets(line, sizeof(line), pf);
    fgets(line, sizeof(line), pf);
    fprintf(pf1,"%s",line);
    fprintf(pf1, " %15.8E %15.8E %15.8E %15.8E %15.8E\n", L2[1][1], L2[2][1], L2[3][1], L2[4][1], L2[5][1]);
    fprintf(pf1, " %15.8E %15.8E %15.8E\n", L2[6][1], L2[7][1], L2[8][1]);
    fgets(line, sizeof(line), pf);
    fgets(line, sizeof(line), pf);
    fgets(line, sizeof(line), pf);
    fprintf(pf1,"%s",line);
    fprintf(pf1, " %15.8E %15.8E %15.8E %15.8E %15.8E\n", L3[1][1], L3[2][1], L3[3][1], L3[4][1], L3[5][1]);
    fprintf(pf1, " %15.8E %15.8E %15.8E\n", L3[6][1], L3[7][1], L3[8][1]);
    fgets(line, sizeof(line), pf);
    fgets(line, sizeof(line), pf);
    fgets(line, sizeof(line), pf);
    fprintf(pf1,"%s",line);
    fprintf(pf1, " %15.8E %15.8E %15.8E %15.8E %15.8E\n", L4[1][1], L4[2][1], L4[3][1], L4[4][1], L4[5][1]);
    fprintf(pf1, " %15.8E %15.8E %15.8E\n", L4[6][1], L4[7][1], L4[8][1]);
    fgets(line, sizeof(line), pf);
    fgets(line, sizeof(line), pf);
    for(i=1; i <= 2*n_OC+4;i++)
    {
       fgets(line, sizeof(line), pf);
       fprintf(pf1,"%s",line);
    }
    for(i=36; i <= 48;i++)
    {
       fgets(line, sizeof(line), pf);
       fprintf(pf1,"%s",line);
    }

/** #### Imprimer des exemples de matrices #### **/

    printf("Matrice E:\n");
    for(i=1;i<= n_OC;i++)
    {
        for (j=1; j <= n_OM; j++)
        {
            printf("%9.7E ", E[i][j]);
        }
            printf("\n");
    }
        printf("\n");
        printf("\n");
        printf("\n");

    printf("Matrice R - Par example R1(Liaison C-H1):\n");
    for(i=1;i<= n_OC;i++)
    {
        for (j=1; j <= n_OC+1; j++)
        {
            printf("%9.7E ", R1[i][j]);
        }
        printf("\n");
    }
    printf("\n");
    printf("\n");
    printf("\n");

    printf("Matrice T - Par example T1(Transposé R1):\n");
    for(i=1;i <=n_OC+1;i++)
    {
        for (j=1; j <= n_OC; j++)
        {
            printf("%9.7E ", T1[i][j]);
        }
       printf("\n");
    }
    printf("\n");
    printf("\n");
    printf("\n");

    printf("Matrice P - Par example P1(R1xT1):\n");
    for(i=1; i <=n_OC; i++)
    {
        for (j=1;j<=n_OC; j++)
        {
            printf("%9.7E ", P1[i][j]);
        }
        printf("\n");
    }
    printf("\n");
    printf("\n");
    printf("\n");


    printf("Matrice X - Par Example X1(Diagnolasisation P1):\n");
    for(i=1; i <=n_OC; i++)
    {
        printf("%9.7E ", X1[i][1]);
        printf("\n");
    }
    printf("\n");
    printf("\n");
    printf("\n");

    printf("### Matrice L ###");
    printf("\n");
    printf("\n");

    printf("Matriz L1:\n");
    for(i=1; i <=n_OM; i++)
    {
        printf("%9.7E ", L1[i][1]);
        printf("\n");
    }
    printf("\n");
    printf("\n");
    printf("\n");

    printf("Matriz L2:\n");
    for(i=1; i <=n_OM; i++)
    {
        printf("%9.7E ", L2[i][1]);
        printf("\n");
    }
    printf("\n");
    printf("\n");
    printf("\n");

    printf("Matriz L3:\n");
    for(i=1; i <=n_OM; i++)
    {
        printf("%9.7E ", L3[i][1]);
        printf("\n");
    }
    printf("\n");
    printf("\n");
    printf("\n");

    printf("Matriz L4:\n");
    for(i=1; i <=n_OM; i++)
    {
        printf("%9.7E ", L4[i][1]);
        printf("\n");
    }

    fclose(pf);
    fclose(pf1);
    fclose(pf2);

    return 0;
}
