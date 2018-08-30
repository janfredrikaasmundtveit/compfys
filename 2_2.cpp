#include<cstdlib>
#include<iostream>
#include<cmath>
int main(int argc, char* argv[]){
double tup;
int k=atoi(argv[1]);
//double* sup[];
double *sup = new double [k];
for(int i=1;i<=k;i++){
tup=1.0/(1.0+(i-1)*1.0);
sup[i-1]=tup;
//printf("%f\n",tup);
}

double tsup=0;
for(int i=1;i<=k;i++){
tsup=tsup+sup[i-1];
}
printf("u %f\n",tsup);

double tdown;

double *sdown = new double [k];
for(int i=k;i>=0;i--){
tdown=1.0/(1.0+(i-1)*1.0);
sdown[i-1]=tdown;
//printf("%f\n",tup);
}

double tsdown=0;
for(int i=1;i<=k;i++){
tsdown=tsdown+sdown[i-1];
}
printf("d %f\n",tsdown);


return 0;

}
