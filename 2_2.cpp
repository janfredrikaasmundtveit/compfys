#include<cstdlib>
#include<iostream>
#include<cmath>
int main(){
double tup;
int k=5;
double tdown;
double sup[k];
for(int i=1;i<=k;i++){
tup=1.0/(1.0+(i-1)*1.0);
sup[i-1]=tup;
//printf("%f\n",tup);
}

double tsup=0;
for(int i=1;i<=k;i++){
tsup=tsup+sup[i-1];
}
printf("%f\n",tsup);

return 0;

}
