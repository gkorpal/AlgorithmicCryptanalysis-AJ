// Program 4.2: Code for low-memory usage Segmented Eratostenes's sieve

#include <stdio.h>
#include <stdlib.h>


typedef unsigned int packedbool;
typedef long long int int64;
packedbool *_IsBPrime;

packedbool *_IsPrime;
int wcount=0;
#define Npck (8*sizeof(packedbool))
#define GetBasePrime(x) (_IsBPrime[(x)/Npck]>>((x)%Npck))&1
#define SetBaseCompo(x) _IsBPrime[(x)/Npck]&=~(1UL<<((x)%Npck))
#define GetIsPrime(x) (_IsPrime[(x)/Npck]>>((x)%Npck))&1
#define SetCompo(x) _IsPrime[(x)/Npck]&=~(1UL<<((x)%Npck))


int64 *Offset;

int InitialSievePrimes(int Limit, int do_print)
{

  int i, j, count, tabLimit;

  count=0;

  tabLimit=(Limit+Npck-1)/Npck;
  for (i=0;i<tabLimit;i++) _IsBPrime[i]=~0;

  for (i=2;i*i<Limit;i++){
    if (GetBasePrime(i)) {

      count++;
      for (j=i*i;j<Limit;j+=i){

        SetBaseCompo(j);
      }
    }
  }
  for (;i<Limit;i++){
    if (GetBasePrime(i)) count++;
  }

  Offset=(int64 *)malloc(count*sizeof(int64));
  count=0;

  for (i=2;i<Limit;i++){
    if (GetBasePrime(i)) {

      if (do_print) printf("%d\n",i);
      j=Limit%i; if (j!=0) j=i-j; Offset[count]=j;

      count++;
    }
  }
}

int SievePrimesInterval(int64 offset, int Length, int do_print)
{

  int i, j, count, tabLimit;

  count=0;

  tabLimit=(Length+Npck-1)/Npck;
  for (i=0;i<tabLimit;i++) _IsPrime[i]=~0;

  for (i=2;i<Length;i++) {
    if (GetBasePrime(i)) {

      for (j=Offset[count];j<Length;j+=i) {

        SetCompo(j);
      }
      Offset[count]=j-Length;
      count++;
    }
  }

  for (i=0;i<Length;i++){
    if (GetIsPrime(i)) {

      wcount++;
      if (do_print)
	printf("%lld\n",offset+i);
    }
  }
}


main()
{
  int i,j;
  int64 Limit, tmp;

  int sqrt;
  
  printf("Enter limit ?");
  scanf("%lld",&Limit);

  for(tmp=0;tmp*tmp<=Limit;tmp++);
  sqrt=tmp;

  _IsBPrime=(packedbool *)malloc((sqrt+Npck-1)/8);

  _IsPrime=(packedbool *)malloc((sqrt+Npck-1)/8);

  InitialSievePrimes(sqrt,0);
  for(tmp=sqrt;tmp<Limit-sqrt;tmp+=sqrt)

    SievePrimesInterval(tmp,sqrt,0);
  SievePrimesInterval(tmp,Limit-tmp,1);

  free(_IsPrime);
  free(_IsBPrime);
}

