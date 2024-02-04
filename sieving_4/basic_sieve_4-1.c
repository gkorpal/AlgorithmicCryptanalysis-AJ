// Program 4.1: Basic code for Eratostenes's sieve

#include <stdio.h>
#include <stdlib.h>


typedef unsigned int packedbool;
packedbool *_IsPrime;

#define NpckBool (8*sizeof(packedbool))
#define GetIsPrime(x) \
   (_IsPrime[(x-1)/NpckBool]>>((x-1)%NpckBool))&1
#define SetCompo(x) \
   _IsPrime[(x-1)/NpckBool]&=~(1UL<<((x-1)%NpckBool))


SievePrimes(int Limit)
{
  int i,j, tabLimit;

  tabLimit=(Limit+NpckBool-1)/NpckBool;
  for (i=0;i<tabLimit;i++) _IsPrime[i]=~0;

  for (i=2;i*i<=Limit;i++){

    if (GetIsPrime(i)) {
      for (j=i*i;j<=Limit;j+=i){

        SetCompo(j); } } }

  printf("List of primes up to %d:\n",Limit);
  for (i=2;i<=Limit;i++){

    if (GetIsPrime(i)) {
      printf("%d\n",i); } } }

main()
{
  int Limit;
  printf("Enter limit ?");
  scanf("%d",&Limit);

  _IsPrime=(packedbool *)malloc((Limit+NpckBool-1)/8);

  SievePrimes(Limit);
  free(_IsPrime);
}
