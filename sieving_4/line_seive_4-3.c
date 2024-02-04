// Program 4.3: Line sieving with two levels of cached memory

#define SIZE 510000     /* Max. Number of prime pairs */
#define ZONE 9000000    /* Size of each zone, mult of CACHE2 */
#define CACHE1 25000    /* Approximate size of L1 CACHE */
#define CACHE2 750000   /* Approximate size of L2 CACHE */
#define START_ZONE -500 /* Number of first sieve zone */
#define END_ZONE 500    /* Number of last sieve zone */
#define THRESHOLD 12    /* Num. of hits per smooth candidate */
#define START 30        /* Prime pairs below START not sieved */
int root[SIZE], mod[SIZE]; /* Arrays for prime pairs */

char Zone[ZONE];           /* Array for sieving */
int Size;         /*True number of pairs */    
int Size1, Size2; /*Limits to fit in CACHE1 and CACHE2 */

main()
{ int i,j,offset;  int limit2,limit1;

  Read();/*Get prime pairs, set Size1/2 to fit CACHE1/2*/
  for(i=START;i<Size;i++) {

    int r,q; /* Loop shifts primes pairs to first zone */
    q=mod[i]; offset=MulMod(START_ZONE,ZONE,q);

    r=root[i]-offset; if (r<0) r+=q; root[i]=r; }

  for (j=START_ZONE;j<END_ZONE;j++){
    for(i=0;i<ZONE;i++) Zone[i]=0; /*Reset sieve zone */

    for(limit1=limit2=0;limit2<=ZONE;limit2+=CACHE2) {

      for(;limit1<=limit2;limit1+=CACHE1) {
        for(i=START;i<Size1;i++) { /* Sieve small*/

          int r,m;  r=root[i]; m=mod[i];

          while (r<limit1) {Zone[r]++; r+=m;}

          root[i]=r;}}
      for(i=Size1;i<Size2;i++) { /* Sieve medium*/

        int r,m; r=root[i]; m=mod[i];

        while (r<limit2) { Zone[r]++; r+=m;}

        root[i]=r;}}
    for(i=START;i<Size2;i++) { /* Shift to next zone */

      root[i]=root[i]-ZONE;}
    for(i=Size2;i<Size;i++) { /* Sieve large */

      int r,m; r=root[i]; m=mod[i];

      while (r<ZONE) {Zone[r]++; r+=m;}

      root[i]=r-ZONE;}
    for(i=0;i<ZONE;i++){ /* Detect and print smooth candidates*/

      if (Zone[i]>=THRESHOLD) {printf("F(%d*SZ+%d);\n",j,i);

	fflush(stdout); 
      }
    }
  }
}
