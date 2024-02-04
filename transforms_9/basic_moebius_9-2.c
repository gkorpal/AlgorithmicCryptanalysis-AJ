// Program 9.2: Basic Moebius transform code

#define TYPE unsigned int

/*In place Moebius transform*/
void Moebius(TYPE *Tab, TYPE size)
{ 
  int Wsize;

  TYPE i,i0,i1;  TYPE step;

  Wsize=size/(8*sizeof(TYPE));

  
  /*Moebius transform for high order bits, using word ops*/
  for (step=1;step<Wsize;step<<=1) {

    for (i1=0;i1<Wsize;i1+=2*step) {

      for (i0=0;i0<step;i0++) {
        i=i1+i0;

        Tab[i+step]^=Tab[i];
      }
    }
  }
  /*Moebius transform for low order bits, within words*/
  /* Assumes 8*sizeof(TYPE)=32 */
  for(i=0;i<Wsize;i++) {

    TYPE tmp;
    tmp=Tab[i];
    tmp^=(tmp<<16);

    tmp^=(tmp&0xff00ff)<<8;
    tmp^=(tmp&0xf0f0f0f)<<4;

    tmp^=(tmp&0x33333333)<<2;
    tmp^=(tmp&0x55555555)<<1;

    Tab[i]=tmp;
  }
}
