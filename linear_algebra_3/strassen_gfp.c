#include <stdio.h>
#include <stdlib.h>
#include <sys/times.h>
#include <unistd.h>

#define TYPE unsigned short
#define CUTOFF VarCutOff
#define MODUP VarModUp
#define MODULO 65521

#define access(M,i,j,size) M[(i)+((j)*(size))]

int VarCutOff=2;
int VarModUp=4;

double runtime ()
{
  struct tms t;
  long clk_tck;
  double ret;

  clk_tck=sysconf(_SC_CLK_TCK);
  times(&t);
  ret=t.tms_utime;
  ret/=clk_tck;
  return(ret);
}


void matmul(TYPE * A, TYPE * B, TYPE *Res, int size) {
  int i,j,k;
  unsigned int tmp;
  
  for(i=0;i<size;i++) 
    for(j=0;j<size;j++) {
      tmp=0;
      for(k=0;k<size;k++) {
	tmp=(tmp+access(A,i,k,size)*access(B,k,j,size))%MODULO;
      }
      access(Res,i,j,size)=tmp;
    }
}

void matmul_fewmod(TYPE * A, TYPE * B, TYPE *Res, int size) {
  int i,j,k;
  unsigned int tmp_lo,tmp_hi,mul;
  
  for(i=0;i<size;i++) 
    for(j=0;j<size;j++) {
      tmp_lo=tmp_hi=0;
      for(k=0;k<size;k++) {
	mul=access(A,i,k,size);
	mul*=access(B,k,j,size);
	tmp_lo+=mul&0xffff;
	tmp_hi+=mul>>16;
      }
      access(Res,i,j,size)=(((tmp_hi%MODULO)<<16)+(tmp_lo%MODULO))%MODULO;
    }
}


int mult_bloc(TYPE * A, TYPE * B, TYPE *Res, int true_size)
{
  int i,j;
  int hsize,size;
  
  TYPE *tmp1, *tmp2, *tmp3;
  int val;

  if (true_size<=CUTOFF){
    matmul_fewmod(A,B,Res,true_size);
  }
  else {
    for(i=0;i<true_size*true_size;i++) Res[i]=0;
    if (true_size==3) {
      hsize=1;
      size=2;
    }
    else if (true_size%(MODUP)==(MODUP-1)) {
      hsize=(true_size+1)/2;
      size=true_size;
    }
    else {
      hsize=true_size/2;
      size=2*hsize;
    }
    tmp1=(TYPE*)malloc((hsize*hsize)*sizeof(*tmp1));
    tmp2=(TYPE*)malloc((hsize*hsize)*sizeof(*tmp2));
    tmp3=(TYPE*)malloc((hsize*hsize)*sizeof(*tmp3));
    
    for(i=0;i<hsize;i++){
      for(j=0;j<hsize;j++){
	val=access(A,i,j,true_size);
	if ((i+hsize)<size) val+=access(A,i+hsize,j,true_size);
	if (val>=MODULO) val-=MODULO;
	access(tmp1,i,j,hsize)=val;
	
	val=access(B,i,j,true_size);
	if ((j+hsize)<size) val+=access(B,i,j+hsize,true_size);
	if (val>=MODULO) val-=MODULO;
	access(tmp2,i,j,hsize)=val;
      }
    }
    mult_bloc(tmp1,tmp2,tmp3,hsize);
    for(i=0;i<hsize;i++){
      for(j=0;j<hsize;j++){
	val=access(Res,i,j,true_size)+access(tmp3,i,j,hsize);
	if (val>=MODULO) val-=MODULO;
	access(Res,i,j,true_size)=val;
      }
    }
    
    
    for(i=0;i<hsize;i++){
      for(j=0;j<hsize;j++){
	if ((j+hsize)<size) {
	  val=access(A,i,j+hsize,true_size);
	  if ((i+hsize)<size) val+=access(A,i+hsize,j+hsize,true_size);
	  if (val>=MODULO) val-=MODULO;
	}
	else val=0;
	access(tmp1,i,j,hsize)=val;
	
	if ((i+hsize)<size) {
	  val=access(B,i+hsize,j,true_size);
	  if ((j+hsize)<size) val+=access(B,i+hsize,j+hsize,true_size);
	  if (val>=MODULO) val-=MODULO;
	}
	else val=0;
	access(tmp2,i,j,hsize)=val;
      }
    }
    mult_bloc(tmp1,tmp2,tmp3,hsize);
    for(i=0;i<size-hsize;i++){
      for(j=0;j<size-hsize;j++){
	val=access(Res,i+hsize,j+hsize,true_size);
	val+=access(tmp3,i,j,hsize);
	if (val>=MODULO) val-=MODULO;
	access(Res,i+hsize,j+hsize,true_size)=val;
      }
    }
    
    for(i=0;i<hsize;i++){
      for(j=0;j<hsize;j++){
	if ((j+hsize)<size) 
	  val=access(A,i,j+hsize,true_size);
	else
	  val=0;
	if ((i+hsize)<size)
	  val+=access(A,i+hsize,j,true_size);
	if (val>=MODULO) val-=MODULO;
	access(tmp1,i,j,hsize)=val;
	
	if ((j+hsize)<size) 
	  val=-access(B,i,j+hsize,true_size);
	else
	  val=0;
	if ((i+hsize)<size) 
	  val+=access(B,i+hsize,j,true_size);
	if (val<0) val+=MODULO;
	access(tmp2,i,j,hsize)=val;
      }
    }
    mult_bloc(tmp1,tmp2,tmp3,hsize);
    for(i=0;i<hsize;i++){
      for(j=0;j<hsize;j++){
	val=access(Res,i,j,true_size);
	val+=access(tmp3,i,j,hsize);
	if (val>=MODULO) val-=MODULO;
	access(Res,i,j,true_size)=val;
      }
    }

    for(i=0;i<size-hsize;i++){
      for(j=0;j<size-hsize;j++){
	val=access(Res,i+hsize,j+hsize,true_size);
	val-=access(tmp3,i,j,hsize);
	if (val<0) val+=MODULO;
	access(Res,i+hsize,j+hsize,true_size)=val;
      }
    }
    
    for(i=0;i<hsize;i++){
      for(j=0;j<hsize;j++){
	val=access(B,i,j,true_size);
	if ((i+hsize)<size) {
	  access(tmp1,i,j,hsize)=access(A,i+hsize,j,true_size);
	  val+=access(B,i+hsize,j,true_size);
	}
	else
	  access(tmp1,i,j,hsize)=0;
	if (val>=MODULO) val-=MODULO;
	access(tmp2,i,j,hsize)=val;
      }
    }
    mult_bloc(tmp1,tmp2,tmp3,hsize);
    for(i=0;i<hsize;i++){
      for(j=0;j<hsize;j++){
	val=access(Res,i,j,true_size);
	val-=access(tmp3,i,j,hsize);
	if (val<0) val+=MODULO;
	access(Res,i,j,true_size)=val;
	if ((i+hsize)<size) {
	  val=access(Res,i+hsize,j,true_size);
	  val+=access(tmp3,i,j,hsize);
	  if (val>=MODULO) val-=MODULO;
	  access(Res,i+hsize,j,true_size)=val;
	}
      }
    }
    
    for(i=0;i<hsize;i++){
      for(j=0;j<hsize;j++){
	val=access(A,i,j,true_size);
	if ((j+hsize)<size)
	  val-=access(A,i,j+hsize,true_size);
	if (val<0) val+=MODULO;
	access(tmp1,i,j,hsize)=val;
	if ((j+hsize)<size)
	  access(tmp2,i,j,hsize)=access(B,i,j+hsize,true_size);
	else
	  access(tmp2,i,j,hsize)=0;
      }
    }
    mult_bloc(tmp1,tmp2,tmp3,hsize);
    for(i=0;i<hsize;i++){
      for(j=0;j<hsize;j++){
	val=access(Res,i,j,true_size);
	val-=access(tmp3,i,j,hsize);
	if (val<0) val+=MODULO;
	access(Res,i,j,true_size)=val;
	if ((j+hsize)<size){
	  val=access(Res,i,j+hsize,true_size);
	  val+=access(tmp3,i,j,hsize);
	  if (val>=MODULO) val-=MODULO;
	  access(Res,i,j+hsize,true_size)=val;
	}
      }
    }
    
    
    for(i=0;i<hsize;i++){
      for(j=0;j<hsize;j++){
	if ((j+hsize)<size) {
	  access(tmp1,i,j,hsize)=access(A,i,j+hsize,true_size);
	  val=access(B,i,j+hsize,true_size);
	  if ((i+hsize)<size) val+=access(B,i+hsize,j+hsize,true_size);
	  if (val>=MODULO) val-=MODULO;
	  access(tmp2,i,j,hsize)=val;	
	}
	else {
	  access(tmp1,i,j,hsize)=0;
	  access(tmp2,i,j,hsize)=0;
	}


      }
    }
    mult_bloc(tmp1,tmp2,tmp3,hsize);
    for(i=0;i<hsize;i++){
      for(j=0;j<size-hsize;j++){
	if ((i+hsize)<size) {
	  val=access(Res,i+hsize,j+hsize,true_size);
	  val-=access(tmp3,i,j,hsize);
	  if (val<0) val+=MODULO;
	  access(Res,i+hsize,j+hsize,true_size)=val;
	}
	val=access(Res,i,j+hsize,true_size);
	val+=access(tmp3,i,j,hsize);
	if (val>=MODULO) val-=MODULO;
	access(Res,i,j+hsize,true_size)=val;
      }
    }
    
    for(i=0;i<hsize;i++){
      for(j=0;j<hsize;j++){
	if ((i+hsize)<size) {
	  val=-access(A,i+hsize,j,true_size);
	  if ((j+hsize)<size)
	    val+=access(A,i+hsize,j+hsize,true_size);
	  if (val<0) val+=MODULO;
	  access(tmp1,i,j,hsize)=val;
	  access(tmp2,i,j,hsize)=access(B,i+hsize,j,true_size);
	}
	else {
	  access(tmp1,i,j,hsize)=0;
	  access(tmp2,i,j,hsize)=0;
	}
      }
    }
    mult_bloc(tmp1,tmp2,tmp3,hsize);
    for(i=0;i<size-hsize;i++){
      for(j=0;j<hsize;j++){
	if ((j+hsize)<size) {
	  val=access(Res,i+hsize,j+hsize,true_size);
	  val-=access(tmp3,i,j,hsize);
	  if (val<0) val+=MODULO;
	  access(Res,i+hsize,j+hsize,true_size)=val;
	}
	val=access(Res,i+hsize,j,true_size);
	val+=access(tmp3,i,j,hsize);
	if (val>=MODULO) val-=MODULO;
	access(Res,i+hsize,j,true_size)=val;
      }
    }
    free(tmp1); free(tmp2); free(tmp3);

    if (size==(true_size-1)) {
      unsigned int mul;

      for (i=0;i<size;i++) {
	for (j=0;j<size;j++) {
	  mul=access(A,i,true_size-1,true_size)*access(B,true_size-1,j,true_size);
	  mul%=MODULO;
	  val=access(Res,i,j,true_size)+mul;
	  if (val>=MODULO) val-=MODULO;
	  access(Res,i,j,true_size)=val;
	}
      }

      for (i=0;i<true_size-1;i++) {
	for (j=0;j<true_size;j++) {
	  mul=access(A,i,j,true_size)*access(B,j,true_size-1,true_size);
	  mul%=MODULO;
	  val=access(Res,i,true_size-1,true_size)+mul;
	  if (val>=MODULO) val-=MODULO;
	  access(Res,i,true_size-1,true_size)=val;

	  mul=access(A,true_size-1,j,true_size)*access(B,j,i,true_size);
	  mul%=MODULO;
	  val=access(Res,true_size-1,i,true_size)+mul;
	  if (val>=MODULO) val-=MODULO;
	  access(Res,true_size-1,i,true_size)=val;
	}
      }
      for (i=0;i<true_size;i++) {
	mul=access(A,true_size-1,i,true_size)*access(B,i,true_size-1,true_size);
	mul%=MODULO;
	val=access(Res,true_size-1,true_size-1,true_size)+mul;
	if (val>=MODULO) val-=MODULO;
	access(Res,true_size-1,true_size-1,true_size)=val;
      }
    }
  }
}

void printmat(TYPE *mat, int size) 
{
  int i,j;
  
  for(i=0;i<size;i++){
    for(j=0;j<size;j++)
      printf("%d ",access(mat,i,j,size));
    printf("\n");
  }
}

#define MAXSIZE 300
main()
{
  int i,j,k,sum;
  double timer;

  int SIZE;
  TYPE *m1,*m2,*m3,*m4;

  for(SIZE=128;SIZE<MAXSIZE;SIZE++) {
    m1=malloc(SIZE*SIZE*sizeof(*m1));
    m2=malloc(SIZE*SIZE*sizeof(*m2));
    m3=malloc(SIZE*SIZE*sizeof(*m3));
    m4=malloc(SIZE*SIZE*sizeof(*m4));
    
    for(i=0;i<SIZE;i++)
      for(j=0;j<SIZE;j++) {
	access(m1,i,j,SIZE)=random()%MODULO;
	access(m2,i,j,SIZE)=random()%MODULO;
      }

#if 1
    for (VarCutOff=32;VarCutOff<=128;VarCutOff<<=1) {
      for (VarModUp=4;VarModUp<=32;VarModUp<<=1) {
	timer=runtime();
	for(i=0;i<100;i++)
	  mult_bloc(m1,m2,m3,SIZE);
	timer=runtime()-timer;
	printf("%d %0.2f %d %d\n",SIZE,timer,VarCutOff,VarModUp);
	fflush(stdout);
      }
    }
#endif
    
    timer=runtime();
    for(i=0;i<100;i++)
      matmul(m1,m2,m4,SIZE);
    timer=runtime()-timer;
    printf("%d %0.2f 0 0\n",SIZE,timer);
    fflush(stdout);

    timer=runtime();
    for(i=0;i<100;i++)
      matmul_fewmod(m1,m2,m4,SIZE);
    timer=runtime()-timer;
    printf("%d %0.2f 0 1\n",SIZE,timer);
    fflush(stdout);

#if 1
    for(i=0;i<SIZE;i++)
      for(j=0;j<SIZE;j++) {
	if (access(m3,i,j,SIZE)!=access(m4,i,j,SIZE)) {
	  printf("ARGGGG (%d %d)-> %d, %d !!\n",i,j,access(m3,i,j,SIZE),access(m4,i,j,SIZE));
	}
      }
#endif
    free(m1);
    free(m2);
    free(m3);
    free(m4);

  }
}
