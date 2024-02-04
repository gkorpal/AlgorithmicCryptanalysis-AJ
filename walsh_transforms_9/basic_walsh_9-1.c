// Program 9.1: Basic Walsh transform code

#define TYPE int

/*In place Walsh transform*/
void Walsh(TYPE *Tab, TYPE size)
{ 
  TYPE i,i0,i1;  TYPE step;

  TYPE sum,diff;
  for (step=1;step<size;step<<=1) {

    for (i1=0;i1<size;i1+=2*step) {

      for (i0=0;i0<step;i0++) {
        i=i1+i0;

        sum=Tab[i]+Tab[i+step];diff=Tab[i]-Tab[i+step];

        Tab[i]=sum;Tab[i+step]=diff;
      }
    }
  }
}

