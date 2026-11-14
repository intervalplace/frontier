#include <stdio.h>
#include <stdlib.h>
int main(){ static long long v[200000]; int n=0; while(scanf("%lld",&v[n])==1) n++;
  __int128 best=0; int arg=-1; __int128 S=0; for(int i=0;i<n;i++) S+=v[i];
  for(int m=0;m<2*n-1;m++){ __int128 c=0; int lo=m-n+1>0?m-n+1:0, hi=m<n-1?m:n-1;
    for(int i=lo;i<=hi;i++) c+=(__int128)v[i]*v[m-i];
    if(c<0)c=-c; if(c>best){best=c;arg=m;} }
  // print 2n*best and S^2 as decimal strings
  char a[64],b[64]; int la=0,lb=0; __int128 x=(__int128)2*n*best, y=S*S;
  while(x>0){a[la++]='0'+(int)(x%10);x/=10;} while(y>0){b[lb++]='0'+(int)(y%10);y/=10;}
  printf("n %d argmax %d num ",n,arg); for(int i=la-1;i>=0;i--)putchar(a[i]); printf(" den "); for(int i=lb-1;i>=0;i--)putchar(b[i]); printf("\n"); }
