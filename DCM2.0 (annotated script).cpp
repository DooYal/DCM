#include"stdio.h"
#include"stdlib.h"
#include"malloc.h"
#define LEN sizeof(s)
struct s
{
	char type[6];
int times;  //repeat time of a SSR
int length; //size (bp) of a SSR
int start; //Genomic locational start of a SSR
int end; //Genomic locational end of a SSR
struct s *next;
	};
int n;
struct s *creat(FILE*fp)
{struct s *head;
struct s *p1,*p2;
n=0;
p1=p2=(struct s*)malloc(LEN);
fscanf(fp,"%s %d %d %d %d %d",&(p1->type),&(p1->times),&(p1->length),&(p1->start),&(p1->end)); //load the data from the SSR extraction file
head=NULL;
while((p1->type!=0)&&(!feof(fp)))
{n=n+1;
if(n==1)
head=p1;
else p2->next=p1;
p2=p1;
p1=(struct s *)malloc(LEN);
fscanf(fp,"%s %d %d %d %d %d",&(p1->type),&(p1->times),&(p1->length),&(p1->start),&(p1->end));
}
p2->next=NULL;
return(head);
}
struct txtname
{
	char name[10];
}filename[100]={
{"human y ssr.txt"}, //read the output of IMEx file, name.txt
}
,*pname=filename;

int main()
{
int j;
int f=1;
int f1=0;
struct s* pshu;
FILE *fp2;

struct filestore
{char c[30];
}filestorename[2]={{"result.txt"}},*filestorep;  //generate output file
filestorep=filestorename;
int n,g,k,z;
int position=0;
int ex_num;
int ex_den;
n=0;
printf("Please input the largest number of genome\n");
scanf("%d",&z);
printf("Please input the window level (solution)(100-4500)\n");  //set the size of differential unit
scanf("%d",&k);
g=z/k;  //differential the sample genomic sequences
if((fp2=fopen(filestorep->c,"a"))==0)
{
	printf("\nThe file result.txt are not established!!!\n");
}
fprintf(fp2,"Type\tposition\tratio_mono\tratio_di\tratio_tri\tratio_tetra\tratio_penta\tratio_hexa\tratio_total\tratio_len_mono\tratio_len_di\tratio_len_tri\tratio_len_tetra\tratio_len_penta\tratio_len_hexa\tratio_len_total\tex_num\tex_den\n");
fclose(fp2);
for(j=0;;)	
{FILE *fp1;
if((fp1=fopen(pname->name,"r"))==0)
{printf("\nThe files can not be read or have all been read off!!!\n");
int e;
scanf("%d",&e);}
int a;
for (a=0;a<2;a++) fscanf(fp1,"%*[^\n]%*c");
pshu=creat(fp1);
fclose(fp1);
struct outputs //set the output data structure, which includes the number and length of mono- to hexanucleotide repeats
{
int mono;
int di;
int tri;
int tetra;
int penta;
int hexa;
int total;
int len;
int mono_len;
int di_len;
int tri_len;
int tetra_len;
int penta_len;
int hexa_len;
int total_len;
int ex_num;
int ex_den;
}
ab[350],*p=ab;
for(n=0;n<g;n++)
{int m=0,d=0,t=0,te=0,pe=0,h=0,len=0,mono_len=0,di_len=0,tri_len=0,tetra_len=0,penta_len=0,hexa_len=0,total_len=0;

	while((pshu->start)>(k*n)&&(pshu->start)<=(k*(n+1))&&pshu->next!=NULL)
{

switch((pshu->length)/(pshu->times))   //calculated the SSR length per differential unit
{
case 1: m++; len=pshu->length; mono_len=len+mono_len; break;
case 2: d++; len=pshu->length; di_len=len+di_len; break;
case 3: t++; len=pshu->length; tri_len=len+tri_len; break;
case 4: te++; len=pshu->length; tetra_len=len+tetra_len; break;
case 5: pe++; len=pshu->length; penta_len=len+penta_len; break;
case 6: h++; len=pshu->length; hexa_len=len+hexa_len; break;
default:printf("default");
}
pshu=pshu->next;
}
p->mono=m;
p->di=d;
p->tri=t;
p->tetra=te;
p->penta=pe;
p->hexa=h;
p->total=(m+d+t+te+pe+h);
p->mono_len=mono_len;
p->di_len=di_len;
p->tri_len=tri_len;
p->tetra_len=tetra_len;
p->penta_len=penta_len;
p->hexa_len=hexa_len;
p->total_len=(mono_len+di_len+tri_len+tetra_len+penta_len+hexa_len);
p++;}
p=ab;
printf("%d*******************************************\n",f);
f1++;
filestorep=filestorename;	
for(n=0;n<g;n++)  //Set the extra outputting columns for the map drawing process next 
{
    position=k+position;
	ex_num=100;
	ex_den=600; 
	printf("position=%8d\tmono=%3d\tdi=%3d\ttri=%3d\ttetra=%3d\tpenta=%3d\thexa=%3d\ttotal=%3d\tmono_len=%4d\tdi_len=%4d\ttri_len=%4d\ttetra_len=%4d\tpenta_len=%4d\thexa_len=%4d\ttotal_len=%4d\tex_num=3%d\tex_den=3%d\n",position,p->mono,p->di,p->tri,p->tetra,p->penta,p->hexa,p->total,p->mono_len,p->di_len,p->tri_len,p->tetra_len,p->penta_len,p->hexa_len,p->total_len,ex_num,ex_den);
	if((fp2=fopen(filestorep->c,"a"))==0)
{
	printf("\nresult.txt are not created!!!\n");
}
	fprintf(fp2,"%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",f1,position,p->mono,p->di,p->tri,p->tetra,p->penta,p->hexa,p->total,p->mono_len,p->di_len,p->tri_len,p->tetra_len,p->penta_len,p->hexa_len,p->total_len,ex_num,ex_den);
	fclose(fp2);	
p++;}
pname++;
f++;
}
char y;
printf("\n the text has been written!!!\n");
scanf("%d",&y);
}
