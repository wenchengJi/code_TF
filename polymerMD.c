//%%%%%%%%%%%%      polymer overdamped 2022_7_17    %%%%%%%%%%%%%
// Developed by Wencheng Ji
//Contact: wencheng.ji_at_weizmann.ac.il



#include <string.h>
#include <time.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>


//#define	N							100
//#define	afactor						sqrt(3.0)
#define    TWO_PI                    6.28318530717958
#define    UNI        ((double)rand()/((double)RAND_MAX + 1.0))
#define    TFd                   25
#define    acs                   1e-10


int N,k,naccount,nbound1,k_flag,nn,F1_label,mm,tupup;
double DOUBLE_N,DOUBLE_TFd,Rcut, Double_d;
double dt,tt,Keff,tem_fmax;

double rho_s,rho_t;
double wds,wd0,Ewca;
double x_mean,y_mean,z_mean;


void calculateForces();
void Binding_initial();
void initial_run();

double *rx;		/*	x component of position	*/
double *ry;		/*	y component of position	*/
double *rz;		/*	y component of position	*/
double *fx;		/*	x component of position	*/
double *fy;		/*	y component of position	*/
double *fz;		/*	y component of position	*/

double *fx_bead, *fy_bead, *fz_bead,*fx_tem, *fy_tem, *fz_tem;

double *rxP;
int ns,Neff,klabel;
double RB,wd,wd0,dt_f,EB;
int d;





void Binding_initial(){
    int j;
    if (ns>0){
    for(j=0;j<ns;j++){
        rxP[j]=(-ns/2.0+0.5+j)*Double_d;}}
    return;
}


void initial_run(){
    double xa,ya,za,la,x1tem,y1tem,z1tem,a2;
    int i,j;
    x1tem=0.0; y1tem=0.0; z1tem=0.0;
    
    
    xa=UNI-0.5;
    ya=UNI-0.5;
    za=UNI-0.5;
    la=sqrt(xa*xa+ya*ya+za*za);
    xa=xa/la;
    ya=ya/la;
    za=za/la;
    
    for (i=0;i<N;i++){
    rx[i]=(i-N/2)*sqrt(DOUBLE_TFd)*xa+0.1*(UNI-0.5)*xa;
    ry[i]=(i-N/2)*sqrt(DOUBLE_TFd)*ya+0.1*(UNI-0.5)*ya;
    rz[i]=(i-N/2)*sqrt(DOUBLE_TFd)*za+0.1*(UNI-0.5)*za;
 //  printf("rx0=%.5f, ry0=%.5f, rz0=%.5f\n", rx[i],ry[i],rz[i]);
        
    }

    a2=sqrt(2*dt/DOUBLE_TFd);
    if (N>1){
        for (j=0;j<Neff*Neff/dt/6.0;j++){
            calculateForces();
            for (i=0;i<N;i++){
             rx[i]=rx[i]+1.0/DOUBLE_TFd*fx[i]*dt+a2*sqrt(-2*log(UNI+acs))*cos(UNI*TWO_PI);
             ry[i]=ry[i]+1.0/DOUBLE_TFd*fy[i]*dt+a2*sqrt(-2*log(UNI+acs))*cos(UNI*TWO_PI);
             rz[i]=rz[i]+1.0/DOUBLE_TFd*fz[i]*dt+a2*sqrt(-2*log(UNI+acs))*cos(UNI*TWO_PI);}
            
        //printf("j=%d, tt=%.2f. fx[0]=%.5f,%.5f\n",j,dt*j,fx[0],sqrt(-2*log(UNI))*cos(UNI*TWO_PI));
            
        }}
    x_mean=0.0;y_mean=0.0;z_mean=0.0;
    for (i=0;i<N;i++){x_mean+=rx[i];y_mean+=ry[i];z_mean+=rz[i];}
    srand(k+10000);
  

    while (1){
        
    x1tem=(UNI-0.5)*2*RB;
    y1tem=(UNI-0.5)*2*RB;
    z1tem=(UNI-0.5)*2*RB;

        if (x1tem*x1tem+y1tem*y1tem+z1tem*z1tem<RB*RB){
            for (i=0;i<N;i++){
                rx[i]=rx[i]-x_mean/DOUBLE_N+x1tem;
                ry[i]=ry[i]-y_mean/DOUBLE_N+y1tem;
                rz[i]=rz[i]-z_mean/DOUBLE_N+z1tem;
            }
                
    printf("rxc=%.5f, ryc=%.5f, rzc=%.5f\n",x1tem,y1tem,z1tem);

                break;}}
    
    
   return;
}
  

 
    
    
    
    


void calculateForces(){
    double rr,rd_low,tem,tem_x,r_dna_min,tem_x_max,tem_x_min,tem_max_rxP;
    double f_tem_intpart,tme_rwca,temp,tem_rq6_inv,tem_rq;
    int i,j;
    
        /////exclude interaction
         for (i=0;i<N;i++){
             fx_bead[i]=0.0;  fy_bead[i]=0.0; fz_bead[i]=0.0;}
        
             if (N>1 && F1_label>0){
               for (i=0;i<N-1;i++){
                   for (j=i+1;j<N;j++){
                   tem_rq=(rx[i]-rx[j])*(rx[i]-rx[j])+(ry[i]-ry[j])*(ry[i]-ry[j])+(rz[i]-rz[j])*(rz[i]-rz[j]);
                   if (tem_rq<wd*wd){
                   tem_rq6_inv=1.0/(tem_rq*tem_rq*tem_rq);
                  
                   tme_rwca=(tem_rq6_inv*tem_rq6_inv-tem_rq6_inv)/tem_rq;
                           
                   temp=Ewca*12*(tme_rwca*(rx[i]-rx[j]));
                   fx_bead[i]+=temp;
                   fx_bead[j]-=temp;
                           
                   temp=Ewca*12*(tme_rwca*(ry[i]-ry[j]));
                   fy_bead[i]+=temp;
                   fy_bead[j]-=temp;
                           
                   temp=Ewca*12*(tme_rwca*(rz[i]-rz[j]));
                   fz_bead[i]+=temp;
                   fz_bead[j]-=temp;
                
                   }}}}
            
        /////%%%%%%%%%%%%%
        
        
        // binding force
        for (i=0;i<N;i++){
            fx_tem[i]=0.0; fy_tem[i]=0.0;  fz_tem[i]=0.0;}
        
        naccount=1;
         
        if (ns>0&& F1_label>0){
         for (i=0;i<N;i++){
             for (j=0;j<ns;j++){
                 
                   temp=0.0;
                   temp=(rx[i]-rxP[j])*(rx[i]-rxP[j])+ry[i]*ry[i]+ rz[i]*rz[i];
                 if (temp<Rcut){
                   f_tem_intpart=0.0;
                   f_tem_intpart=-EB*exp(-temp/(2*wds))/wds;
                   fx_tem[i]=f_tem_intpart*(rx[i]-rxP[j]);
                   fy_tem[i]=f_tem_intpart*ry[i];
                   fz_tem[i]=f_tem_intpart*rz[i];
                 
                     if ((rx[i]-rxP[j])*(rx[i]-rxP[j])+ry[i]*ry[i]+rz[i]*rz[i]<wds){
                         naccount=naccount+1;}}}}}
            
        nbound1=naccount-1;
          
    for (i=0;i<N;i++){
        fx[i]=0.0;
        fy[i]=0.0;
        fz[i]=0.0;}
       
    
    
    if (N>1){
    for (i=0;i<N;i++){
        if (i==0){
            fx[i]=Keff*(rx[i+1]-rx[i])+fx_tem[i]+fx_bead[i];
            fy[i]=Keff*(ry[i+1]-ry[i])+fy_tem[i]+fy_bead[i];
            fz[i]=Keff*(rz[i+1]-rz[i])+fz_tem[i]+fz_bead[i];}
        
        if (i==N-1){
            fx[i]=Keff*(-rx[i]+rx[i-1])+fx_tem[i]+fx_bead[i];
            fy[i]=Keff*(-ry[i]+ry[i-1])+fy_tem[i]+fy_bead[i];
            fz[i]=Keff*(-rz[i]+rz[i-1])+fz_tem[i]+fz_bead[i];}
        
        if (i<N-1&i>0){
            fx[i]=Keff*(rx[i+1]+rx[i-1]-2*rx[i])+fx_tem[i]+fx_bead[i];
            fy[i]=Keff*(ry[i+1]+ry[i-1]-2*ry[i])+fy_tem[i]+fy_bead[i];
            fz[i]=Keff*(rz[i+1]+rz[i-1]-2*rz[i])+fz_tem[i]+fz_bead[i];}
       }}
    else{i=0;
        fx[i]=fx_tem[i];
        fy[i]=fy_tem[i];
        fz[i]=fz_tem[i];}
        
        
    return;

        
}




void run_MD(){
    double a2,tem,tem_x_max,tem_x_min,r_dna_min,tem_max_rxP,f_tem_intpart;
    double NRc,rand_x,rand_y,rand_z,tem_fmax_bead,dt_f_tem;
    int i;
    r_dna_min=ry[0]*ry[0]+rz[0]*rz[0];
    tem_x_max=fabs(rx[0]);
    tem_x_min=fabs(rx[0]);
    
    for (i=0;i<N;i++){
        tem=ry[i]*ry[i]+rz[i]*rz[i];
        if (tem<r_dna_min) {r_dna_min=tem;}
        tem=fabs(rx[i]);
        if (tem<tem_x_min) {tem_x_min=tem;}
        if (tem>tem_x_max) {tem_x_max=tem;}}
        r_dna_min=sqrt(r_dna_min);
    
    if (ns>1) {tem_max_rxP=rxP[ns-1];}
    else{tem_max_rxP=0.0;}
    
    if ( r_dna_min<10.0 && (tem_x_min<tem_max_rxP+10.0 | tem_x_max< tem_max_rxP+10.0)){
     
        calculateForces();
        tem_fmax=fx[0]*fx[0]+fy[0]*fy[0]+fz[0]*fz[0];
        for (i=0;i<N;i++){tem=fx[i]*fx[i]+fy[i]*fy[i]+fz[i]*fz[i];
            if (tem> tem_fmax){tem_fmax=tem;}}
        tem_fmax=sqrt(tem_fmax);
      
       // tem_fmax_bead=0;
        //for (i=0;i<N;i++){tem=fx_bead[i]*fx_bead[i]+fy_bead[i]*fy_bead[i]+fz_bead[i]*fz_bead[i];
         //   if (tem>tem_fmax_bead){tem_fmax_bead=tem;}}
        //tem_fmax_bead=sqrt(tem_fmax_bead);
        
        
        if (r_dna_min<10.0 && (tem_x_min<10.0 | tem_x_max<10.0)){
             dt_f_tem=0.001;}
        else {dt_f_tem=dt_f;}
         
        
        if (tem_fmax>1){dt=dt_f_tem/tem_fmax*DOUBLE_TFd;}
        else {dt=dt_f_tem*DOUBLE_TFd;}
        
        a2=sqrt(2*dt/DOUBLE_TFd);
         
        for (i=0;i<N;i++){
           rx[i]+=1.0/DOUBLE_TFd*fx[i]*dt+a2*sqrt(-2*log(UNI+acs))*cos(UNI*TWO_PI);
           ry[i]+=1.0/DOUBLE_TFd*fy[i]*dt+a2*sqrt(-2*log(UNI+acs))*cos(UNI*TWO_PI);
           rz[i]+=1.0/DOUBLE_TFd*fz[i]*dt+a2*sqrt(-2*log(UNI+acs))*cos(UNI*TWO_PI);}
            //if (nn%1000000==0){
            //printf("%.5f %.5f %.5f, f_max=%.10f,f_max_bead=%.10f\n",rx[0],ry[0],rz[0],tem_fmax,tem_fmax_bead);}
    }
    else
    {if (r_dna_min<200 && tem_x_min>-tem_max_rxP-200.0 && tem_x_max<tem_max_rxP+200.0) {dt=0.1*Neff;}
    else {dt=10*Neff;}
    
        rand_x=sqrt(-2*log(UNI+acs))*cos(UNI*TWO_PI);
        rand_y=sqrt(-2*log(UNI+acs))*cos(UNI*TWO_PI);
        rand_z=sqrt(-2*log(UNI+acs))*cos(UNI*TWO_PI);
       
        for (i=0;i<N;i++){
            rx[i]+=sqrt(2*dt/Neff)*rand_x;
            ry[i]+=sqrt(2*dt/Neff)*rand_y;
            rz[i]+=sqrt(2*dt/Neff)*rand_z;}
                
        // center of mass
        x_mean=0.0;y_mean=0.0;z_mean=0.0;
        for (i=0;i<N;i++){x_mean+=rx[i];y_mean+=ry[i];z_mean+=rz[i];}
        x_mean=x_mean/DOUBLE_N;
        y_mean=y_mean/DOUBLE_N;
        z_mean=z_mean/DOUBLE_N;
        // end enter of mass

        NRc=sqrt(x_mean*x_mean+y_mean*y_mean+z_mean*z_mean);
        
        if (NRc>RB){
            for (i=0;i<N;i++){
              rx[i]=2*x_mean/NRc*RB-rx[i];
              ry[i]=2*y_mean/NRc*RB-ry[i];
              rz[i]=2*z_mean/NRc*RB-rz[i];}
            
          //  printf("%.5f\n",sqrt(rx[0]*rx[0]+ry[0]*ry[0]+rz[0]*rz[0]));
         }
        
        
    }
    
    // center of mass
    x_mean=0.0;y_mean=0.0;z_mean=0.0;
    for (i=0;i<N;i++){x_mean+=rx[i];y_mean+=ry[i];z_mean+=rz[i];}
    x_mean=x_mean/DOUBLE_N;
    y_mean=y_mean/DOUBLE_N;
    z_mean=z_mean/DOUBLE_N;
    // end enter of mass
    

    return;


}



void saveCon(char *snapFileName){
    int i;
    FILE *snapFile;
    snapFile = fopen(snapFileName,"wb");
    fprintf(snapFile,"%d\t%.15g\t%d\n",nn,tt,mm);
    for (i=0; i<N; i++)
    fprintf(snapFile,"%.15g\t%.15g\t%.15g\n",rx[i],ry[i],rz[i]);
    fclose(snapFile);
    return;
}



int readCon(char *snapFileName){
    int i;
    FILE *snapFile;
    snapFile = fopen(snapFileName,"rb");
    if ( snapFile ){
    fscanf(snapFile,"%d  %lf %d\n",&nn,&tt,&mm);
    for (i=0; i<N; i++)
    fscanf(snapFile,"%lf %lf %lf",&(rx[i]),&(ry[i]),&(rz[i]));
    fclose(snapFile);
    return 1;}
    return 0;
}




















void allocate_memory(){
 
    rx = (double *)malloc(sizeof(double)*N);
    ry = (double *)malloc(sizeof(double)*N);
    rz = (double *)malloc(sizeof(double)*N);
    
    rxP = (double *)malloc(sizeof(double)*ns);

    fx = (double *)malloc(sizeof(double)*N);
    fy = (double *)malloc(sizeof(double)*N);
    fz = (double *)malloc(sizeof(double)*N);
    fx_tem = (double *)malloc(sizeof(double)*N);
    fy_tem = (double *)malloc(sizeof(double)*N);
    fz_tem = (double *)malloc(sizeof(double)*N);
    
    fx_bead = (double *)malloc(sizeof(double)*N);
    fy_bead = (double *)malloc(sizeof(double)*N);
    fz_bead = (double *)malloc(sizeof(double)*N);

    return;
}




void free_everything(){
    int i;
    free(rx);       free(ry);       free(rz);
    free(fx);      free(fy);      free(fz);
    free(fx_tem);      free(fy_tem);      free(fz_tem);
    free(fx_bead);      free(fy_bead);      free(fz_bead);
    free(rxP);
   
    return;
}




int main(int argc,char *argv[]){
    //int main(int n,char **inputStrings){
    int q,numOfSnaps,t0,t1,i;
    char snapFileName[1024];
    double rxC,tem0,tem,r_to_DNA;
    FILE *File1;
    t0 = time(0);
   
    EB=atof(argv[1]);
    N=atoi(argv[2]);
    double LL=atof(argv[3]);
   
    k=atoi(argv[4]);
    
    
    Neff=N*TFd;//atoi(argv[1]);
    d=10;//atoi(argv[2]);
    wd=1;//atof(argv[3]);
    


    ns=(int) LL/d;//atoi(argv[6]);
    RB=500;//atof(argv[7]);
    wd0=0.5;//atof(argv[8]);
    klabel=0;//atof(argv[9]);
    dt_f=0.01;//atof(argv[10]);
    
    
    DOUBLE_N=(double) N;
    DOUBLE_TFd=(double) TFd;
    Double_d=(double) d;
    Keff=3.0/DOUBLE_TFd;
    
    ///%%%%%%potential radius
    wds=wd*wd;
    Rcut=16*wds;
    Ewca=1.0;
    ////%%%%%%%%%%%%%%%%%%%%
 
    dt=0.05;
    allocate_memory();
    srand(k);
    Binding_initial();
    
    F1_label=0;
    srand(k);
    initial_run();
    nn=0; tt=0.0; mm=0;
    F1_label=1;
    
    if (F1_label ==1){
    sprintf(snapFileName,"Search_n=%d_EB=%.2f_L=%.1f_%d.dat",N,EB,LL,k);
    File1 = fopen(snapFileName,"ab");
        x_mean=0.0;y_mean=0.0;z_mean=0.0;
        for (i=0;i<N;i++){x_mean+=rx[i];y_mean+=ry[i];z_mean+=rz[i];}
        x_mean=x_mean/DOUBLE_N;
        y_mean=y_mean/DOUBLE_N;
        z_mean=z_mean/DOUBLE_N;
    
        printf("Xc=%.3f,Yc=%.3f,Zc=%.3f\n", x_mean, y_mean,z_mean);
        for (i=0;i<N;i++){printf("i=%d, rx0=%.5f, ry0=%.5f, rz0=%.5f\n", i,rx[i],ry[i],rz[i]);}
         
        
    while(1){
        /// please be careful of stopping condition!
        //for (i=0;i<N;i++){
       //  if (rx[i]*rx[i]+ry[i]*ry[i]+rz[i]*rz[i]<wd0*wd0){
        if (rx[0]*rx[0]+ry[0]*ry[0]+rz[0]*rz[0]<wd0*wd0){

            tem=0.0;tem0=0.0;
            for (i=0;i<N;i++){
            tem+=ry[i]*ry[i]+rz[i]*rz[i];
            tem0+=rx[i];}
            r_to_DNA=sqrt(tem);
            rxC=tem0/DOUBLE_N;
            fprintf(File1,"%2.3e %d %.3f\n",tt,nbound1,rxC);
            printf("%2.3e %d %.3f\n",tt,nbound1,rxC);
             break;}
        run_MD();
        nn=nn+1;
        tt=tt+dt;
        
        if (tt*1e-5>mm){
            mm=mm+1;
        tem=0.0;tem0=0.0;
        for (i=0;i<N;i++){
        tem+=ry[i]*ry[i]+rz[i]*rz[i];
        tem0+=rx[i];}
        r_to_DNA=sqrt(tem);
        rxC=tem0/DOUBLE_N;
        fprintf(File1,"%2.3e %d %.3f\n",tt,nbound1,rxC);
        printf("%2.3e %d %.3f\n",tt,nbound1,rxC);
}}
    fclose(File1);
    }
    
    

    
    free_everything();
    t1 = time(0);
    printf("time in seconds: %d\n\n\n", t1-t0);
    return 0;
}


