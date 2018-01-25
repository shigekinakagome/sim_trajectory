// This program generates frequency trajectories of an allele that was observed to have frequency = nder/nsam.
// It assumes that the mutation arose once, at a random time in the past, and no further mutation occurs.
// The command line arguments are "nreps ageln h nsam nder ft seed".
// nreps: number of trajectories to generate.
// ageln: the parameter of prior distribution for onset of selection on standing variation which was present at ft.
// h: dominance coefficient.
// nsam: sample size.
// nder: number of copies of the derived allele found in the sample of size nsam.
// ft: frequency of standing mutation at t when selection happened.
// The function popsize(j), returns the diploid population size at time j generations in the past.
// The user must supply this function. Examples are provided.
// A Wright-Fisher model is assumed.
// The first line of the output is nreps.
// The following line contains "npoints: n1 (where n1 is the number of time points until present)."
// Following that are n1 lines where each line contains two numbers, the time point(=generation/(4*N)),and the frequency.
// Time is measured in units of 4N generations.
// Compilation: "gcc -o traj -lm Trajectory_SSV.c popsize.c binomial_mt.c mt.c"
// usage: "./traj nreps ageln h nsam nder ft seed | ./stepftn > XXX.out"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

extern int npres;
extern double r;

#define MAX 960000
#define MAXLEN 1000

#define SFSsamp 1000000
#define PERIODS 4
#define t1 500
#define t2 1100
#define t3 1600

main(int argc,char **argv){

    FILE *fp;

    int nreps,gen,ngen,genfwd,genbwd,i,popmax,nsam,nder,count,nsuccess,seed;
    int flag_fwd,flag_bwd,ran,flag;
    double s,s_tmp,h,p,pfwd,pbwd,wbar,pprime,*freq,*freqsel,*freqneu,pmax,prob,ft;
    double ngentmp,ngenmu,ageln,alpha;

    int popsize(int);
    int GetRandom(int min,int max);
    float bnldev(float,int);
    double gasdev(void);
    double genrand_real3(void);

    if(argc<7){
        fprintf(stderr,"./traj nreps ageln h nsam nder ft seed\n");
        exit(1);
    }
    printf("// ");
    for(i=0;i<argc;i++) printf(" %s",argv[i]);
    printf("\n");

    nreps=strtol(argv[1],NULL,10);          // number of replicates.
    ageln=strtod(argv[2],NULL)-log(2);      // a parameter of log-normal distribution (t).
    h=strtod(argv[3],NULL);                 // a dominance effect.
    nsam=strtol(argv[4],NULL,10);           // number of samples.
    nder=strtol(argv[5],NULL,10);           // number of derived alleles.
    ft=strtod(argv[6],NULL);                // maximum frequency of ft (initial frequency of selected alleles).
    sscanf(argv[7],"%x",&seed);

    void init_by_array(unsigned long init_key[], int key_length);
    unsigned long init[4]={seed, 0x265, 0x367, 0x569}, length=4;
    init_by_array(init, length);
    srand((unsigned) time(NULL));

    freqsel=(double *)malloc((unsigned)((MAX/2.0)+1)*sizeof(double));     // trajectory of selected allele
    freqneu=(double *)malloc((unsigned)((MAX/2.0)+1)*sizeof(double));     // trajectory of neutral allele
    freq=(double *)malloc((unsigned)(MAX+1)*sizeof(double));
    pmax=nder/(double)nsam;

    popmax=0;
    for(i=0;i<MAX;i++) if(popsize(i)>popmax)popmax=popsize(i);
    printf("%d N0: %d\n",nreps,popsize(0));
    nsuccess=0;
    count=0;

    flag=0;
    flag_fwd=0;
    flag_bwd=0;

    fp=fopen("NewSFS.csv","r");
    char line[MAXLEN];
    double **SFS;
    SFS=(double **)malloc(SFSsamp*sizeof(double*));
    for(i=0;i<SFSsamp;i++){
        *(SFS+i)=(double *)malloc(PERIODS*sizeof(double));
        if(*(SFS+i)==NULL){
            printf("BUG\n");
            exit(0);
        }
    }
    for(i=0;i<SFSsamp;i++){
        fgets(line,MAXLEN,fp);
        sscanf(line,"%lf,%lf,%lf,%lf",&SFS[i][0],&SFS[i][1],&SFS[i][2],&SFS[i][3]);
    }
    fclose(fp);

    int tflag;
    while(nsuccess<nreps){

        count++;
        
    // Log-normal distribution (variance=mean^2)
        ngentmp=gasdev();
        ngenmu=log(ageln+sqrt(2))-(log(2)/2);
        ngen=(int)(exp(ngenmu+ngentmp*sqrt(log(2))));

    // Sample selection coefficient
        s_tmp=(genrand_real3()*(-2.5))-0.5;
        s=pow(10.0,s_tmp);

    // Sample initial frequency of focal allele
        if(ngen<=t1){tflag=0;
        }else if(ngen<=t2){tflag=1;
        }else if(ngen<=t3){tflag=2;
        }else{tflag=3;}
        while(1){
            p=0.0;
            ran=GetRandom(0,SFSsamp-1);
            if(ft>SFS[ran][tflag]){
                p=SFS[ran][tflag];
                break;
            }
        }

        pfwd=p;
        genfwd=0;
        freqsel[genfwd]=pfwd;
        while((genfwd<ngen) && (pfwd>0) && (pfwd<1.0)){
            genfwd++;
            wbar=pfwd*pfwd*(1.+s)+2.*pfwd*(1.-pfwd)*(1.+s*h)+(1.-pfwd)*(1.-pfwd);
            pprime=(pfwd*pfwd*(1.+s)+pfwd*(1.-pfwd)*(1.+s*h))/wbar;
            pfwd=bnldev(pprime,2*popsize(ngen-genfwd))/(2.0*popsize(ngen-genfwd));
            freqsel[genfwd]=pfwd;
        }

        prob=pow(pfwd/pmax,(double)nder)*pow((1.0-pfwd)/(1.0-pmax),(double)(nsam-nder));
        if(genrand_real3()<prob && genfwd==ngen){

            while(flag_bwd==0){

                pbwd=p;
                genbwd=0;
                freqneu[genbwd]=pbwd;

                while((ngen+genbwd<MAX) && (pbwd>0) && (pbwd<1.0)){
                    genbwd++;
                    pprime=pbwd;
                    pbwd=bnldev(pprime,2*popsize(ngen+genbwd))/(2.0*popsize(ngen+genbwd));
                    freqneu[genbwd]=pbwd;
                }
                if(pbwd==0 && (genrand_real3()<((double)popsize(genbwd)/(double)popmax))){
                    flag_fwd=1;
                    flag_bwd=1;
                    break;
                }
            }
        }

        if((flag_fwd==1) && (flag_bwd==1)){

            flag_fwd=0;
            flag_bwd=0;
            gen=0;
            for(i=genbwd;i>0;i--){
                freq[gen]=freqneu[i];
                gen++;
            }
            for(i=0;i<genfwd;i++){
                freq[gen]=freqsel[i];
                gen++;
            }
            nsuccess++;
            if((nreps>9) && (nsuccess %(int)(0.1*nreps)==0))fprintf(stderr,". ");
            printf("# s: %lf age: %lf freq: %lf\n",s,ngen/(4.0*popsize(0)),p);
            for(i=0;i<gen;i++)printf("%lf\t%lf\n",i/(4.0*popsize(0)),freq[i]);
        }
    }

    free(freq);
        freq=NULL;    
    free(freqsel);
        freqsel=NULL;    
    free(freqneu);
        freqneu=NULL;    
    free(SFS);
        SFS=NULL;

    fprintf(stderr," !%lf \n",nreps/(double)count);
    return 0;

}

double gasdev(void){

    double genrand_real3(void);
    static int iset=0;
    static double gset;
    double fac,rsq,v1,v2;

    if(iset==0){                            // We don't have an extra deviate handy, so

        do{
            v1=2.0*genrand_real3()-1.0;     // pick two uniform numbers in the square extending
            v2=2.0*genrand_real3()-1.0;     // fromb -1 to +1 in each direction.
            rsq=v1*v1+v2*v2;                // see if they are in the unit circle.
        }while(rsq>=1.0 || rsq==0.0);       // and if they are not, try again.

        fac=sqrt(-2.0*log(rsq)/rsq);        // Now make the Box-Muller transformation to get
        gset=v1*fac;                        // two normal deviates.
        iset=1;
        return v2*fac;                      // Return one and save the other for next time.
    }else{                                  // We have an extra deviate handy,
        iset=0;                             // so unset the flag,
        return gset;                        // and return it.
    }
}

int GetRandom(int min,int max){
    return min+(int)(rand()*(max-min+1.0)/(1.0+RAND_MAX));
}
