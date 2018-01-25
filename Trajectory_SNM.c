// This program generates frequency trajectories of an allele that was observed to have frequency = nder/nsam.   
// It assumes that the mutation arose once, at a random time in the past, and no further mutation occurs.
// The command line arguments are "nreps genmax h nsam nder."
// nreps: number of trajectories to generate.
// ngen: the parameter of prior distribution for the time of selection on standing variation which was present at ft.
// h: dominance coefficient.
// nsam: sample size.
// nder: number of copies of the derived allele found in the sample of size nsam.
// The function popsize(j), returns the diploid population size at time j generations in the past.
// The user must supply this function. Examples are provided.
// A Wright-Fisher model is assumed.  
// The first line of the output is nreps.
// The following line contains "npoints: n1 (where n1 is the number of time points until present)."
// Following that are n1 lines where each line contains two numbers, the time point(=generation/(4*N)),and the frequency.
// Time is measured in units of 4N generations.
// Compilation: "gcc -o traj -lm Trajectory_SNM.c popsize.c binomial_mt.c mt.c"
// usage: "./traj nreps ageln h nsam nder seed | ./stepftn > XXX.traj"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

extern int npres;
extern double r;

#define MAX 960000

main(int argc,char **argv){
    
    int nreps,gen,ngen,i,popmax,nsam,nder,count,nsuccess,seed;
    double s,s_tmp,h,p,pfwd,wbar,pprime,*freq,pmax,prob;
    double ngentmp,ngenmu,ageln,alpha;

    int popsize(int);
    float bnldev(float,int);
    double gasdev(void);
    double genrand_real3(void);

    if(argc<6){
        fprintf(stderr,"./traj nreps ageln h nsam nder seed\n");
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
    sscanf(argv[6],"%x",&seed);

    void init_by_array(unsigned long init_key[], int key_length);
    unsigned long init[4]={seed, 0x265, 0x367, 0x569}, length=4;
    init_by_array(init, length);
    
    freq=(double *)malloc((unsigned)(MAX+1)*sizeof(double));
    pmax=nder/(double)nsam;

    popmax=0;
    for(i=0;i<MAX;i++) if(popsize(i)>popmax)popmax=popsize(i);
    printf("%d N0: %d\n",nreps,popsize(0));
    nsuccess=0;
    count=0;
                    
    while(nsuccess<nreps){
        count++;

    // Log-normal distribution (variance=mean^2)
        ngentmp=gasdev();
        ngenmu=log(ageln+sqrt(2))-(log(2)/2);
        ngen=(int)(exp(ngenmu+ngentmp*sqrt(log(2))));
        
    // ft (fixed value)        
        p=1.0/(2.0*popsize(ngen));

    // Sample selection coefficient
        s_tmp=(genrand_real3()*(-2.5))-0.5;
        s=pow(10.0,s_tmp);

        freq[0]=0.0;
        pfwd=p;
        gen=1;
        freq[gen]=pfwd;
        while((gen<ngen) && (pfwd>0) && (pfwd<1.0)){
            gen++;
            wbar=pfwd*pfwd*(1.+s)+2.*pfwd*(1.-pfwd)*(1.+s*h)+(1.-pfwd)*(1.-pfwd);
            pprime=(pfwd*pfwd*(1.+s)+pfwd*(1.-pfwd)*(1.+s*h))/wbar;
            pfwd=bnldev(pprime,2*popsize(ngen-gen))/(2.0*popsize(ngen-gen));
            freq[gen]=pfwd;
        }
        
        prob=pow(pfwd/pmax,(double)nder)*pow((1.0-pfwd)/(1.0-pmax),(double)(nsam-nder));
        if(genrand_real3()<prob && gen==ngen){
            
            nsuccess++;
            printf("# s: %lf age: %.10f freq: %lf\n",s,ngen/(4.0*popsize(0)),p);
            for(i=0;i<ngen;i++)printf("%lf\t%lf\n",i/(4.0*popsize(0)),freq[i]);
        }
    }
    
    free(freq);
        freq=NULL;    
    fprintf(stderr," !%.10f \n",nreps/(double)count);

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

