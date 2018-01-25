// This program generates frequency trajectories of an allele that was observed to have frequency = nder/nsam.   
// It assumes that the mutation arose once, at a random time in the past, and no further mutation occurs.
// The command line arguments are "nreps nsam nder seed".
// nreps: number of trajectories to generate.
// nsam: sample size.
// nder: number of copies of the derived allele found in the sample of size nsam.
// seed: random number.
// The function popsize(j), returns the diploid population size at time j generations in the past.
// The user must supply this function. Examples are provided.
// A Wright-Fisher model is assumed.  
// The first line of the output is nreps.
// The following line contains "npoints: n1 (where n1 is the number of time points until present)."
// Following that are n1 lines where each line contains two numbers, the time point(=generation/(4*N)),and the frequency.
// Time is measured in units of 4N generations.
// Compilation: "gcc -o traj -lm Trajectory_NTL.c popsize_stv.c binomial_mt.c mt.c"
// usage: "./traj nreps nsam nder seed | ./stepftn > XXX.out"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

extern int npres;
extern double r;

#define MAX 960000

main(int argc,char **argv){

    int nreps,gen,i,popmax,nsam,nder,count,nsuccess,seed;
    double p,pbwd,wbar,pprime,*freq,pmax,prob;

    int popsize(int);
    float bnldev(float,int);
    double genrand_real3(void);    // uniform dev.

    if(argc<4){
        fprintf(stderr,"./traj nreps nsam nder seed\n");
        exit(1);
    }
    printf("// ");
    for(i=0;i<argc;i++) printf(" %s",argv[i]);
    printf("\n");

    nreps=strtol(argv[1],NULL,10);          // number of replicates.
    nsam=strtol(argv[2],NULL,10);           // number of samples.
    nder=strtol(argv[3],NULL,10);           // number of derived alleles.
    sscanf(argv[4],"%x",&seed);             // random seed.

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

        gen=0;
        pbwd=pmax;
        freq[gen]=pbwd;

        while((gen<MAX) && (pbwd>0.0) && (pbwd<1.0)){
            gen++;
            pprime=pbwd;
            pbwd=bnldev(pprime,2*popsize(gen))/(2.0*popsize(gen));
            freq[gen]=pbwd;
        }

        if(pbwd==0.0 && (genrand_real3()<((double)popsize(gen)/(double)popmax))){
            nsuccess++;
            if((nreps>9) && (nsuccess %(int)(0.1*nreps)==0))fprintf(stderr,". ");		
            printf("# s: 0.0 age: %lf freq: %lf\n",gen/(4.0*popsize(0)),pbwd);
            for(i=0;i<gen;i++)printf("%lf\t%lf\n",i/(4.0*popsize(0)),freq[gen-i]);
        }
    }
    
    free(freq);    
        freq=NULL;    

    fprintf(stderr," !%.10f \n",nreps/(double)count);
    return 0;    
}

