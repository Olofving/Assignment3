#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <sys/time.h>
#include "graphics.h"

//Particle object
typedef struct particle{
  double velx;
  double vely;
  double posx;
  double posy;
  double forcex;
  double forcey;
	double mass;
	double brightness;
} particle;

// function for timetaking
static double get_wall_seconds() {
  struct timeval tv;
  gettimeofday(&tv, NULL);
  double seconds = tv.tv_sec + (double)tv.tv_usec / 1000000;
  return seconds;
}

//Calculation of the forces between the particles
void force(particle *particle_list, int N, double G){

	double epsilon = 0.001;
	double r;
	
	
	for(int i=0; i<N; i++){
		particle_list[i].forcex=0.0;
		particle_list[i].forcey=0.0;
		
				
			for(int j=0;j<i; j++){
						r = sqrtf( (particle_list[i].posx - particle_list[j].posx)*(particle_list[i].posx - particle_list[j].posx)  + (particle_list[i].posy - particle_list[j].posy)* (particle_list[i].posy - particle_list[j].posy) ) + epsilon;
								
					particle_list[i].forcex += (particle_list[j].mass)/( r*r*r )*(particle_list[i].posx-particle_list[j].posx);				
					particle_list[i].forcey += (particle_list[j].mass)/( r*r*r )*(particle_list[i].posy-particle_list[j].posy);
			}
			
			for(int j=i+1;j<N; j++){
					r = sqrtf(   (particle_list[i].posx - particle_list[j].posx)*(particle_list[i].posx - particle_list[j].posx)  + (particle_list[i].posy - particle_list[j].posy)* (particle_list[i].posy - particle_list[j].posy) ) + epsilon;				
					particle_list[i].forcex += (particle_list[j].mass)/(r*r*r)*(particle_list[i].posx-particle_list[j].posx);				
					particle_list[i].forcey += (particle_list[j].mass)/(r*r*r)*(particle_list[i].posy-particle_list[j].posy);
			}
			particle_list[i].forcex = -G*particle_list[i].forcex;
			particle_list[i].forcey = -G*particle_list[i].forcey;
	}
}



// uppdates all the velocities and positions, with timestep delta_t
void timestep(particle *particle_list, int N, double delta_t){
	double accx;
	double accy;
	
	for(int i=0;i<N;i++){
		accx = particle_list[i].forcex;
		accy = particle_list[i].forcey;		
		particle_list[i].velx += delta_t*accx;
		particle_list[i].vely += delta_t*accy;		
		particle_list[i].posx += delta_t*particle_list[i].velx;
		particle_list[i].posy += delta_t*particle_list[i].vely;
	}
}
 

int main(int argc, const char** argv){

// Timer
double time1 = get_wall_seconds();

//Checking for correct ammount of arguments
	if(argc != 6) {
		printf("ERROR! Expected five input arguments: N filename nsteps delta_t graphics\n");
  	return -1;
	}

//declaring variables for input args
	const double N = atoi(argv[1]);
	const char* filename = argv[2];
	const int nsteps = atoi(argv[3]);
	const double delta_t = atof(argv[4]);
	const int graphics = atoi(argv[5]);
	const double G = 100/N;

// array of all the particles
	particle *particle_list = (particle*)malloc(N*sizeof(particle));

//opening and reading of file

	FILE* fptr;
	fptr = fopen(filename,"rb");

	int size_double=sizeof(double);
	double totalmass;

	for(int i =0;i<N;i++){  
		fread(&particle_list[i].posx,size_double,1,fptr);
		fread(&particle_list[i].posy,size_double,1,fptr);
		fread(&particle_list[i].mass,size_double,1,fptr);
		fread(&particle_list[i].velx,size_double,1,fptr);
		fread(&particle_list[i].vely,size_double,1,fptr);
		fread(&particle_list[i].brightness,size_double,1,fptr);
		printf("mass: %lf possx: %lf possy: %lf\n",particle_list[i].mass,particle_list[i].posx,particle_list[i].posy);
		totalmass += particle_list[i].mass;
	}
	printf("Total mass: %lf\n",totalmass);
	fclose(fptr);

//Simulation
	int i;

// graphics related stuff

	const float circleRadius=0.004, circleColor=0;
	float L=1, W=1;
	if(graphics == 1){
		InitializeGraphics("FATAL ERROR!",1400,800);
		SetCAxes(0,1);
	}
  
	for(i=0; i<nsteps; i++){

	force(particle_list,N,G);
	timestep(particle_list,N,delta_t);


    /* Call graphics routines. */
  if(graphics == 1){
	ClearScreen();
    for(int j=0;j<N;j++){
    	DrawCircle(particle_list[j].posx, particle_list[j].posy, L, W, circleRadius, circleColor);
    }   
       Refresh();
    /* Sleep a short while to avoid screen flickering. */
    usleep(3000);   
	}     
}

// Writing data to .gal-file	
	FILE *fp;
	fp = fopen("result.gal","wb");

	for(int i =0;i<N;i++){
		fwrite(&particle_list[i].posx,size_double,1,fp);
		fwrite(&particle_list[i].posy,size_double,1,fp);
		fwrite(&particle_list[i].mass,size_double,1,fp);
		fwrite(&particle_list[i].velx,size_double,1,fp);
		fwrite(&particle_list[i].vely,size_double,1,fp);
		fwrite(&particle_list[i].brightness,size_double,1,fp);
	}
	fclose(fp);
	free(particle_list);
	
	// printing of computation time in wall seconds
	printf("computation time: %7.3f wall seconds\n",get_wall_seconds()-time1);
	return 0;
}
