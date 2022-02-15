/* 
class_tsunami.c
Program to simulate tsunami waves on a 2D Cartesian grid. The program uses
a 4th-order finite difference solution of the equation
Ptt = div * GH grad P
where P is the height of the Tsunami wave above sea level, G is the
acceleration due to gravity (a constant), and H is the ocean depth.
The speed of the wave is sqrt(GH).

R. Clayton, Caltech, Jan 2005

COMPILE:
  gcc -o class_tsunami class_tsunami.c -lm

RUN:
  ./class_tsunami
*/
#include	<stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <fcntl.h> // open function
#include <unistd.h> // close function

/* WE - added a few extra arrays - V now just holds velocity where as Vth holds V*((dt/h)^2); Vgradx and Vgrady hold gradients of x,y   */
float *v, *vth, *vgradx, *vgrady, *p1, *p2;
float *f1, *f2;
int	ord	= 4;
int	nx	=1000;
int	ny	=800;
int	nt	=10000;
float	h	=10.0;
float	dt	=1.0;
char	vmodel[]	="bathymetry.in";
char	output[]	="hetero4th.out";
int	itprint	=100;   /* time steps between print messages */
int	itslice	=100;  /* time steps between slice outputs */
float	latref	=-40.0;
float	lonref	=35;
float	slat	=3.30;
float	slon	=95.87;

/*WE*/
int hetero = 1;

#define V(ix,iy)		v[(ix) +nx*(iy)]
#define Vth(ix,iy)	    vth[(ix) +nx*(iy)]
#define Vgradx(ix,iy)	vgradx[(ix) +nx*(iy)]
#define Vgrady(ix,iy)	vgrady[(ix) +nx*(iy)]

#define P1(ix,iy)		p1[(ix) +nx*(iy)]
#define P2(ix,iy)		p2[(ix) +nx*(iy)]

/* 2nd order second-derviative coefficients */
#define B1	-2.0	
#define B2	 1.0	

/* WE 4th order second-derviative coefficients */
#define C0   0.6666666666
#define C1	 0.0833333333
#define C2   16.000000000
#define C3   30.000000000
#define C4   0.0069444444
#define C5   8.0000000000



float	mindepth	= 10.0; /* minimum depth to consider still ocean (m) */
#define G 0.00985	/* acc of gravity in km/sec**2 */
int ixref	=0;
int iyref	=0;



main(int ac, char **av)
   {
	void zap(), output_slice(), xcoord_convert();

	int it, ix, iy, fd;
	int ixs, iys;
	float *tmp, vel, f, val, velmax, dtdh2;
	double norm(), sqrt();

	/*WE temp variables */
	double lhs_term, homo_termx, homo_termy, homogenous, heterox, heteroy;


	fprintf(stderr,"order= %d\n",ord);
	printf("Heterogenous (1) or Homogenous (0): %d\n", hetero);
	v= (float *)(malloc(4*nx*ny));
	vth= (float *)(malloc(4*nx*ny));
	vgradx = (float *)(malloc(4*nx*ny));
	vgrady = (float *)(malloc(4*nx*ny));

	f1= (float *)(malloc(4*nx*ny));
	f2= (float *)(malloc(4*nx*ny));


	if(v == NULL || f1 == NULL || f2 == NULL)
	   {
		fprintf(stderr,"cannot alloc memory\n");
		exit(-1);
	   }
	if( (fd= open(vmodel,0)) < 0)
	   {
		fprintf(stderr,"cannot open velocity file=%s\n",vmodel);
		exit(-1);
	   }
	if( read(fd,v,4*nx*ny) != 4*nx*ny )
	   {
		fprintf(stderr,"read error in velocity file=%s\n",vmodel);
		exit(-1);
	   }
	close(fd);
	output_slice(v,nx,ny,-1.0);

	/* convert depth to velocity v= sqrt(g*depth).
	 * set values for land (pos. depths) to negative as flag
	 */
	
	/* Calculate dt/dh squared */
	dtdh2 = (dt*dt)/(h*h);

	velmax= 0.0;
	for(iy=0; iy<ny; iy++)
	for(ix=0; ix<nx; ix++)
	   {
		val= -V(ix,iy); /* make depth positive */
		/* note 0.001 to convert depth (m) to km */
		if(val > mindepth) vel= sqrt(G*val*0.001);
		 else		   vel= -0.001;
		
		if(vel > velmax) velmax= vel;
		/* WE - define V as velocity squared where Vth is velocity squared with dt/dx factor */
		if(vel > 0.0) 
		  {
			V(ix,iy) = vel*vel;
			Vth(ix,iy) = vel*vel*dtdh2;
		  }	
		else	 
		  {
            V(ix,iy)= -0.001;
			Vth(ix,iy)= -0.001;
	      }
	   }  


	/*WE edit - pre-calculate the gradient of the velocity for computation efficiency - note that V no longer is multiplied by dt^2/h^2 */
	for(iy=1; iy < ny-1; iy++)
		for(ix=1; ix < nx-1; ix++)
		   {
			   Vgradx(ix, iy) = (1/h) * (-C1*V(ix+2, iy) + C0*V(ix+1, iy) - C0*V(ix-1, iy) + C1*V(ix-2, iy));
			   Vgrady(ix, iy) = (1/h) * (-C1*V(ix, iy+2) + C0*V(ix, iy+1) - C0*V(ix, iy-1) + C1*V(ix, iy-2));
		   } 


	for(iy=0; iy<ny; iy++)
	for(ix=0; ix<nx; ix++)

	/*fprintf(stdout,"maximum velocity = %8.4f (km/s)\n",velmax); 
	fprintf(stdout,"nx= %d ny=%d nt=%d h=%8.4f dt=%8.4f\n",nx,ny,nt,h,dt); */

	/* point the memory planes to real memory and zero it */
	p1= f1;
	p2= f2;

	zap(p1,nx*ny);
	zap(p2,nx*ny);


	/* add source */
	xcoord_convert(slat,slon,&ixs,&iys);
	fprintf(stderr,"source %8.3f %9.3f %4d %4d\n",slat,slon,ixs,iys);
	/*  source is placed on a grid:
	 *    1/4  1/2  1/4
	 *    1/2   1   1/2
	 *    1/4  1/2  1/4
	 */
	P2(ixs,iys)= 1.0;
	P2(ixs +1,iys  )= P2(ixs -1,iys  )= P2(ixs   ,iys+1)= P2(ixs   ,iys-1)= 0.5;
	P2(ixs+1,iys+1)= P2(ixs-1,iys+1)= P2(ixs+1,iys-1)= P2(ixs-1,iys-1)= 0.25;

	/* loop over time steps */
	for(it= 0; it<nt; it++)
	   {
		/* loop over x-y plane */
		for(iy=1; iy < ny-1; iy++)
		for(ix=1; ix < nx-1; ix++)
		   {
			/* ignore points on land */
			if(V(ix,iy) < 0.0)
			   {
				P1(ix,iy)= 0.0;
				continue;
			   }
			if(ord == 2 || ix == 1 || ix == nx-2 || iy ==1 || iy == ny-2)
			   {
				/* 2nd-order */
				P1(ix,iy)= 2.0*(1.0+B1*V(ix,iy))*P2(ix,iy) - P1(ix,iy)
					 + B2*V(ix,iy)*(P2(ix+1,iy)+P2(ix-1,iy)
						       +P2(ix,iy+1)+P2(ix,iy-1));
			   }
			 else
			   {

                 /* 4th-order 
                  * P1 = t-1 timestep  
                  * P2 = current timestep  
                  * The newest time step is calculated and stored into P1 in the above 2nd order system
                  * C1 = 1/12
                  * C2 = 16
                  * C3 = 30
                  * There is an implicit assumption that dx = dy for these calculations otherwise using V as defined above does not work 
                  *  */
                       
                   lhs_term   = 2*P2(ix,iy) - P1(ix,iy);
        
                   homo_termx =  C1*Vth(ix, iy)*(-P2(ix+2,iy) + C2*P2(ix+1,iy) - C3*P2(ix,iy) + C2*P2(ix-1,iy) - P2(ix-2,iy)) ;
                   homo_termy =  C1*Vth(ix, iy)*(-P2(ix,iy+2) + C2*P2(ix,iy+1) - C3*P2(ix,iy) + C2*P2(ix,iy-1) - P2(ix,iy-2)) ;
                       
                   homogenous = lhs_term + homo_termx + homo_termy;
                    

                 /* Heterogenous terms: 
                  * C0 = 2/3
				  * C1 = 1/12
				  * C4 = 1/144
                  * C5 = 8.000
                  * Note below that hetero is an integer that is 0 for homogenous case and 1 for heterogenous case
                  *  */


                 if(hetero==1){                     

				   
                   heterox = dtdh2*C4*(-P2(ix+2,iy) + C5*(P2(ix+1,iy)- P2(ix-1,iy)) + P2(ix-2,iy))
                                       *(- V(ix+2,iy) + C5*(V(ix+1,iy) - V(ix-1,iy)) + V(ix-2,iy));
                   

                   heteroy = dtdh2*C4*(P2(ix,iy-2)-P2(ix,iy+2) + C5*(P2(ix,iy+1) - P2(ix,iy-1)))
                               *(V(ix,iy-2) - V(ix,iy+2) + C5*(V(ix,iy+1) -  V(ix,iy-1)) ); 
				   
				   
				   /* 
					THIS METHOD BLOWS UP FASTER
				   heterox = (dt*dt) * (-C1*P2(ix+2,iy) + C0*(P2(ix+1,iy) - C1*P2(ix-1,iy)) + C0*P2(ix-2,iy)) * Vgradx(ix,iy)/h;
				   heteroy = (dt*dt) * (-C1*P2(ix,iy+2) + C0*(P2(ix,iy+1) - C1*P2(ix,iy-1)) + C0*P2(ix,iy-2)) * Vgrady(ix,iy)/h; 
				   */

                   }

                 P1(ix,iy) = homogenous + (heterox + heteroy)*hetero ; 


		   		}
		   }
	


		/* Dirichlet boundary conditions */
		for(ix=0,    iy=0;    ix<nx; ix++) P1(ix,iy)= 0.0;
		for(ix=0,    iy=ny-1; ix<nx; ix++) P1(ix,iy)= 0.0;
		for(ix=0,    iy=0;    iy<ny; iy++) P1(ix,iy)= 0.0;
		for(ix=nx-1, iy=0;    iy<ny; iy++) P1(ix,iy)= 0.0;

		if(it%itprint == 0)
			fprintf(stderr,"done it=%3d norm=%14.3e \n",it,
				norm(&P1(0,0),nx*ny));
		if(it%itslice == 0)
			output_slice(p1,nx,ny,(double)(it*dt));
		/* rotate the memory pointers */
		tmp= p1; p1= p2; p2= tmp;
	   }
	}




#define DEG2KM	111.195
#define DEG2R	  0.0174532
void xcoord_convert(double slat, double slon,int *ixs,int *iys)
   {
	/* convert (lat,lon) to grid coords (ixs,iys) */
	double cos();
	double x, y;
	y= (slat-latref)*DEG2KM;
	x= (slon-lonref)*DEG2KM*cos(slat*DEG2R);
	*ixs = x/h + ixref;
	*iys = y/h + iyref;
   }

double norm(float *x, int n)
   {
	/* compute the norm of a vector or plane */
	double sum, val, sqrt();
	int i;

	sum= 0.0;
	for(i=0; i<n; i++) sum += x[i]*x[i];
	val= sqrt( sum/(double)(n) );
	return(val);
   }

void zap(float *x, int n)
   {
	/* zero out a field */
	int i;
	for(i=0; i<n; i++) x[i]= 0.0;
   }

double getmax(float *x, int n)
   {
	/* find the absolute max of a plane */
	int i;
	double fabs(), max;
	max= fabs(x[0]);
	for(i=1; i<n; i++) if(fabs(x[i]) > max) max= fabs(x[i]);
	return(max);
   }

/* output a slice of the field.
 * Note: the field is output at every istep samples
 *       the field is reversed in y
 */
float line[1000];
int outfd	= -1;
int	istep	= 2;
void output_slice(float *x,int nx,int ny,double t)
   {
	double max, getmax(), val;
	int ix, iy, i, ival;

	if(outfd < 0) outfd= creat(output,0664);
	if(outfd < 0)
	   {
		fprintf(stderr,"cannot create plot file= %s\n",output);
		exit(-1);
	   }
	for(iy=ny-1; iy >= 0; iy -= istep)
	   {
		for(ix=0, i=0; ix <nx; ix += istep, i++)
			line[i]= x[ix+iy*nx];
		write(outfd,line,4*i);
	   }
   }
