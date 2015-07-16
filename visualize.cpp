#include "headers.h"
#include "init_2.h"
#include "declarations_2.h"

///////////////////Function overloaded of Writing to a file////////////////////
//Writing data at frequency specified by file_freq in the "init" file

void write_to_file(vertex * node, fval * fvar, int t)
{
	ofstream out_put, vel_put0, vel_put1; //0->x-velocity along y-direction at geometric center and 1->y-velocity along x-direction at geom center  	
	
	//ofstream outpmtv, outvmtv, outvortmtv, outtmtv; 
	
	int p,p1,p2;
	double udy=0.0, vdx=0.0; 

        tdma1y(fvar,0);
	tdma1x(fvar,1);  

	string file("4bdata_2D");	
	
	/*string filep("4bpressure_mtv"); 
	string filev("4bvector_mtv"); 
	string filevort("4bvorticity_mtv"); 
	string filetemp("4btemp_mtv");*/
	
	stringstream tag; 

	tag<<t;

	file = file + "_" + tag.str() + ".vtk";   
	
	/*filep = filep + "_" + tag.str() + ".dat";   
	filev = filev + "_" + tag.str() + ".dat";   
	filevort = filevort + "_" + tag.str() + ".dat";   
	filetemp = filetemp + "_" + tag.str() + ".dat";*/

	out_put.open(file.c_str());
	/*outpmtv.open(filep.c_str());
	outvmtv.open(filev.c_str());
	outvortmtv.open(filevort.c_str());
	outtmtv.open(filetemp.c_str());*/

	out_put<<"# vtk DataFile Version 2.0"<<"\n"; 
	out_put<<"2DGrid"<<"\n";
	out_put<<"ASCII"<<"\n";
	out_put<<"DATASET"<<" "<<"STRUCTURED_GRID"<<"\n";
	out_put<<"DIMENSIONS "<<nx+1<<" "<<ny+1<<" 1"<<"\n";
	out_put<<"POINTS"<<" "<<tot_p<<" float"<<"\n";

	for (int j=0;j<=ny;j++)
	{ 
		for (int i=0;i<=nx;i++)
		{
			p = i + j*str_x; 
			out_put<<node[p].x[0]<<" "<<node[p].x[1]<<" 0.0"<<"\n"; 
		}
	}

	out_put<<"POINT_DATA "<<tot_p<<"\n";
	out_put<<"SCALARS"<<" Pressure"<<" double"<<"\n"; 
	out_put<<"LOOKUP_TABLE"<<" Pressure_lp"<<"\n";

	for (int j=0;j<=ny;j++)
	{ 
		for (int i=0;i<=nx;i++)
		{
			p = i + j*str_x; 
			out_put<<fvar[p].u[2]<<"\n"; 
		}
	}
	
	out_put<<"SCALARS"<<" Temperature"<<" double"<<"\n"; 
	out_put<<"LOOKUP_TABLE"<<" res_lp"<<"\n";

	for (int j=0;j<=ny;j++)
	{ 
		for (int i=0;i<=nx;i++)
		{
			p = i + j*str_x; 
			out_put<<fvar[p].u[5]<<"\n"; 
		}
	}
	
	out_put<<"SCALARS"<<" divergence"<<" double"<<"\n"; 
	out_put<<"LOOKUP_TABLE"<<" div_lp"<<"\n";

	for (int j=0;j<=ny;j++)
	{ 
		for (int i=0;i<=nx;i++)
		{
			p = i + j*str_x; 
			out_put<<fvar[p].div<<"\n"; 
		}
	}
	out_put<<"SCALARS"<<" residual"<<" double"<<"\n"; 
	out_put<<"LOOKUP_TABLE"<<" res_lp"<<"\n";

	for(int j=0;j<=ny;j++)
	{ 
		for (int i=0;i<=nx;i++)
		{
			p = i + j*str_x; 
			out_put<<fvar[p].res<<"\n"; 
		}
	}
	
	out_put<<"SCALARS"<<" Vorticity"<<" double"<<"\n";
        out_put<<"LOOKUP_TABLE"<<" Vorticity_lp"<<"\n";

        for(int j=0;j<=ny;j++)
        {
                for (int i=0;i<=nx;i++)
                {
                        p = i + j*str_x;
                        out_put<<fvar[p].uy[0]-fvar[p].ux[1]<<"\n";
                }
        }


	out_put<<"VECTORS"<<" Velocity"<<" double"<<"\n"; 

	for (int j=0;j<=ny;j++)
	{ 
		for (int i=0;i<=nx;i++)
		{
			p = i + j*str_x; 
			out_put<<fvar[p].u[0]<<" "<<fvar[p].u[1]<<" "<<"0.0"<<"\n";
		}
	}
        	
	/****************************************************************/
	/************Writing data in the PLOTMTV format******************/	

	/*outpmtv<<"$DATA=CONTOUR\n";
	outpmtv<<"%xmin= 0"<<"\n";
	outpmtv<<"%ymin= 0"<<"\n";
	outpmtv<<"%xmax= "<<l_x<<"\n";
	outpmtv<<"%ymax= "<<l_y<<"\n";
	outpmtv<<"%xlabel=\"x\""<<"\n";
	outpmtv<<"%ylabel=\"y\""<<"\n";
	outpmtv<<"%contstyle=2"<<"\n";
	outpmtv<<"%nx="<<nx+1<<"\n";
	outpmtv<<"%ny="<<ny+1<<"\n";
	outpmtv<<"%cmin=-17000"<<"\n";
	outpmtv<<"%cmax=17000"<<"\n";
	outpmtv<<"%nsteps=50"<<"\n";
	outpmtv<<"%INTERPOLATE=2"<<"\n";
	
	for (int j=0;j<=ny;j++)
	{ 
		for (int i=0;i<=nx;i++)
		{
			p = i + j*str_x; 
			outpmtv<<fvar[p].u[2]<<"\n"; 
		}
	}
	
	outvortmtv<<"$DATA=CONTOUR\n";
	outvortmtv<<"%xmin= 0"<<"\n";
	outvortmtv<<"%ymin= 0"<<"\n";
	outvortmtv<<"%xmax= "<<l_x<<"\n";
	outvortmtv<<"%ymax= "<<l_y<<"\n";
	outvortmtv<<"%xlabel=\"x\""<<"\n";
	outvortmtv<<"%ylabel=\"y\""<<"\n";
	outvortmtv<<"%contstyle=2"<<"\n";
	outvortmtv<<"%nx="<<nx+1<<"\n";
	outvortmtv<<"%ny="<<ny+1<<"\n";
	outvortmtv<<"%cmin=-17000"<<"\n";
	outvortmtv<<"%cmax=17000"<<"\n";
	outvortmtv<<"%nsteps=50"<<"\n";
	outvortmtv<<"%INTERPOLATE=2"<<"\n";
	
	for(int j=0;j<=ny;j++)
        {
                for (int i=0;i<=nx;i++)
                {
                        p = i + j*str_x;
                        outvortmtv<<fvar[p].uy[0]-fvar[p].ux[1]<<"\n";
                }
        }

	outtmtv<<"$DATA=CONTOUR\n";
	outtmtv<<"%xmin= 0"<<"\n";
	outtmtv<<"%ymin= 0"<<"\n";
	outtmtv<<"%xmax= "<<l_x<<"\n";
	outtmtv<<"%ymax= "<<l_y<<"\n";
	outtmtv<<"%xlabel=\"x\""<<"\n";
	outtmtv<<"%ylabel=\"y\""<<"\n";
	outtmtv<<"%contstyle=2"<<"\n";
	outtmtv<<"%nx="<<nx+1<<"\n";
	outtmtv<<"%ny="<<ny+1<<"\n";
	outtmtv<<"%cmin=-17000"<<"\n";
	outtmtv<<"%cmax=17000"<<"\n";
	outtmtv<<"%nsteps=50"<<"\n";
	outtmtv<<"%INTERPOLATE=2"<<"\n";
	
	for(int j=0;j<=ny;j++)
        {
                for (int i=0;i<=nx;i++)
                {
                        p = i + j*str_x;
                        outtmtv<<fvar[p].u[5]<<"\n";
                }
        }	
	
	outvmtv<<"$DATA=VECTOR\n";	
	outvmtv<<"%xmin=0"<<"	"<<"xmax = "<<l_x<<"\n";
	outvmtv<<"%ymin=0"<<"	"<<"ymax = "<<l_y<<"\n";
	outvmtv<<"%zmin=0"<<"	"<<"zmax = 0"<<"\n";	
	outvmtv<<"%vscale = 0.015"<<"\n";	
	
	for (int j=0;j<=ny;j++)
	{ 
		for (int i=0;i<=nx;i++)
		{
			p = i + j*str_x; 
			outvmtv<<node[p].x[0]<<" "<<node[p].x[1]<<" 0.0"<<" "<<fvar[p].u[0]<<" "<<fvar[p].u[1]<<" "<<"0.0"<<"\n";
		}
	}		       
		
	out_put.close(); 	 	
	outpmtv.close(); 
	outvmtv.close();*/
}

///////////////////Writing to a file////////////////////

void write_to_file(vertex * node, fval * fvar)
{

	ofstream out_put;  
	int p; 

	out_put.open("plane_jet.vtk");

	out_put<<"# vtk DataFile Version 2.0"<<"\n"; 
	out_put<<"2DGrid"<<"\n";
	out_put<<"ASCII"<<"\n";
	out_put<<"DATASET"<<" "<<"STRUCTURED_GRID"<<"\n";
	out_put<<"DIMENSIONS "<<nx+1<<" "<<ny+1<<" 1"<<"\n";
	out_put<<"POINTS"<<" "<<tot_p<<" float"<<"\n";

	for (int j=0;j<=ny;j++)
	{ 
		for (int i=0;i<=nx;i++)
		{
			p = i + j*str_x; 
			out_put<<node[p].x[0]<<" "<<node[p].x[1]<<" 0.0"<<"\n"; 
		}
	}

	out_put<<"POINT_DATA "<<tot_p<<"\n";
	out_put<<"SCALARS"<<" Pressure"<<" double"<<"\n"; 
	out_put<<"LOOKUP_TABLE"<<" Pressure_lp"<<"\n";

	for (int j=0;j<=ny;j++)
	{ 
		for (int i=0;i<=nx;i++)
		{
			p = i + j*str_x; 
			out_put<<fvar[p].u[2]<<"\n"; 			
						
		}
	}

	out_put<<"SCALARS"<<" rhs"<<" double"<<"\n"; 
	out_put<<"LOOKUP_TABLE"<<" rhs"<<"\n";

	for (int j=0;j<=ny;j++)
	{ 
		for (int i=0;i<=nx;i++)
		{
			p = i + j*str_x; 
			out_put<<fvar[p].rhs<<"\n"; 
		}
	}
	
	out_put<<"SCALARS"<<" F"<<" double"<<"\n"; 
	out_put<<"LOOKUP_TABLE"<<" F"<<"\n";

	for (int j=0;j<=ny;j++)
	{ 
		for (int i=0;i<=nx;i++)
		{
			p = i + j*str_x; 
			out_put<<fvar[p].F<<"\n"; 
		}
	}
	
	out_put<<"SCALARS"<<" disc_error"<<" double"<<"\n"; 
	out_put<<"LOOKUP_TABLE"<<" disc_error"<<"\n";

	for (int j=0;j<=ny;j++)
	{ 
		for (int i=0;i<=nx;i++)
		{
			p = i + j*str_x; 
			out_put<<fvar[p].u[2]- cos(M_PI*i*dx)*cos(M_PI*j*dy)<<"\n";			
			//out_put<<fvar[p].u[2]- sin(M_PI*i*dx)<<"\n";  
			//out_put<<fvar[p].u[2]- cos(M_PI*j*dy)<<"\n";		
		}
	}	
	
	
	out_put<<"VECTORS"<<" Velocity"<<" double"<<"\n"; 

	for (int j=0;j<=ny;j++)
	{ 
		for (int i=0;i<=nx;i++)
		{
			p = i + j*str_x; 
			out_put<<fvar[p].u[0]<<" "<<fvar[p].u[1]<<" "<<"0.0"<<"\n";
		}
	}

	out_put.close(); 

}

/*********************************************************************/
//Writing the unsteady time data to a file
//Written for only Bousinessq problem at hand

void write_time_file(fval * fvar, int t, ofstream & tp)
{
	double umax, umin, vmax, vmin, umax1, umin1; 
	
	//Writing the max and min of the X-velocity along the vertical line passing through geometric center 
        
        umax = fvar[nx/2].u[0]; 
        umin = fvar[nx/2].u[0];

	umax1 = fvar[nx/4].u[0]; 
	umin1 = fvar[nx/4].u[0];   
        
	for(int j=1;j<=ny;j++)
	{		
		if(fvar[nx/2+j*str_x].u[0]>=umax) umax = fvar[nx/2+j*str_x].u[0];
		if(fvar[nx/4+j*str_x].u[0]>=umax1) umax1 = fvar[nx/4+j*str_x].u[0]; 
 
		if(fvar[nx/2+j*str_x].u[0]<=umin) umin = fvar[nx/2+j*str_x].u[0]; 		
		if(fvar[nx/4+j*str_x].u[0]<=umin1) umin1 = fvar[nx/4+j*str_x].u[0]; 
	}		
	
	//Writing the max and min of the Y-velocity along the horizontal line passing through geometric center 

	vmax = fvar[(ny/2)*str_x].u[1]; 
	vmin = fvar[(ny/2)*str_x].u[1]; 
	
	for(int i=1;i<=nx;i++)
	{
		if(fvar[i+(ny/2)*str_x].u[1]>=vmax) vmax = fvar[i+(ny/2)*str_x].u[1]; 
		if(fvar[i+(ny/2)*str_x].u[1]<=vmin) vmin = fvar[i+(ny/2)*str_x].u[1];	
	}
	
	tp<<t<<"	"<<t*dt<<"	"<<umin<<"	"<<umax<<"	"<<vmin<<"	"<<vmax<<"	"<<umin1<<"	"<<umax1<<"	"<<"\n";
}

/***************************************************************************************************************/
