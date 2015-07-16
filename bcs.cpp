#include "headers.h"
#include "init_2.h"
#include "declarations_2.h"

//////////////////////Boundary conditions on the top, bottom and side walls 

void bcs(vertex * node, fval * fvar, double Q)
{	
	int ind; 
	double ux,vx,px; 
	
	tdma1y(fvar,0); //Computing du/dy for the boundary

	tdma2x(fvar,0);  //d2u/dx2  //Computing the second derivatives needed for pressure BCs
	tdma2y(fvar,0);  //d2u/dy2  //Computing the second derivatives needed for pressure BCs
	
	tdma2x(fvar,1);  //d2v/dx2  //Computing the second derivatives needed for pressure BCs
	tdma2y(fvar,1);  //d2v/dy2 
	
	tdma1x(fvar,0);  //du/dx //Computing du/dx required for the outflow pressure boundary condition
	tdma1x(fvar,1);  //dv/dx //Computing dv/dx required for the outflow pressure boundary condition

	//Top boundary 

	for(int i=1;i<nx;i++)
	{
		ind = i + ny*str_x; 					
		
		//BFS BC
		
		fvar[ind].u[0] = 0.0; 
		fvar[ind].u[1] = 0.0; 	
		fvar[ind].u[2] = -(1./B1)*dy*( (1./Re)*fvar[ind].uyy[1] )- B2*fvar[ind-str_x].u[2] - B3*fvar[ind-2*str_x].u[2] - B4*fvar[ind-3*str_x].u[2];
		
		//Kovasznay BC
		
		/*fvar[ind].u[0] = 1. - exp(lambda*(node[ind].x[0]-0.5))*cos(2.*M_PI*(node[ind].x[1]-0.5));  		
		fvar[ind].u[1] = (lambda/(2.*M_PI))*exp(lambda*(node[ind].x[0]-0.5))*sin(2.*M_PI*(node[ind].x[1]-0.5));
		fvar[ind].u[2] = - B2*fvar[ind-str_x].u[2] - B3*fvar[ind-2*str_x].u[2] - B4*fvar[ind-3*str_x].u[2];*/		
	}

	//Right boundary 

	for(int j=0;j<=ny;j++)
	{
		ind = nx + j*str_x; 				
		
		//Kovasznay BC
		 
		/*ux = -lambda*exp(lambda*(node[ind].x[0]-0.5))*cos(2.*M_PI*(node[ind].x[1]-0.5)); //For Kovasznay flow
		vx = (lambda*lambda/(2.*M_PI))*exp(lambda*(node[ind].x[0]-0.5))*sin(2.*M_PI*(node[ind].x[1]-0.5));
		
		fvar[ind].u[0] = -(1./B1)*dx*(ux) - ( B2*fvar[ind-1].u[0] + B3*fvar[ind-2].u[0] + B4*fvar[ind-3].u[0] ) ; 
		fvar[ind].u[1] = -(1./B1)*dx*(vx) - ( B2*fvar[ind-1].u[1] + B3*fvar[ind-2].u[1] + B4*fvar[ind-3].u[1] ) ;			
		fvar[ind].u[2] = (1./Re)*fvar[ind].ux[0] - (1./4.)*( (fvar[ind].u[0]*fvar[ind].u[0]) + (fvar[ind].u[1]*fvar[ind].u[1]) )*(1. - tanh(fvar[ind].u[0]/delta)); //Dirichlet boundary condition from the paper by Dong et al  //OBC -C type 		*/
		
		/*BFS BC*/
		
		fvar[ind].u[0] = fvar[ind].u[0] - fvar[ind].ux[0]*(dt/Q); 
		fvar[ind].u[1] = fvar[ind].u[1] - fvar[ind].ux[1]*(dt/Q);
		fvar[ind].u[2] = (1./Re)*fvar[ind].ux[0];		
	}	

	//Bottom boundary

	for(int i=1;i<=nx;i++)
	{			
		//BFS BC
		
		fvar[i].u[0] = 0.0;
		fvar[i].u[1] = 0.0;
		fvar[i].u[2] = (1./B1)*dy*( (1./Re)*fvar[i].uyy[1]) - ( B2*fvar[i+str_x].u[2] + B3*fvar[i+2*str_x].u[2] + B4*fvar[i+3*str_x].u[2] );
		
		//Kovasznay BC
		
		/*fvar[i].u[0] = 1. - exp(lambda*(node[i].x[0]-0.5))*cos(2.*M_PI*(node[i].x[1]-0.5));  //For Kovasznay flow
		fvar[i].u[1] = (lambda/(2.*M_PI))*exp(lambda*(node[i].x[0]-0.5))*sin(2.*M_PI*(node[i].x[1]-0.5));
		fvar[i].u[2] = - ( B2*fvar[i+str_x].u[2] + B3*fvar[i+2*str_x].u[2] + B4*fvar[i+3*str_x].u[2] );*/		
	}

	//Left boundary

	for(int j=0;j<=ny;j++)
	{
	 	ind = j*str_x;
	 	
	 	//BFS BC  	 		

		if(node[ind].x[1]<0.5)
		{
			fvar[ind].u[0] = 0.0;
			fvar[ind].u[1] = 0.0; 			
	                fvar[ind].u[2] = (1./B1)*dx*(1./Re)*(fvar[ind].uxx[0]) - ( B2*fvar[ind+1].u[2] + B3*fvar[ind+2].u[2] + B4*fvar[ind+3].u[2]) ;
		}	
		else
		{
		  	fvar[ind].u[0] = 24.*(1.-node[ind].x[1])*(node[ind].x[1]-0.5);
                        fvar[ind].u[1] = 0.0;
                        fvar[ind].u[2] = (1./B1)*dx*(1./Re)*(fvar[ind].uxx[0] + fvar[ind].uyy[0]) - ( B2*fvar[ind+1].u[2] + B3*fvar[ind+2].u[2] + B4*fvar[ind+3].u[2]) ;			

		}
		
		//Kovasznay BC
		
		/*px = -lambda*exp(2.*lambda*(node[ind].x[0]-0.5)); 
		
		fvar[ind].u[0] = 1. - exp(lambda*(node[ind].x[0]-0.5))*cos(2.*M_PI*(node[ind].x[1]-0.5));  //For Kovasznay flow
		fvar[ind].u[1] = (lambda/(2.*M_PI))*exp(lambda*(node[ind].x[0]-0.5))*sin(2.*M_PI*(node[ind].x[1]-0.5));		
		fvar[ind].u[2] = (1./B1)*dx*px - (B2*fvar[ind+1].u[2] + B3*fvar[ind+2].u[2] + B4*fvar[ind+3].u[2]);*/		
	}
}
