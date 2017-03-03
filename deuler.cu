//  2016 Juncheng E at PIMS.

#define _USE_MATH_DEFINES
#include <math.h>
#include <stdio.h>
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>  
#include <ctime>
#include <stdlib.h>
#include <string.h>

using namespace std;

int count = 0;
int nmax; // maximum atom number
int natom; // actual atom number
int numgrain; // the number of grains
int nx,ny,nz; // number of primitive cells in each direction
float a; // lattice constant
float a2; // lattice constant^2
float lx,ly,lz; // the size of the simulation cell
float *alpha, *beta, *gama;
float3 *gr_centerp; // the centers of each of the grains
float ratio = 1.0;
float3 *r; // atom positon
float temp,mass; // temperature, mass
int *atom_grain, *atom_neigh;
bool *atom_id;
//float *atom_mini;
//float *drlist;
int3 DIM;
__constant__ float d_lx, d_ly, d_lz, d_a;
__constant__ int d_nmax, d_numgrain, d_nx[2], d_ny[2], d_nz[2], d_natom, d_GBatoms;


int read_config(char* ifn)
{
	int i;
	ifstream ifile;
	string linebuffer;
	stringstream ss;
	ifile.open(ifn);
	cout << "Read " << ifn << "..." << endl;
	getline(ifile, linebuffer);
	ss << linebuffer;
	ss >> a;
	a2 = a*a;
	ss.str(""); // Clean up ss
	ss.clear(); // Clean up ss
	getline(ifile, linebuffer);
	ss << linebuffer;
	ss >> numgrain >> lx >> ly >> lz;
	ss.str(""); // Clean up ss
	ss.clear(); // Clean up ss
	printf("Number of grains: %d\n",numgrain);
	gr_centerp = new float3 [numgrain];
	alpha = new float [numgrain];
	beta = new float [numgrain];
	gama = new float [numgrain];
	getline(ifile, linebuffer);
	ss << linebuffer;
	ss >> mass >> temp;
	ss.str(""); // Clean up ss
	ss.clear(); // Clean up ss
	for ( i = 0; i < numgrain; ++i )
	{
		getline(ifile, linebuffer);
		ss << linebuffer;
		ss >> gr_centerp[i].x >> gr_centerp[i].y >> gr_centerp[i].z >> alpha[i] >> beta[i] >> gama[i];
//		alpha[i] = alpha[i]*M_PI/180.0;
//		beta[i] = beta[i]*M_PI/180.0;
//		gama[i] = gama[i]*M_PI/180.0;
		ss.str(""); // Clean up ss
		ss.clear(); // Clean up ss
	}
	ifile.close();
	return 0;
}

__device__ float pos_PBC(float pos, float l)
	// This function calculates and returns the positions of the 
	// atoms with periodic boundary conditions used.  
{
	float pos_PBC;
	if (pos < (0.0))
		pos_PBC = pos + l;
	else if (pos > (l)) 
		pos_PBC = pos - l;
	else
		pos_PBC = pos;
	return pos_PBC;
}

__device__ float separation_PBC(float ds, float l)
{
	float s_PBC;
	if (ds > (0.5*l)) 
		s_PBC = ds - l;
	else if (ds < (-0.5*l)) 
		s_PBC = ds + l;
	else
		s_PBC = ds;
	return s_PBC;
}

__device__ int getGlobalIdx_3D_3D_l(int l)
{
	int blockId = blockIdx.x + blockIdx.y * gridDim.x + gridDim.x * gridDim.y * blockIdx.z;
	int threadId = blockId * (blockDim.x * blockDim.y * blockDim.z * 4) + (threadIdx.z * (blockDim.x * blockDim.y * 4)) + (threadIdx.y * blockDim.x * 4 )+ threadIdx.x * 4 + l;
	return threadId;
}

__device__ int getGlobalIdx_3D_1D() {
int blockId = blockIdx.x + blockIdx.y * gridDim.x + gridDim.x * gridDim.y * blockIdx.z;
int threadId = blockId * blockDim.x + threadIdx.x;
return threadId;
}

__device__ int getGlobalIdx_3D_3D()
{
	int blockId = blockIdx.x + blockIdx.y * gridDim.x + gridDim.x * gridDim.y * blockIdx.z;
	int threadId = blockId * (blockDim.x * blockDim.y * blockDim.z ) + (threadIdx.z * (blockDim.x * blockDim.y )) + (threadIdx.y * blockDim.x )+ threadIdx.x ;
	return threadId;
}

__device__ int getGlobalIdx_1D_1D()
{
	return blockIdx.x *blockDim.x + threadIdx.x;
}

__device__ int check_position (float3 *d_gr_centerp, float x, float y, float z, int * grain)
	// This function checks to see if the atom's position is
	// closer to the center of the current grain than it is to 
	// any other grain. IF so, check is assigned 1. If not, check 
	// is assigned 0.
{
	int i, check;
	float r12,r22,dx,dy,dz;

	check = 1;

	//check if atom is outside the outer periodic image cells
	if (x >= 2.0*d_lx || x <= -d_lx) 
		check = 0;
	else if (y >= 2.0*d_ly || y <= -d_ly)
		check = 0;
	else if (z >= 2.0*d_lz || z <= -d_lz) 
		check = 0;
	if (check == 0) return 0;

	dx = d_gr_centerp[*grain].x - x;
	dy = d_gr_centerp[*grain].y - y;
	dz = d_gr_centerp[*grain].z - z;

	//check if atom is nearest to the actual grain center (and not it's image)
	if (abs(dx) > 0.5*d_lx) 
		check = 0;
	else if (abs(dy) > 0.5*d_ly) 
		check = 0;
	else if (abs(dz) > 0.5*d_lz) 
		check = 0;
	if (check == 0) return 0;

	//check if atom is closest to current grain center 
	r12 = dx*dx+dy*dy+dz*dz;
	for ( i = 0; i < d_numgrain; ++i)
	{
		if (i == *grain) continue;
		dx = d_gr_centerp[i].x - x;
		dy = d_gr_centerp[i].y - y;
		dz = d_gr_centerp[i].z - z;
		dx = separation_PBC(dx,d_lx);
		dy = separation_PBC(dy,d_ly);
		dz = separation_PBC(dz,d_lz);
		r22 = dx*dx+dy*dy+dz*dz;
		if (r22 <= r12) 
		{	check = 0;
		break;
		}
	}

	return check;
}


// Heavy calculation
__global__ void assign_initial_positions(float3 *d_gr_centerp, float3 *d_r, float *d_alpha, float *d_beta, float *d_gama, bool *d_atom_id,int *d_grain, int *d_l1)
{
	int check,l1;
	int i = d_nx[0] + threadIdx.x + blockIdx.x * blockDim.x;
	int j = d_ny[0] + threadIdx.y + blockIdx.y * blockDim.y;
	int k = d_nz[0] + threadIdx.z + blockIdx.z * blockDim.z;
	int n1;
	float x1,y1,z1,x_rot,y_rot,z_rot;
    float h11,h12,h13;
    float h21,h22,h23;
    float h31,h32,h33;
	float basis[4][3];
	float phi1, phi2, phi3;


	// Distribution threads here
	if ( i < d_nx[1] && j < d_ny[1] && k < d_nz[1] )
	{
		basis[0][0]=0.00;
		basis[0][1]=0.00;
		basis[0][2]=0.00;
		basis[1][0]=0.50;
		basis[1][1]=0.50;
		basis[1][2]=0.00;
		basis[2][0]=0.00;
		basis[2][1]=0.50;
		basis[2][2]=0.50;
		basis[3][0]=0.50;
		basis[3][1]=0.00;
		basis[3][2]=0.50;

		// Tilt the grains.
		phi1 = d_alpha[*d_grain]*M_PI/180.0;
		phi2 = d_beta[*d_grain]*M_PI/180.0;
		phi3 = d_gama[*d_grain]*M_PI/180.0;

	    h11=cos(phi1)*cos(phi3)-sin(phi1)*sin(phi3)*cos(phi2);
	    h12=sin(phi1)*cos(phi3)+cos(phi1)*sin(phi3)*cos(phi2);
	    h13=sin(phi3)*sin(phi2);
	    h21=-(cos(phi1)*sin(phi3)+sin(phi1)*cos(phi3)*cos(phi2));
	    h22=-sin(phi1)*sin(phi3)+cos(phi1)*cos(phi3)*cos(phi2);
	    h23=cos(phi3)*sin(phi2);
	    h31=sin(phi1)*sin(phi2);
	    h32=-cos(phi1)*sin(phi2);
	    h33=cos(phi2);

			l1=*d_l1;
			n1 = getGlobalIdx_3D_3D();
			x1 = i*d_a + basis[l1][0]*d_a;
			y1 = j*d_a + basis[l1][1]*d_a;
			z1 = k*d_a + basis[l1][2]*d_a;
			x_rot = (x1*h11 + y1*h21 + z1*h31)+ d_gr_centerp[*d_grain].x;
			y_rot = (x1*h12 + y1*h22 + z1*h32)+ d_gr_centerp[*d_grain].y;
			z_rot = (x1*h13 + y1*h23 + z1*h33)+ d_gr_centerp[*d_grain].z;

			check = check_position(d_gr_centerp, x_rot,y_rot,z_rot,d_grain);
			if (check == 1)	
			{
				d_r[n1].x = pos_PBC(x_rot,d_lx);
				d_r[n1].y = pos_PBC(y_rot,d_ly);
				d_r[n1].z = pos_PBC(z_rot,d_lz);
				d_atom_id[n1]  = 1;
			}

	}
	//__syncthreads();


}

__global__ void get_GBlist(float3 *d_gr_centerp, float3 *d_r, int *d_atom_grain,int *d_atom_neigh, bool *d_tag)
{
	int i, j, mygrain;
	float dx1, dx2, dx, dy1, dy2, dy, dz1, dz2, dz, r12, r22, r32, r1, r3;
	float co, projec, dis;
	float d_mini;
	i = getGlobalIdx_3D_3D();
	if ( i < d_natom)
	{
		mygrain = d_atom_grain[i];
		dx1 = separation_PBC(d_r[i].x - d_gr_centerp[mygrain].x,d_lx);
		dy1 = separation_PBC(d_r[i].y - d_gr_centerp[mygrain].y,d_ly);
		dz1 = separation_PBC(d_r[i].z - d_gr_centerp[mygrain].z,d_lz);
		r12 = dx1*dx1+dy1*dy1+dz1*dz1;
		r1 = sqrt(r12);
		d_mini = d_a;
		for ( j = 0; j < d_numgrain; ++j)
		{
			if ( j == mygrain ) continue;
			dx = separation_PBC(d_r[i].x - d_gr_centerp[j].x,d_lx);
			dy = separation_PBC(d_r[i].y - d_gr_centerp[j].y,d_ly);
			dz = separation_PBC(d_r[i].z - d_gr_centerp[j].z,d_lz);
			r22 = dx*dx+dy*dy+dz*dz;
			dx2 = separation_PBC(d_gr_centerp[mygrain].x - d_gr_centerp[j].x,d_lx);
			dy2 = separation_PBC(d_gr_centerp[mygrain].y - d_gr_centerp[j].y,d_ly);
			dz2 = separation_PBC(d_gr_centerp[mygrain].z - d_gr_centerp[j].z,d_lz);
			r32 = dx2*dx2+dy2*dy2+dz2*dz2;
			r3 = sqrt(r32);
			// What's this?
			co = (r12+r32-r22)/2.0/r1/r3;
			projec = r1*co;
			dis = r3/2.0 - projec;
			if (i == 0 && j == 0)
				printf("d_natom = %d, dis = %f\n",d_natom,dis);
//			if (dis < 0.22*d_a)
			if (dis < 0.27*d_a)
			{			
				d_tag[i] = 1;
				if (dis < d_mini)
				{
					d_mini = dis;
					d_atom_neigh[i] = j;
				}
			}
		}
	}
}

__global__ void clean_grain_boundaries(float3 *d_r, int *d_atom_grain,int *d_atom_neigh, int *d_GBlist, bool *d_tag)
{

	int i, j, ii, jj;
	float dx, dy, dz, dr2;
	float a2 = d_a*d_a;
	ii = getGlobalIdx_3D_1D();
//	d_tag[d_GBlist[getGlobalIdx_1D_1D()]] = getGlobalIdx_1D_1D();
//	if ( ii < d_GBatoms/10 && ii > d_GBatoms/100)
		if (ii < d_GBatoms)
	{	
		i = d_GBlist[ii];
		for ( jj = ii+1; jj < d_GBatoms; ++jj)
		{
			j = d_GBlist[jj];
			if (d_atom_neigh[i] != d_atom_grain[j] || d_atom_neigh[j] != d_atom_grain[i])
				continue;
			dx = d_r[i].x - d_r[j].x;
			dy = d_r[i].y - d_r[j].y;
			dz = d_r[i].z - d_r[j].z;
			dx = separation_PBC(dx,d_lx);
			dy = separation_PBC(dy,d_ly);
			dz = separation_PBC(dz,d_lz);
			dr2 = dx*dx+dy*dy+dz*dz;
//				if (dr2 <= 0.17*a2)
			if (dr2 <= 0.215*a2)
			{
				d_tag[i] = 1;
//				d_drlist[ii] = dr2;
				break;
			}
		}
	}

}

void  create_sample()
{
	int i, grain,l;
	int nx[2],ny[2],nz[2];
	int aindex;
	float3 *d_gr_centerp, *d_r;
	float3 *rr;
//	float *d_atom_mini;

	int *d_atom_grain, *d_atom_neigh, *d_grain, *d_l1; 
	float *d_alpha ,*d_beta ,*d_gama;
	bool *d_atom_id;
	nx[0] = int(-0.88*(lx/a));
	nx[1] = int(0.88*(lx/a));
	ny[0] = int(-0.88*(ly/a)*ratio);
	ny[1] = int(0.88*(ly/a)*ratio);
	nz[0] = int(-0.88*(lz/a)*ratio);
	nz[1] = int(0.88*(lz/a)*ratio);
	DIM.x = nx[1]-nx[0]+1;
	DIM.y = ny[1]-ny[0]+1;
	DIM.z = nz[1]-nz[0]+1;
	printf("DIM.x = %d, DIM.y = %d, DIM.z = %d\n", DIM.x, DIM.y, DIM.z);
	if (DIM.x <10) DIM.x = 10;
	if (DIM.y <10) DIM.y = 10;
	if (DIM.z <10) DIM.z = 10;

	nmax = DIM.x * DIM.y * DIM.z;
	printf("Maximum atom number: %d\n", nmax);
	r = new float3[nmax/5*4];
	rr = new float3[nmax];
	atom_id = new bool[nmax];
	memset(atom_id, 0, nmax*sizeof(bool));
	atom_grain = new int[nmax/5*4];
//	atom_neigh = new int[nmax];
//	atom_mini = new float[nmax];
	

	clock_t begin = clock();
	// CUDA memoray allocation
	cudaMalloc(&d_r, nmax*sizeof(float3));
	cudaMalloc(&d_atom_id, nmax*sizeof(bool));
	cudaMalloc(&d_grain, 1*sizeof(int));
	cudaMalloc(&d_gr_centerp, numgrain*sizeof(float3));
	cudaMalloc(&d_alpha, numgrain*sizeof(float));
	cudaMalloc(&d_beta, numgrain*sizeof(float));
	cudaMalloc(&d_gama, numgrain*sizeof(float));
	cudaMalloc(&d_l1, 1*sizeof(int));




	// Device initiallization
	cudaMemcpyToSymbol(d_lx, &lx, sizeof(float));
	cudaMemcpyToSymbol(d_ly, &ly, sizeof(float));
	cudaMemcpyToSymbol(d_lz, &lz, sizeof(float));
	cudaMemcpyToSymbol(d_numgrain, &numgrain, sizeof(int));
	cudaMemcpyToSymbol(d_nmax, &nmax, sizeof(int));
	cudaMemcpyToSymbol(d_a, &a, sizeof(float));
	cudaMemcpyToSymbol(d_nx, nx, 2*sizeof(int));
	cudaMemcpyToSymbol(d_ny, ny, 2*sizeof(int));
	cudaMemcpyToSymbol(d_nz, nz, 2*sizeof(int));
	cudaMemcpy(d_gr_centerp, gr_centerp, numgrain*sizeof(float3), cudaMemcpyHostToDevice);
	cudaMemcpy(d_alpha, alpha, numgrain*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_beta, beta, numgrain*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_gama, gama, numgrain*sizeof(float), cudaMemcpyHostToDevice);
//	cudaMemcpy(d_atom_id, atom_id, nmax*sizeof(bool), cudaMemcpyHostToDevice);


	//if (DIM.z > 640)
	//{
	//	printf("Warning: DIM.z is larger than 640, set DIM.z to 640");
	//	DIM.z = 640;
	//}
	dim3 blocks((DIM.x+8-1)/8, (DIM.y+8-1)/8, (DIM.z+8-1)/8);
	dim3 threads(8, 8, 8);

	// Initlal positions
	for ( grain=0; grain < numgrain; ++grain )
	{
		printf ("%d\n",grain);
		cudaMemcpy(d_grain, &grain, 1*sizeof(int), cudaMemcpyHostToDevice);
		//CUDA//
	for (l=0;l<4;++l){
		cudaMemset(d_atom_id, 0, nmax*sizeof(bool));
		cudaMemcpy(d_l1, &l, 1*sizeof(int), cudaMemcpyHostToDevice);
		assign_initial_positions<<< blocks, threads >>>(d_gr_centerp, d_r, d_alpha, d_beta, d_gama, d_atom_id, d_grain,d_l1);
		//CUDA END//
		cudaMemcpy(rr, d_r, nmax*sizeof(float3), cudaMemcpyDeviceToHost);
		cudaMemcpy(atom_id, d_atom_id, nmax*sizeof(bool), cudaMemcpyDeviceToHost);
		for (i =0; i < nmax; ++i)
		{
			if (atom_id[i] == 1)	
			{
				count = count+1;
				aindex = count-1;
				r[aindex].x = rr[i].x;
				r[aindex].y = rr[i].y;
				r[aindex].z = rr[i].z;
				atom_grain[aindex] = grain;			
			}
		}
		}
	}
	free(rr);



	natom = count;
	printf ("Initial atom number: %d\n",natom);
	
	clock_t end = clock();
	float elapsed_secs = float(end - begin) / CLOCKS_PER_SEC;
	cout << "Time elapsed: " << elapsed_secs << " s" << endl;
	
	clock_t begin1 = clock();
	//Clean grain boundaries 
	cudaFree(d_r);
	cudaFree(d_atom_id);
	cudaFree(d_alpha);
	cudaFree(d_beta);
	cudaFree(d_gama);

	int * d_GBlist, * GBlist;
	int counter2,counter3;
	bool *d_tag, *tag;
//	float * d_drlist;

	int GBatoms = 0;
	GBlist = new int[natom];
	tag = new bool[natom];
	atom_neigh = new int[natom];
//	atom_mini = new float[natom];


	cudaMalloc(&d_r, natom*sizeof(float3));
	cudaMalloc(&d_atom_grain, natom*sizeof(int));
	cudaMalloc(&d_atom_neigh, natom*sizeof(int));
//	cudaMalloc(&d_atom_mini, natom*sizeof(float));
	cudaMalloc(&d_tag, natom*sizeof(bool));

	cudaMemset(d_tag, 0, natom*sizeof(bool));
	cudaMemcpyToSymbol(d_natom, &natom, sizeof(int));
	cudaMemcpy(d_r, r, natom*sizeof(float3), cudaMemcpyHostToDevice);
	cudaMemcpy(d_atom_grain, atom_grain, natom*sizeof(int), cudaMemcpyHostToDevice);

	dim3 blocks2((natom+32768-1)/32768, 8, 8);
	dim3 threads2(8, 8, 8);

	get_GBlist <<< blocks2, threads2 >>> (d_gr_centerp, d_r, d_atom_grain, d_atom_neigh, d_tag);
	cudaMemcpy(tag, d_tag, natom*sizeof(bool), cudaMemcpyDeviceToHost);
	for ( i = 0; i < natom; ++i)
	{
if (tag[i] == 1)
		{
			++GBatoms;
			GBlist[GBatoms-1] = i;
	//cout << i << " " << tag[i] << endl;
		}
	}
	printf ("GBatoms: %d\n",GBatoms);

	clock_t end1 = clock();
	float elapsed_secs1 = float(end1 - begin1) / CLOCKS_PER_SEC;
	cout << "Time elapsed: " << elapsed_secs1 << " s" << endl;
	
//	drlist = new float[GBatoms];
	clock_t begin2 = clock();

//	cudaFree(d_atom_mini);
	cudaFree(d_gr_centerp);

	cudaMemset(d_tag, 0, natom*sizeof(bool));
	cudaMalloc(&d_GBlist, GBatoms*sizeof(int));
	cudaMemcpy(d_GBlist, GBlist, GBatoms*sizeof(int), cudaMemcpyHostToDevice);

//	cudaMalloc(&d_drlist, GBatoms*sizeof(float));
//	cudaMemset(d_drlist, 100 , GBatoms*sizeof(float));
	for (i=0;i<GBatoms;++i){tag[GBlist[i]]=0;}
//	cudaMemcpy(d_drlist, drlist, GBatoms*sizeof(float), cudaMemcpyHostToDevice);
//	cudaMemcpy(drlist, d_drlist, GBatoms*sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(GBlist, d_GBlist, GBatoms*sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpyToSymbol(d_GBatoms, &GBatoms, sizeof(int));
//	cudaMemcpy(tag, d_tag, natom*sizeof(int), cudaMemcpyDeviceToHost);
	printf("%d %d\n",natom,GBatoms);
//	for (i=0;i<GBatoms;++i){cout << tag[GBlist[i]] << " " << drlist[i] << " " << GBlist[i] << endl;}
	dim3 blocks3((GBatoms+16384-1)/16384,8,8);
	dim3 threads3(256);
	clean_grain_boundaries <<< blocks3, threads3 >>> (d_r, d_atom_grain, d_atom_neigh, d_GBlist, d_tag);
//cudaError_t error = cudaGetLastError();
//printf("CUDA error: %s\n", cudaGetErrorString(error));
	counter2 = -1;
	counter3 = 0;
//	cudaMemcpy(drlist, d_drlist, GBatoms*sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(tag, d_tag, natom*sizeof(bool), cudaMemcpyDeviceToHost);
//	for (i=0;i<GBatoms;++i){cout << tag[GBlist[i]] << " " << drlist[i] << " " << GBlist[i] << endl;}

//	cout << natom << endl;
	for (i = 0; i < natom; ++i)
	{
		if (tag[i] == 0) 
		{
			counter2 = counter2 + 1;
			r[counter2].x = r[i].x;
			r[counter2].y = r[i].y;
			r[counter2].z = r[i].z;
		        atom_grain[counter2] = atom_grain[i];			
//			cout << i << " " << tag[i] << " " << counter2 << endl;
		}
//			else {counter3++;cout << i << " " << tag[i] << " " << counter2 << endl;}
				else {counter3++;}
	}
	natom = counter2+1;
	printf ("Atom number (after cleaning): %d %d\n",natom,counter3);
	
	clock_t end2 = clock();
	float elapsed_secs2 = float(end2 - begin2) / CLOCKS_PER_SEC;
	cout << "Time elapsed: " << elapsed_secs2 << " s" << endl;

}

int write_output_files(char* ofn)
{
	int i;
	ofstream ofile;
	ofile.open(ofn);
	ofile << "# Position data for Cu system" << endl;
	ofile << "" << endl;
	ofile << natom << " atoms" << endl;
	ofile << "1 atom types" << endl;
	ofile << "" << endl;
	ofile << "0 " << lx << " xlo xhi" << endl;
	ofile << "0 " << ly << " ylo yhi" << endl;
	ofile << "0 " << lz << " zlo zhi" << endl;
	ofile << "" << endl;
	ofile << "Masses" << endl;
	ofile << "" << endl;
	ofile << "1 63.55" << endl;
	ofile << "" << endl;
	ofile << "Atoms" << endl;
	ofile << "" << endl;
	for ( i = 0; i < natom; ++i)
	{
		ofile << i+1 << " 1 ";
		ofile << setprecision(6) << r[i].x <<  " " << r[i].y << " " << r[i].z << endl;
	}
	ofile.close();
	return 0;
}

void write_output_cfg(char* ofn)
{
	int i;
	ofstream ofile;
	ofile.open(ofn);
	ofile << "Number of particles = " << natom << endl;
	ofile << "A = 1 Angstrom (basic length-scale)" << endl;
	ofile << "H0(1,1) = " << lx << " A" << endl;
	ofile << "H0(1,2) = 0 A" << endl;
	ofile << "H0(1,3) = 0 A" << endl;
	ofile << "H0(2,1) = 0 A" << endl;
	ofile << "H0(2,2) = " << ly << " A" << endl;
	ofile << "H0(2,3) = 0 A" << endl;
	ofile << "H0(3,1) = 0 A" << endl;
	ofile << "H0(3,2) = 0 A" << endl;
	ofile << "H0(3,3) = " << lz << " A" << endl;
	ofile << ".NO_VELOCITY." << endl;
	ofile << "entry_count = 4" << endl;
	ofile << "auxiliary[0] = grain" << endl;
	ofile << "63.55" << endl;
	ofile << "Cu" << endl;
	for ( i = 0; i < natom; ++i)
	{
		ofile << setprecision(5) << r[i].x/lx <<  " " << r[i].y/ly << " " << r[i].z/lz << " " << atom_grain[i] << endl;
	}
	ofile.close();
}

int main(int argc, char* argv[])
{
	int deviceCount;

	cudaGetDeviceCount(&deviceCount);
	printf("Number of GPU devices: %d\n", deviceCount);


	clock_t begin3 = clock();


	char* ofn; //output filename
	char* ifn; //input filename
	if (argc < 3 || strncmp(argv[1],"-h",2) == 0 || strncmp(argv[1],"--help",6) == 0)
	{cout << "./ggp input output" << endl;cout << "Example: ./ggp input.txt a.out" << endl;return 0;}
	ifn = argv[1];
	ofn = argv[2];
	read_config(ifn);
	create_sample();
	
	clock_t end3 = clock();
	
	float elapsed_secs3 = float(end3 - begin3) / CLOCKS_PER_SEC;

	cout << "Total time elapsed: " << elapsed_secs3 << " s" << endl;
	
	cout << "Writing file..." << endl;
//	write_output_files(ofn);
	write_output_cfg(ofn);
	cout << "Done" << endl;
	
	

	
	return 0;

}

