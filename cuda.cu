#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/generate.h>
#include <thrust/sort.h>
#include <thrust/copy.h>
#include <bits/stdc++.h>
#include <curand.h>
#include <curand_kernel.h>


using namespace std;

const int MAX_FES = 50000;
int NL = 10;
int LS = 4;
int dim1 = 20;
int dim2 = 1000000;
const int num_points = NL*LS;
const double phi = 0.1;
const double tol = 0.000001;

struct coor {
    double a0,a1,a2,a3,a4,a5,a6;
};

struct coor1 {
    double* a1;
    double a2;
    double a3;
    double* a4;
    double a5;
    double a0;
};

__device__ void func_eval_input_data(coor1* d_coo, double* y, double* sum)
{
    // int tId = threadIdx.x + (blockIdx.x * blockDim.x);    
    double temp = (*y);
    *sum *= temp;
    
}

__host__ __device__ void func_eval_input_data(coor1* d_coo, double* x, double* y, double* total) {
    double* d_sum, *sum;
    
    sum = (double *) malloc(sizeof(double));
    *sum = 1;
    
    
    // cudaMalloc((void **)&sum,sizeof(double));
    
    // int tId = threadIdx.x + (blockIdx.x * blockDim.x);
    
    // func_eval_input_data<<<(*x)/512,512>>>(d_coo, y, d_sum);
    
    
    
    *total += *sum;
    
    free(sum);
  
}

__host__ __device__ double eval(double x, double y) {
    // Function to be evaluated. Can be changed
    return sin(x)*cos(y);
}

/*
bool comp(coor temp1, coor temp2) {
    // Comparator function for std::sort()
    return temp1.a2 < temp2.a2;
}
*/

__host__ __device__ bool operator < (const coor& a,const coor& b){return a.a2 < b.a2;}



__global__ void noicetoit(coor* h_coo, int* levels, double* best_x, double* best_y, double* best_fx, int* i, int* LS, int* NL) {

    int tId = threadIdx.x + (blockIdx.x * blockDim.x);
    curandState state;
    curand_init((unsigned long long)clock() + tId, 0, 0, &state);
    
    int index = threadIdx.x + blockIdx.x * blockDim.x;

    double rand1 = curand_uniform_double(&state);
    double rand2 = curand_uniform_double(&state);
    double rand3 = curand_uniform_double(&state);
    double rand4 = curand_uniform_double(&state);
    double rand5 = curand_uniform_double(&state);
    double rand6 = curand_uniform_double(&state);
    double rand7 = curand_uniform_double(&state);
    
    int pt = levels[(*i-1)*(*NL) + index];
    // Choosing a random level
    int lev1 = (int)rand1 % (*i);
    int lev2 = (int)rand2 % (*i);

    if(lev2 < lev1) {
        int temp = lev2;
        lev2 = lev1;
        lev1 = temp;
    }
    // Choosing random points from those levels
    int pt1 = (int)rand3 % (*LS);
    int pt2 = (int)rand4 % (*LS);
    int temp1 = levels[lev1*(*NL) + pt1];
    int temp2 = levels[lev2*(*NL) + pt2];

    int r1 = ((double) rand5 / (RAND_MAX));
    int r2 = ((double) rand6 / (RAND_MAX));
    int r3 = ((double) rand7 / (RAND_MAX));

    // Update Functions
    
    h_coo[pt].a4 = r1*(h_coo[pt].a4) + r2*((h_coo[temp1].a0) - (h_coo[pt].a0)) + phi*r3*((h_coo[temp2].a0) - (h_coo[pt].a0));
    h_coo[pt].a5 = r1*h_coo[pt].a5 + r2*(h_coo[temp1].a1 - h_coo[pt].a1) + phi*r3*(h_coo[temp2].a1 - h_coo[pt].a1);
    h_coo[pt].a0 = h_coo[pt].a0 + h_coo[pt].a4;
    h_coo[pt].a1 = h_coo[pt].a1 + h_coo[pt].a5;

    double fx = eval(h_coo[pt].a0, h_coo[pt].a1);
    if(abs(fx - *best_fx) < tol) {
        *best_x = h_coo[pt].a0;
        *best_y = h_coo[pt].a1;
        *best_fx = fx;
    }
    h_coo[pt].a2 = fx;


}

int main() {
    
    clock_t tStart = clock();
    
    int fes = 0;
    // Domain for coordinates
    int coor_low_lim = -10;
    int coor_high_lim = 50;
    // Domain for velocities
    double vel_low_lim = -0.1;
    double vel_high_lim = 0.1;

    double best_x;
    double best_y;
    double best_fx;

    
    //parallel programming
    coor* h_coo, *d_coo;
    
    h_coo = (coor *)malloc(num_points*sizeof(coor));
    
    coor1* h_coo1, *d_coo1;
    
    h_coo1 = (coor1 *)malloc(num_points*sizeof(coor));
    
    for(int i=0; i<num_points; i++)
    {
        h_coo1[i].a1 = (double *)malloc(dim1*sizeof(double));
        h_coo1[i].a4 = (double *)malloc(dim1*sizeof(double));
    }
    
    cudaMalloc((void **)&d_coo, num_points*sizeof(coor));
    
    int *levels, *d_levels;
    levels = (int *)malloc(sizeof(int)*NL*LS);
    
    cudaMalloc((void **)&d_levels, sizeof(int)*NL*LS);
    
//    vector<vector<int> > levels(NL, vector<int> (LS, 0));
    
    // Seeding random
    srand(static_cast <unsigned> (time(0)));

    
    for(int i=0; i< num_points; i++) {
        /*
        coor[i][...] contains property of each point.
        coor[i][0] -> x-coordinate
        coor[i][1] -> y-coordinate
        coor[i][2] -> function evaluation
        coor[i][3] -> level number
        coor[i][4] -> x-velocity
        coor[i][5] -> y-velocity
        */
        h_coo[i].a0 = (coor_low_lim + static_cast <double> (rand()) /( static_cast <double> (RAND_MAX/(coor_high_lim - coor_low_lim))));
        h_coo[i].a1 = (coor_low_lim + static_cast <double> (rand()) /( static_cast <double> (RAND_MAX/(coor_high_lim - coor_low_lim))));
        h_coo[i].a2 = eval(h_coo[i].a0, h_coo[i].a1);
        h_coo[i].a4 = (vel_low_lim + static_cast <double> (rand()) /( static_cast <double> (RAND_MAX/(vel_high_lim - vel_low_lim))));
        h_coo[i].a5 = (vel_low_lim + static_cast <double> (rand()) /( static_cast <double> (RAND_MAX/(vel_high_lim - vel_low_lim))));
    }
    fes += num_points;
    
    cudaMemcpy(d_coo, h_coo, num_points*sizeof(coor), cudaMemcpyHostToDevice);
  
    
    while(fes < MAX_FES) {
        
        
        thrust::sort(h_coo, h_coo + num_points);
 
        
        
        
        best_fx = h_coo[0].a2;
        best_x = h_coo[0].a0;
        best_y = h_coo[0].a1;
        

        // Segregating points into levels
        for(int i=0; i<num_points; i++) {
            /*
            Levels basically acts as a lookup for coor array.
            Dimensions of levels: levels[NL][LS]
            levels[i] denotes (i+1)th level and each element in levels[i][...] denotes the number of
            the point in the coor array.

            For instance,
            levels[1][2] = 5 denotes that the 3rd point in 2nd level corresponds to point 5 in coor array,
            i.e. coor[5][...]
            */
            
            h_coo[i].a3 = i/LS;
            levels[(i/LS)*NL + i%LS] = i;
        }

        cudaMemcpy(d_coo, h_coo, num_points*sizeof(coor), cudaMemcpyHostToDevice);
        double *d_best_x, *d_best_y,  *d_best_fx;
        int* d_LS, *d_i, *d_NL;

        cudaMalloc((void **)&d_best_x, sizeof(double));
        cudaMalloc((void **)&d_best_y, sizeof(double));
        cudaMalloc((void **)&d_best_fx, sizeof(double));
        cudaMalloc((void **)&d_i, sizeof(int));
        cudaMalloc((void **)&d_LS, sizeof(int));
        cudaMalloc((void **)&d_NL, sizeof(int));

        cudaMemcpy(d_best_x, &best_x, sizeof(double), cudaMemcpyHostToDevice);
        cudaMemcpy(d_best_y, &best_y, sizeof(double), cudaMemcpyHostToDevice);
        cudaMemcpy(d_best_fx, &best_fx, sizeof(double), cudaMemcpyHostToDevice);
        
        cudaMemcpy(d_LS, &LS, sizeof(int), cudaMemcpyHostToDevice);
        cudaMemcpy(d_NL, &NL, sizeof(int), cudaMemcpyHostToDevice);
        cudaMemcpy(d_levels, levels, sizeof(int)*NL*LS, cudaMemcpyHostToDevice);
        
        for(int i=NL; i>=3; i--) {
            cudaMemcpy(d_i, &i, sizeof(int), cudaMemcpyHostToDevice);
            noicetoit<<<LS/512,512>>>(d_coo,d_levels, d_best_x, d_best_y,  d_best_fx, d_i, d_LS, d_NL);          
            
            fes+= LS;
        }
        
        // cudaMemcpy(h_coo, d_coo, num_points*sizeof(coor), cudaMemcpyDeviceToHost);

        for(int i=0; i<LS; i++) {
            int pt1 = 0 + (rand() % static_cast<int>(LS));
            int pt2 = 0 + (rand() % static_cast<int>(LS));
            int pt = levels[1*NL + i];

            int temp1 = levels[0*NL + pt1];
            int temp2 = levels[0*NL + pt2];
            if(abs(eval(h_coo[temp2].a0, h_coo[temp2].a1) - eval(h_coo[temp1].a0, h_coo[temp1].a1)) < tol) {
                swap(temp1, temp2);
            }

            int r1 = ((double) rand() / (RAND_MAX));
            int r2 = ((double) rand() / (RAND_MAX));
            int r3 = ((double) rand() / (RAND_MAX));
         // Update Functions
            h_coo[pt].a4 = r1*h_coo[pt].a4 + r2*(h_coo[temp1].a0 - h_coo[pt].a0) + phi*r3*(h_coo[temp2].a0 - h_coo[pt].a0);
            h_coo[pt].a5 = r1*h_coo[pt].a5 + r2*(h_coo[temp1].a1 - h_coo[pt].a1) + phi*r3*(h_coo[temp2].a1 - h_coo[pt].a1);
            h_coo[pt].a0 = h_coo[pt].a0 + h_coo[pt].a4;
            h_coo[pt].a1 = h_coo[pt].a1 + h_coo[pt].a5;

            double fx = eval(h_coo[pt].a0, h_coo[pt].a1);
            if(abs(fx - best_fx) < tol) {
                best_x = h_coo[pt].a0;
                best_y = h_coo[pt].a1;
                best_fx = fx;
            }
            h_coo[pt].a2 = fx;


        }
        fes+= LS;

    }

    cout << "Time: " << (double)(clock() - tStart)/CLOCKS_PER_SEC << endl;
    
    cout << "FINAL RESULTS: " << endl;
    cout << "Best x: " << best_x << endl;
    cout << "Best y: " << best_y << endl;
    // cout << "Best evaluation: " <<  best_fx << endl;
    return 0;
}