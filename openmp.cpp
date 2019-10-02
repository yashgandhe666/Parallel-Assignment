#include <omp.h>
#include <bits/stdc++.h>

using namespace std;

const int MAX_FES = 500000;
const int NL = 10;
const int LS = 4;
const int num_points = NL*LS;
const double phi = 0.1;
const double tol = 0.000001;

struct coor {
    double a1;
    double a2;
    double a3;
    double a4;
    double a5;
    double a0;
};


double eval(double x, double y) {
    // Function to be evaluated. Can be changed
    return sin(x)*cos(y);
}


bool comp(coor temp1, coor temp2) {
    // Comparator function for std::sort()
    return temp1.a2 < temp2.a2;
}


int main() {
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
    // thrust::host_vector<thrust::host_vector<double> > coor(num_points, thrust::host_vector<double> (6, 0));
 //   thrust::host_vector<double> coor[num_points][6];
   // thrust::device_vector<double> d_coor = coor;    
    
    vector<vector<int> > levels(NL, vector<int> (LS, 0));
    
    // Seeding random
    srand(static_cast <unsigned> (time(0)));

    vector<coor> coo(num_points);
    #pragma omp parallel for default(shared)
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
        coo[i].a0 = (coor_low_lim + static_cast <double> (rand()) /( static_cast <double> (RAND_MAX/(coor_high_lim - coor_low_lim))));
        coo[i].a1 = (coor_low_lim + static_cast <double> (rand()) /( static_cast <double> (RAND_MAX/(coor_high_lim - coor_low_lim))));
        coo[i].a4 = (vel_low_lim + static_cast <double> (rand()) /( static_cast <double> (RAND_MAX/(vel_high_lim - vel_low_lim))));
        coo[i].a5 = (vel_low_lim + static_cast <double> (rand()) /( static_cast <double> (RAND_MAX/(vel_high_lim - vel_low_lim))));
    }
    for(int i=0; i<num_points; i++) {
        coo[i].a2 = eval(coo[i].a0, coo[i].a1);
    }
    fes += num_points;
    
    
    
    while(fes < MAX_FES) {
        sort(coo.begin(), coo.end(), comp);
        best_fx = coo[0].a2;
        best_x = coo[0].a0;
        best_y = coo[0].a1;
        
        // Segregating points into levels
        #pragma omp parallel for default(shared)
        for(int i=0; i<coo.size(); i++) {
            /*
            Levels basically acts as a lookup for coor array.
            Dimensions of levels: levels[NL][LS]
            levels[i] denotes (i+1)th level and each element in levels[i][...] denotes the number of
            the point in the coor array.

            For instance,
            levels[1][2] = 5 denotes that the 3rd point in 2nd level corresponds to point 5 in coor array,
            i.e. coor[5][...]
            */
            
            coo[i].a3 = i/LS;
            levels[i/LS][i%LS] = i;
        }

        for(int i=NL; i>=3; i--) {
            for(int k=0; k<LS; k++) {
                #pragma omp parallel sections {
                    #pragma omp section {
                        int pt = levels[i-1][k];
                    }
                    // Choosing a random level
                    #pragma omp section {
                        int lev1 = rand() % static_cast<int>(i);
                        int lev2 = rand() % static_cast<int>(i);
                    }
                    #pragma omp section {
                        if(lev2 < lev1) {
                            swap(lev2, lev1);
                        }
                    }
                    // Choosing random points from those levels
                    #pragma omp section {
                        int pt1 = rand() % static_cast<int>(LS);
                        int pt2 = rand() % static_cast<int>(LS);
                    }
                    #pragma omp section {
                        int temp1 = levels[lev1][pt1];
                        int temp2 = levels[lev2][pt2];

                        int r1 = ((double) rand() / (RAND_MAX));
                        int r2 = ((double) rand() / (RAND_MAX));
                        int r3 = ((double) rand() / (RAND_MAX));
                    }
                    // Update Functions
                    #pragma omp section {
                        coo[pt].a4 = r1*coo[pt].a4 + r2*(coo[temp1].a0 - coo[pt].a0) + phi*r3*(coo[temp2].a0 - coo[pt].a0);
                        coo[pt].a5 = r1*coo[pt].a5 + r2*(coo[temp1].a1 - coo[pt].a1) + phi*r3*(coo[temp2].a1 - coo[pt].a1);
                    }
                    #pragma omp section {
                        coo[pt].a0 = coo[pt].a0 + coo[pt].a4;
                        coo[pt].a1 = coo[pt].a1 + coo[pt].a5;
                    }
                    #pragma omp section {
                        double fx = eval(coo[pt].a0, coo[pt].a1);
                    }
                    if(abs(fx - best_fx) < tol) {
                        best_x = coo[pt].a0;
                        best_y = coo[pt].a1;
                        best_fx = fx;
                    }
                    coo[pt].a2 = fx;
                    
                }
                fes+= LS;
            }
        }

        for(int i=0; i<LS; i++) {
            #pragma omp parallel sections {
            #pragma omp section {
            int pt1 = 0 + (rand() % static_cast<int>(LS));
            int pt2 = 0 + (rand() % static_cast<int>(LS));
            }
            int pt = levels[1][i];

            int temp1 = levels[0][pt1];
            int temp2 = levels[0][pt2];
            if(abs(eval(coo[temp2].a0, coo[temp2].a1) - eval(coo[temp1].a0, coo[temp1].a1)) < tol) {
                swap(temp1, temp2);
            }
            #pragma omp section {
                int r1 = ((double) rand() / (RAND_MAX));
                int r2 = ((double) rand() / (RAND_MAX));
                int r3 = ((double) rand() / (RAND_MAX));
            }
            // Update Functions
                #pragma omp section {
                    coo[pt].a4 = r1*coo[pt].a4 + r2*(coo[temp1].a0 - coo[pt].a0) + phi*r3*(coo[temp2].a0 - coo[pt].a0);
                    coo[pt].a5 = r1*coo[pt].a5 + r2*(coo[temp1].a1 - coo[pt].a1) + phi*r3*(coo[temp2].a1 - coo[pt].a1);
                }
                #pragma omp section {
                    coo[pt].a0 = coo[pt].a0 + coo[pt].a4;
                    coo[pt].a1 = coo[pt].a1 + coo[pt].a5;
                }
                #pragma omp section {
                    double fx = eval(coo[pt].a0, coo[pt].a1);
                }
                if(abs(fx - best_fx) < tol) {
                    best_x = coo[pt].a0;
                    best_y = coo[pt].a1;
                    best_fx = fx;
                }
                coo[pt].a2 = fx;
                
            }
            fes+= LS;
        }
    }
    cout << "FINAL RESULTS: " << endl;
    cout << "Best x: " << best_x << endl;
    cout << "Best y: " << best_y << endl;
    cout << "Best evaluation: " <<  best_fx << endl;
    return 0;
}
    