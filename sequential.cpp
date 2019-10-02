#include <bits/stdc++.h>

using namespace std;

const int MAX_FES = 500000;
const int NL = 10;
const int LS = 4;
const int num_points = NL*LS;
const double phi = 0.1;
const double tol = 0.000001;

double eval(double x, double y) {
    // Function to be evaluated. Can be changed
    return sin(x)*cos(y);
}

bool comp(vector<double> temp1, vector<double> temp2) {
    // Comparator function for std::sort()
    return temp1[2] < temp2[2];
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

    vector<vector<double> > coor(num_points, vector<double> (6, 0));
    vector<vector<int> > levels(NL, vector<int> (LS, 0));
    
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
        coor[i][0] = (coor_low_lim + static_cast <double> (rand()) /( static_cast <double> (RAND_MAX/(coor_high_lim - coor_low_lim))));
        coor[i][1] = (coor_low_lim + static_cast <double> (rand()) /( static_cast <double> (RAND_MAX/(coor_high_lim - coor_low_lim))));
        coor[i][2] = eval(coor[i][0], coor[i][1]);
        coor[i][4] = (vel_low_lim + static_cast <double> (rand()) /( static_cast <double> (RAND_MAX/(vel_high_lim - vel_low_lim))));
        coor[i][5] = (vel_low_lim + static_cast <double> (rand()) /( static_cast <double> (RAND_MAX/(vel_high_lim - vel_low_lim))));
    }
    fes += num_points;
    
    while(fes < MAX_FES) {
        sort(coor.begin(), coor.end(), comp);
        best_fx = coor[0][2];
        best_x = coor[0][0];
        best_y = coor[0][1];
        
        // Segregating points into levels
        for(int i=0; i<coor.size(); i++) {
            /*
            Levels basically acts as a lookup for coor array.
            Dimensions of levels: levels[NL][LS]
            levels[i] denotes (i+1)th level and each element in levels[i][...] denotes the number of
            the point in the coor array.

            For instance,
            levels[1][2] = 5 denotes that the 3rd point in 2nd level corresponds to point 5 in coor array,
            i.e. coor[5][...]
            */
            coor[i][3] = i/LS;
            levels[i/LS][i%LS] = i;
        }

        for(int i=NL; i>=3; i--) {
            for(int k=0; k<LS; k++) {
                int pt = levels[i-1][k];
                // Choosing a random level
                int lev1 = rand() % static_cast<int>(i);
                int lev2 = rand() % static_cast<int>(i);

                if(lev2 < lev1) {
                    swap(lev2, lev1);
                }
                // Choosing random points from those levels
                int pt1 = rand() % static_cast<int>(LS);
                int pt2 = rand() % static_cast<int>(LS);
                int temp1 = levels[lev1][pt1];
                int temp2 = levels[lev2][pt2];

                int r1 = ((double) rand() / (RAND_MAX));
                int r2 = ((double) rand() / (RAND_MAX));
                int r3 = ((double) rand() / (RAND_MAX));

                // Update Functions
                coor[pt][4] = r1*coor[pt][4] + r2*(coor[temp1][0] - coor[pt][0]) + phi*r3*(coor[temp2][0] - coor[pt][0]);
                coor[pt][5] = r1*coor[pt][5] + r2*(coor[temp1][1] - coor[pt][1]) + phi*r3*(coor[temp2][1] - coor[pt][1]);
                coor[pt][0] = coor[pt][0] + coor[pt][4];
                coor[pt][1] = coor[pt][1] + coor[pt][5];

                double fx = eval(coor[pt][0], coor[pt][1]);
                if(abs(fx - best_fx) < tol) {
                    best_x = coor[pt][0];
                    best_y = coor[pt][1];
                    best_fx = fx;
                }
                coor[pt][2] = fx;
                
            }
            fes+= LS;
        }

        for(int i=0; i<LS; i++) {
            int pt1 = 0 + (rand() % static_cast<int>(LS));
            int pt2 = 0 + (rand() % static_cast<int>(LS));
            int pt = levels[1][i];

            int temp1 = levels[0][pt1];
            int temp2 = levels[0][pt2];
            if(abs(eval(coor[temp2][0], coor[temp2][1]) - eval(coor[temp1][0], coor[temp1][1])) < tol) {
                swap(temp1, temp2);
            }

            int r1 = ((double) rand() / (RAND_MAX));
            int r2 = ((double) rand() / (RAND_MAX));
            int r3 = ((double) rand() / (RAND_MAX));
            coor[pt][4] = r1*coor[pt][4] + r2*(coor[temp1][0] - coor[pt][0]) + phi*r3*(coor[temp2][0] - coor[pt][0]);
            coor[pt][5] = r1*coor[pt][5] + r2*(coor[temp1][1] - coor[pt][1]) + phi*r3*(coor[temp2][1] - coor[pt][1]);
            coor[pt][0] = coor[pt][0] + coor[pt][4];
            coor[pt][1] = coor[pt][1] + coor[pt][5];

            double fx = eval(coor[pt][0], coor[pt][1]);
            if(abs(fx - best_fx) < tol) {
                best_x = coor[pt][0];
                best_y = coor[pt][1];
                best_fx = fx;
            }
            coor[pt][2] = fx;
        
        }
        fes += LS;

    }
    cout << "FINAL RESULTS: " << endl;
    cout << "Best x: " << best_x << endl;
    cout << "Best y: " << best_y << endl;
    cout << "Best evaluation: " <<  best_fx << endl;
    return 0;
}