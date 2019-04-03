#include <bits/stdc++.h>

using namespace std;

const int MAX_FES = 50000;
const int NL = 10;
const int LS = 4;
const int num_points = NL*LS;
const double phi = 0.5;
const double tol = 0.0001;

double eval(double x, double y) {
    return sin(x)*cos(y);
}

bool comp(vector<double> temp1, vector<double> temp2) {
    return temp1[2] < temp2[2];
}

int main() {
    int fes = 0;
    int coor_low_lim = -2;
    int coor_high_lim = 2;
    double vel_low_lim = -0.1;
    double vel_high_lim = 0.1;

    double best_x;
    double best_y;
    double best_fx;

    vector<vector<double> > coor (num_points, vector<double> (6, 0));
    // MAKE LEVELS INTO A MAP
    vector<vector<int> > levels(NL, vector<int> (LS, 0));
    
    // Seeding random
    srand (static_cast <unsigned> (time(0)));

    for(int i=0; i< num_points; i++) {
        // range from [-2, 1]
        coor[i][0] = (coor_low_lim + static_cast <double> (rand()) /( static_cast <double> (RAND_MAX/(coor_high_lim - coor_low_lim))));
        coor[i][1] = (coor_low_lim + static_cast <double> (rand()) /( static_cast <double> (RAND_MAX/(coor_high_lim - coor_low_lim))));
        coor[i][2] = eval(coor[i][0], coor[i][1]);
        coor[i][4] = (vel_low_lim + static_cast <double> (rand()) /( static_cast <double> (RAND_MAX/(vel_high_lim - vel_low_lim))));
        coor[i][5] = (vel_low_lim + static_cast <double> (rand()) /( static_cast <double> (RAND_MAX/(vel_high_lim - vel_low_lim))));
    }
    // cout << "BEST: "<< best_fx << endl; 
    fes += num_points;
    
    while(fes < MAX_FES) {
        cout << "ITER " << fes<< endl;
        sort(coor.begin(), coor.end(), comp);
        best_fx = coor[0][2];
        best_x = coor[0][0];
        best_y = coor[0][1];
        // Segregating into levels with different level sizes
        // for(int i=0; i<coor.size(); i++) {
        //     for(int k=0; k<NL; k++) {
        //         double factor = (double)(coor_high_lim - coor_low_lim - 1)/NL;
        //         if(coor[i][2] >= coor_low_lim + (k*factor) && coor[i][2] <= coor_low_lim + ((k+1)*factor)) {
        //             coor[i][3] = k+1;
        //             break;
        //         }
        //     }
        // }
        for(int i=0; i<coor.size(); i++) {
            cout << coor[i][0] << endl;
            cout << coor[i][1] << endl;
        }
        // Segregating into levels with same level size
        for(int i=0; i<coor.size(); i++) {
            // With 0-indexed: i/LS
            coor[i][3] = i/LS;
            levels[i/LS][i%LS] = i;
        }
        cout << "TEMP" << endl;
        for(int i=0; i<levels.size(); i++) {
            for(int j=0; j<levels[0].size(); j++) {
                cout << levels[i][j] << " ";
            }
            cout << endl;
        }

        for(int i=NL; i>=3; i--) {
            for(int k=0; k<LS; k++) {
                // Check if the level is the same
                // levels[i-1][k]
                    int pt = levels[i-1][k];
                    int lev1 = 0 + (rand() % static_cast<int>(i));
                    int lev2 = 0 + (rand() % static_cast<int>(i));
                    // cout << i << " "<< lev1 << "  "<<  lev2 <<endl; 
                    if(lev2 < lev1) {
                        swap(lev2, lev1);
                    }
                    int pt1 = 0 + (rand() % static_cast<int>(LS));
                    int pt2 = 0 + (rand() % static_cast<int>(LS));
                    // cout << " TEST1: " << pt1 << "  "<< pt2 << endl; 
                    int temp1 = levels[lev1][pt1];
                    int temp2 = levels[lev2][pt2];
                    // cout << " TEST2: " << temp1 << "  "<< temp2 << endl; 


                    // check

                    int r1 = ((double) rand() / (RAND_MAX));
                    int r2 = ((double) rand() / (RAND_MAX));
                    int r3 = ((double) rand() / (RAND_MAX));


                    coor[pt][4] = r1*coor[pt][4] + r2*(coor[temp1][0] - coor[pt][0]) + phi*r3*(coor[temp2][0] - coor[pt][0]);
                    coor[pt][5] = r1*coor[pt][5] + r2*(coor[temp1][1] - coor[pt][1]) + phi*r3*(coor[temp2][1] - coor[pt][1]);
                    // cout << "VELOCITY:  "<< coor[pt][4] << "  "<< coor[pt][5] << endl;
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


        // cout << "RESULTS: " << endl;
        // for(int i=0; i<coor.size(); i++) {
        //     cout << coor[i][2] << endl;
        // }

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
            // cout << "VELOCITY:  "<< coor[pt][4] << "  "<< coor[pt][5] << endl;
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
        
        // cout << "RESULTS: "<< endl;
        // for(int i=0; i< coor.size(); i++) {
        //     cout << coor[i][2] << " " <<coor[i][4] << " " <<coor[i][5] << endl;
        // }
    }
    cout << "FINAL RESULTS: " << endl;
    cout << best_fx;
    return 0;
}