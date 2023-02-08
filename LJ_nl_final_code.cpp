#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <random>
#include <time.h>

double r_mi(const std::vector<double>& r_input, double L){
    int m = r_input.size();
    std::vector<double> r_res(m), n(m);

    for (size_t i=0; i!=m; ++i){

        if (r_input[i] >= L/2.0){
            n[i] = -1;
        }
        else if(r_input[i] <= - L/2.0 ){
            n[i] = 1;
        }
        else {
            n[i] = 0;
        }
    
    r_res[i] = r_input[i] + n[i]*L;

     if (r_res[i] >= L/2.0){  
            while (r_res[i] >= L/2.0)
            {
                r_res[i] = r_res[i] - L;
            }
            
        }
    
    if (r_res[i] <=- L/2.0){  
            while (r_res[i] <=- L/2.0)
            {
                r_res[i] = r_res[i] + L;
            }
            
        }
   
    }
    return  r_res[0]*r_res[0] + r_res[1]*r_res[1] + r_res[2]*r_res[2];

}
int main(){
    int N=729, c;
    int accept, update, y, N_accum = 0, sweeps_count = 0, sum = 0,  total_steps = 2500, equilibration= 2000;
    int M = round(std::cbrt(N))+1; // M is first int after cubic root of N
    size_t num_of_points = M*M*M, n=3; 
    double r_cutoff = 3.0, temp = 1.2, delta_r_max = 0.2, delta_r = 0.01;
    auto L=10.14;
    double A, B, delta_V, alpha, density = N/(pow(L, 3));
    double r_skin = 3*delta_r_max;
    double r_n = r_cutoff + r_skin; 
    double intervals = L/(2*delta_r);
    int n_max = round(2*3.14*pow(r_n,3)*4/3); 

    clock_t tStart = clock();

    std::vector<int> num_of_neighbours(N); //num_of_neighbours_upd(N); // current_num_nei(N);
    std::vector<std::vector<int> > neighb_indices(N, std::vector<int>(n_max)); //neighb_indices_upd(N, std::vector<int>(n_max)); // current_ind(N, std::vector<int>(n_max)); // here we collect indices of neighbour atoms for each atom
    std::vector<std::vector<double> > nei_list_struct(N, std::vector<double>(n)); //nei_list_struct_upd(N, std::vector<double>(n)); //current structure for nei_list comparison
    std::vector<double> changed_atom(n), r_ij(n),  r_trial(n), r_init(n), histogram(intervals), g(intervals);
    std::vector<std::vector<double> > structure(num_of_points, std::vector<double>(n)); //generate a matrix of M**3 atoms with 3 coordinates
    std::vector<std::vector<double> > initial_struct(N, std::vector<double>(n)); //a matrix to compute succes rate
    
      for (size_t i = 0; i != num_of_points; ++i) { //fill the matrix
        
          structure[i][2]=L/M*(i%M + 1); // z-coordinate is growing fastest, units place in decimal notation
          structure[i][1]=L/M*((i%(M*M))/M + 1); //analog of tens place in decimal notation
          structure[i][0]=L/M*(i/(M*M) + 1); // analog of hundreds place in decimal notation 
      }

//for i in range MMM-N randomly erase atoms from structure vector
std::random_device r; // Seed with a real random value, if available
std::default_random_engine e1(r());  

while (structure.size()!=N){
    std::uniform_int_distribution<int> uniform_dist(0, structure.size()-1); // Choose a random mean between 1 and current size
    int pos = uniform_dist(e1);
    structure.erase(structure.begin()+pos); //-1 at the end not to get past the vector end position
}

initial_struct  = structure;
nei_list_struct = structure; //1st time assign nei list str;

//fill in the neighb_indices for the first time

for(int i=0; i!=N; ++i){
    c = num_of_neighbours[i];
	for(int j=0; j!=N; ++j){
        int b = num_of_neighbours[j] ;
        if(j>i){
            for (size_t z=0; z!=r_ij.size(); ++z){
                            r_ij[z] = nei_list_struct[i][z] - nei_list_struct[j][z];   
                        }  
            double dist_r_ij_n = (r_mi(r_ij, L));  //dist_r_ij_n is to differentiate from dist_r_ij further
            
            if (dist_r_ij_n<=pow(r_n, 2)){
                
                neighb_indices[i][c] = j; 
                c=c+1; 
                
                neighb_indices[j][b] = i;
                num_of_neighbours[j] +=1;
            }
        }
    }
num_of_neighbours[i]=c;
}

std::uniform_int_distribution<int> uniform_dist(0, structure.size()-1); // structure.size have changed, so we call it again
std::uniform_real_distribution<> dis(0.0, 1.0);
for (int step=0; step!=total_steps; ++step){

    // lets run 1 sweep and calculate success - proportion of atoms that have changed position
    for (int l=0; l!=N; ++l){

        int k = uniform_dist(e1);  //generate a trial move
        double wx = dis(e1);
        double wy = dis(e1);
        double wz = dis(e1);

        changed_atom[0] = structure[k][0]; // old atom, before change
        changed_atom[1] = structure[k][1];
        changed_atom[2] = structure[k][2];

        structure[k][0] = structure[k][0] + 2*delta_r_max*(wx-0.5);
        structure[k][1] = structure[k][1] + 2*delta_r_max*(wy-0.5);
        structure[k][2] = structure[k][2] + 2*delta_r_max*(wz-0.5);

        //Calculate delta V
        for (size_t i = 0; i != num_of_neighbours[k]; ++i) {  //iterate over the neighbour list of displaced atom
            auto index = neighb_indices[k][i];                 //get the index of a neigh atom 
            //calc A
            for (size_t j=0; j!=r_trial.size(); ++j){       //calculate vector from displaced atom to a neigh atom             
                r_trial[j] = structure[index][j] - structure[k][j];
            }  
            double r_trial_dist = r_mi(r_trial, L);             //calculate  r_trial[i,k] in min im conv
            
            if (r_trial_dist<=pow(r_cutoff, 2)){                 //calculate A, first part of delta V
                A = pow((1.0/r_trial_dist), 6) - pow((1.0/r_trial_dist), 3);
            }
            else {
                A =0;
                }
            //calc B
            for (size_t j=0; j!=r_init.size(); ++j){            //calculate vector r_init[i,k]
                r_init[j] = structure[index][j] - changed_atom[j];
            } 
            double r_init_dist = r_mi(r_init, L);               //calculate  r_init[i,k] in min im conv
            
            if (r_init_dist<=pow(r_cutoff, 2)){                 //calculate B, second part of delta V
                B = pow((1.0/r_init_dist), 6) - pow((1.0/r_init_dist), 3);
            }
            else {
                B =0;
                }
        
        delta_V = delta_V + 4*(A - B);

        }
        
        
        //std :: cout << delta_V << "\n";
        //rule for accepting 
        accept = 0;
        if (delta_V <= 0){
            accept = 1;
        }
        else {
            alpha = exp(-delta_V/temp);
            double w = dis(e1);
            if (alpha>=w){
                accept = 1;
        }
        }
        delta_V  = 0; 
    if (accept==0){
        structure[k][0] = changed_atom[0];
        structure[k][1] = changed_atom[1];
        structure[k][2] = changed_atom[2];
        }
    else{
        //should we update the neighbor lists?
            for (size_t j=0; j!=3; ++j){
                r_trial[j] = nei_list_struct[k][j] - structure[k][j];
            } 
            double r_trial_dist = r_mi(r_trial, L);    //calculate  r_trial[i,k] in min im conv
            
            update = 0;
            if (r_trial_dist >=pow(r_skin/2, 2)){       // update the nei lists!!
                update = 1;
                nei_list_struct = structure;
                std::vector<int> num_of_neighbours_upd(N); 
                std::vector<std::vector<int> > neighb_indices_upd(N, std::vector<int>(n_max)); 
                for(int i=0; i!=N; ++i){
                    c = num_of_neighbours_upd[i];
                    for(int j=0; j!=N; ++j){
                        int b = num_of_neighbours_upd[j] ;
                         if(j>i){
                            for (size_t z=0; z!=r_ij.size(); ++z){
                                r_ij[z] = structure[i][z] - structure[j][z];   
                            }  
                            double dist_r_ij_n = (r_mi(r_ij, L));  
                            if (dist_r_ij_n<=pow(r_n, 2)){
                                neighb_indices_upd[i][c] = j; 
                                c=c+1; 
                
                                neighb_indices_upd[j][b] = i;
                                num_of_neighbours_upd[j] +=1;
                            }
                        }
                    }   
                num_of_neighbours_upd[i]=c;          
                }
            num_of_neighbours = num_of_neighbours_upd;
            neighb_indices = neighb_indices_upd;
            }
        }
    }
        double count = 0.0;    //to count acceptance rate
        for (size_t i=0; i!=structure.size(); ++i){
            if (initial_struct[i]!=structure[i]){
                count = count +1.0;
            }
        } 

    

    //std :: cout << count/N << "\n"; -  to calculate success


    sweeps_count = sweeps_count +1;
    std::cout << sweeps_count << "\n";

        if (sweeps_count>equilibration){
            for (int i=0; i!=structure.size(); ++i){
                for (int j=0; j!=structure.size(); ++j){
                    if (j>i){
                        for (size_t z=0; z!=r_ij.size(); ++z){
                            r_ij[z] = structure[i][z] - structure[j][z];
                            
                        } 
                    double dist_r_ij = sqrt(r_mi(r_ij, L));
                    if (dist_r_ij<=L/2){
                    y = dist_r_ij/delta_r; 
                    histogram[y] = histogram[y] + 1;
                   // sum = sum +1;
                   }
                    }
                }
            }
        N_accum = N_accum + 1;
        }
    }

for (size_t i=0; i!=histogram.size(); ++i){
        g[i] = histogram[i]/(N_accum*N/2*4/3*density*3.1415*(pow(((i+1)*delta_r), 3) - pow(i*delta_r, 3)));
         sum = sum + histogram[i];
}
std::ofstream myfile ("output_g.txt");
if (myfile.is_open()){
    for (size_t i=0;i!=g.size(); ++i){
        myfile << i*delta_r << "    " << g[i] << "\n";
    }
    
    myfile.close();
}
 
  

  printf("Time taken: %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);
}

