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

    if (r_res[i] >= L/2.0){  //again to cover many cases and avoid while loop :)
            n[i] = -1;
        }
        else if(r_res[i] <= - L/2.0 ){
            n[i] = 1;
        }
        else {
            n[i] = 0;
        }
     r_res[i] = r_res[i] + n[i]*L;

     if (r_res[i] >= L/2.0){  
            while (r_res[i] >= L/2.0)
            {
                r_res[i] -= L;
            }
        }
    
    if (r_res[i] <=- L/2.0){  
            while (r_res[i] <=- L/2.0)
            {
                r_res[i] += L;
            }
        }
   
    }
    return  r_res[0]*r_res[0] + r_res[1]*r_res[1] + r_res[2]*r_res[2];
}

double cos(const std::vector<double>& r_input_1, const std::vector<double>& r_input_2, double L){
    int m = r_input_1.size();
    double result;
    std::vector<double> r_res_1(m),r_res_2(m), n(m);

    r_res_1 = r_input_1;
    r_res_2 = r_input_2;

    for (size_t i=0; i!=m; ++i){
     if (r_res_1[i] >= L/2.0){  
            while (r_res_1[i] >= L/2.0)
            {
                r_res_1[i] -= L;
            }   
        }

    if (r_res_2[i] >= L/2.0){  
            while (r_res_2[i] >= L/2.0)
            {
                r_res_2[i] -= L;
            }   
        }

    if (r_res_1[i] <=- L/2.0){  
            while (r_res_1[i] <=- L/2.0)
            {
                r_res_1[i] += L;
            }
        }

    if (r_res_2[i] <=- L/2.0){  
            while (r_res_2[i] <=- L/2.0)
            {
                r_res_2[i] += L;
            }
        }

    }
    result = r_res_1[0]*r_res_2[0] + r_res_1[1]*r_res_2[1] + r_res_1[2]*r_res_2[2];
    result = result/(sqrt(r_res_1[0]*r_res_1[0] + r_res_1[1]*r_res_1[1] + r_res_1[2]*r_res_1[2]));
    result = result/(sqrt(r_res_2[0]*r_res_2[0] + r_res_2[1]*r_res_2[1] + r_res_2[2]*r_res_2[2]));
    
 return  result;

}



int main(){
    int N,  n=3, c;
    double L=7.6461;
    int accept, update, o, N_accum = 0, sweeps_count = 0, sum = 0,  total_steps = 4500, equilibration= 2500;
    double A=7.049556277, B=0.6022245584, a= 1.80, lambda = 21.0, gamma = 1.20;
    double  temp = 0.1, delta_r_max = 0.15, delta_r = 0.025, r_skin = 0.4;
    double r_n = a + r_skin; 
    double intervals = L/(2*delta_r);
    int n_max = round(2*3.14*pow(r_n,3)*4/3);  
    double x, y, z;
    std::ifstream file ("prod_struc.xyz");
    std::string word, line, atom;
    std::getline(file, line);
    N = std::stoi(line);
    std::cout << N << "\n";
    std::vector<std::vector<double> > structure(N, std::vector<double>(n)), initial_struct(N, std::vector<double>(n));// keep structure here, and initial structure to calc acceptance
    std::getline(file, line);

    for (int i=0;i!=N;++i){
        file >> word >> x >> y >> z ;
        structure[i][0] = x/2.0951; // divide by sigma to turn into reduced units
        structure[i][1] = y/2.0951;
        structure[i][2] = z/2.0951;
    }

    double delta_V, delta_V_3_2=0.0, delta_V_3_1=0.0, delta_V_2=0.0, trial_f2, init_f2, trial_v3_first, init_v3_first, trial_v3_second, init_v3_second, angle, alpha, density = N/(pow(L, 3));
    
    clock_t tStart = clock(); 

    std::vector<int> num_of_neighbours(N); 
    std::vector<std::vector<int> > neighb_indices(N, std::vector<int>(n_max));   // here we collect indices of neighbour atoms for each atom
    std::vector<std::vector<double> > nei_list_struct(N, std::vector<double>(n));  //current structure for nei_list comparison
    std::vector<double> changed_atom(n), r_ij(n),  r_trial(n), r_init(n), histogram(intervals), g(intervals);

    std::random_device r; // Seed with a real random value, if available
    std::default_random_engine e1(r()); 

   initial_struct  = structure;
    nei_list_struct = structure; //1st time assign nei list str;

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

std::uniform_int_distribution<int> uniform_dist(0, N-1); 
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
        structure[k][2] = structure[k][2] + 2*delta_r_max*(wy-0.5);

        //Calculate delta V

        //first calc f_2 (two bodies potential)

        for (size_t i = 0; i != num_of_neighbours[k]; ++i) {
        if (i!=k) {
            auto index = neighb_indices[k][i];
           // calculate vector r_trial[i,k]
            for (size_t j=0; j!=r_trial.size(); ++j){
                r_trial[j] = structure[index][j] - structure[k][j]; 
         } 
            double r_trial_dist = r_mi(r_trial, L);       //calculate  r_trial[i,k] in min im conv
            
            
           trial_f2 = 0;                                   //calculate trial_f2, first part of delta V
            if (r_trial_dist<=pow(a, 2)){  
               trial_f2 = A*(B*pow((1.0/r_trial_dist), 2) - 1)*exp(1.0/(sqrt(r_trial_dist) - a));     
            }
            
         for (size_t j=0; j!=r_init.size(); ++j){           //calculate vector r_init[i,k]
                r_init[j] = structure[index][j] - changed_atom[j];
            } 
        double r_init_dist = r_mi(r_init, L);

            init_f2 = 0;                                    //calculate init_f2, second part of delta V
            if (r_init_dist<pow(a, 2)){
                    
                init_f2 = A*(B*pow((1.0/r_init_dist), 2) - 1)*exp(1.0/(sqrt(r_init_dist) - a));
                
            }
        }
        delta_V_2 = delta_V_2 + trial_f2 - init_f2; 
       
        }
        
        for (size_t j=0; j!=num_of_neighbours[k]; ++j){              //calc of V_3_first
            for (size_t i=0; i!=num_of_neighbours[k]; ++i){
                if(i>j){
                    //trial 
                    auto index_j = neighb_indices[k][j];
                    auto index_i = neighb_indices[k][i];
                    std::vector<double>  r_kj(n);
                    std::vector<double>  r_ki(n);
                    for (size_t z=0; z!=3; ++z){
                        r_kj[z] = structure[k][z] - structure[index_j][z]; 
                        r_ki[z] = structure[k][z] - structure[index_i][z]; 
                    }
                    double r_trial_kj = r_mi(r_kj, L);
                    double r_trial_ki = r_mi(r_ki, L);
                    
                    trial_v3_first = 0;
                    if ((r_trial_kj<pow(a, 2)) and (r_trial_ki<pow(a, 2))){
                            angle = cos(r_ki, r_kj, L);   
                           
                            trial_v3_first = lambda*exp((gamma/((sqrt(r_trial_kj) - a)) + (gamma/( sqrt(r_trial_ki)- a))))*pow((angle +1./3.),2); 
                    }   

                    //init 
                    for (size_t z=0; z!=3; ++z){
                        r_kj[z] = changed_atom[z] - structure[index_j][z] ; 
                        r_ki[z] = changed_atom[z] - structure[index_i][z] ; 
                    }
                    double r_init_kj = r_mi(r_kj, L);
                    double r_init_ki = r_mi(r_ki, L);
                    
                    init_v3_first = 0;
                    if ((r_init_ki<pow(a, 2)) and (r_init_kj<pow(a, 2))){
                            angle = cos(r_ki, r_kj, L);     
                            
                            init_v3_first = lambda*exp((gamma/((sqrt(r_init_kj) - a)) + (gamma/( sqrt(r_init_ki)- a))))*pow((angle +1./3.),2); 
                    }   
            
            delta_V_3_1 = delta_V_3_1  + trial_v3_first - init_v3_first; //(trial_v3_first - init_v3_first)
            
                }
            }
            
        }
        
        
        
        //calculate v3_second
        // k - selected atom; difference with pdf is that k:=i, i:=k, j=j :)
        for (size_t j=0; j!=num_of_neighbours[k]; ++j){ // go over neighbours of sel atom
            for (size_t i=0; i!=num_of_neighbours[neighb_indices[k][j]]; ++i){ // go over neighbours of neighbours of sel atom
                    
                    auto index_j = neighb_indices[k][j];
                    std::vector<double>  r_jk(n);
                    for (size_t z=0; z!=r_jk.size(); ++z){
                        r_jk[z] = structure[index_j][z] - structure[k][z] ; 
                    }
                    double r_trial_jk = r_mi(r_jk, L);
    
                    auto index_i = neighb_indices[index_j][i]; // go over neighbours of neighbours of sel atom
                    std::vector<double>  r_ji(n);
                    for (size_t z=0; z!=r_ji.size(); ++z){
                        r_ji[z] = structure[index_j][z] - structure[index_i][z] ; 
                    }
                    double r_trial_ji = r_mi(r_ji, L);
                    
                    trial_v3_second = 0;
                    if (index_i != k){
                    if ((r_trial_jk<pow(a, 2)) and (r_trial_ji<pow(a, 2))){
                        
                            angle = cos(r_ji, r_jk, L );
                            trial_v3_second = lambda*exp((gamma/( sqrt(r_trial_jk) - a)) + (gamma/(sqrt(r_trial_ji) - a)))*(angle +1./3.)*(angle +1./3.); 
                    
                    }  
                    }
                    
                    for (size_t z=0; z!=r_jk.size(); ++z){
                        r_jk[z] = structure[index_j][z] - changed_atom[z]  ; 
                    }
                    double r_init_jk = sqrt(r_mi(r_jk, L));
                
                    
                    double r_init_ji = sqrt(r_mi(r_ji, L));
                    
                    init_v3_second = 0;
                    if (index_i != k){
                    if ((r_init_jk<a) and (r_init_ji<a)){
                        
                            angle = cos(r_ji, r_jk, L );
                            init_v3_second = lambda*exp((gamma/(r_init_jk - a)) + (gamma/(r_init_ji - a)))*(angle +1./3.)*(angle +1./3.);   
                    } 
                    }
                delta_V_3_2 = delta_V_3_2 + trial_v3_second  - init_v3_second  ;  
                
            }
        }//v3_second
       
       
        delta_V = delta_V_2 + delta_V_3_1 +  delta_V_3_2 ; 
    

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
        delta_V_3_2 = 0;
        delta_V_3_1 = 0;
        delta_V_2 =0;

        if (accept==0){
        structure[k][0] = changed_atom[0];
        structure[k][1] = changed_atom[1];
        structure[k][2] = changed_atom[2];
        }
        else{                       //should we update the neighbor lists?
        
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
            }// should we update or not end 
        } 
     

    } // 1 sweep end

    double count = 0.0;    //to count acceptance rate
        for (size_t i=0; i!=structure.size(); ++i){
           if (initial_struct[i]!=structure[i]){
                count = count +1.0;
           }

     } 
    initial_struct = structure;
    //std::cout << count/N << "\n"; // - uncomment to show acceptance rate of this sweep
    
sweeps_count = sweeps_count +1;
  // std::cout << sweeps_count << "\n";  // - uncomment to show progress

        if (sweeps_count>equilibration){
            for (int i=0; i!=structure.size(); ++i){
                for (int j=0; j!=structure.size(); ++j){
                    if (j>i){
                        for (size_t z=0; z!=r_ij.size(); ++z){
                            r_ij[z] = structure[i][z] - structure[j][z];
                            
                        } 
                    double dist_r_ij = sqrt(r_mi(r_ij, L));
                    if (dist_r_ij<=L/2){
                    o = dist_r_ij/delta_r; 
                    histogram[o] = histogram[o] + 1;
                  
                   }
                    }
                }
            }
        N_accum = N_accum + 1;
        }
} // N steps end
 

for (size_t i=0; i!=histogram.size(); ++i){
        g[i] = histogram[i]/(N_accum*N/2*4/3*density*3.1415*(pow(((i+1)*delta_r), 3) - pow(i*delta_r, 3)));
         sum = sum + histogram[i];
}
std::ofstream myfile ("output_still_nl_eff.txt");
if (myfile.is_open()){
    for (size_t i=0;i!=g.size(); ++i){
        myfile << i*delta_r << "    " << g[i] << "\n";
    }
    
    myfile.close();
}
    std::cout << sum << "\n";
    std::cout << sweeps_count << "\n";
    std::cout << N_accum << "\n";

  printf("Time taken: %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);
} //int main end