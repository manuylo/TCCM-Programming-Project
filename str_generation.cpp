#include <iostream>
#include <cmath>
#include <vector>
#include <random>

int main(){
    int N=125;
    auto L=6.14;
    int M = round(std::cbrt(N))+1; // M is first int after cubic root of N
    size_t num_of_atoms = M*M*M, n=3; 
    
    std::vector<std::vector<double> > structure(num_of_atoms, std::vector<double>(n)); //generate a matrix of M**3 atoms with 3 coordinates
  

      for (size_t i = 0; i != num_of_atoms; ++i) { //fill the matrix
        
          structure[i][2]=L/M*(i%M + 1); // z-coordinate is growing fastest, units place in decimal notation
          structure[i][1]=L/M*((i%(M*M))/M + 1); //analog of tens place in decimal notation
          structure[i][0]=L/M*(i/(M*M) + 1); // analog of hundreds place in decimal notation 
    
      }
      
//for i in range MMM-N randomly erase atoms from structure vector
std::random_device r; // Seed with a real random value, if available
std::default_random_engine e1(r());  
while (structure.size()!=N){
    std::uniform_int_distribution<int> uniform_dist(1, structure.size()); // Choose a random mean between 1 and current size
    int pos = uniform_dist(e1);
    structure.erase(structure.begin()+pos-1); //-1 at the end not to get past the vector end position
}
}
