#include <Rcpp.h>
using namespace Rcpp;

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar)

// For more on using Rcpp click the Help button on the editor toolbar

// [[Rcpp::export]]
List DLA_func_radial(int size, int boundary) {
  //Set Seed
  IntegerMatrix f_dla_tree(size,2);
  f_dla_tree(0,0) = 1; f_dla_tree(0,1) = 1;
  f_dla_tree(1,0) = 1; f_dla_tree(1,1) = 2;
  IntegerMatrix link_tree(size-1,2);
  link_tree(0,0) = 1; link_tree(0,1) = 2;
  int num_node = 2;
  
  //Find which f_dla_tree spaces are not equal to zero
//  int num_node = 0;
//  for (int i=0;i<size; i++) {
//    if (f_dla_tree(i,0) > 0) {num_node = num_node + 1;}
//  }
  
  while (num_node < size) {
    
    //Subset the matrix for smaller workload
    //This may not be necessary in C++?
    IntegerMatrix dla_tree(num_node,2);
    for (int i=0;i<num_node;i++) {
      for (int j=0;j<2;j++) {
        dla_tree(i,j) = f_dla_tree(i,j);
      }
    }
    
    //Define maximum boundary
    //by concatonating both vectors of dla_tree
    IntegerVector dla_vector(2*num_node);
    int tic = 0;
    for (int i=0;i<num_node;i++) {
      for (int j=0;j<2;j++) {
        dla_vector(tic) = dla_tree(i,j);
        tic = tic+1;
      }
    }
    //Define maximum boundary with respect to the current dla tree
    //int max_bound = max(dla_vector) + boundary;
    
    //Take 2
    //Define maximum boundary of current dla tree using euclidean distance
    //Calculate distance btw the origin and each node in the dla tree
    NumericVector euclid_dist(num_node);
    for (int i=0;i<num_node;i++) {
      euclid_dist(i) = sqrt(pow(0 - dla_tree(i,0),2.0) + pow(0 - dla_tree(i,1),2.0));
    }
    //Find maximum euclid_dist
    double max_dist = max(euclid_dist);
    //max_bound is the maximum euclidean distance 
    double max_bound_t = max_dist + boundary;
    //turn it into integer
    int max_bound = (int)max_bound_t;
    
//    //Find Random starting point that is not part of DLA tree but within boundaries
//    int seed_in_tree = 1; // This is just to start the while loop
//    IntegerVector set_seed(2);
//    while(seed_in_tree >= 1) {
//      seed_in_tree = 0; // Set to zero... if seed is in tree, this will increase
//      //Draw 2 random numbers between 1 and boundary
//      for (int i=0;i<2;i++){
//        set_seed(i) = rand() % (max_bound + 1 - 1) + 1;
//      }
//      //Determine whether this vector is in dla_tree
//      for (int i=0;i<num_node;i++) {
//        if ((dla_tree(i,0) == set_seed(0)) && (dla_tree(i,1) == set_seed(1))) {
//          seed_in_tree = seed_in_tree + 1;
//        }
//      }
//    } // end seed while loop
    
    //Take 2... initiate seed at random point on the boundary
    //1) Define boundary
    int boundary_tics = 1000;
    double pi = 3.141592654;
    double t = 0.0;
    IntegerMatrix boundary_points(boundary_tics,2);
    for (int i=0;i<max_bound;i++) {
      double t = t + (2*pi)/boundary_tics;
      double x_coord = max_bound * cos(t);
      double y_coord = max_bound * sin(t);
      boundary_points(i,0) = x_coord;
      boundary_points(i,1) = y_coord;
    }
    //2) Draw random number and sample from boundary
    int set_seed = rand() % (boundary_tics + 1 - 1) + 1;
    //3) Define seed
    IntegerVector seed(2);
    seed(0) = boundary_points(set_seed,0); seed(1) = boundary_points(set_seed,1);
    
    //Rcout << "seed0 =  " << seed(0) << std::endl;
    //Rcout << "seed1 =  " << seed(1) << std::endl;
    
    //Begin Random Walk until it hits dla.tree or leaves the Zone
    int check = 1;
    int rr_tic = 0;
    while (check == 1) {
      rr_tic = rr_tic + 1;
      
      IntegerVector seed_new(2);
      for (int i=0;i<2;i++){
        seed_new(i) = seed(i) + rand() % (1 + 1 - (-1)) + (-1);
      }
      
      //Calculate distance btw the new seed and the origin
      double euclid_dist_seed = sqrt(pow(seed_new(0) - 0,2.0) + pow(seed_new(1) - 0,2.0));
      
      //Calculte euclidean distance between seed and of the tree (for later)
      NumericVector euclid_dist(num_node);
      for (int i=0;i<num_node;i++) {
        euclid_dist(i) = sqrt(pow(seed_new(0) - dla_tree(i,0),2.0) + pow(seed_new(1) - dla_tree(i,1),2.0));
      }
      
      //Find minimum euclid_dist
      //double min_dist = euclid_dist(which_min(euclid_dist));
      
      // Is the random walk within the radius of the seeding circle? (with a little lee-way)
      if (euclid_dist_seed < (max_bound + 10)) {
        
        // Is the NEW seed part of the dla tree?
        //Determine whether this vector is in dla_tree
        //4-14-2014: I changed i=1 in the for loop;; bc if seed hits confluence, I don't want it to stick
        int seed_in_tree = 0;
        for (int i=1;i<num_node;i++) { 
          if ((dla_tree(i,0) == seed_new(0)) && (dla_tree(i,1) == seed_new(1))) {
            seed_in_tree = seed_in_tree + 1;
          }
        }
        
        if (seed_in_tree == 0) {
          //Find any nodes with  0 < dist < sqrt(2)
          IntegerVector hit_test(num_node);
          hit_test(0) = 0; //just ensures that nothing sticks to first node (confluence)
          for (int i=1;i<num_node;i++) {
            if ((euclid_dist(i) <= sqrt(2)) && (euclid_dist(i) > 0)) {
              hit_test(i) = 1;
              }
          }
          int hit;
          int hit_node;
          if (sum(hit_test) > 0) {
            hit_node = which_max(hit_test);
            hit = 1;
          } else {hit = 0;}
          
          //If you hit the tree
          if (hit == 1) {
            //Update the full dla tree (f_dla_tree)
            num_node = num_node + 1;
            f_dla_tree(num_node-1,0) = seed_new(0); //-1 bc C++ starts at 0
            f_dla_tree(num_node-1,1) = seed_new(1); //-1 bc C++ starts at 0
            //Note the link
            link_tree((num_node-2),0) = hit_node+1; //-2 bc C++ starts at 0 // +1 for same reason
            link_tree((num_node-2),1) = num_node; //-2 bc C++ starts at 0
            check = 0;
            Rcout << "Number of nodes =  " << num_node << std::endl;
            //PRINT?
          } else {
            seed = seed_new;
          }
        } else {check = 0;}
      } else {check = 0;} // end random walk loop
    
    } // end check loop
  
    // Update num_node!
    //Find which f_dla_tree spaces are not equal to zero
//    int num_node = 0;
//    for (int i=0;i<size; i++) {
//      if (f_dla_tree(i,0) > 0) {num_node = num_node + 1;}
//    }
  
  } // end tree while loop
  
  // Return the edge matrix
  List out(2);
  out(0) = link_tree;
  out(1) = f_dla_tree;
  return out;
}
