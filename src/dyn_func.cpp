#include <Rcpp.h>
using namespace Rcpp;

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar)

// For more on using Rcpp click the Help button on the editor toolbar

// [[Rcpp::export]]
List dyn_func(int t_term, int reps, List nn_list, List downstream_list, NumericVector ext_vector, double pr_col, double pr_fme, double pr_tc) {
  //Setup
  int l = ext_vector.size(); //Number of extinction ratio iterations
  int n = nn_list.size(); //Number of nodes
  double t_mid_d = t_term/2;
  int t_mid = (int)t_mid_d;
  
  
  
  NumericVector ratio(l);
  for (int i=0; i<l; i++) { //Define extinction/colonization ratio
    ratio[i] = ext_vector[i]/pr_col;
  }

  //Begin replication iterations
  NumericVector mean_occupancy(l); //Mean occupancy across reps for each extinction iteration
  NumericVector sd_occupancy(l); //St Dev occupancy accross reps for each extinction iteration
  NumericMatrix mean_per_node_occupancy(n,l);
  NumericMatrix sd_per_node_occupancy(n,l);
  
  for (int p= 0; p < l; p++) { //Begin extinction vector iterations
    
    double pr_ext = ext_vector(p);
    Rcout << "Pr(Extinction) =  " << pr_ext << std::endl;

    NumericVector rep_occupancy(reps);
    NumericMatrix per_node_occupancy(n,reps);
    
    for (int r = 0; r < reps; r++) { //Begin replication iterations
      
      NumericVector node_occupancy(t_term);
      IntegerMatrix node_cond(n,t_term); //Build State Matrix
      
      for (int i = 0; i < n; i++) { //Set the initial state
        node_cond(i,0) = 1;
      }
      //Begin Time Series
      for (int t = 0; t < t_term-1; t++) {
        
        //Build the extinction vector... which nodes go extinct?
        IntegerVector ext(n);
        NumericVector ext_rand = runif(n);
        for (int i=0; i<n; i++) {
          if (ext_rand(i) < pr_ext) {ext(i) = 1;} else{ext(i) = 0;}
        }
        
        //Build the colonization information for each node's nn's, separately
        //i.e. is node 'x' colonized? Doesn't matter who from...
        //Only those nodes alive can colonize!
        IntegerVector col(n);
        for (int x=0; x<n; x++) {
          //If population is present at node x, then colonize to neighboring nodes
          if (node_cond(x,t) == 1) {
            NumericVector node_nn = as<NumericVector>(nn_list[x]);
            int num_node_nn = node_nn.size();
            NumericVector col_rand = runif(num_node_nn);
            for (int i=0; i<num_node_nn; i++) { //Across each nearest neighbor
              int node_id = node_nn(i);
              if (col_rand(i) < pr_col) {col(node_id) = 1;}
              //else {col(node_id-1) = 0;} // -1 to account for R:1; C++:0
            }
          }
        }
        
        //Determine flow-mediated extinctions
        IntegerVector fme(n);
        //Will a flow-mediated extinction occur? Draw probability
        NumericVector fme_rand_temp = runif(1);
        double fme_rand = as<double>(fme_rand_temp);
        if (fme_rand < pr_fme) {
          //If it does occur choose RANDOM starting node
          int max_num = n-1; // -1 because 0 is counted
          int min_num = 0;
          int fme_start = rand() % (max_num + 1 - min_num) + min_num;
          //Find list of downstream nodes
          NumericVector downstream_nodes = as<NumericVector>(downstream_list[fme_start]);
          int ds_size = downstream_nodes.size();
          int tic;
          for (int x=0;x<n;x++) {
            tic = 0;
            for (int y=0;y<ds_size;y++) {
              //If node x is on the list of downstream nodes, the pr_fme for node x is 1
              if (x == downstream_nodes(y)) {
                tic = tic + 1;
              }
            }
            if (tic > 0) {fme(x) = 1;}
          }
        } //else, pr_fme is a vector of zeros
        
        //Determine probabilistic transect colonizations... for each node
        IntegerVector tc(n);
        //For each node in the system...
        for (int x=0;x<n;x++) {
          //Draw probability to determine if the node x (if present) colonizes another random node
          if (node_cond(x,t) == 1) {
            NumericVector tc_rand_temp = runif(1);
            double tc_rand = as<double>(tc_rand_temp);
            //If transecting colonization is successful
            if (tc_rand < pr_tc) {
              //Randomly choose which node is colonized
              int max_num = n-1; // -1 because 0 is counted
              int min_num = 0;
              int tc_node = rand() % (max_num + 1 - min_num) + min_num;
              //Now note that the node tc_node will be colonized
              tc(tc_node) = 1;
            }
          }
        }
        
        //if (t==10) {Rcout << "The value is " << col(10) << std::endl;}
        
        //First, set the next time step equal to current time step
        for (int x=0; x<n; x++) {
          node_cond(x,t+1) = node_cond(x,t);
        }
        //Modifications of next time step based on extinction and colonization:::
        //First apply colonizations, and then apply extinctions...
        //This means that colonizations happen at the end of the interval...
        for (int x=0; x<n; x++) {
          if (col(x) == 1) {node_cond(x,t+1) = 1;} //Normal colonization
          if (tc(x) == 1) {node_cond(x,t+1) = 1;} //Transecting colonization
          if (ext(x) == 1) {node_cond(x,t+1) = 0;} //Normal extinctions
          if (fme(x) == 1) {node_cond(x,t+1) = 0;} //Flow-mediated extinctions
        }
        
        //Calculate the proportion of nodes at the current time t that are occupied
        int count=0;
        double tot_nodes = n;
        for (int x=0; x<n; x++) {
          if (node_cond(x,t) == 1) {count = count + 1;}
        }
        node_occupancy(t) = count/tot_nodes;
      } // end time steps
      
      
      //Average ocupancy from time steps: middle:end
      int burn = t_term-t_mid;
      NumericVector occupancy(burn);
      for (int i=0; i<burn; i++) {
        occupancy(i) = node_occupancy(i+burn);
      }
      rep_occupancy(r) = mean(occupancy);
      
      //Average occupancy PER NODE from time steps: middle:end
      double burn_double = burn;
      for (int x=0;x<n;x++) {
        double tic = 0;
        for (int t=t_mid;t<t_term;t++) {
          tic = tic + node_cond(x,t); //tic advances if node is occupied
        }
        per_node_occupancy(x,r) = tic/burn_double; //Number of occupations / Number of possible occupations
      }
      
    } // End reps
  
  //Across Node Occupancy calculation
  mean_occupancy(p) = mean(rep_occupancy);
  sd_occupancy(p) = sd(rep_occupancy);
  
  //Per Node Occupancy calculation
  for (int x=0;x<n;x++) {
    NumericVector per_node_occupancy_vector(reps);
    for (int r=0;r<reps;r++) {
      per_node_occupancy_vector(r) = per_node_occupancy(x,r);
    }
    //Calculate means and sds across replicates
    mean_per_node_occupancy(x,p) = mean(per_node_occupancy_vector);
    sd_per_node_occupancy(x,p) = sd(per_node_occupancy_vector);
  }
  
  
  } // End ext.vector (p) iterations
  
  List out(4);
  out(0) = mean_occupancy;
  out(1) = sd_occupancy;
  out(2) = mean_per_node_occupancy;
  out(3) = sd_per_node_occupancy;
  //out(2) = tree_state;
  return out;
} //End function
