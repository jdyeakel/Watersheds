#include <Rcpp.h>
using namespace Rcpp;

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar)

// For more on using Rcpp click the Help button on the editor toolbar

// [[Rcpp::export]]
List ESL_meta2(int n, int t_term, int repetitions, IntegerMatrix X, IntegerMatrix Meta, List nn, double c, double r, double m, NumericVector ext_seq) {
  
  //Cypher
  
  //E = 0
  //S = 1
  //L = 2
  int extl = ext_seq.size(); 
  List Meta_ext(extl);
  
  for (int extic=0;extic<extl;extic++) {
    
    //Define extinction values
    double el = ext_seq(extic);
    double es = 2*el;
    
    NumericMatrix Meta_mean(3,repetitions);
    
    //Begin repitition iterations
    for (int rep=0;rep<repetitions;rep++) {
      
      //Probability matrices
      NumericMatrix pr_Colonize(n,t_term);
      NumericMatrix pr_Grow(n,t_term);
      NumericMatrix pr_NoGrowRescue(n,t_term);
      NumericMatrix pr_NoGrowNoRescueStay(n,t_term);
      NumericMatrix pr_NoGrowNoRescueExtinct(n,t_term);
      
      //Begin time iterations (t)
      for (int t=0;t<(t_term-1);t++) {
        
        //Record the states of patches at time t
        IntegerVector states(n);
        for (int j=0;j<n;j++) {
          states(j) = X(j,t);
        }
        
        //Set the 'next' counters at 0
        //These counters will count over iterations i
        int XE_next = 0;
        int XS_next = 0;
        int XL_next = 0;
        
        
        //Begin individual patch iterations (i)
        for (int i=0;i<n;i++) {
          
          //Large nodes counter :: how many large nodes are connected to node 'i'?
          int num_L = 0;
          
          //Find the states of nearest neighbors for 'i'
          IntegerVector nn_patches = as<IntegerVector>(nn[i]);
          int l_patches = nn_patches.size();
          IntegerVector nn_states(l_patches);
          
          //Loop through nearest neighbor patches
          for (int j=0;j<l_patches;j++) {
            int j_patch = nn_patches(j);
            nn_states(j) = states(j_patch);
            
            //How many Large nodes are nearest neighbors to node i?
            if (nn_states(j) == 2) {num_L = num_L + 1;}
          } 
          
          //Dynamics of the current patch
          
          //If the current state is Empty
          if (states(i) == 0) {
            
            double pr_ES = 1 - pow((1-c),num_L);
            NumericVector draw_ES_temp = runif(1);
            double draw_ES = as<double>(draw_ES_temp);
            
            if (draw_ES < pr_ES) {
              X(i,t+1) = 1;
            } else {
              X(i,t+1) = 0;
            }
            
            pr_Colonize(i,t) = pr_ES;
            
          }
          
          //If the current state is Small
          if (states(i) == 1) {
            
            double pr_SL_grow = r;
            double pr_SL_rescue = 1 - pow((1-m),num_L);
            double pr_S_stay = (1 - es);
            
            NumericVector draw_SL_grow_temp = runif(1);
            double draw_SL_grow = as<double>(draw_SL_grow_temp);
            
            if (draw_SL_grow < pr_SL_grow) {
              X(i,t+1) = 2;
            } else {
              NumericVector draw_SL_rescue_temp = runif(1);
              double draw_SL_rescue = as<double>(draw_SL_rescue_temp);
              
              if (draw_SL_rescue < pr_SL_rescue) {
                X(i,t+1) = 2;
              } else {
                NumericVector draw_S_stay_temp = runif(1);
                double draw_S_stay = as<double>(draw_S_stay_temp);
                
                if (draw_S_stay < pr_S_stay) {
                  X(i,t+1) = 1;
                } else {X(i,t+1) = 0;}
              }
            }
            
            pr_Grow(i,t) = pr_SL_grow;
            pr_NoGrowRescue(i,t) = (1-pr_SL_grow)*pr_SL_rescue;
            pr_NoGrowNoRescueStay(i,t) = (1-pr_SL_grow)*(1-pr_SL_rescue)*(1-es);
            pr_NoGrowNoRescueExtinct(i,t) = (1-pr_SL_grow)*(1-pr_SL_rescue)*es;
            
          }
          
          //If the current state is Large
          if (states(i) == 2) {
            
            double pr_LS_shrink = el;
            NumericVector draw_LS_shrink_temp = runif(1);
            double draw_LS_shrink = as<double>(draw_LS_shrink_temp);
            
            if (draw_LS_shrink < pr_LS_shrink) {
              X(i,t+1) = 1;
            } else {
              X(i,t+1) = 2;
            }
            
          }
          
          if (X(i,t+1) == 0) {
            XE_next = XE_next + 1;
          }
          if (X(i,t+1) == 1) {
            XS_next = XS_next + 1;
          }
          if (X(i,t+1) == 2) {
            XL_next = XL_next + 1;
          }
          
        } // end i
        
        Meta(0,t+1) = XE_next;
        Meta(1,t+1) = XS_next;
        Meta(2,t+1) = XL_next;
        
      }// end t
      
      
      //Burn off first half of simulation
      double t_mid_d = t_term/2;
      int t_mid = (int)t_mid_d;
      int interval = t_term - t_mid;
      NumericVector Meta_burnE(interval);
      NumericVector Meta_burnS(interval);
      NumericVector Meta_burnL(interval);
      for (int i=0;i<interval;i++) {
        Meta_burnE(i) = Meta(0,t_mid+i);
        Meta_burnS(i) = Meta(1,t_mid+i);
        Meta_burnL(i) = Meta(2,t_mid+i);
      }
      //Record Summary Statistics
      Meta_mean(0,rep) = mean(Meta_burnE);
      Meta_mean(1,rep) = mean(Meta_burnS);
      Meta_mean(2,rep) = mean(Meta_burnL);
      
    }// end reps
    
    Meta_ext(extic) = Meta_mean;
    
  } // end extinction sequence
  
  List esl_out(7);
  
  esl_out(0) = Meta_ext;
  //esl_out(1) = X;
  //esl_out(2) = pr_Colonize;
  //esl_out(3) = pr_Grow;
  //esl_out(4) = pr_NoGrowRescue;
  //esl_out(5) = pr_NoGrowNoRescueStay;
  //esl_out(6) = pr_NoGrowNoRescueExtinct;
  
  return esl_out;
  
}

