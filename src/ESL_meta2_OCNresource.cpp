#include <Rcpp.h>
using namespace Rcpp;

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar)

// Importing distinct up and downstream colonization and rescue rates

// [[Rcpp::export]]
List ESL_meta2_OCNresource(int n, int t_term, int repetitions, List nn_up, List nn_down, double r, double cup, double cdown, double mup, double mdown, NumericVector ext_seq, NumericVector rarea) {
  
  //Cypher
  
  //E = 0
  //S = 1
  //L = 2
  int extl = ext_seq.size(); 
  List Meta_ext(extl);
  
  List Meta_nodeS_ext(extl);
  List Meta_nodeL_ext(extl);
  
  int erl = extl*repetitions;
  List Xstates(erl);
  
  double t_mid_d = t_term/2;
  int t_mid = (int)t_mid_d;
  int interval = t_term - t_mid;
  
  int exreptic = -1;
  
  for (int extic=0;extic<extl;extic++) {
    
    //Define extinction values
    double el = ext_seq(extic);
    double es = 2*el;
    
    NumericMatrix Meta_mean(3,repetitions);
    
    // List Meta_nodeS(repetitions,);
    // List Meta_nodeL(repetitions);
    
    // NumericMatrix Meta_burnE(interval,repetitions);
    // NumericMatrix Meta_burnS(interval,repetitions);
    // NumericMatrix Meta_burnL(interval,repetitions);
    
    //Begin repitition iterations
    for (int rep=0;rep<repetitions;rep++) {
      
      exreptic = exreptic + 1;
      
      //Probability matrices
      NumericMatrix pr_Colonize(n,t_term);
      NumericMatrix pr_Grow(n,t_term);
      NumericMatrix pr_NoGrowRescue(n,t_term);
      NumericMatrix pr_NoGrowNoRescueStay(n,t_term);
      NumericMatrix pr_NoGrowNoRescueExtinct(n,t_term);
      
      //Initlialize Res, X, and Meta
      IntegerVector Res(n);
      IntegerMatrix X(n,t_term);
      IntegerMatrix Meta(3,t_term);
      for (int i=0;i<n;i++) {
        NumericVector rdraw_temp = runif(1);
        double rdraw = as<double>(rdraw_temp);
        // double rdraw; rdraw = rand();
        if (rdraw < 1/3) {
          X(i,0) = 0;
          Meta(0,0) = Meta(0,0) + 1;
          //Probability that resource is + 1
          // NumericVector resdraw_temp = runif(1);
          // double resdraw = as<double>(resdraw_temp);
          NumericVector resdraw_temp = runif(1);
          double resdraw = as<double>(resdraw_temp);
          // double resdraw; resdraw = rand();
          //Population not present, so Resource may be +1 or -1
          if (resdraw < 0.5) {
            Res(i) = 1;
          } else {
            Res(i) = -1;
          }
        }
        if (rdraw > 1/3) {
          if (rdraw < 2/3) {
            X(i,0) = 1;
            Meta(1,0) = Meta(1,0) + 1;
            //Resource =1 because a population is present
            Res(i) = 1;
          }
          if (rdraw > 2/3) {
            X(i,0) = 2;
            Meta(2,0) = Meta(2,0) + 1;
            //Resource =1 because a population is present
            Res(i) = 1;
          }
        }
      }
      
      //Begin time iterations (t)
      for (int t=0;t<(t_term-1);t++) {
        
        //Update resource spins by the collective influence of their neighbors
        //Currently does not depend on the state of its own site
        //Modify this by int influence = Res(i);
        for (int i=0;i<n;i++) {
          IntegerVector nn_patchesup = as<IntegerVector>(nn_up[i]);
          IntegerVector nn_patchesdown = as<IntegerVector>(nn_down[i]);
          int l_patchesup = nn_patchesup.size();
          int l_patchesdown = nn_patchesdown.size();
          
          //Loop through nearest UPSTREAM neighbor patches
          int influence = 0;
          for (int j=0;j<l_patchesup;j++) {
            int j_patch = nn_patchesup(j);
            influence = influence + (Res(j_patch)*rarea(j_patch)); //*rarea(j_patch)
          } 
          for (int j=0;j<l_patchesdown;j++) {
            int j_patch = nn_patchesdown(j);
            influence = influence + (Res(j_patch)*rarea(j_patch)); //*rarea(j_patch)
          } 
          if (influence > 0) {
            Res(i) = 1;
          } else {
            Res(i) = -1;
          }
        }
        
        
        
        //ESL UPDATING
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
          int num_Lup = 0;
          int num_Ldown = 0;
          
          //Find the states of nearest neighbors for 'i'
          IntegerVector nn_patchesup = as<IntegerVector>(nn_up[i]);
          IntegerVector nn_patchesdown = as<IntegerVector>(nn_down[i]);
          int l_patchesup = nn_patchesup.size();
          int l_patchesdown = nn_patchesdown.size();
          IntegerVector nn_statesup(l_patchesup);
          IntegerVector nn_statesdown(l_patchesdown);
          
          //Loop through nearest UPSTREAM neighbor patches
          for (int j=0;j<l_patchesup;j++) {
            int j_patch = nn_patchesup(j);
            nn_statesup(j) = states(j_patch);
            //How many Large nodes are nearest neighbors to node i?
            if (nn_statesup(j) == 2) {num_Lup = num_Lup + 1;}
          } 
          for (int j=0;j<l_patchesdown;j++) {
            int j_patch = nn_patchesdown(j);
            nn_statesdown(j) = states(j_patch);
            //How many Large nodes are nearest neighbors to node i?
            if (nn_statesdown(j) == 2) {num_Ldown = num_Ldown + 1;}
          } 
          
          //Dynamics of the current patch
          
          //If the current state is Empty
          if (states(i) == 0) {
            // What is the probability that empty nodes is colonized by a neighboring large node?
            double pr_ES = 1 - (pow((1-cup),num_Lup)*pow((1-cdown),num_Ldown));
            NumericVector draw_ES_temp = runif(1);
            double draw_ES = as<double>(draw_ES_temp);
            // double draw_ES; draw_ES = rand();
            
            if (draw_ES < pr_ES) {
              //Species move to node to colonize, but are there resources?
              if (Res(i) > 0) {
                //If it is colonized, 0 -> 1
                X(i,t+1) = 1;
                //Don't need to update Res because it is already +1
              } else {
                X(i,t+1) = 0;
                //Don't need to update Res because it is already -1
              }
            } else {
              //If it is not colonized, 0 -> 0
              X(i,t+1) = 0;
              //Res stays the same
            }
            //save the probability or colonization
            pr_Colonize(i,t) = pr_ES;
            
          }
          
          //If the current state is Small
          if (states(i) == 1) {
            //Probability that S -> L from growth
            double pr_SL_grow = r;
            //Probability that S -> L from rescue
            double pr_SL_rescue = 1 - (pow((1-mup),num_Lup)*pow((1-mdown),num_Ldown));
            //Probability that S -> S... no growth, no rescue, no extinction
            double pr_S_stay = (1 - es);
            NumericVector draw_SL_grow_temp = runif(1);
            double draw_SL_grow = as<double>(draw_SL_grow_temp);
            // double draw_SL_grow; draw_SL_grow = rand();
            
            //Does S grow to L?
            if (draw_SL_grow < pr_SL_grow) {
              X(i,t+1) = 2;
              //Resources are maintained by their population
              Res(i) = 1;
            } else {
              //If S does not grow, it can 1) get rescued, stay small, or go extinct
              NumericVector draw_SL_rescue_temp = runif(1);
              double draw_SL_rescue = as<double>(draw_SL_rescue_temp);
              // double draw_SL_rescue; draw_SL_rescue = rand();
              
              //Does S get rescued to L?
              if (draw_SL_rescue < pr_SL_rescue) {
                X(i,t+1) = 2;
                //Resources are maintained by their population
                Res(i) = 1;
              } else {
                //If not it can stay or go extinct
                NumericVector draw_S_stay_temp = runif(1);
                double draw_S_stay = as<double>(draw_S_stay_temp);
                // double draw_S_stay; draw_S_stay = rand();
                //Does S -> S?
                if (draw_S_stay < pr_S_stay) {
                  X(i,t+1) = 1;
                  //Resources are maintained by their population
                  Res(i) = 1;
                  //If not, S -> 0 and extinction occurs
                } else {
                  X(i,t+1) = 0;
                  //Resources flip from +1 to -1
                  Res(i) = -1;
                }
              }
            }
            
            pr_Grow(i,t) = pr_SL_grow;
            pr_NoGrowRescue(i,t) = (1-pr_SL_grow)*pr_SL_rescue;
            pr_NoGrowNoRescueStay(i,t) = (1-pr_SL_grow)*(1-pr_SL_rescue)*(1-es);
            pr_NoGrowNoRescueExtinct(i,t) = (1-pr_SL_grow)*(1-pr_SL_rescue)*es;
            
          }
          
          //If the current state is Large
          if (states(i) == 2) {
            //L->S 
            double pr_LS_shrink = el;
            NumericVector draw_LS_shrink_temp = runif(1);
            double draw_LS_shrink = as<double>(draw_LS_shrink_temp);
            // double draw_LS_shrink; draw_LS_shrink = rand();
            
            if (draw_LS_shrink < pr_LS_shrink) {
              X(i,t+1) = 1;
              //Resources are maintained by their population
              Res(i) = 1;
            } else {
              X(i,t+1) = 2;
              //Resources are maintained by their population
              Res(i) = 1;
            }
            
          }
          
          //Record how many of ea
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
        
        //COUNTS of E S L
        Meta(0,t+1) = XE_next;
        Meta(1,t+1) = XS_next;
        Meta(2,t+1) = XL_next;
        
      }// end t
      
      Xstates(exreptic) = X;
      
      //Burn off first half of simulation
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
      
      // Meta_nodeS(rep) = Meta_burnS;
      // Meta_nodeL(rep) = Meta_burnL;
      
    }// end reps
    
    Meta_ext(extic) = Meta_mean;
    // Meta_nodeS_ext(extic) = Meta_burnS;
    // Meta_nodeL_ext(extic) = Meta_burnL;
    
  } // end extinction sequence
  
  List esl_out(2);
  
  esl_out(0) = Meta_ext; // [extl][3,rep]
  // esl_out(1) = Meta_nodeS_ext;
  // esl_out(2) = Meta_nodeL_ext;
  esl_out(1) = Xstates; // [extl*repititions][numnodes,t_term]
  //esl_out(1) = X;
  //esl_out(2) = pr_Colonize;
  //esl_out(3) = pr_Grow;
  //esl_out(4) = pr_NoGrowRescue;
  //esl_out(5) = pr_NoGrowNoRescueStay;
  //esl_out(6) = pr_NoGrowNoRescueExtinct;
  
  return esl_out;
  
}

