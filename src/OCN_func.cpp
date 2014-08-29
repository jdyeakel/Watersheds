#include <Rcpp.h>
using namespace Rcpp;

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar)

// For more on using Rcpp click the Help button on the editor toolbar

// [[Rcpp::export]]
List OCN_func(int L, int iterations) {
   //Define coordinates of Nearest Neighbors
   //Should be pretty easy to change to diagonals now
   int nn_num = 4;
   int tic = 0;
   List nn(L*L);
   IntegerMatrix neighbors(nn_num,2);
   neighbors(0,0) = 0;neighbors(0,1) = 1;
   neighbors(1,0) = 1;neighbors(1,1) = 0;
   neighbors(2,0) = 0;neighbors(2,1) = -1;
   neighbors(3,0) = -1;neighbors(3,1) = 0;
   int toc = 0;
   IntegerMatrix coords(L*L,2);
   IntegerMatrix coords_matrix(L,L);
   for(int i=0;i<L;i++) {
     for(int j=0;j<L;j++) {
       coords(toc,0) = i;
       coords(toc,1) = j;
       coords_matrix(i,j) = toc;
       IntegerMatrix nn_ij(nn_num,2);
       IntegerVector keep(nn_num); for (int k=0;k<nn_num;k++) {keep(k)= -10;}
       int tic = 0;
       for (int k=0;k<nn_num;k++) {
         int i_new = i + neighbors(k,0);
         int j_new = j + neighbors(k,1);
         if ((i_new >= 0) && (i_new < L) && (j_new >= 0) && (j_new < L)) {
           keep(tic) = k;
           tic = tic + 1;
         }
         nn_ij(k,0) = i_new;
         nn_ij(k,1) = j_new;
       }
       int cnt = 0;
       for (int k=0;k<nn_num;k++) {
         if (keep(k)>-10) {cnt = cnt + 1;}
       }
       int nn_size = cnt;
       IntegerMatrix nn_ij_final(nn_size,2);
       for (int k=0;k<cnt;k++) {
         nn_ij_final(k,0) = nn_ij(keep(k),0);
         nn_ij_final(k,1) = nn_ij(keep(k),1);
       }
       nn(toc) = nn_ij_final;
       toc = toc + 1;
     }
   }
   
   
   
   
   //Build early river from random walk
   int tree_size = 1;
   int max_size = L*L-1;
   //Link Matrix
   IntegerMatrix tree_link(max_size,2);
   //List of nodes in tree
   IntegerVector nodes_in_tree(L*L);
   //Starting node values... places confluence in the lower left of graph
   int node = 0;
   int step = 0;
   while (step < max_size) {
     
     //Find nearest neighbors of node
     IntegerMatrix nn_node = nn(node);
     int nn_num = nn_node.nrow();
     IntegerVector nn_id(nn_num);
     //Obtain nearest neighbor identities
     for (int i=0;i<nn_num;i++) {
       nn_id(i) = coords_matrix(nn_node(i,0),nn_node(i,1));
     }
     
     //Select nearest neighbor not part of nodes_in_tree
     //Which nodes are part of the tree?
     IntegerVector open_nn(nn_num);
     //Look across each nearest neighbor
     for (int i=0;i<nn_num;i++) {
       int tic = 0;
       //Does the nearest neighbor match a node in the tree?
       for (int j=0;j<tree_size;j++) {
         if (nn_id(i) == nodes_in_tree(j)) {
           tic = tic + 1; //If it does, advance tic
         }
       }
     if (tic == 0) {open_nn(i) = 1;} //labels the nodes that are NOT in the tree
     }
     
     int open_nn_num = sum(open_nn); //Number of open nearest neighbors (not in tree)
     
     if (open_nn_num >= 1) { //IF THERE ARE ANY OPEN NEAREST NEIGHBORS
        //Create list of open nearest neighbors NOT in the tree
        IntegerVector nn_open(open_nn_num);
        int tic = 0;
        for (int i=0;i<nn_num;i++) {
          if (open_nn(i) == 1) {
            nn_open(tic) = nn_id(i); 
            tic = tic + 1;
          }
        }
     
        //Draw random integer between 0 and length of nn_open
        int max_num = open_nn_num-1; // -1 because 0 is counted
        int min_num = 0;
        int r_num = rand() % (max_num + 1 - min_num) + min_num;
        //Define the linking nearest neighbor
        //This will be the new node for the next iteration
        int link_node = nn_open(r_num);
        
        //Update bookkeeping
        //Establish link between node and link_node & update attributes
        tree_link(step,0) = node;
        tree_link(step,1) = link_node;
        tree_size = tree_size + 1;
        
        //Set new node to link_node
        node = link_node;
        nodes_in_tree(step) = node;
        //Advance counter
        //step = step + 1;
      } 
      else { //IF THERE ARE NOT ANY OPEN NEAREST NEIGHBORS
        
        //Find a node on the spanning tree with open nearest neighbors to re-start RR
        
        //Find nodes on tree with open nearest neighbors (this is all we need to do here)
        //Sample among these nodes for new starting location
        //List open_nn_tree(tree_size);
        IntegerVector open_tree_size(tree_size);
        for (int k=0;k<tree_size;k++) {
          //Look across each node in tree
          int node_tree = nodes_in_tree(k);
          //Find nearest neighbors of node
          IntegerMatrix nn_node_tree = nn(node_tree);
          int nn_num = nn_node_tree.nrow();
          IntegerVector nn_tree_id(nn_num);
          //Obtain nearest neighbor identities
          int tic = 0;
          for (int j=0;j<nn_num;j++) {
            nn_tree_id(j) = coords_matrix(nn_node_tree(j,0),nn_node_tree(j,1));
            //Is this nearest neighbor on the spanning tree?
            for (int i=0;i<tree_size;i++) {
              if(nn_tree_id(j) != nodes_in_tree(i)) {tic = tic + 1;}
            }
          }
          //If tic is greater than zero, there are nn at this tree node that are also not part of tree!
          if (tic > 0) {open_tree_size(k) = 1;}
        }  
        
        //Compile the nodes where open_tree_size == 1 from which to sample next random node  
        int num_open = sum(open_tree_size);
        IntegerVector open_nodes(num_open);
        int tic = 0;
        for (int i=0;i<num_open;i++) {
          if (open_tree_size(i) == 1) {
            open_nodes(tic) = nodes_in_tree(i);
          }
          tic = tic + 1;
        }
        
        //Choose a new starting location ON THE TREE
        int max_num = num_open-1; // -1 because 0 is counted
        int min_num = 1; //Do not choose the confluence node
        int r_num = rand() % (max_num + 1 - min_num) + min_num;
        node = open_nodes(r_num);
        
        //DO NOT increase tree size
        //Do NOT form a new link
        //This will be done on the next step.
        step = step - 1; //this ensure the counter (next step) does not advance
                        //Because no NEW node has been chosen... only an old node.
      } //End else

      //Advance counter
      step = step + 1;
    
   } //end while loop for the early random walk river net
   
   //Assign River Area
   
   //Get the degree for each node
   NumericVector link_vector((L*L-1)*2);
   int cnt = 0;
   for (int i=0;i<(L*L-1);i++) {
     for (int j=0;j<2;j++) {
       link_vector(cnt) = tree_link(i,j);
       cnt = cnt + 1;
     }
   }
   NumericVector degree(L*L);
   int tribs = 0;
   for (int i=0;i<L*L;i++) {
     NumericVector node_present((L*L-1)*2);
     for (int j=0;j<((L*L-1)*2);j++) {
       if (i == link_vector(j)) {node_present(j)=1;}
     }
     degree(i) = sum(node_present);
     //Find number of tributaries that do not count the confluence
     if (degree(i) == 1) {tribs = tribs + 1;}
   }
   tribs = tribs - 1; //Eliminate the confluence node
   
   //Make a list of the tributaries
   IntegerVector trib_node(tribs);
   cnt = 0;
   for (int i=0;i<L*L;i++) {
     if ((degree(i)==1) && (i != 0)) {
       trib_node(cnt) = i;
       cnt = cnt + 1;
     }
   }
   
   //Find graph distance between each tributary and the confluence.
   for (int i=0;i<tribs;i++) {
     //Define starting point
     int trib_i = trib_node(i);
     //Define ending point
     int conf = 0;
   }
   
   
   
   
   
   
   
   List data(5);
   data(0) = coords;
   data(1) = coords_matrix;
   data(2) = tree_link;
   data(3) = degree;
   data(4) = trib_node;
   return data;
}
