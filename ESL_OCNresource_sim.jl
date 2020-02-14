using Distributed
@everywhere using RCall
@everywhere using SharedArrays
@everywhere using JLD2
@everywhere using Distributions
@everywhere include("$(homedir())/Dropbox/Postdoc/2014_Empirical_Watersheds/Optimal_Channel_Nets/src/smartpath.jl")
#Set up R environment
R"""
options(warn = -1)
library(OCNet)
library(igraph)
library(RColorBrewer)
"""

numnodes = 30^2;
numedges = numnodes -1; #only true for trees!
ngraphs = 1;
# cd = 0.2;
udscalevec = [1,4];
ludscalevec = length(udscalevec);
# maxlinks = 30*30;

#Growthrate
r = 0.2

ext_seq = collect(0:0.02:0.4);
extl = length(ext_seq);
t_term = 1000;
repetitions = 10;

its = ngraphs*ludscalevec;

#downstream edgelist
dedges = SharedArray{Int64}(numedges,2*ngraphs);
#upstream edgelist
uedges = SharedArray{Int64}(numedges,2*ngraphs);
#River area
area = SharedArray{Float64}(numnodes,ngraphs);
#make the graphs beforehand
# R"graphs <- list(); dedges <- list(); uedges <- list()";
@time @sync @distributed for i=1:ngraphs
    R"""
    library(OCNet)
    library(igraph)
    OCN <- aggregate_OCN(landscape_OCN(create_OCN(sqrt($numnodes), sqrt($numnodes))), thrA = 4)
    a <- OCN$FD$A
    # draw_simple_OCN(OCN)
    g <- OCN_to_igraph(OCN, level = "FD")
    de <- get.edgelist(g);
    ue <- de[,c(2,1)];
    # dedges[[$i]] <- de;
    # uedges[[$i]] <- ue;
    # graphs[[$i]] <- g;
    # numv = V(g);
    # nume = E(g);
    """;
    @rget de;
    @rget ue;
    @rget a;
    # @rget numv;
    # @rget nume;
    pos = i+(i-1)
    dedges[1:numedges,pos:(pos+1)] = de;
    uedges[1:numedges,pos:(pos+1)] = ue;
    area[:,i] = a;
end

# M_e = SharedArray{Float64}(ngraphs,ludscalevec,length(ext_seq),repetitions);
# M_s = SharedArray{Float64}(ngraphs,ludscalevec,length(ext_seq),repetitions);
# M_l = SharedArray{Float64}(ngraphs,ludscalevec,length(ext_seq),repetitions);

numout = SharedArray{Float64}(ngraphs,ludscalevec,extl,3,repetitions);
statesout = SharedArray{Int64}(ngraphs,ludscalevec,extl,repetitions,numnodes,t_term);

filename = "data/ESLresource/sim_settings.jld";
namespace = smartpath(filename);
@save namespace udscalevec ext_seq t_term repetitions numnodes dedges uedges area r;


@time @sync @distributed for i=0:(its-1)
    #Across graphs
    a = mod(i,ngraphs) + 1;
    #Across udscale
    b = Int64(floor(i/ngraphs)) + 1;
    
    pos = a + (a-1);
    
    #import edgelists
    ue = uedges[:,pos:(pos+1)];
    de = dedges[:,pos:(pos+1)];
    
    udscale = udscalevec[b];
    
    #define colonization upstream and downstream
    #upstream colonization/rescue is greater than downstream colonization/rescue
    #NOTE: 2/12/20 should we scale it so sum(cdown,cup) and sum(mdown,mup) is constant regardless of the udscale that is applied? Because right now, with increased udscale there is a increased overall c and m... which confounds the signal apart from up/down structure.
    cdown = copy(r);
    cup = cdown*udscale;
    mdown = copy(r);
    mup = mdown*udscale;
    
    # initial = rand([0,1,2],numnodes);
    #get nearest neighbors for up/downstream
    #Make list of nearest neighbors for each node
    R"""
    library(igraph)
    library(Rcpp)
    sourceCpp('/Users/jdyeakel/Dropbox/Postdoc/2014_Empirical_Watersheds/Optimal_Channel_Nets/src/ESL_meta2_OCN.cpp')
    sourceCpp('/Users/jdyeakel/Dropbox/Postdoc/2014_Empirical_Watersheds/Optimal_Channel_Nets/src/ESL_meta2_OCNresource.cpp')

    adj.m <- get.adjacency(graph_from_edgelist($ue))
    nn_up <- list()
    for (j in 1:$numnodes) {
        arow <- adj.m[j,]
        nn_up[[j]] <- which(arow==1)
    }
    nn_upc <- lapply(nn_up,function(x){x-1})
    
    adj.m <- get.adjacency(graph_from_edgelist($de))
    nn_down <- list()
    for (j in 1:$numnodes) {
        arow <- adj.m[j,]
        nn_down[[j]] <- which(arow==1)
    }
    nn_downc <- lapply(nn_down,function(x){x-1})
    
    
    #Initiate the simulation
    # X <- matrix(0,$numnodes,$t_term)
    # Meta <- matrix(0,3,$t_term)
    
    #Initiate starting vector
    # initial <- c($numnodes - 2*round($numnodes/3,0),round($numnodes/3,0),round($numnodes/3,0))
    # X[,1] <- c(rep(0,initial[1]),rep(1,initial[2]),rep(2,initial[3]))
    # X[,1] <- $initial
    # Meta[,1] <- c(length(which($initial==0)),length(which($initial==1)),length(which($initial==2)))

    #Perform Simulation
    #without resource lattice
    ESL_out <- ESL_meta2_OCN($numnodes, $t_term, $repetitions, nn_upc, nn_downc, $r, $cup, $cdown, $mup, $mdown, $ext_seq)
    #with resource lattice
    ESL_outres <- ESL_meta2_OCNresource($numnodes, $t_term, $repetitions, nn_upc, nn_downc, $r, $cup, $cdown, $mup, $mdown, $ext_seq)
    
    # m_e <- matrix(0,length($ext_seq),$repetitions)
    # m_s <- matrix(0,length($ext_seq),$repetitions)
    # m_l <- matrix(0,length($ext_seq),$repetitions)
    # for (j in 1:length($ext_seq)) {
    #     m_e[j,] <- ESL_out[[1]][[j]][1,]/$numnodes
    #     m_s[j,] <- ESL_out[[1]][[j]][2,]/$numnodes
    #     m_l[j,] <- ESL_out[[1]][[j]][3,]/$numnodes
    # }  
    """;
    # @rget m_e;
    # @rget m_s;
    # @rget m_l;
    # 
    @rget ESL_out;
    @rget ESL_outres;
    
    #save data
    indices = [a,b];
    filename = "data/ESLresource/sim.jld";
    namespace = smartpath(filename,indices);
    @save namespace ESL_out;
    indices = [a,b];
    filename = "data/ESLresource/sim_res.jld";
    namespace = smartpath(filename,indices);
    @save namespace ESL_out;
    
    # let tic=0
    #     for i=1:extl
    #         numout[a,b,i,:,:] = ESL_out[1][i];
    #         for r=1:repetitions
    #             tic=tic+1;
    #             statesout[a,b,i,r,:,:] = ESL_out[2][tic];
    #             println(i)
    #         end
    #     end
    # end
    # 
    # #save data
    # M_e[a,b,:,:] = m_e;
    # M_s[a,b,:,:] = m_s;
    # M_l[a,b,:,:] = m_l;
end
# 
# @save "/Users/jdyeakel/Dropbox/Postdoc/2014_Empirical_Watersheds/Optimal_Channel_Nets/data/ESL_meta2_OCN.jld" numnodes ngraphs uedges dedges area udscalevec t_term repetitions r numout statesout;

# @save "/Users/jdyeakel/Dropbox/Postdoc/2014_Empirical_Watersheds/Optimal_Channel_Nets/data/ESL_meta2_OCN.jld" numnodes ngraphs uedges dedges area udscalevec t_term repetitions r  M_e M_s M_l;
