using Distributed
using LightGraphs
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


filename = "data/ESLresource/sim_settings.jld";
namespace = smartpath(filename);
@load namespace udscalevec ext_seq t_term repetitions numnodes dedges uedges area r;

lext_seq = length(ext_seq);
ludscalevec = length(udscalevec);

#for graph 1
a=1;
for b=1:ludscalevec
    indices = [a,b];
    filename = "data/ESLresource/sim.jld";
    namespace = smartpath(filename,indices);
    @load namespace ESL_out;
    
    #Network data
    pos = a + (a-1);
    #import edgelists
    ue = uedges[:,pos:(pos+1)];
    de = dedges[:,pos:(pos+1)];
    R"""
    g <- graph_from_edgelist($ue);
    adjm <- as.matrix(get.adjacency(g))
    deg <- degree(g);
    #NOTE: deg from ue and de is the same
    """;
    @rget deg;
    @rget adjm;
    
    rarea = area[:,a];
    confnode = findmax(rarea)[2];
    g = DiGraph(adjm);
    #Calculate this from ue bc it starts from the confnode and explores out
    cdist = gdistances(g, confnode; sort_alg=QuickSort)
    
    
    #Metapopulation data
    numout = Array{Float64}(undef,lext_seq,3,repetitions);
    statesout = Array{Int64}(undef,lext_seq,repetitions,numnodes,t_term);
    let tic=0
        for i=1:lext_seq
            numout[i,:,:] = ESL_out[1][i];
            for r=1:repetitions
                tic=tic+1;
                statesout[i,r,:,:] = ESL_out[2][tic];
            end
        end
    end
    
    # lineplot(vec(sum(statesout[1,3,:,1:1000],dims=1)))
    # lineplot(vec(sum(numout[1,2:3,:],dims=1)))
    
    #8
    #11
    
    # scatterplot(repeat(cdist,inner=repetitions), msd ./ mm)
    
    
    m_e = Array{Float64}(undef,length(ext_seq),repetitions);
    m_s = Array{Float64}(undef,length(ext_seq),repetitions);
    m_l = Array{Float64}(undef,length(ext_seq),repetitions);
    for j=1:length(ext_seq)
        m_e[j,:] = ESL_out[1][j][1,:] ./ numnodes;
        m_s[j,:] = ESL_out[1][j][2,:] ./ numnodes;
        m_l[j,:] = ESL_out[1][j][3,:] ./ numnodes;
    end
    filename="ResultsRes/extplot_noflow.pdf"
    namespace=smartpath(filename);
    R"""
    pdf($namespace,width=5,height=4)
    plot($ext_seq/$r,$(mean(m_s,dims=2)),type='l',ylim=c(0,1),xlab='e/c',ylab='Proportion')
    lines($ext_seq/$r,$(mean(m_l,dims=2)),lty=2)
    dev.off()
    """
    
    # findmax(mean(m_s,dims=2))
    #look at results at the extinction value where S is maximized
    extindex = findmax(vec(mean(m_s,dims=2)))[2];
    
    filename="ResultsRes/nodevtime_noflow.pdf"
    namespace=smartpath(filename);
    R"""
    pdf($namespace,width=12,height=8)
    image(x=seq(1,$t_term),y=seq(1,$numnodes),z=t($(statesout[extindex,1,:,1:1000])),col=c('white','gray','black'),xlab='Time',ylab='Sites')
    # image(x=seq(1,$t_term),y=seq(1,$numnodes),z=t($(ESL_out[2][100])),col=c('white','gray','black'))
    dev.off()
    """
    
    msd = vec(std(statesout[extindex,:,:,1:1000],dims=3));
    mm = vec(mean(statesout[extindex,:,:,1:1000],dims=3));
    
    
    filename="ResultsRes/degreevCV_noflow.pdf"
    namespace=smartpath(filename);
    R"""
    pdf($namespace,width=5,height=4)
    plot(x=jitter($(repeat(deg,inner=repetitions))),y=$(msd ./ mm),pch='.',xlab='Site degree',ylab='Coefficient of Variation',ylim=c(0,max($(msd ./ mm))))
    # image(x=seq(1,$t_term),y=seq(1,$numnodes),z=t($(ESL_out[2][100])),col=c('white','gray','black'))
    dev.off()
    """
    
    filename="ResultsRes/distvCV_noflow.pdf"
    namespace=smartpath(filename);
    R"""
    pdf($namespace,width=5,height=4)
    plot(x=jitter($(repeat(cdist,inner=repetitions))),y=$(msd ./ mm),pch='.',xlab='Distance to confluence',ylab='Coefficient of Variation',ylim=c(0,max($(msd ./ mm))))
    # image(x=seq(1,$t_term),y=seq(1,$numnodes),z=t($(ESL_out[2][100])),col=c('white','gray','black'))
    dev.off()
    """
    filename="ResultsRes/areavCV_noflow.pdf"
    namespace=smartpath(filename);
    R"""
    pdf($namespace,width=5,height=4)
    plot(x=jitter($(repeat(rarea,inner=repetitions))),y=$(msd ./ mm),pch='.',xlab='River area',ylab='Coefficient of Variation',log='x',ylim=c(0,max($(msd ./ mm))))
    # image(x=seq(1,$t_term),y=seq(1,$numnodes),z=t($(ESL_out[2][100])),col=c('white','gray','black'))
    dev.off()
    """
    
    
    
    
    for i=1:lext_seq
        #mean across reps
        #mean across t=500:1000
        meanstatebysite = mean(mean(statesout[i,:,:,:],dims=1)[1,:,500:1000],dims=2);









@load "/Users/jdyeakel/Dropbox/Postdoc/2014_Empirical_Watersheds/Optimal_Channel_Nets/data/ESL_meta2_OCN.jld" numnodes ngraphs uedges dedges area udscalevec t_term repetitions r numout statesout;

# @load "/Users/jdyeakel/Dropbox/Postdoc/2014_Empirical_Watersheds/Optimal_Channel_Nets/data/ESL_meta2_OCN.jld" numnodes ngraphs uedges dedges area udscalevec t_term repetitions r  M_e M_s M_l;

ludscalevec = length(udscalevec);
lext_seq = length(ext_seq);

me_g_ud = Array{Float64}(undef,ngraphs,ludscalevec,lext_seq);
ms_g_ud = Array{Float64}(undef,ngraphs,ludscalevec,lext_seq);
ml_g_ud = Array{Float64}(undef,ngraphs,ludscalevec,lext_seq);
for i=1:ngraphs
    for j=1:ludscalevec
        #mean across repetitions
        me_g_ud[i,j,:] = mean(M_e[i,j,:,:],dims=2);
        ms_g_ud[i,j,:] = mean(M_s[i,j,:,:],dims=2);
        ml_g_ud[i,j,:] = mean(M_l[i,j,:,:],dims=2);
    end
end
#Take mean across graphs
me_ud = Array{Float64}(undef,ludscalevec,lext_seq);
ms_ud = Array{Float64}(undef,ludscalevec,lext_seq);
ml_ud = Array{Float64}(undef,ludscalevec,lext_seq);
for j=1:ludscalevec
    #mean across repetitions
    me_ud[j,:] = mean(me_g_ud[:,j,:],dims=1);
    ms_ud[j,:] = mean(ms_g_ud[:,j,:],dims=1);
    ml_ud[j,:] = mean(ml_g_ud[:,j,:],dims=1);
end

#Examine differences for increasing upstream effects over the course
R"""
pdf("/Users/jdyeakel/Dropbox/Postdoc/2014_Empirical_Watersheds/Optimal_Channel_Nets/Results2/fig_OCN.pdf",width=8,height=8)
layout(matrix(seq(1,16),4,4,byrow=TRUE), widths=rep(1,16), heights=rep(0.5,16))
par(oma = c(4, 4, 1, 1), mar = c(2, 2, 1, 1)) #,mai=c(0.6,0.6,0,0.1)
"""
for i=1:ludscalevec
    R"""
    plot($ext_seq/$r,$(ms_ud[i,:]),type='l',lty=2,xlab='e/c',ylab='Proportion nodes',ylim=c(0,1))
    lines($ext_seq/$r,$(ml_ud[i,:]),lty=1)
    """
end

R"""
title(xlab='e/c',ylab='Proportion filled',outer=TRUE,line=1,cex.lab=2)
dev.off()
"""

#Examine differences for increasing upstream effects over the course
R"""
pdf("/Users/jdyeakel/Dropbox/Postdoc/2014_Empirical_Watersheds/Optimal_Channel_Nets/Results2/fig_OCN2.pdf",width=5,height=4)
library(RColorBrewer)
pal=brewer.pal(3,'Set1')
plot($ext_seq/$r,$(ms_ud[1,:]),type='l',lty=2,xlab='e/c',ylab='Proportion nodes',ylim=c(0,1),lwd=2)
lines($ext_seq/$r,$(ml_ud[1,:]),lty=1,lwd=2)

lines($ext_seq/$r,$(ms_ud[16,:]),type='l',lty=2,col=pal[1],lwd=2)
lines($ext_seq/$r,$(ml_ud[16,:]),lty=1,col=pal[1],lwd=2)
dev.off()
"""

#Examine differences for increasing upstream effects over the course
R"""
pdf("/Users/jdyeakel/Dropbox/Postdoc/2014_Empirical_Watersheds/Optimal_Channel_Nets/Results2/fig_OCN3.pdf",width=5,height=4)
library(RColorBrewer)
pal=brewer.pal(3,'Set1')
plot($ext_seq/$r,$(ms_ud[1,:])+$(ml_ud[1,:]),type='l',lty=1,xlab='e/c',ylab='Proportion nodes',ylim=c(0,1),lwd=2)
# lines($ext_seq/$r,$(ml_ud[1,:]),lty=1,lwd=2)

lines($ext_seq/$r,$(ms_ud[16,:])+$(ml_ud[16,:]),type='l',lty=1,col=pal[1],lwd=2)
# lines($ext_seq/$r,$(ml_ud[16,:]),lty=1,col=pal[1],lwd=2)
dev.off()
"""

