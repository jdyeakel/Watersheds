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

flowvec = ["noflow","flow"];

msd = Array{Float64}(undef,ludscalevec,numnodes*repetitions);
mm = Array{Float64}(undef,ludscalevec,numnodes*repetitions);
msdres = Array{Float64}(undef,ludscalevec,numnodes*repetitions);
mmres = Array{Float64}(undef,ludscalevec,numnodes*repetitions);

#for graph 1
a=1;
for b=1:ludscalevec
    indices = [a,b];
    filename = "data/ESLresource/sim.jld";
    namespace = smartpath(filename,indices);
    @load namespace ESL_out;
    
    filename = "data/ESLresource/sim_res.jld";
    namespace = smartpath(filename,indices);
    @load namespace ESL_outres;
    
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
    numout_res = Array{Float64}(undef,lext_seq,3,repetitions);
    statesout_res = Array{Int64}(undef,lext_seq,repetitions,numnodes,t_term);
    let tic=0
        for i=1:lext_seq
            numout[i,:,:] = ESL_out[1][i];
            numout_res[i,:,:] = ESL_outres[1][i];
            for r=1:repetitions
                tic=tic+1;
                statesout[i,r,:,:] = ESL_out[2][tic];
                statesout_res[i,r,:,:] = ESL_outres[2][tic];
            end
        end
    end
    
    # lineplot(vec(sum(statesout[1,3,:,1:1000],dims=1)))
    # lineplot(vec(sum(numout[1,2:3,:],dims=1)))
    
    #8
    #11
    
    # scatterplot(repeat(cdist,inner=repetitions), msd[b,:] ./ mm[b,:])
    
    
    m_e = Array{Float64}(undef,length(ext_seq),repetitions);
    m_s = Array{Float64}(undef,length(ext_seq),repetitions);
    m_l = Array{Float64}(undef,length(ext_seq),repetitions);
    m_eres = Array{Float64}(undef,length(ext_seq),repetitions);
    m_sres = Array{Float64}(undef,length(ext_seq),repetitions);
    m_lres = Array{Float64}(undef,length(ext_seq),repetitions);
    for j=1:length(ext_seq)
        m_e[j,:] = ESL_out[1][j][1,:] ./ numnodes;
        m_s[j,:] = ESL_out[1][j][2,:] ./ numnodes;
        m_l[j,:] = ESL_out[1][j][3,:] ./ numnodes;
        m_eres[j,:] = ESL_outres[1][j][1,:] ./ numnodes;
        m_sres[j,:] = ESL_outres[1][j][2,:] ./ numnodes;
        m_lres[j,:] = ESL_outres[1][j][3,:] ./ numnodes;
    end
    extindex = findmax(vec(mean(m_s,dims=2)))[2];
    extindexres = findmax(vec(mean(m_sres,dims=2)))[2];
    
    filename=string("ResultsRes/extplot_",flowvec[b],".pdf");
    namespace=smartpath(filename);
    R"""
    pdf($namespace,width=5,height=4)
    plot($ext_seq/$r,$(mean(m_s,dims=2)),type='l',ylim=c(0,1),xlab='e/c',ylab='Proportion')
    lines($ext_seq/$r,$(mean(m_l,dims=2)),lty=2)
    dev.off()
    """
    filename=string("ResultsRes/extplot_",flowvec[b],"res.pdf");
    namespace=smartpath(filename);
    R"""
    pdf($namespace,width=5,height=4)
    plot($ext_seq/$r,$(mean(m_sres,dims=2)),type='l',ylim=c(0,1),xlab='e/c',ylab='Proportion')
    lines($ext_seq/$r,$(mean(m_lres,dims=2)),lty=2)
    dev.off()
    """
    filename=string("ResultsRes/extplot_",flowvec[b],"both.pdf");
    namespace=smartpath(filename);
    R"""
    library(RColorBrewer)
    pal=brewer.pal(3,'Set1')
    pdf($namespace,width=5,height=4)
    plot($ext_seq/$r,$(mean(m_s,dims=2)),type='l',ylim=c(0,1),xlab='e/c',ylab='Proportion')
    lines($ext_seq/$r,$(mean(m_l,dims=2)),lty=2)
    lines($ext_seq/$r,$(mean(m_sres,dims=2)),lty=1,col=pal[2])
    lines($ext_seq/$r,$(mean(m_lres,dims=2)),lty=2,col=pal[2])
    dev.off()
    """
    
    # findmax(mean(m_s,dims=2))
    #look at results at the extinction value where S is maximized
    
    
    filename=string("ResultsRes/nodevtime_",flowvec[b],".pdf");
    namespace=smartpath(filename);
    R"""
    pdf($namespace,width=12,height=8)
    image(x=seq(1,$t_term),y=seq(1,$numnodes),z=t($(statesout[extindex,1,:,1:1000])),col=c('white','gray','black'),xlab='Time',ylab='Sites')
    dev.off()
    """
    filename=string("ResultsRes/nodevtime_",flowvec[b],"res.pdf");
    namespace=smartpath(filename);
    R"""
    pdf($namespace,width=12,height=8)
    image(x=seq(1,$t_term),y=seq(1,$numnodes),z=t($(statesout[extindexres,1,:,1:1000])),col=c('white','gray','black'),xlab='Time',ylab='Sites')
    dev.off()
    """
    
    msd[b,:] = vec(std(statesout[extindex,:,:,1:1000],dims=3));
    mm[b,:] = vec(mean(statesout[extindex,:,:,1:1000],dims=3));
    
    msdres[b,:] = vec(std(statesout_res[extindexres,:,:,1:1000],dims=3));
    mmres[b,:] = vec(mean(statesout_res[extindexres,:,:,1:1000],dims=3));
    
    
    filename=string("ResultsRes/degreevCV_",flowvec[b],".pdf");
    namespace=smartpath(filename);
    R"""
    pdf($namespace,width=5,height=4)
    plot(x=jitter($(repeat(deg,inner=repetitions))),y=$(msd[b,:] ./ mm[b,:]),pch='.',xlab='Site degree',ylab='Coefficient of Variation',ylim=c(0,max($(msd[b,:] ./ mm[b,:]))))
    dev.off()
    """
    filename=string("ResultsRes/degreevCV_",flowvec[b],"res.pdf");
    namespace=smartpath(filename);
    R"""
    pdf($namespace,width=5,height=4)
    plot(x=jitter($(repeat(deg,inner=repetitions))),y=$(msdres[b,:] ./ mmres[b,:]),pch='.',xlab='Site degree',ylab='Coefficient of Variation',ylim=c(0,max($(msdres[b,:] ./ mmres[b,:]))))
    dev.off()
    """
    
    filename=string("ResultsRes/distvCV_",flowvec[b],".pdf");
    namespace=smartpath(filename);
    R"""
    pdf($namespace,width=5,height=4)
    plot(x=jitter($(repeat(cdist,inner=repetitions))),y=$(msd[b,:] ./ mm[b,:]),pch='.',xlab='Distance to confluence',ylab='Coefficient of Variation',ylim=c(0,max($(msd[b,:] ./ mm[b,:]))))
    dev.off()
    """
    filename=string("ResultsRes/distvCV_",flowvec[b],"res.pdf");
    namespace=smartpath(filename);
    R"""
    pdf($namespace,width=5,height=4)
    plot(x=jitter($(repeat(cdist,inner=repetitions))),y=$(msdres[b,:] ./ mmres[b,:]),pch='.',xlab='Distance to confluence',ylab='Coefficient of Variation',ylim=c(0,max($(msdres[b,:] ./ mmres[b,:]))))
    dev.off()
    """
    
    filename=string("ResultsRes/areavCV_",flowvec[b],".pdf");
    namespace=smartpath(filename);
    R"""
    pdf($namespace,width=5,height=4)
    plot(x=jitter($(repeat(rarea,inner=repetitions))),y=$(msd[b,:] ./ mm[b,:]),pch='.',xlab='River area',ylab='Coefficient of Variation',log='x',ylim=c(0,max($(msd[b,:] ./ mm[b,:]))))
    dev.off()
    """
    filename=string("ResultsRes/areavCV_",flowvec[b],"res.pdf");
    namespace=smartpath(filename);
    R"""
    pdf($namespace,width=5,height=4)
    plot(x=jitter($(repeat(rarea,inner=repetitions))),y=$(msdres[b,:] ./ mmres[b,:]),pch='.',xlab='River area',ylab='Coefficient of Variation',log='x',ylim=c(0,max($(msdres[b,:] ./ mmres[b,:]))))
    dev.off()
    """
end    
    
    
filename=string("ResultsRes/areavCVratio.pdf");
namespace=smartpath(filename);
R"""
pdf($namespace,width=5,height=4)
plot(x=jitter($(repeat(rarea,inner=repetitions))),y=$(msd[2,:] ./ mm[2,:])/$(msd[1,:] ./ mm[1,:]),pch='.',xlab='River area',ylab='CV ratio: HiFlow/NoFlow',log='x',ylim=c(0,max($(msd[2,:] ./ mm[2,:])/$(msd[1,:] ./ mm[1,:]))))
dev.off()
"""    

filename=string("ResultsRes/areavCVratio_res.pdf");
namespace=smartpath(filename);
R"""
pdf($namespace,width=5,height=4)
plot(x=jitter($(repeat(rarea,inner=repetitions))),y=$(msdres[2,:] ./ mmres[2,:])/$(msdres[1,:] ./ mmres[1,:]),pch='.',xlab='River area',ylab='CV ratio: HiFlow/NoFlow',log='x',ylim=c(0,max($(msdres[2,:] ./ mmres[2,:])/$(msdres[1,:] ./ mmres[1,:]))))
dev.off()
"""    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
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

