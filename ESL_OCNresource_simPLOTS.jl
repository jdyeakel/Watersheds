using LightGraphs
using RCall
using JLD2
using Distributions
include("$(homedir())/Dropbox/Postdoc/2014_Empirical_Watersheds/Optimal_Channel_Nets/src/smartpath.jl")
#Set up R environment
R"""
options(warn = -1)
library(OCNet)
library(igraph)
library(RColorBrewer)
"""


filename = "data/ESLresource/sim_settings.jld";
namespace = smartpath(filename);
@load namespace udscalevec ext_seq t_term repetitions numnodes dedges uedges area r LM;

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
    g2 <- graph_from_edgelist($de);
    adjm <- as.matrix(get.adjacency(g))
    deg <- degree(g);
    #NOTE: deg from ue and de is the same
    """;
    @rget deg;
    @rget adjm;
    
    layoutmatrix = LM[a];
    
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
    
    filename=string("ResultsResArea/extplot_",flowvec[b],".pdf");
    namespace=smartpath(filename);
    R"""
    pdf($namespace,width=5,height=4)
    plot($ext_seq/$r,$(mean(m_s,dims=2)),type='l',ylim=c(0,1),xlab='e/c',ylab='Proportion')
    lines($ext_seq/$r,$(mean(m_l,dims=2)),lty=2)
    dev.off()
    """
    filename=string("ResultsResArea/extplot_",flowvec[b],"res.pdf");
    namespace=smartpath(filename);
    R"""
    pdf($namespace,width=5,height=4)
    plot($ext_seq/$r,$(mean(m_sres,dims=2)),type='l',ylim=c(0,1),xlab='e/c',ylab='Proportion')
    lines($ext_seq/$r,$(mean(m_lres,dims=2)),lty=2)
    dev.off()
    """
    
    extindex = findmax(vec(mean(m_s,dims=2)))[2];
    extindexres = findmax(vec(mean(m_sres,dims=2)))[2];
    
    
    filename=string("ResultsResArea/extplot_",flowvec[b],"both.pdf");
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
    
    
    filename=string("ResultsResArea/nodevtime_",flowvec[b],".pdf");
    namespace=smartpath(filename);
    R"""
    pdf($namespace,width=12,height=8)
    image(x=seq(1,$t_term),y=seq(1,$numnodes),z=t($(statesout[extindex,1,:,1:1000])),col=c('white','gray','black'),xlab='Time',ylab='Sites')
    dev.off()
    """
    filename=string("ResultsResArea/nodevtime_",flowvec[b],"res.pdf");
    namespace=smartpath(filename);
    R"""
    pdf($namespace,width=12,height=8)
    image(x=seq(1,$t_term),y=seq(1,$numnodes),z=t($(statesout_res[extindexres,1,:,1:1000])),col=c('white','gray','black'),xlab='Time',ylab='Sites')
    dev.off()
    """
    

    
    msd[b,:] = vec(std(statesout[extindex,:,:,800:1000],dims=3));
    mm[b,:] = vec(mean(statesout[extindex,:,:,800:1000],dims=3));
    
    msdres[b,:] = vec(std(statesout_res[extindexres,:,:,800:1000],dims=3));
    mmres[b,:] = vec(mean(statesout_res[extindexres,:,:,800:1000],dims=3));
    
    #River plot
    filename=string("ResultsResArea/river_",flowvec[b],".pdf");
    namespace=smartpath(filename);
    R"""
    library(RColorBrewer)
    pal = colorRampPalette(brewer.pal(11,'Greys'))(100)
    colorindex <- floor((2 - $(mm[b,:]))/(2-0)*100)
    colorindex[which(colorindex==0)] = 1;
    pdf($namespace,width=10,height=10)
    plot(g2, layout = $layoutmatrix,vertex.size=log($rarea+1),vertex.color=pal[colorindex],edge.color="darkgrey",edge.arrow.size=0.3,vertex.label=NA)
    dev.off()
    """
    filename=string("ResultsResArea/river_",flowvec[b],"res.pdf");
    namespace=smartpath(filename);
    R"""
    library(RColorBrewer)
    pal = colorRampPalette(brewer.pal(11,'Greys'))(100)
    colorindex <- floor((2 - $(mmres[b,:]))/(2-0)*100)
    colorindex[which(colorindex==0)] = 1;
    pdf($namespace,width=10,height=10)
    plot(g2, layout = $layoutmatrix, vertex.size=log($rarea+1), vertex.color=pal[colorindex], edge.color="darkgrey", edge.arrow.size=0.3, vertex.label=NA)
    dev.off()
    """
    
    filename=string("ResultsResArea/degreevCV_",flowvec[b],".pdf");
    namespace=smartpath(filename);
    R"""
    pdf($namespace,width=5,height=4)
    plot(x=jitter($(repeat(deg,inner=repetitions))),y=$(msd[b,:] ./ mm[b,:]),pch='.',xlab='Site degree',ylab='Coefficient of Variation',ylim=c(0,max($(msd[b,:] ./ mm[b,:]))))
    dev.off()
    """
    filename=string("ResultsResArea/degreevCV_",flowvec[b],"res.pdf");
    namespace=smartpath(filename);
    R"""
    pdf($namespace,width=5,height=4)
    plot(x=jitter($(repeat(deg,inner=repetitions))),y=$(msdres[b,:] ./ mmres[b,:]),pch='.',xlab='Site degree',ylab='Coefficient of Variation',ylim=c(0,max($(msdres[b,:] ./ mmres[b,:]))))
    dev.off()
    """
    
    filename=string("ResultsResArea/distvCV_",flowvec[b],".pdf");
    namespace=smartpath(filename);
    R"""
    pdf($namespace,width=5,height=4)
    plot(x=jitter($(repeat(cdist,inner=repetitions))),y=$(msd[b,:] ./ mm[b,:]),pch='.',xlab='Distance to confluence',ylab='Coefficient of Variation',ylim=c(0,max($(msd[b,:] ./ mm[b,:]))))
    dev.off()
    """
    filename=string("ResultsResArea/distvCV_",flowvec[b],"res.pdf");
    namespace=smartpath(filename);
    R"""
    pdf($namespace,width=5,height=4)
    plot(x=jitter($(repeat(cdist,inner=repetitions))),y=$(msdres[b,:] ./ mmres[b,:]),pch='.',xlab='Distance to confluence',ylab='Coefficient of Variation',ylim=c(0,max($(msdres[b,:] ./ mmres[b,:]))))
    dev.off()
    """
    
    filename=string("ResultsResArea/areavCV_",flowvec[b],".pdf");
    namespace=smartpath(filename);
    R"""
    pdf($namespace,width=5,height=4)
    plot(x=jitter($(repeat(rarea,inner=repetitions))),y=$(msd[b,:] ./ mm[b,:]),pch='.',xlab='River area',ylab='Coefficient of Variation',log='x',ylim=c(0,max($(msd[b,:] ./ mm[b,:]))))
    dev.off()
    """
    filename=string("ResultsResArea/areavCV_",flowvec[b],"res.pdf");
    namespace=smartpath(filename);
    R"""
    pdf($namespace,width=5,height=4)
    plot(x=jitter($(repeat(rarea,inner=repetitions))),y=$(msdres[b,:] ./ mmres[b,:]),pch='.',xlab='River area',ylab='Coefficient of Variation',log='x',ylim=c(0,max($(msdres[b,:] ./ mmres[b,:]))))
    dev.off()
    """
    
    filename=string("ResultsResArea/degreevSS_",flowvec[b],".pdf");
    namespace=smartpath(filename);
    R"""
    pdf($namespace,width=5,height=4)
    plot(x=jitter($(repeat(deg,inner=repetitions))),y=$(mm[b,:]),pch='.',xlab='Site degree',ylab='Steady state',ylim=c(0,max($(mm[b,:]))))
    points
    dev.off()
    """
    filename=string("ResultsResArea/degreevSS_",flowvec[b],"res.pdf");
    namespace=smartpath(filename);
    R"""
    pdf($namespace,width=5,height=4)
    plot(x=jitter($(repeat(deg,inner=repetitions))),y=$(mmres[b,:]),pch='.',xlab='Site degree',ylab='Steady state (res)',ylim=c(0,max($(mmres[b,:]))))
    dev.off()
    """
    
    filename=string("ResultsResArea/SSvSSres",flowvec[b],".pdf");
    namespace=smartpath(filename);
    R"""
    pdf($namespace,width=5,height=4)
    plot(x=$(mm[b,:]),y=$(mmres[b,:]),pch='.',xlab='Mean SS',ylab='Mean SS (res)',xlim=c(0,2),ylim=c(0,2))
    dev.off()
    """
    
    
    
end    

rarea = area[:,1];    
    
filename=string("ResultsResArea/areavCVratio.pdf");
namespace=smartpath(filename);
R"""
pdf($namespace,width=5,height=4)
plot(x=jitter($(repeat(rarea,inner=repetitions))),y=$(msd[2,:] ./ mm[2,:])/$(msd[1,:] ./ mm[1,:]),pch='.',xlab='River area',ylab='CV ratio: HiFlow/NoFlow',log='x',ylim=c(0,max($(msd[2,:] ./ mm[2,:])/$(msd[1,:] ./ mm[1,:]))))
dev.off()
"""    

filename=string("ResultsResArea/areavCVratio_res.pdf");
namespace=smartpath(filename);
R"""
pdf($namespace,width=5,height=4)
plot(x=jitter($(repeat(rarea,inner=repetitions))),y=$(msdres[2,:] ./ mmres[2,:])/$(msdres[1,:] ./ mmres[1,:]),pch='.',xlab='River area',ylab='CV ratio: HiFlow/NoFlow',log='x',ylim=c(0,max($(msdres[2,:] ./ mmres[2,:])/$(msdres[1,:] ./ mmres[1,:]))))
dev.off()
"""    
    
    
    
    
    
    
    
    
# 
# mean(vec((msd[2,:] ./ mm[2,:]))[vec(findall(x->x==900,vec(repeat(rarea,inner=repetitions))))]) / mean(vec((msd[1,:] ./ mm[1,:]))[vec(findall(x->x==900,vec(repeat(rarea,inner=repetitions))))])
# 
# 
    
    
    
    
    
    
    
    
    
    
    
