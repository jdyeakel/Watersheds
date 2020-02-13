
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

