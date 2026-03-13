function [r,p] = pvalfromsim(ur,uct)
N = length(ur);
counts = sum(ur);

xr = sqrt(ur);
xct = sqrt(uct);
mur = mean(xr,1);
muct = mean(xct,1);
sigr = std(xr,0,1);
sigct = std(xct,0,1);
xr = (xr-mur)./sigr;
xct = (xct-muct)./sigct;
r = xr'*xct/(N-1);

numSamples = 10^5;
hits = randi(N,counts,numSamples);
edges = (1:N)-.5;
us = histc(hits,edges,1);
xs = sqrt(us);
mus = mean(xs,1);
sigs = std(xs,0,1);
xs = (xs-mus)./sigs;
rs = xs'*xct/(N-1);


p = sum(rs>r)/length(rs);


end