function [r,p,R] = getcostheta(u1,u2,numSamples)
% inputs:
%   u1    a length-N vector of non-negative integers which sum to C
%   u2    a non-negative vector of the same length as u1
% outputs;
%   r     the cosine of the angle between u1 and u2
%   p     the probability of obtaining a cosine equal or larger than r by
%         randomly putting C balls into N buckets

% compute the cos of the angle between the two vectors using
%  v1'*v2 = |v1||v2|cos(theta)

% assert(sum(floor(u1)==u1)==length(u1),'first vector must be integer valued')
assert(sum(u1 >=0)==length(u1),'first vector must be non-negative')
assert(sum(u2 >=0)==length(u2),'second vector must be non-negative')
assert(length(u1)==length(u2),'vectors must have equal lengths')

u2 = u2/vecnorm(u2);

r = u1'*u2/vecnorm(u1);

N = length(u1);
C = round(sum(u1));   % edit: added round 
hits = randi(N,C,numSamples);
edges = (1:N)-.5;
U = histc(hits,edges,1);
R = (U'*u2)./(vecnorm(U)');


p = sum(R >= r)/numSamples;

end

