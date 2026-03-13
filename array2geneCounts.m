function A = array2geneCounts(C,S)
A = C;
A.counts = zeros(size(C.names)); % total counts summed over spatial locations

A.vec = zeros(length(S.x),size(C.names,1),size(C.names,2));
for m = 1:numel(C.names)
    A.vec(:,m) = geneCounts(A.names{m},S);
    A.counts(m) = sum(A.vec(:,m));
end
clear m


end

