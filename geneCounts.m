function v = geneCounts(geneName,S)
% given a gene name and a data structure containg 10x spatial
% transcriptomics data, return a vector of counts for that gene for
% each spatial point labeled by a barcode

gid = find(strcmp(S.features(:,1),geneName));
v = zeros(length(S.x),1);
    
if ~isempty(gid)
    tmp = S.matrix(S.matrix(:,1)==gid,2:3);
    v(tmp(:,1)) = tmp(:,2);
end

end

