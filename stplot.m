function stplot(geneName,S,sz)

% 5/11/2022
% Michael G. Moore, Michigan State University

% spatial transcriptomics plot
v = geneCounts(geneName,S);
colormap jet
scatter(S.x,S.y,sz,'k')
hold on
scatter(S.x(v~=0),S.y(v~=0),sz,v(v~=0),'filled')
hold off
axis image
colorbar
title(geneName)

end

