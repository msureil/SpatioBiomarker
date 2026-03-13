function stvplot(v,S,sz)

% 5/11/2022
% Michael G. Moore, Michigan State University

% spatial transcriptomics plot given a vector of counts, v

colormap jet
scatter(S.x,S.y,sz,'k')
hold on
scatter(S.x(v>=0),S.y(v>=0),sz,v(v>=0),'filled')
hold off
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'ytick',[])
set(gca,'yticklabel',[])
% axis image
axis square
box on
% caxis([1 5])
colorbar

end

