function [S] = readSampleData(sampleDir,type)
% 5/10/2022
% Michael Moore
% Read Spatial Transscriptomics data for a single sample

% 10.31.2022  updated to read the SCT normalized only

switch type
    case 'filtered'
        fdir = 'filtered_feature_bc_matrix';
    case 'raw'
        fdir = 'raw_feature_bc_matrix';
    case 'sct'
        fdir = 'filtered_feature_bc_matrix';
end

%% basic data structure

% S.features 
%   3 columns: Ensemble Gene ID, Gene Name, Property
%   column 3 is just 'Gene Expression', not important
%   number of rows = number of genes (e.g. 25499)
%   the row index is used as a unique gene identifier (gid)

% S.barcodes
%   1 column: a dna sequence called a barcode
%   each row is a spatial location (e.g. 1379)
%   row index can be thought of as a unique location identifer (lid)

% S.matrix
%   a sparse matrix encoding of the levels
%       col 1 refers to a gene by gid
%       col 2 refers to a location by lid
%       col 3 gives the level of the gene at the location

% S.coords
%   a 5 column matrix which collates barcode to x-y coordinates
%   to be used for plots

% Note there are more barcodes in S.coords than in S.barcodes


%%

featuresTSV = fullfile(sampleDir,fdir,'features.tsv');
barcodesTSV = fullfile(sampleDir,fdir,'barcodes.tsv');
matrixMTX = fullfile(sampleDir,fdir,'matrix.mtx');
tissueCSV = fullfile(sampleDir,fdir,'spatial','tissue_positions_list.csv');
SCTnormalizedCSV = fullfile(sampleDir,fdir,'SCTnormalized.csv');

% S = struct;

% opts = detectImportOptions(featuresTSV, 'NumHeaderLines', 0,'ReadVariableNames',false,'FileType','text');
% opts.DataLines = [1 Inf];  
% opts.Delimiter = {'\t'};
% opts.VariableNames = {'ID','Name','Property'};
% S.features = readtable(featuresTSV,opts);
% S.features = table2cell(S.features(:,1:2));
% S.featuresColNames = opts.VariableNames(1:2);
% 
% opts = detectImportOptions(matrixMTX, 'NumHeaderLines', 0,'ReadVariableNames',false,'FileType','text');
% opts.DataLines = [4 Inf];  
% opts.VariableNames = {'Gene','Location','Level'};
% S.matrix = readtable(matrixMTX,opts);
% S.matrix = table2array(S.matrix);
% S.matrixColNames = opts.VariableNames;

opts = detectImportOptions(barcodesTSV, 'NumHeaderLines', 0,'ReadVariableNames',false,'FileType','text');
opts.DataLines = [1 Inf];  
opts.VariableNames = {'Barcode'};
barcodes = readtable(barcodesTSV,opts);
barcodes = table2cell(barcodes);

opts = detectImportOptions(tissueCSV, 'NumHeaderLines', 0,'ReadVariableNames',false);
opts.DataLines = [1 Inf];  
opts.VariableNames = {'barcode','inTissue','row','col','rowFullRes','colFullRes'};
coords=readtable(tissueCSV,opts);
coords = table2cell(coords);
coordsColNames = opts.VariableNames;

opts = detectImportOptions(SCTnormalizedCSV, 'NumHeaderLines',0,'ReadVariableNames',false,'FileType','text');
opts.DataLines = [2 Inf];  
opts.VariableNames = {'row','gene','barcode','expression'};
sct = readtable(SCTnormalizedCSV,opts);
sct = table2cell(sct);
sctColNames = opts.VariableNames;

%% The sct matrix is sorted by barcode

%%
S = struct;
S.features(:,1) = unique(sct(:,2));
S.matrix = zeros(size(sct,1),3);
S.matrixColNames = {'Gene','Location','Level'};
% put the normalized values into a matrix col 3
S.matrix(:,3) = [sct{:,4}];

%% compute the location index from barcode
% -------------------------------------------------------------------------
% should take 20 seconds to run this:
for m = 1:length(barcodes)
    test = strcmp(sct(:,3),barcodes{m});
    S.matrix(test,2) = m;
end
clear m
%% Similarly convert the gene name into a number
% this should take several minutes to run because there are 12674 gene
% names
%      might consider finding a faster algorithm
%      faster approach might involve sorting
for m = 1:length(S.features)
    test = strcmp(sct(:,2),S.features{m,1});
    S.matrix(test,1) = m;
end
clear m



%% we want a fullres x-y pair for each row in barcodes

[C,ia,ib] = intersect(barcodes,coords(:,1),'stable');
S.x0 = cell2mat(coords(:,5));
S.y0 = cell2mat(coords(:,6));
S.x = cell2mat(coords(ib,5));
S.y = cell2mat(coords(ib,6));

end

