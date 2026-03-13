function C = readCelltypeGeneCSV(gfapFlag)
if gfapFlag
    fname = '/home/anirban/ubssd/codeforfigure_gfap_bhavnaupdated/CellTypeGenes_Gfap.csv';
else
    fname = '/home/anirban/ubssd/codeforfigure_gfap_bhavnaupdated/CellTypeGenes.csv';
end
opts = detectImportOptions(fname);
opts.DataLines = [1 Inf];  
opts.VariableNamesLine = 1;
C = struct;
C.names =readtable(fname,opts);
C.names = table2cell(C.names);
end

