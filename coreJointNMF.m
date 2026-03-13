% Joint NMF analysis across grouped datasets
% Pools all spatial spots within each time point group, runs a single
% group-level NMF factorization, and computes Pearson correlation and
% cosine similarity between NMF factors and reactivity gene distributions.

% --- Configuration ---
sampleDir = '/home/anirban/ubssd/BhanvaPaper_NMF';

% Sample groups
datasets_1week = {'1weekB', '1weekC', '1weekD', '1weekE', '1weekF', '1weekG', '1weekH'};
datasets_6week = {'6weekB', '6weekC', '6weekF', '6weekG', '6weekH', '6weekI', '6weekJ'};

% Reactivity genes per group
%genes_1week = {'Cxcl3', 'Cxcl2', 'Ccl2', 'Il1b', 'Slpi', 'S100a9', ...
%               'Hmox1', 'Lgals3', 'Gpnmb', 'Spp1', 'Ighm', 'Ctsk'};
genes_1week = {'S100a4', 'S100a9', 'Cxcl3', 'Cxcl2', 'Ccl2', 'Il1b',...
                   'Hmox1', 'Lgals3', 'Gpnmb', 'Spp1','Ighm', 'Ctsk','Slpi'};
genes_6week = {'Spp1', 'A2m', 'Lgals3', 'Gpnmb', 'S100a4', 'Serping1', ...
               'Vim', 'Lyz2', 'S100a6', 'S100a9', 'Slpi', 'Ccl6', ...
               'Hmox1', 'Il1b', 'Cxcl3', 'Il1rn', 'Srgn', 'Ccl2', 'Cxcl2'};

% NMF rank per cell type: [Astrocytes, Neurons, OPCs, Microglia]
% Change values here — everything downstream adapts automatically
Klist = [1, 2, 1, 1];

% Cell type names (must match order of Klist)
cellTypeNames = {'Astrocytes', 'Neurons', 'OPCs', 'Microglia'};

% Groups to process
groups     = {'1week',        '6week'       };
allDatasets = {datasets_1week, datasets_6week};
allGenes    = {genes_1week,    genes_6week   };
% --- End Configuration ---

% Make profile names automatically from Klist and cellTypeNames
% e.g. Klist=[1,2,1,1] -> {'Astrocyte','Neuron 1','Neuron 2','OPC','Microglia'}
profileNames = {};
for ct = 1:length(Klist)
    baseName = cellTypeNames{ct};
    % Remove trailing 's' for singular form (e.g. 'Astrocytes' -> 'Astrocyte')
    if baseName(end) == 's'
        baseName = baseName(1:end-1);
    end
    if Klist(ct) == 1
        profileNames{end+1} = baseName;
    else
        for f = 1:Klist(ct)
            profileNames{end+1} = [baseName, ' ', num2str(f)];
        end
    end
end
nFactors = sum(Klist); % Total NMF factors across all cell types

% Map each factor index back to its cell type (used for gene labels in plots)
% e.g. Klist=[1,2,1,1] -> factorCellTypeIdx = [1, 2, 2, 3, 4]
factorCellTypeIdx = [];
for ct = 1:length(Klist)
    factorCellTypeIdx = [factorCellTypeIdx, repmat(ct, 1, Klist(ct))];
end

% --- Main Loop: one iteration per group (1week, 6week) ---
for g = 1:length(groups)
    groupName       = groups{g};
    current_datasets = allDatasets{g};
    current_genes    = allGenes{g};

    disp(['=== Starting Joint Analysis for Group: ', groupName, ' ===']);

    % Create output directory 
    outputDir = fullfile(sampleDir, ['Joint_Analysis_', groupName]);
    if ~exist(outputDir, 'dir')
        mkdir(outputDir);
    end

 
    disp('Step 1: Pooling data from all samples...');

    pooled_Y    = cell(1, length(Klist));  % One cell per cell type: [TotalSpots x 40 genes]
    pooled_R_vec = [];                     % Reactivity gene counts: [TotalSpots x nReactivityGenes]
    spot_origin  = [];                     % Tracks which sample each spot came from

    for k = 1:length(current_datasets)
        fdir     = current_datasets{k};
        mat_file = fullfile(sampleDir, fdir, [fdir, '_normalized.mat']);

        if ~exist(mat_file, 'file')
            disp(['  Skipping ', fdir, ' (file not found)']);
            continue;
        end

        fprintf('  Loading %s... ', fdir);
        load(mat_file, 'S');


        C = readCelltypeGeneCSV(1);
        C = array2geneCounts(C, S);


        R        = struct;
        R.names  = current_genes;
        R        = array2geneCounts(R, S);

        for ct = 1:length(Klist)
            if k == 1
                pooled_Y{ct} = C.vec(:, :, ct);
            else
                pooled_Y{ct} = [pooled_Y{ct}; C.vec(:, :, ct)];
            end
        end


        pooled_R_vec = [pooled_R_vec; R.vec];
        num_spots   = size(C.vec, 1);
        spot_origin = [spot_origin; repmat({fdir}, num_spots, 1)];

        fprintf('Done. (%d spots added)\n', num_spots);
        clear S C R
    end


    disp('Step 2: Running joint NMF models...');

  
    C_ref = readCelltypeGeneCSV(1);

    U_joint        = cell(1, length(Klist));  % Spatial loadings per cell type
    V_joint        = cell(1, length(Klist));  % Gene profiles per cell type
    U_total_matrix = [];                      % All factors concatenated: [TotalSpots x nFactors]

    for ct = 1:length(Klist)
        disp(['  Processing cell type: ', cellTypeNames{ct}]);

        Y      = pooled_Y{ct};
        K      = Klist(ct);
        nruns  = 500;   
        niters = 1000;
        rho    = 1e2;
        flag   = [0, 0];

        [U, V] = RINMF(Y, K, nruns, niters, rho, flag);

        % Round U if binary flag is set
        if isequal(flag, [1, 0])
            U = round(U);
        end

        U_joint{ct} = U;
        V_joint{ct} = V;


        U_total_matrix = [U_total_matrix, U];


        for k_nmf = 1:K
            h_profile = figure('Visible', 'off', 'Position', [100, 100, 800, 1000]);
            barh(V(:, k_nmf))
            xlabel('Gene Contribution (Global)')
            set(gca, 'YDir', 'reverse')
            set(gca, 'YTick', 1:40);


            ct_idx = factorCellTypeIdx(sum(Klist(1:ct-1)) + k_nmf);
            try
                set(gca, 'YTickLabel', C_ref.names(:, ct_idx));
            catch
                disp('  Warning: Could not set gene labels.');
            end

            title({['Global Joint Factor ', num2str(k_nmf)], ...
                   [cellTypeNames{ct}, ' (', groupName, ')']});

            fname = sprintf('Joint_Profile_%s_%s_Factor%d.png', groupName, cellTypeNames{ct}, k_nmf);
            saveas(h_profile, fullfile(outputDir, fname));
            close(h_profile);
        end
    end
    clear ct k_nmf

  
    disp('Step 4: Calculating joint statistics (Pearson & Cosine)...');

    numSamples = 10^5;  % Random permutations for cosine similarity p-value

    temp = zeros(nFactors, length(current_genes));


    pearson_r_table = array2table(temp, 'RowNames', profileNames, 'VariableNames', current_genes);
    pearson_p_table = array2table(temp, 'RowNames', profileNames, 'VariableNames', current_genes);


    cosine_theta_table = array2table(temp, 'RowNames', profileNames, 'VariableNames', current_genes);
    cosine_p_table     = array2table(temp, 'RowNames', profileNames, 'VariableNames', current_genes);

    for mr = 1:length(current_genes)
        for mct = 1:nFactors

            % Pearson: sqrt transform             
            R_vec_pearson = sqrt(pooled_R_vec(:, mr));
            U_vec_pearson = sqrt(U_total_matrix(:, mct));

            % Cosine: raw counts (no transform )
            R_vec_cosine = pooled_R_vec(:, mr);
            U_vec_cosine = U_total_matrix(:, mct);


            min_len       = min(length(R_vec_pearson), length(U_vec_pearson));
            R_vec_pearson = R_vec_pearson(1:min_len);
            U_vec_pearson = U_vec_pearson(1:min_len);
            R_vec_cosine  = R_vec_cosine(1:min_len);
            U_vec_cosine  = U_vec_cosine(1:min_len);

            % Pearson correlation
            [r, p] = corrcoef(R_vec_pearson, U_vec_pearson);
            pearson_r_table{mct, mr} = r(1, 2);
            pearson_p_table{mct, mr} = p(1, 2);

            % Cosine similarity with empirical p-value
            [cos_r, cos_p, ~] = getcostheta(R_vec_cosine, U_vec_cosine, numSamples);
            cosine_theta_table{mct, mr} = cos_r;
            cosine_p_table{mct, mr}     = cos_p;
        end
    end
    clear mr mct

    
    pearson_r_table.Properties.VariableNames    = strcat('Pearson_r_',     pearson_r_table.Properties.VariableNames);
    pearson_p_table.Properties.VariableNames    = strcat('Pearson_p_',     pearson_p_table.Properties.VariableNames);
    cosine_theta_table.Properties.VariableNames = strcat('Cosine_theta_',  cosine_theta_table.Properties.VariableNames);
    cosine_p_table.Properties.VariableNames     = strcat('Cosine_p_',      cosine_p_table.Properties.VariableNames);

 
    combined_joint_table = [pearson_r_table, pearson_p_table, cosine_theta_table, cosine_p_table];


    csvName = fullfile(outputDir, ['Joint_Stats_', groupName, '.csv']);
    writetable(combined_joint_table, csvName, 'WriteRowNames', true);
    disp(['  Saved joint statistics to: ', csvName]);


    disp('Step 5: Saving joint results...');
    saveFile = fullfile(outputDir, ['Joint_Results_', groupName, '.mat']);
    save(saveFile, 'U_joint', 'V_joint', 'pooled_Y', 'spot_origin', 'cellTypeNames', 'profileNames', 'combined_joint_table');
    disp(['  Saved joint MAT file to: ', saveFile]);

    disp(['=== Finished Joint Analysis for ', groupName, ' ===']);
    disp(' ');
end
clear g

disp('All joint analyses completed.');