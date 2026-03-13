% Final plots: NMF analysis across all datasets
% Runs NMF per cell type, saves spatial maps, gene profiles, and correlations

% --- Configuration ---
sampleDir = '/home/anirban/ubssd/BhanvaPaper_NMF';
datasets = {'1weekB', '1weekC', '1weekD', '1weekE', '1weekF', '1weekG', '1weekH', ...
            '6weekB', '6weekC', '6weekF', '6weekG', '6weekH', '6weekI', '6weekJ'};

% NMF rank per cell type: [Astrocytes, Neurons, OPCs, Microglia]
% Change values here — everything downstream adapts automatically
Klist = [1, 2, 1, 1];

% Cell type names (must match order of Klist)
cellTypeNames = {'Astrocytes', 'Neurons', 'OPCs', 'Microglia'};
% --- End Configuration ---

% Build profile names automatically from Klist and cellTypeNames
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
nFactors = sum(Klist); % Total number of NMF factors across all cell types

% Map each factor index back to its cell type index (used for gene labels)
% e.g. Klist=[1,2,1,1] -> factorCellTypeIdx = [1, 2, 2, 3, 4]
factorCellTypeIdx = [];
for ct = 1:length(Klist)
    factorCellTypeIdx = [factorCellTypeIdx, repmat(ct, 1, Klist(ct))];
end

% Tables to accumulate correlation and cosine similarity results across all datasets
compiled_1week = table();
compiled_6week = table();
compiled_sim_1week = table();
compiled_sim_6week = table();

for ds = 1:length(datasets)
    fdir = datasets{ds};
    disp(['--- Processing dataset: ', fdir, ' ---']);

    % Create output directory for figures if it doesn't exist
    outputFigDir = fullfile(sampleDir, fdir, 'figures');
    if ~exist(outputFigDir, 'dir')
        mkdir(outputFigDir);
    end

    % Load normalized data
    mat_file = fullfile(sampleDir, fdir, [fdir, '_normalized.mat']);

    if ~exist(mat_file, 'file')
        disp(['!!! Skipping ', fdir, ': Normalized .mat file not found at: ', mat_file]);
        continue;
    end

    disp(['Loading: ', mat_file]);
    load(mat_file); % Loads the 'S' structure

    % Load cell type gene reference matrix
    C = readCelltypeGeneCSV(1);
    C = array2geneCounts(C, S);
    C.colNames = cellTypeNames;

    R = struct;

    % Select reactivity genes based on time point (1week vs 6week)
    if contains(fdir, '1week', 'IgnoreCase', true)
        R.names = {'S100a4', 'S100a9', 'Cxcl3', 'Cxcl2', 'Ccl2', 'Il1b',...
                   'Hmox1', 'Lgals3', 'Gpnmb', 'Spp1','Ighm', 'Ctsk','Slpi'};
    elseif contains(fdir, '6week', 'IgnoreCase', true)
        R.names = {'Spp1', 'A2m', 'Lgals3', 'Gpnmb', 'S100a4', 'Serping1', ...
                   'Vim', 'Lyz2', 'S100a6', 'S100a9', 'Slpi', 'Ccl6', ...
                   'Hmox1', 'Il1b', 'Cxcl3', 'Il1rn', 'Srgn', 'Ccl2', 'Cxcl2'};
    else
        % Fallback for unrecognized time points
        disp('Warning: Dataset does not match 1week or 6week. Using default genes.');
        R.names = {'Ccl12', 'C3', 'Itgax', 'Gfap'};
    end

    disp(['Using Reactivity Genes: ', strjoin(R.names, ', ')]);
    R = array2geneCounts(R, S);

    % --- NMF Models ---
    Uct = [];   % Cell spatial loadings (all cell types concatenated)
    Vct = [];   % Gene profiles (all cell types concatenated)

    disp('Running NMF models...');
    for ct = 1:length(Klist)
        Y = C.vec(:, :, ct);
        K = Klist(ct);
        nruns = 500;
        niters = 1000;
        rho = 1e2;
        flag = [0, 0];
        [U, V] = RINMF(Y, K, nruns, niters, rho, flag);

        % Round U if binary flag is set
        if isequal(flag, [1, 0])
            U = round(U);
        end

        Uct = cat(2, Uct, U);
        Vct = cat(2, Vct, V);

        % Save a figure for each NMF factor of this cell type
        for k_nmf = 1:K
            % Cell count spatial map
            h_counts = figure('Visible', 'off');
            stvplot(U(:, k_nmf), S, 20/K)
            title(['Factor ', num2str(k_nmf), ' Cell Counts - ', C.colNames{ct}])
            counts_filename = fullfile(outputFigDir, sprintf('NMF_Fig_%d_%s_Factor%d_CellCounts.png', ct, C.colNames{ct}, k_nmf));
            saveas(h_counts, counts_filename);
            close(h_counts);

            % Gene profile bar chart
            h_profile = figure('Visible', 'off');
            barh(V(:, k_nmf))
            xlabel('Gene counts')
            set(gca, 'YDir', 'reverse')
            set(gca, 'YTick', 1:40);
            set(gca, 'YTickLabel', {C.names{:, ct}});
            title(['Factor ', num2str(k_nmf), ' Gene Profile - ', C.colNames{ct}])

            % Add overall title on the last factor's plot for this cell type
            if k_nmf == K
                sgtitle([C.colNames{ct}, ': Rank ', num2str(K), ' Decomposition'], 'FontSize', 14, 'FontWeight', 'bold', 'Interpreter', 'none');
                set(gca, 'Position', get(gca, 'Position') + [0, -0.05, 0, 0])
            end

            profile_filename = fullfile(outputFigDir, sprintf('NMF_Fig_%d_%s_Factor%d_GeneProfile.png', ct, C.colNames{ct}, k_nmf));
            saveas(h_profile, profile_filename);
            close(h_profile);
        end
    end
    clear ct k_nmf

    % Save canonical NMF factors for this dataset
    canonical_filename = fullfile(sampleDir, fdir, ['canonicalFactors_', fdir, '_new_GFP.mat']);
    save(canonical_filename, 'Uct', 'Vct', 'C', 'profileNames');

    % --- Spatial Maps ---
    disp('Generating spatial maps...');
    sz = 10;
    nTopDisplay = min(5, length(R.names));   % Cap composite top row at 5 genes
    nBotDisplay = min(5, nFactors);          % Cap composite bottom row at 5 factors

    % Composite figure: reactivity genes (top row) and NMF factors (bottom row)
    h_spatial_composite = figure('Units', 'normalized', 'OuterPosition', [0 0 1 1]);

    % Top row: plot up to nTopDisplay reactivity genes; save all individually
    for m = 1:length(R.names)
        if m <= nTopDisplay
            subplot(2, nTopDisplay, m);
            stplot(R.names{m}, S, sz);
            set(gca, 'xtick', [], 'ytick', []);
            axis image; box on;
            title(R.names{m}, 'FontSize', 14);
        end

        % Save individual gene plot
        newFig = figure('Visible', 'off');
        stplot(R.names{m}, S, sz);
        set(gca, 'xtick', [], 'ytick', []);
        axis image; box on;
        title(R.names{m}, 'FontSize', 14);
        saveas(newFig, fullfile(outputFigDir, sprintf('TopRow_Plot_%d_%s.png', m, R.names{m})));
        close(newFig);
    end

    % Bottom row: plot NMF factor spatial maps
    for m = 1:nFactors
        % Normalize factor vector and threshold at 10% of max
        v = Uct(:, m);
        v = v ./ vecnorm(v);
        above_thresh = v >= 0.1 * max(v(:));  % Explicit variable — avoids leaking 'test'

        % Add to composite figure (only first nBotDisplay factors)
        if m <= nBotDisplay
            subplot(2, nTopDisplay, nTopDisplay + m);
            colormap jet;
            scatter(S.x, S.y, sz, 'k');        % Background spots in black
            hold on;
            scatter(S.x(above_thresh), S.y(above_thresh), sz, v(above_thresh), 'filled');
            hold off;
            set(gca, 'xtick', [], 'ytick', []);
            axis image; box on; colorbar;
            title(profileNames{m}, 'FontSize', 14);
        end

        % Save individual NMF spatial map (all factors)
        newFig = figure('Visible', 'off');
        colormap jet;
        scatter(S.x, S.y, sz, 'k');
        hold on;
        scatter(S.x(above_thresh), S.y(above_thresh), sz, v(above_thresh), 'filled');
        hold off;
        set(gca, 'xtick', [], 'ytick', []);
        axis image; box on; colorbar;
        title(profileNames{m}, 'FontSize', 14);
        saveas(newFig, fullfile(outputFigDir, sprintf('BottomRow_Plot_%d_%s.png', m, profileNames{m})));
        close(newFig);
    end

    saveas(h_spatial_composite, fullfile(outputFigDir, 'Spatial_Composite_Map.png'));
    close(h_spatial_composite);

    % --- Gene Profiles (Wide composite figure) ---
    disp('Generating gene profile maps...');
    h_profiles_composite = figure('Visible', 'off', 'Position', [100, 100, 1800, 600]);
    for m = 1:nFactors
        subplot(1, nFactors, m)
        barh(Vct(:, m))
        xlabel('Gene counts')
        set(gca, 'YDir', 'reverse')
        set(gca, 'YTick', 1:40);

        % Use factorCellTypeIdx to look up the correct gene labels for this factor
        % This generalizes correctly regardless of how many factors each cell type has
        ct_idx = factorCellTypeIdx(m);
        set(gca, 'YTickLabel', {C.names{:, ct_idx}});

        title(profileNames{m})
    end
    clear m
    sgtitle('NMF Cell type factor gene profiles')

    saveas(h_profiles_composite, fullfile(outputFigDir, 'GeneProfile_Composite.png'));
    close(h_profiles_composite);

    % --- Correlations: reactivity genes vs NMF factors ---
    disp('Calculating correlations...');
    temp = zeros(nFactors, length(R.names));
    corr_table = array2table(temp, 'RowNames', profileNames, 'VariableNames', R.names);
    pval       = array2table(temp, 'RowNames', profileNames, 'VariableNames', R.names);

    for mr = 1:length(R.names)
        for mct = 1:nFactors
            % Square-root transform before correlating
            R_vec_mr = sqrt(R.vec(:, mr));
            Uct_mct  = sqrt(Uct(:, mct));

            % Align lengths in case of mismatch
            min_len  = min(length(R_vec_mr), length(Uct_mct));
            R_vec_mr = R_vec_mr(1:min_len);
            Uct_mct  = Uct_mct(1:min_len);

            [r, p] = corrcoef(R_vec_mr, Uct_mct);
            corr_table{mct, mr} = r(1, 2);
            pval{mct, mr}       = p(1, 2);
        end
    end
    clear mr mct

    % Prefix p-value column names and merge with correlation table
    pval.Properties.VariableNames = strcat('pval_', pval.Properties.VariableNames);
    combined_table = [corr_table, pval];

    % Tag each row with sample name for later compilation
    combined_table.SampleName = repmat({fdir}, height(combined_table), 1);

    % Save per-dataset results
    corr_filename = fullfile(sampleDir, fdir, ['combined_results_', fdir, '_new.csv']);
    writetable(combined_table, corr_filename, 'WriteRowNames', true);
    disp(['Correlation results saved to ', corr_filename]);

    % Prepare table for cross-dataset compilation
    temp_tbl = combined_table;
    temp_tbl.CellType = combined_table.Properties.RowNames;  % Preserve row labels as a column
    temp_tbl.Properties.RowNames = {};                        % Clear row names to allow stacking

    if contains(fdir, '1week', 'IgnoreCase', true)
        compiled_1week = [compiled_1week; temp_tbl];
    elseif contains(fdir, '6week', 'IgnoreCase', true)
        compiled_6week = [compiled_6week; temp_tbl];
    end

    % --- Cosine Similarity: reactivity genes vs NMF factors ---
    disp('Calculating cosine similarity...');
    numSamples = 10^5;  % Number of random samples for empirical p-value estimation

    % Initialize similarity result tables (nFactors x nGenes)
    temp_sim = zeros(nFactors, length(R.names));
    sim_costheta = array2table(temp_sim, 'RowNames', profileNames, 'VariableNames', R.names);
    sim_pval     = array2table(temp_sim, 'RowNames', profileNames, 'VariableNames', R.names);

    for mr = 1:length(R.names)
        for mct = 1:nFactors
            R_vec_mr = R.vec(:, mr);
            Uct_mct  = Uct(:, mct);

            % Align lengths in case of mismatch
            min_len  = min(length(R_vec_mr), length(Uct_mct));
            R_vec_mr = R_vec_mr(1:min_len);
            Uct_mct  = Uct_mct(1:min_len);

            % getcostheta returns cos(theta), empirical p-value, and null distribution
            [r, p, ~] = getcostheta(R_vec_mr, Uct_mct, numSamples);
            sim_costheta{mct, mr} = r;
            sim_pval{mct, mr}     = p;
        end
    end
    clear mr mct

    % Display results (p=0 means < 1/numSamples; p=1 means all nulls had larger overlap)
    disp('Cosine Similarity (cos theta):');
    disp(sim_costheta);
    disp('Cosine Similarity p-values:');
    disp(sim_pval);
    disp(['Note: p=0 means < ', num2str(1/numSamples), '; p=1 means all null simulations had larger overlap']);

    % Prefix p-value column names and merge cosine similarity tables
    sim_pval.Properties.VariableNames = strcat('pval_', sim_pval.Properties.VariableNames);
    sim_combined = [sim_costheta, sim_pval];
    sim_combined.SampleName = repmat({fdir}, height(sim_combined), 1);

    % Save cosine similarity results for this dataset
    sim_filename = fullfile(sampleDir, fdir, ['cosine_similarity_', fdir, '.csv']);
    writetable(sim_combined, sim_filename, 'WriteRowNames', true);
    disp(['Cosine similarity results saved to ', sim_filename]);

    % Accumulate cosine similarity for cross-dataset compilation
    temp_sim_tbl = sim_combined;
    temp_sim_tbl.CellType = sim_combined.Properties.RowNames;
    temp_sim_tbl.Properties.RowNames = {};

    if contains(fdir, '1week', 'IgnoreCase', true)
        compiled_sim_1week = [compiled_sim_1week; temp_sim_tbl];
    elseif contains(fdir, '6week', 'IgnoreCase', true)
        compiled_sim_6week = [compiled_sim_6week; temp_sim_tbl];
    end

    disp(['--- Finished: ', fdir, ' ---']);
    disp(' ');

    clear S C R Uct Vct h_counts h_profile h_spatial_composite h_profiles_composite newFig
end
clear ds

% --- Save Compiled Results Across All Datasets ---
disp('--- Compiling All Results ---');

if ~isempty(compiled_1week)
    outfile_1wk = fullfile(sampleDir, 'Compiled_1Week_Results.csv');
    writetable(compiled_1week, outfile_1wk);
    disp(['Saved compiled 1-Week results to: ', outfile_1wk]);
end

if ~isempty(compiled_6week)
    outfile_6wk = fullfile(sampleDir, 'Compiled_6Week_Results.csv');
    writetable(compiled_6week, outfile_6wk);
    disp(['Saved compiled 6-Week results to: ', outfile_6wk]);
end

% Save compiled cosine similarity results
if ~isempty(compiled_sim_1week)
    outfile_sim_1wk = fullfile(sampleDir, 'Compiled_1Week_CosineSimilarity.csv');
    writetable(compiled_sim_1week, outfile_sim_1wk);
    disp(['Saved compiled 1-Week cosine similarity to: ', outfile_sim_1wk]);
end

if ~isempty(compiled_sim_6week)
    outfile_sim_6wk = fullfile(sampleDir, 'Compiled_6Week_CosineSimilarity.csv');
    writetable(compiled_sim_6week, outfile_sim_6wk);
    disp(['Saved compiled 6-Week cosine similarity to: ', outfile_sim_6wk]);
end

disp('All datasets processed and compiled.');