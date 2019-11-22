function [] = Zonation_main_with_cell_groups_p56_new(leave_out, leave_log, prof_6)
% Algorithm for spatial reconstruction of mouse cortical astrocyte layer profiles
% Main steps:
% (1) Estimate the distributions of cellular expression of each of the landmark genes in each layer, based on smFISH data.
% (2) For each cell, obtain the probabilities to have the given landmark gene expression vector at each layer (sampling distribution).
% (3) Calculate for each cell, the posterior probability as the product of the sampling distribution and 
%     the prior probability of the cell to belong to that layer (prior).
% (4) Multiply the UMI table (the expression of each gene in each cell) by the (posterior) probability matrix (the probability for each cell to be at each zone)

% Input: load the following variables:
% 'seq_data': UMI counts per gene and cell (after filtering and normalisation).
%             Rows = genes, cols = cells.
% skip row and column one (gene names and cell ids)
seq_data = csvread('../data/Holt_astrocytes/normalised_expression_mat_UMI_TRUE.csv'); 
% 'all_genes': the names of the genes (not including ERCC spike-in).
T = readtable('../data/Holt_astrocytes/normalised_expression_genes_UMI_TRUE.csv', ...
    'ReadVariableNames',false);
all_genes = table2array(T(:,1));
% load spatial profiles
if prof_6 == true
    spatial = tdfread('../data/Omer_astrocytes/20190122_layerAst_P56cor.tsv'); 
    spatial_int = spatial.x0x22spotcounts0x22;
    spatial_genes = spatial.x0x22genes0x22;
    spatial_depth = spatial.x0x22normalisedDepth0x22;
    spatial_layer = spatial.x0x22ctxdepthInterval0x22;
end

if prof_6 ~= true
    spatial = tdfread('../data/Omer_astrocytes/p56_mix.csv', ','); 
    spatial_int = spatial.spotcounts;
    spatial_genes = string(cellstr(spatial.genes));
    spatial_depth = spatial.normalisedDepth;
    spatial_layer = string(spatial.ctxdepthInterval);
end


% Prior probability of sampling a cell at each  layer - reflect the
% number of cells in these zones - compute using proportion of cells in
% each zone
P_Z=[0.1767786 0.1869367    0.2338591    0.2233959    0.1790297]';

% set parameters
analysis_name = 'astro_16mark_UMI_Ntotal180000_1000k_out_norm'
NUMZONES = size(unique(spatial_layer), 1); % We divide the cortex into 'NUMZONES' zones
Bins_num_UMIs = 8; % 8 bins representing cells with similar sampling levels 
Ntotal = 280000; % an estimate of the total number of mRNA molecules in a tetraploid hepatocyte / 28 to produce the number expected for astrocyte
normalize = true; % normalize so that the profiles (genes*layers) represent the total fraction of UMI
GENES = unique(spatial_genes); % set of landmark genes

if leave_log == true
    leave_out = leave_out; % leave out this landmark gene
    GENES = GENES(GENES ~= leave_out);
    analysis_name = strcat(analysis_name, '_', leave_out)
end 

% find proportion of cells expressing each landmark gene
gene_means = mean(seq_data > 0, 2);
all_genes(ismember(all_genes, GENES))
gene_means(ismember(all_genes, GENES))

reps = 1000; % number of booststrapping runs
N = 1000000; % number of sampled cells
if leave_log == true % use less cells for leave out experiment for speed
    N = 50000;
end
Vf = 1 / 0.1; % a factor used to correct for the fact that the smFISH measurements do not capture the entire hepatocyte volume
n_cores = 4; % number of cores for parallel evaluation of bootstrapping
rescale_by_scrna = true; % re-scale smFISH counts so that the mean equals fraction of reads for that gene in scRNA-seq data (add nonscaled to dataset name if false) 

% *** load 'seq_data' & 'all_genes' ***
number_of_cells = length(seq_data(1,:));

% calculate the fraction of each LM gene in the sequenced cells
fraction_gene = nan(length(GENES), 1);
totalUMIs = sum(sum(seq_data));
for g=1:length(GENES)
    geneI = strcmpi(GENES{g}, all_genes);
    fraction_gene(g) = sum(seq_data(geneI,:)) / totalUMIs;
end

% ----------- step (1) -----------
% Generate the Gamma distributions
 for g=1:length(GENES)
     p = empirical_distribution_per_gene(GENES{g}, spatial_genes, spatial_layer, spatial_depth, spatial_int, fraction_gene(g), rescale_by_scrna); 
     eval(['LM.' GENES{g} '.Gamma = p ;'])
 end
 clear p

% -- estimate a cell-specific sampling distribution --
% % discretize sampling levels into 8 bins representing cells with similar sampling levels 
 cellGroups = cell(1,Bins_num_UMIs);
 [~,edges] = histcounts(sum(seq_data),Bins_num_UMIs);
 for k=1:Bins_num_UMIs
     ind = find((sum(seq_data) > edges(k)) .* (sum(seq_data) <= edges(k+1)));
     cellGroups{k} = ind;
 end
% clear edges

% compute beta for each cell
 beta_i = sum(seq_data) / Ntotal;
% ---------
 for g=1:length(GENES) % (1) for each LM gene
     geneName = GENES{g};
     GammaParams = LM.(geneName).Gamma;
     
     for i=1:Bins_num_UMIs % (2) for each group of cells (bin)
         Q = nan(N,NUMZONES);
         cellsInGroup = cellGroups{i}; % extract the indices of the cells in group i
         beta = median(beta_i(cellsInGroup)); % compute a median sampling for each group
         
         for z=1:NUMZONES % (3) for each zone
             % sample N cells from the Gamma distribution of gene g in zone z
             % down-sample the values with Poisson and save it in Q
             data = gamrnd(GammaParams(z,1), GammaParams(z,2), [N,1]);
             if rescale_by_scrna == true
                 %multiply by Ntotal to counteract normalisation specified by `rescale_by_scrna`
                 data = data * Ntotal;
             end
                 Q(:,z) = poissrnd(data * Vf * beta) ./ Vf ;
         end % (3)
         
         % save the information of the cell group:
         % Q : [N virtual cells x NUMZONES] matrix
         % ind : the indices of the cells in the group
         % beta : the median beta_i of the cells in the group
         eval(['CellGroups.group', num2str(i), '.Q.' geneName ' = Q;'])
         eval(['CellGroups.group', num2str(i), '.ind = cellsInGroup;'])
         eval(['CellGroups.group', num2str(i), '.beta = beta;'])
     end % (2)
 end % (1)



%%
% ----------- step (2) -----------
% Construct a matrix P (cells x layers) which will hold the probability of each cell to belong to each layer layer, 
% given the vector of expression of the 6 landmark genes 

beta_factor = nan(number_of_cells,1);
for g=1:length(GENES) % (1) for each LM gene
    geneName = GENES{g};
    geneIndex = strcmpi(geneName,all_genes);
    prob_in_zone_g = nan(number_of_cells, NUMZONES);
    
    for i=1:Bins_num_UMIs % (2) for each group of cells (with similar total #UMIs)
        eval(['cellsInGroup = CellGroups.group' num2str(i) '.ind ;' ])
        eval(['Q = CellGroups.group' num2str(i) '.Q.' geneName ';'])
        beta = median(beta_i(cellsInGroup));
        
        for z=1:NUMZONES % (3) for each zone
            binsH=0:max(Q(:,z)); % define the bins for the hist
            if length(binsH)== 1
                % in case most values are zeros
                % there is high prob to be 0 and nearly zero prob to be in 1
                Counts = [1-eps; eps];
                X = [0 1];
            else
                [Counts, X] = hist(Q(:,z), binsH);
                Counts = Counts/sum(Counts);
            end
            
            for c=1:length(cellsInGroup) % (4) for each cell in the group
                % compute the prob of cell i to belong to each zone, based on the
                % profile of the current landmark gene g
                beta_factor(cellsInGroup(c)) = (beta / beta_i(cellsInGroup(c)));
                gene_UMIs = round(beta_factor(cellsInGroup(c)) * seq_data(geneIndex, cellsInGroup(c)));
                 
                [~, indexMin] = min(abs(X - gene_UMIs)); % closest value in the hist
                prob_in_zone_g(cellsInGroup(c), z) =  Counts(indexMin);
            end % (4)
        end % (3)
    end % (2)
    LM.(geneName).cellsProbZones = prob_in_zone_g;
end % (1)

% ----------- step (3) -----------
P = ones(size( LM.(GENES{1}).cellsProbZones));
for g=1:length(GENES)
    P = P.*LM.(GENES{g}).cellsProbZones;
end

P = P .* repmat(P_Z',size(P,1),1); % mult by the probability to see a cell in each zone (smaller for inner zones)
[~, zones] = max(P,[],2);
P(:,NUMZONES+1) = zones; % add to P a column representing the most likely layer for each cell

% ----- Bootstrapping to obtain standard errors for the mean zonation profiles -----
tic
[Fuzzy, KW] = bootStrap(seq_data, reps, P, NUMZONES, n_cores, normalize);
toc

% *** save the results ***
'saving predicted profiles'
res_folder = strcat('./', analysis_name)
status = mkdir(res_folder)
% mean profile and sd
save(strcat(res_folder, '/zonation_profile_and_sd.mat'), 'Fuzzy')
csvwrite(strcat(res_folder, '/zonation_profile_genes_x_zones.csv'), Fuzzy.Mat)
csvwrite(strcat(res_folder, '/zonation_profile_sd_genes_x_zones.csv'), Fuzzy.SE_bootstrap)
% p- and q-value
save(strcat(res_folder, '/zonation_profile_p_q_val.mat'), 'KW')
csvwrite(strcat(res_folder, '/zonation_profile_p_val_genes_x_1.csv'), KW.pval)
csvwrite(strcat(res_folder, '/zonation_profile_q_val_genes_x_1.csv'), KW.qval)

csvwrite(strcat(res_folder, '/prob_cell_assignment_cells_x_zones_plus1.csv'), Fuzzy.P)
end

