function [] = Zonation_main()
% Algorithm for spatial reconstruction of liver zonation profiles
% Main steps:
% (1) Estimate the distributions of cellular expression of each of the landmark genes in each layer, based on smFISH data.
% (2) For each cell, obtain the probabilities to have the given landmark gene expression vector at each layer (sampling distribution).
% (3) Calculate for each cell, the posterior probability as the product of the sampling distribution and 
%     the prior probability of the cell to belong to that layer (prior).
% (4) Multiply the UMI table (the expression of each gene in each cell) by the (posterior) probability matrix (the probability for each cell to be at each zone)

% Input: load the following variables:
% 'seq_data': UMI counts per gene and cell (after non-hepatocytes and cells with low fraction of Albumin were filtered-out).
%             Rows = genes, cols = cells.
% 'all_genes': the names of the genes (not including ERCC spike-in).
% 'P_Z': an array with prior probabilities for sampling a cell from each lobule layer.

NUMZONES = 9; % We divide the lobule into 'NUMZONES' zones
Ntotal = 787054; % an estimate of the total number of mRNA molecules in a tetraploid hepatocyte
GENES = {'Glul','Cyp2f2','Alb','Cyp2e1','Ass1','Asl'}; % set of landmark genes

% *** load 'seq_data' & 'all_genes' ***
number_of_cells = length(seq_data(1,:));

% calculate the fraction of each LM gene in the sequenced cells
fraction_gene = nan(length(GENES), 1);
totalUMIs = sum(sum(seq_data));
for g=1:length(GENES)
    geneI = strcmp(GENES{g}, all_genes);
    fraction_gene(g) = sum(seq_data(geneI,:)) / totalUMIs;
end

% ----------- step (1) -----------
% Generate the Gamma distributions
for g=1:length(GENES)
    p = empirical_distribution_per_gene(GENES{g}, NUMZONES, fraction_gene(g)); 
    eval(['LM.' GENES{g} '.Gamma = p ;'])
end
clear p


% -- estimate a cell-specific sampling distribution --
% compute beta for each cell
beta_i = sum(seq_data) / Ntotal;
% discretize sampling levels into 8 bins representing cells with similar sampling levels 
Bins_num_UMIs = 8; 
cellGroups = cell(1,Bins_num_UMIs);
[~,edges] = histcounts(sum(seq_data),Bins_num_UMIs);
for k=1:Bins_num_UMIs
    ind = find((sum(seq_data) > edges(k)) .* (sum(seq_data) <= edges(k+1)));
    cellGroups{k} = ind;
end
clear edges

% ---------
N = 500000; % number of sampled cells
Vf = 10; % a factor used to correct for the fact that the smFISH measurements do not capture the entire hepatocyte volume

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
            data = Ntotal * gamrnd(GammaParams(z,1), GammaParams(z,2), [N,1]);
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
% Construct a matrix P (cells x lobules) which will hold the probability of each cell to belong to each lobule layer, 
% given the vector of expression of the 6 landmark genes 

beta_factor = nan(number_of_cells,1);
for g=1:length(GENES) % (1) for each LM gene
    geneName = GENES{g};
    geneIndex = strcmp(geneName,all_genes);
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
                
                % For Glul in layer 1. In case the cells has more than the maximal expected #UMIs of Glul,
                % the cell should have high probability to belong to zone 1
                if (strcmp('Glul',geneName) && (z == 1) && gene_UMIs > max(X)) 
                    prob_in_zone_g(cellsInGroup(c), z) = 1;
                else 
                    [~, indexMin] = min(abs(X - gene_UMIs)); % closest value in the hist
                    prob_in_zone_g(cellsInGroup(c), z) =  Counts(indexMin);
                end
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

% Prior probability of sampling a cell at each lobule layer
%  *** load 'P_Z' ***
P = P .* repmat(P_Z',size(P,1),1); % mult by the probability to see a cell in each zone (smaller for inner zones)
[~, zones] = max(P,[],2);
P(:,NUMZONES+1) = zones; % add to P a column representing the most likely layer for each cell

% ----- Bootstrapping to obtain standard errors for the mean zonation profiles -----
reps = 1000; % number of booststrapping runs
[Fuzzy, KW] = bootStrap(seq_data, reps, P, NUMZONES);

% *** save the data ***

end

