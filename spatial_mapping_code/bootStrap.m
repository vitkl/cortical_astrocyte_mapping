function [Fuzzy, KW] = bootStrap(seq_data, reps, P, NUMZONES, n_cores, normalize)
% Obtain standard errors for the mean zonation profiles
% Input: 
% 'seq_data': UMI counts per gene and cell (after non-hepatocytes and cells with low fraction of Albumin were filtered-out).
%             Rows = genes, cols = cells.
% 'reps': number of booststrapping runs
% 'P': cells x (zones+1) matrix. 1-9 columns are the propabilities for a cell to belong to each layer.
%      The last column holds the maximum-likelihood layer.
% 'NUMZONES': the number of zones/layers (9)
% 
% Output: 
% 'Fuzzy': struct which holds: (1) the gene zonation matrix (genes x zones matrix with the average expression level in fraction of 
%                                  total cellular mRNA, attributed to every gene at lobule layer)
%                              (2) the standard errors for the mean zonation profiles (SE)
% 'KW': struct which holds the p-values and q-values of the kruskalwallis-test for each gene.

numCells = length(seq_data(1,:));
bootGenes = zeros(size(seq_data,1),NUMZONES,reps); %  (genes x zones x iterations) matrix 

parfor (i=1:reps, n_cores)
    samples = randsample(1:numCells, numCells, 'true'); % sample 'numCells' with replacment
    sampled_data = seq_data(:,samples);
    if normalize == true
        sampled_data = sampled_data./repmat(sum(sampled_data),size(sampled_data,1),1);
    end
    P_boot = P(samples, 1:NUMZONES);
    P_boot = P_boot ./ repmat(sum(P_boot,2), 1, length(P_boot(1,:)));
    P_boot(isnan(P_boot)) = 0;
    P_boot = P_boot ./ repmat(sum(P_boot),size(P_boot,1),1);
    
    bootGenes(:,:,i)= sampled_data * P_boot;
end

% calculate SE
for i=1:size(bootGenes,1), for j=1:size(bootGenes,2), SE(i,j)=std(bootGenes(i,j,:)); end; end

[Fuzzy, KW] = generate_fuzzy_matrix(seq_data, P, SE, normalize);

end


