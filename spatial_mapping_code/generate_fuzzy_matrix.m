function [Fuzzy, KW] = generate_fuzzy_matrix(seq_data, P, SE, normalize)
% Construct a weighted probability matrix.
% Multiply the UMI table by the weighted probability matrix to obtain the gene zonation matrix
% 
% Input:
% 'seq_data': UMI counts per gene and cell (after filtering).
%            (genes x cells) matrix. Its elements are the UMI count.
% 'P': cells x (zones+1) matrix. 1-9 columns are the propabilities for a cell to belong to each layer.
%      The last column holds the maximum-likelihood layer.
% 'SE': genes x zones matrix. Holds the standard errors for the mean zonation profiles
% 
% Output: 
% 'Fuzzy': struct which holds: (1) the gene zonation matrix (genes x zones matrix with the
%                                  average expression level in the same units as seq_data, 
%                                  or normalised to represent fraction of total UMI)
%                              (2) the standard errors for the mean zonation profiles (SE)
% 'KW': struct which holds the p-values and q-values of the kruskalwallis-test for each gene.

% normalize by the sum
if normalize == true
    seq_data = seq_data ./ repmat(nansum(seq_data), length(seq_data),1);
end
    
numZones = size(P,2)-1;

weights = P(:, 1:numZones) ;
% normalizeover the different zones
% normalize the posterior matrix by the column sums to obtain a weight matrix
% ensures that the number of cells in each lobule layer does not affect the average layer expression
weights = weights ./ repmat(sum(weights,2), 1, length(weights(1,:)));
weights(isnan(weights)) = 0;
weights = weights ./ repmat(sum(weights),size(weights,1),1);
pval = nan(length(seq_data), 1);
% [genes x cells] * [cells x probs] = [sum genes over all cells]
fuzzyMat = seq_data * weights;
Fuzzy.Mat = fuzzyMat;
Fuzzy.SE_bootstrap = SE;
Fuzzy.P = weights;


%  ------------- kruskalwallis -------------
% How significant is the zonation of each gene - are the distributions of
% reads different across cells assigned to different spatial bins
% 
for i=1:size(seq_data,1),
    pval(i) = kruskalwallis(seq_data(i,:),P(:,size(P,2))','off');
end
[~, qval] = mafdr(pval);
KW.pval = pval;
KW.qval = qval;

end

