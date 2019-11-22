% This a master script that runs mapping with cross-validation

% load spatial profiles
spatial = tdfread('../data/Omer_astrocytes/p56_mix.csv', ','); % 16 markers
genes = unique(string(cellstr(spatial.genes)));

% run using all genes 16
Zonation_main_with_cell_groups_p56_new('', false, false)

% leave one out
for ind = 1:size(genes, 1)
   Zonation_main_with_cell_groups_p56_new(genes(ind), true, false)
end
