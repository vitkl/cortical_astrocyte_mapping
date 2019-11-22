function [zoneParams] = empirical_distribution_per_gene(geneName, spatial_genes, spatial_layer, spatial_depth, spatial_int, meanExp, rescale_by_scrna)

% For each landmark gene:
% (1) Pool together segmented cells from smFISH images and divide them into 5 cortical areas. 
% (2) Construct the histogram of expression levels at every layer (sampling distribution of the gene in units of fluorescence intensity concentrations).
% (3) Fit the sampling distributions with gamma functions. 

% Input:
% geneName - a string representing the name of the specific landmark gene
% spatial_genes - an array specifying which genes measured in cells 
% spatial_layer - an array specifying the layer
% spatial_depth - an array specifying the continous depth
% spatial_int - an array containing molecule counts for spatial_genes in
% cells of varying spatial_layer & spatial_depth 
% meanExp - the mean expression of the gene across all sequenced cells 
% Output: 
% 'zoneParams' - parameters of the Gamma distribution at every layer

gene_sel = spatial_genes == geneName;

numZones = length(unique(spatial_layer));
layers = sort(unique(spatial_layer));
spatial_layer = spatial_layer(gene_sel);

Dist_all = spatial_depth(gene_sel);
Int_all = spatial_int(gene_sel);

% remove cells/dots that are not in the region of interest
% i.e. cells that are not on the radial axis between the central to the portal vein
toRemove = isnan(Dist_all);
All = [Dist_all, Int_all];
All(toRemove,:) = [];

DIST_CV = 1;
INTENSITY = 2;

% --- divide into zones ---
lims = linspace(0,1,numZones+1);
Sections = nan(numZones,length(All));
zonesInd = nan(length(All),1);

for index=1:numZones
    Sections(index,:) = spatial_layer == layers(index);
    zonesInd(find(Sections(index,:))) = index;
end


if rescale_by_scrna == true
% -- Divide by the mean and mult by 'meanExp' (the factor calculated from sequenced cell data) --
% (we use mean and not median, because the median for Glul is 0)
    m = mean(All(:,INTENSITY));
    All(:,INTENSITY) = All(:,INTENSITY) / m;
    All(:,INTENSITY) = meanExp * All(:,INTENSITY);
end
    
% --- create Gamma distribution ---
zoneParams = nan(numZones,2);
for i=1:numZones
    zoneParams(i,:) = gamfit(All((zonesInd == i),INTENSITY));
end
zoneParams(isnan(zoneParams)) = 0;

end
