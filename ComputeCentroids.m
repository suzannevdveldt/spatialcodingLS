function centroids = ComputeCentroids(ms,SFP)
%PLOTCELLSCOLORINDEX Summary of this function goes here
%   Detailed explanation goes here

centroids = zeros(ms.numNeurons,2);

for cell_i = 1:ms.numNeurons;
    current_cell = SFP(cell_i,:,:); % get SFP of particular neuron
    current_cell = squeeze(current_cell(1, :, :));
    current_cell_binarized =logical(current_cell);  % binarize this SFP 
    properties = regionprops(current_cell_binarized,'centroid');
    
    centroids(cell_i,:) = properties.Centroid(1,:);
end

end
