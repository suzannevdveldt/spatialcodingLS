function [mean_dispersion, centroid] = computeCentroidDispersion(binarized_trace, interp_behav_vec, inclusion_vector)
%% Extracts spatial coordinates corresponding to active periods. Derives median centroid and dispersion.
% binary_trace: a binary vector of a given neuron
% interp_behav_vec: interpolated, 2D position of the mouse
% inclusion_vector: a logical vector indicating frames to be included


%% Extract coordinates corresponding to active periods
coordinates = [];
for frame_i = 1:length(binarized_trace)  
	if binarized_trace(frame_i) && inclusion_vector(frame_i)
		coordinates(end+1,:) = interp_behav_vec(frame_i,:);
	end
end

if ~isempty(coordinates); 
centroid(1) = median(coordinates(:,1));
centroid(2) = median(coordinates(:,2));

dispersion = []; 
for i = 1:size(coordinates,1)
		dispersion(end+1) = sqrt((coordinates(i,1)-centroid(1)).^2 + (coordinates(i,2)-centroid(2)).^2);
end

mean_dispersion = mean(dispersion);

%% Plotting (optional)
% plot(interp_behav_vec(:,1),interp_behav_vec(:,2),'color', [0.5 0.5 0.5])
% hold on
% scatter(coordinates(:,1), coordinates(:,2),[], '.', 'MarkerFaceColor', 'k')
% scatter(centroid(1), centroid(2),'MarkerFaceColor', 'r')
% %viscircles(centroid(1), centroid(2),mean(dispersion))
else
    mean_dispersion = [NaN]; 
    centroid = [NaN]; 
end
    
end