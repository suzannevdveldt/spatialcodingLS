function [extended_calcium_analysis_output] = extended_calcium_analysis_plotCell(cell_id);

% Potentially run this by running select_folders_for_analysis.m
%%                             
 load 'behav.mat';
  load 'ms.mat';

plotting = 1; 
save_figure = 0;
bin_size = 3;
min_speed_threshold = 5; % 2 cm.s-1
params.smoothing = true;
numShuffles = 10;

if ~isfield(ms,'Binary');
    disp('Extracting binary information');
    ms = msExtractBinary(ms);
end

%% determine behav and ms timestamps and convert to seconds
behav_time = behav.time/1000; %convert to seconds 
behav_vec = behav.position;
[behav_time, idx]=unique(behav_time); %getting rid of duplicate timestamps
behav_vec = behav_vec(idx,:); %excluding the positions associated with these double timestamps
ca_time = ms.time/1000; % convert to seconds 
[ca_time, idx]=unique(ca_time); %getting rid of duplicate timestamps

cell_ids = zeros(ms.numNeurons,1);
is_place_cell = zeros(ms.numNeurons,1);

%% Binarize calcium trace
sampling_frequency = 30; % This data set has been sampled at 30 images per second
%z_threshold = 4; % A 2 standard-deviation threshold is usually optimal to differentiate calcium ativity from background noise
%[binarized_trace] = extract_binary(ca_trace, sampling_frequency, z_threshold);
%ms.BinarizedTraceExtended(:,cell_id) = binarized_trace; 

%% Interpolate behavior
% In most cases, behavior data has to be interpolated to match neural temporal
% activity assuming it has a lower sampling frequency than neural activity
interp_behav_vec(:,1) = interpolate_behavior(behav_vec(:,1), behav_time, ca_time); % in the X dimension
interp_behav_vec(:,2) = interpolate_behavior(behav_vec(:,2), behav_time, ca_time); % in the Y dimension
interp_behav_vec(end,:) = interp_behav_vec(end-1,:); % Make sure to interpolate or remove every NaN so that every timepoint has a corresponding behavioral state

%% Compute velocities
% In this example, we will ignore moments when the mouse is immobile
[velocity] = extract_velocity(interp_behav_vec, ca_time);

%% Find periods of immobility
% This will be usefull later to exclude immobility periods from analysis
running_ts = velocity > min_speed_threshold;

%% Create an inclusion vector to isolate only a specific behavior
inclusion_vector = running_ts; % Only include periods of running

%% Compute occupancy and joint probabilities
% You can use 'min(behav_vec)' and 'max(behav_vec)'  to estimate the
% boundaries of the behavior vector (in our case, location in space)
%% for paper figs
% min_X = 0;
% min_Y = 0;
% max_X = 49; % Set these to the borders of the maze used. !!! Has to be the same for every mouse for consistent MI computation !!!
% max_Y = 49;
% 
% X_bin_vector = min_X:bin_size:max_X+bin_size;
% X_bin_centers_vector = X_bin_vector + bin_size/2;
% X_bin_centers_vector(end) = [];
% Y_bin_vector = min_Y:bin_size:max_Y+bin_size;
% Y_bin_centers_vector = Y_bin_vector + bin_size/2;
% Y_bin_centers_vector(end) = [];

min_X = 0;
min_Y = 0;
max_X = max(behav.position(:,1)); % Set these to the borders of the maze used. !!! Has to be the same for every mouse for consistent MI computation !!!
max_Y = max(behav.position(:,2));

X_bin_vector = min_X:bin_size:max_X+bin_size;
X_bin_centers_vector = X_bin_vector + bin_size/2;
X_bin_centers_vector(end) = [];
Y_bin_vector = min_Y:bin_size:max_Y+bin_size;
Y_bin_centers_vector = Y_bin_vector + bin_size/2;
Y_bin_centers_vector(end) = [];

%% for dreadd exps 

%% 
ca_trace = ms.RawTraces(:,cell_id);

ca_trace = ca_trace(idx); %excluding the positions associated with these double timestamps

%for singel cell activity
binarized_trace = ms.Binary(:, cell_id); 

% For population activity
%binarized_trace = sum(ms.Binary,2)

[MI, posterior, occupancy_map, prob_being_active, likelihood] = extract_2D_information(binarized_trace, interp_behav_vec, X_bin_vector, Y_bin_vector, inclusion_vector);
peak_joint_probability = max(likelihood(:));
% 
% if plotting
% plotting_fig = figure;
% subplot(3,1,1)
% imagesc(X_bin_centers_vector,Y_bin_centers_vector,likelihood)
% daspect([1 1 1])
% title 'Probability of firing in location'
% xlabel 'Position (cm)'
% ylabel 'Position (cm)'
% colorbar
% subplot(3,1,2)
% imagesc(X_bin_centers_vector,Y_bin_centers_vector,occupancy_map)
% daspect([1 1 1])
% title 'Relative occupancy'
% xlabel 'Position (cm)'
% ylabel 'Position (cm)'
% colorbar
% subplot(3,1,3)
% surf(X_bin_centers_vector,Y_bin_centers_vector,posterior)
% daspect([1 1 0.01])
% title 'Posterior probability density function'
% xlabel 'Position (cm)'
% ylabel 'Position (cm)'
% colorbar
% end

if save_figure;
    saveas(plotting_fig, sprintf('%d.jpg',cell_id));
end

%% Shuffle data
shuffled_data = zeros(length(Y_bin_centers_vector), length(X_bin_centers_vector), numShuffles);
shuffled_peak_joint_probability = [];

for k = 1:numShuffles
    random_ts = ceil(rand*length(binarized_trace));
    shuffled_binarized = zeros(length(binarized_trace),1);
    
    % Permute the trace
    shuffled_binarized(1:random_ts) = binarized_trace(end-random_ts+1:end);
    shuffled_binarized(random_ts+1:end) = binarized_trace(1:end-random_ts);
    
    [shuffled_MI(k), shuffled_posterior(:,:,k), ~ , ~, shuffled_likelihood(:,:,k)] = extract_2D_information(shuffled_binarized, interp_behav_vec, X_bin_vector, Y_bin_vector, running_ts);
    shuffled_peak_joint_probability(k) = max(shuffled_likelihood(:));
end

%% Compute significance and determine if this cell is a place cell

pN = sum(shuffled_likelihood > likelihood,3)/numShuffles; %  pN, supra-threshold tests

significant_likelihood = likelihood;
significant_likelihood(pN > 0.01) = NaN; %Determine whether this cell fires significantly different from chance in these particular bins 

place_cell = ~isempty(significant_likelihood); 

is_place_cell(cell_id,1) = place_cell;

% if plotting;
% plotting_fig = figure;
% subplot(3,1,1)
% imagesc(X_bin_centers_vector,Y_bin_centers_vector, pN)
% daspect([1 1 1])
% title 'p-value map'
% colorbar
% subplot(3,1,2)
% imagesc(X_bin_centers_vector,Y_bin_centers_vector,significant_likelihood);
% daspect([1 1 1])
% title 'significant tuning map'
% colorbar
% end

%% Compute split half stability index 
binarized_trace_early = [];
binarized_trace_late = [];
interp_behav_vec_early = [];
interp_behav_vec_late = [];

binarized_trace = ms.Binary(:, cell_id); 
binarized_trace(running_ts == 0) = 0; %exclude periods in immobility
nBinarizedEvents = sum(binarized_trace); % find total number of binarized events in inclusion vector 
half_recording_idx = nBinarizedEvents/2;

for frame = 1:length(binarized_trace);;
     if sum(binarized_trace(1:frame,1)) < half_recording_idx;
     binarized_trace_early(end+1) = binarized_trace(frame);
     
     elseif sum(binarized_trace(1:frame,1)) >= half_recording_idx;
     binarized_trace_late(end+1) = binarized_trace(frame);
     end 
end

interp_behav_vec_early = interp_behav_vec(1:length(binarized_trace_early),:);
interp_behav_vec_late = interp_behav_vec((length(binarized_trace_early)+1):length(interp_behav_vec),:);

inclusion_vector_early = running_ts(1:length(binarized_trace_early),:); 
inclusion_vector_late = running_ts((length(binarized_trace_early)+1):length(interp_behav_vec),:); 

[MI_early, posterior_early, occupancy_map_early, prob_being_active_early, likelihood_early] = extract_2D_information(binarized_trace_early, interp_behav_vec_early, X_bin_vector, Y_bin_vector, inclusion_vector_early);
peak_joint_probability = max(likelihood_early(:));

[MI_late, posterior_late, occupancy_map_late, prob_being_active_late, likelihood_late] = extract_2D_information(binarized_trace_late, interp_behav_vec_late, X_bin_vector, Y_bin_vector, inclusion_vector_late);
peak_joint_probability = max(likelihood_late(:));

if params.smoothing
    kernel_size = [bin_size bin_size];
    occupancy_std = 2;

[Xgrid,Ygrid]=meshgrid(-kernel_size(1)/2: kernel_size(1)/2, -kernel_size(2)/2:kernel_size(2)/2);
    Rgrid=sqrt((Xgrid.^2+Ygrid.^2));
    kernel = pdf('Normal', Rgrid, 0, occupancy_std);
    kernel = kernel./sum(sum(kernel));
    
    likelihood_early = conv2(likelihood_early, kernel, 'same'); % smoothing
    likelihood_late = conv2(likelihood_late, kernel, 'same'); % smoothing
    smoothed_likelihood = conv2(likelihood, kernel, 'same'); % smoothing
end

%correlate the tuning maps 
place_field_stability = corr2(likelihood_early,likelihood_late);

if plotting;
plotting_fig = figure;
subplot(3,1,1)
imagesc(X_bin_centers_vector,Y_bin_centers_vector, smoothed_likelihood)
caxis([0 .15]) %.15 max used for paper
daspect([1 1 1])
title 'tuning map smoothed'
colorbar
subplot(3,1,2)
imagesc(X_bin_centers_vector,Y_bin_centers_vector, likelihood_early)
daspect([1 1 1])
title 'tuning map early'
colorbar
subplot(3,1,3)
imagesc(X_bin_centers_vector,Y_bin_centers_vector,likelihood_late);
daspect([1 1 1])
title 'tuning map late'
colorbar
title(['Stability: ' num2str(place_field_stability)])
end

stability(cell_id,1) = place_field_stability; 

%% Output the data
 cell_ids(cell_id) = cell_id;
 extended_calcium_analysis_output(cell_id).MI = MI;
 extended_calcium_analysis_output(cell_id).prob_being_active = prob_being_active;
 extended_calcium_analysis_output(cell_id).place_field_stability = place_field_stability;
 
 extended_calcium_analysis_output(cell_id).PDF = posterior;
 extended_calcium_analysis_output(cell_id).pN = pN;
 extended_calcium_analysis_output(cell_id).occupancy_map = occupancy_map;
 extended_calcium_analysis_output(cell_id).tuning_map = likelihood;
 extended_calcium_analysis_output(cell_id).stability = stability;

% save('extended_calcium_analysis_output.mat', 'extended_calcium_analysis_output');

end 


