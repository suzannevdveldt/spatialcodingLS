function [extended_calcium_analysis_output] = extended_calcium_analysis_tmaze(ms);

% Potentially run this by running select_folders_for_analysis.m
%%                             
 clear all; 
  load 'behav.mat';
  load 'ms.mat';
   load 'behavDLC.mat'; 
  %load('extended_calcium_analysis_output.mat')
 
poolobj = gcp('nocreate'); % If no pool, do not create new one.
if isempty(poolobj)
    poolsize = 0;
else
    poolsize = poolobj.NumWorkers
end
 
plotting = 0; 
save_figure = 0;
bin_size = 3;
min_speed_threshold = 0; % 2 cm.s-1
params.smoothing = true;
numShuffles = 1000;
centroid_analysis = 0; 
place_field.centroid = zeros(1,2);

if ~isfield(ms,'Binary');
    disp('Extracting binary information');
    ms = msExtractBinary(ms);
end

%% determine behav and ms timestamps and convert to seconds
behav_time = behav.time/1000; %convert to seconds 
behav_vec = behav.position;
[behav_time, idx]=unique(behav_time); %getting rid of duplicate timestamps
behav_vec = behav_vec(idx,:); %excluding the positions associated with these double timestamps

heading_vec = behavDLC.headDirection;
heading_vec = heading_vec(idx);
heading_vec = heading_vec';

ca_time = ms.time/1000; % convert to seconds 
[ca_time, idx]=unique(ca_time); %getting rid of duplicate timestamps

%% Binarize calcium trace
sampling_frequency = 30; % This data set has been sampled at 30 images per second

%% Interpolate behavior
% In most cases, behavior data has to be interpolated to match neural temporal
% activity assuming it has a lower sampling frequency than neural activity
interp_behav_vec(:,1) = interpolate_behavior(behav_vec(:,1), behav_time, ca_time); % in the X dimension
interp_behav_vec(:,2) = interpolate_behavior(behav_vec(:,2), behav_time, ca_time); % in the Y dimension
interp_behav_vec(end,:) = interp_behav_vec(end-1,:); % Make sure to interpolate or remove every NaN so that every timepoint has a corresponding behavioral state

 %Interpolate head direction
[interp_heading_vec] = interpolate_HD(heading_vec, behav_time, ca_time);
interp_heading_vec(end) = interp_heading_vec(end-1); % Make sure to interpolate or remove every NaN so that every timepoint has a corresponding behavioral state

%% Compute velocities
% In this example, we will ignore moments when the mouse is immobile
[velocity] = extract_velocity(interp_behav_vec, ca_time);

%% Compute acceleration
acceleration = diff(velocity);
acceleration(end+1) = 0;

%% Find periods of immobility
% This will be usefull later to exclude immobility periods from analysis
running_ts = velocity > min_speed_threshold;

%% Create an inclusion vector to isolate only a specific behavior
inclusion_vector = running_ts; % Only include periods of running

min_speed_threshold_velocity = 2.5; % 2 cm.s-1
running_ts = velocity > min_speed_threshold_velocity;
inclusion_vector_velocity = running_ts; % Only include periods of running


%% Compute occupancy and joint probabilities
% You can use 'min(behav_vec)' and 'max(behav_vec)'  to estimate the
% boundaries of the behavior vector (in our case, location in space)

min_X = 0;
min_Y = 0;
%max_X = 105; % Set these to the borders of the maze used. !!! Has to be the same for every mouse for consistent MI computation !!!
%max_Y = 55;
max_X = behav.trackLength; % Set these to the borders of the maze used. !!! Has to be the same for every mouse for consistent MI computation !!!
max_Y = behav.trackLength;

X_bin_vector = min_X:bin_size:max_X+bin_size;
X_bin_centers_vector = X_bin_vector + bin_size/2;
X_bin_centers_vector(end) = [];
Y_bin_vector = min_Y:bin_size:max_Y+bin_size;
Y_bin_centers_vector = Y_bin_vector + bin_size/2;
Y_bin_centers_vector(end) = [];

%% Velocity and Acceleration binning 

min_velocity = 2.5;
min_acceleration  = -2; % I am exluding some of the extrimities as these bins are likely empty
max_velocity = 30; 
max_acceleration = 2;

velocity_binsize = (max_velocity-min_velocity)/20;
acceleration_binsize = (-(min_acceleration-max_acceleration))/20; 

velocity_bin_vector = min_velocity:velocity_binsize:max_velocity;
velocity_bin_centers_vector = velocity_bin_vector + velocity_binsize/2;
velocity_bin_centers_vector(end) = [];

acceleration_bin_vector = min_acceleration:acceleration_binsize:max_acceleration;
acceleration_bin_centers_vector = acceleration_bin_vector + acceleration_binsize/2;
velocity_bin_centers_vector(end) = [];

%%  Binning polar space
bin_size = 9; %open field with 3 cm bins in 45 cm environment has 256 bins total 
bin_vector = 0:bin_size:360;

for cell_id = 1:ms.numNeurons %size(binarized_trace,2);
 
    cell_id
 ca_trace = ms.RawTraces(:,cell_id);
 ca_trace = ca_trace(idx); %excluding the positions associated with these double timestamps
 
 binarized_trace = ms.Binary(:, cell_id); 
 binarized_trace = binarized_trace(idx);

[spatial_MI, posterior, occupancy_map, prob_being_active, likelihood] = extract_2D_information(binarized_trace, interp_behav_vec, X_bin_vector, Y_bin_vector, inclusion_vector);
[~, velocity_MI, ~, ~, ~, ~] = extract_1D_information(binarized_trace, velocity, velocity_bin_vector, inclusion_vector_velocity);
[~, acceleration_MI, ~, ~, ~, ~] = extract_1D_information(binarized_trace, acceleration, acceleration_bin_vector, inclusion_vector_velocity);
[~, MI_heading, posterior_heading, occupancy_vector_heading, ~, likelihood_heading ] = extract_1D_information(binarized_trace, interp_heading_vec, bin_vector, inclusion_vector);

for shuffle_i = 1:numShuffles;
    random_ts = ceil(rand*length(ca_time));
    shuffled_binarized = circshift(binarized_trace,random_ts);
 
    [shuffled_spatial_MI(shuffle_i), ~, ~, ~, ~] = extract_2D_information(shuffled_binarized, interp_behav_vec,  X_bin_vector, Y_bin_vector, inclusion_vector);
    [~, shuffled_velocity_MI(shuffle_i), ~, ~, ~, ~] = extract_1D_information(shuffled_binarized, velocity,  velocity_bin_vector, inclusion_vector_velocity);
    [~, shuffled_acceleration_MI(shuffle_i), ~, ~, ~, ~] = extract_1D_information(shuffled_binarized, acceleration,  acceleration_bin_vector, inclusion_vector_velocity);
    [~, shuffled_heading_MI(shuffle_i), ~, ~ , ~, ~] = extract_1D_information(shuffled_binarized, interp_heading_vec, bin_vector, inclusion_vector);

end 

spatial_MI_pvalue = sum(shuffled_spatial_MI>spatial_MI)/numShuffles;
spatial_zMI = (spatial_MI-mean(shuffled_spatial_MI))/std(shuffled_spatial_MI);

velocity_MI_pvalue = sum(shuffled_velocity_MI>velocity_MI)/numShuffles;
velocity_zMI = (velocity_MI-mean(shuffled_velocity_MI))/std(shuffled_velocity_MI);

acceleration_MI_pvalue = sum(shuffled_acceleration_MI>acceleration_MI)/numShuffles;
acceleration_zMI = (acceleration_MI-mean(shuffled_acceleration_MI))/std(shuffled_acceleration_MI);

heading_MI_pvalue = sum(shuffled_heading_MI>MI_heading)/numShuffles;
heading_zMI = (MI_heading-mean(shuffled_heading_MI))/std(shuffled_heading_MI);

%% Centroid analysis (to update) 
    
%  if centroid_analysis;
%  [mean_dispersion, centroid] = computeCentroidDispersion(binarized_trace, interp_behav_vec, inclusion_vector);
%  end 

% for center of the place field take the bin with peak firing rate
[row, col] = find(ismember(likelihood, max(likelihood(:))));

centroid = [row, col];

%% Compute split half stability index 
binarized_trace_early = [];
binarized_trace_late = [];
interp_behav_vec_early = [];
interp_behav_vec_late = [];

%binarized_trace = ms.Binary(:, cell_id); 
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

[~, ~, ~, ~, likelihood_early] = extract_2D_information(binarized_trace_early, interp_behav_vec_early, X_bin_vector, Y_bin_vector, inclusion_vector_early);

[~, ~, ~, ~, likelihood_late] = extract_2D_information(binarized_trace_late, interp_behav_vec_late, X_bin_vector, Y_bin_vector, inclusion_vector_late);

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

%% stability for head direction

binarized_trace_early = [];
binarized_trace_late = [];

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

inclusion_vector_early = running_ts(1:length(binarized_trace_early),:); 
inclusion_vector_late = running_ts((length(binarized_trace_early)+1):length(interp_heading_vec),:); 

interp_heading_early = interp_heading_vec(1:length(binarized_trace_early),:);
interp_heading_late = interp_heading_vec((length(binarized_trace_early)+1):length(interp_heading_vec),:);

[~, ~, ~, ~, ~, likelihood_heading_early] = extract_1D_information(binarized_trace_early, interp_heading_early, bin_vector, inclusion_vector_early);
[~, ~, ~, ~, ~, likelihood_heading_late] = extract_1D_information(binarized_trace_late, interp_heading_late, bin_vector, inclusion_vector_late);

if params.smoothing
    kernel_size = [bin_size bin_size];
    occupancy_std = 2;

[Xgrid,Ygrid]=meshgrid(-kernel_size(1)/2: kernel_size(1)/2, -kernel_size(2)/2:kernel_size(2)/2);
    Rgrid=sqrt((Xgrid.^2+Ygrid.^2));
    kernel = pdf('Normal', Rgrid, 0, occupancy_std);
    kernel = kernel./sum(sum(kernel));
end

likelihood_heading_early = conv2(likelihood_heading_early, kernel, 'same'); % smoothing
likelihood_heading_late = conv2(likelihood_heading_late, kernel, 'same'); % smoothing

heading_stability = corr2(likelihood_heading_early, likelihood_heading_late);

%% Output the data
 extended_calcium_analysis_output(cell_id).cell_id = cell_id;
 extended_calcium_analysis_output(cell_id).MI = spatial_MI;
 extended_calcium_analysis_output(cell_id).prob_being_active = prob_being_active;
 extended_calcium_analysis_output(cell_id).place_field_stability = place_field_stability;
 extended_calcium_analysis_output(cell_id).spatial_MI_pvalue = spatial_MI_pvalue;
 extended_calcium_analysis_output(cell_id).spatial_zMI = spatial_zMI;
 
 extended_calcium_analysis_output(cell_id).PDF = posterior;
 extended_calcium_analysis_output(cell_id).occupancy_map = occupancy_map;
 extended_calcium_analysis_output(cell_id).tuning_map = likelihood;
 
 split_half_stability(cell_id,1) = place_field_stability;
 tuning_map(cell_id).likelihood = likelihood;
  
 %extended_calcium_analysis_output(cell_id).mean_dispersion = mean_dispersion;
 extended_calcium_analysis_output(cell_id).centroid = centroid; 
  
 %place_field.mean_dispersion(cell_id,1) = mean_dispersion;
% place_field.centroid(cell_id,:) = centroid; 

  extended_calcium_analysis_output(cell_id).velocity_MI = velocity_MI;
  extended_calcium_analysis_output(cell_id).velocity_MI_pvalue = velocity_MI_pvalue;
  extended_calcium_analysis_output(cell_id).velocity_zMI = velocity_zMI;
  
  extended_calcium_analysis_output(cell_id).acceleration_MI = acceleration_MI;
  extended_calcium_analysis_output(cell_id).acceleration_MI_pvalue = acceleration_MI_pvalue;
  extended_calcium_analysis_output(cell_id).acceleration_zMI = acceleration_zMI;
  
  extended_calcium_analysis_output(cell_id).MI_HeadDirection = MI_heading;
extended_calcium_analysis_output(cell_id).headdirection_stability = heading_stability;
 extended_calcium_analysis_output(cell_id).headdirection_MI_pvalue = heading_MI_pvalue;
 extended_calcium_analysis_output(cell_id).headdirection_zMI = heading_zMI;
 
extended_calcium_analysis_output(cell_id).PDF_HeadDirection = posterior_heading;
extended_calcium_analysis_output(cell_id).occupancy_vector_HeadDirection = occupancy_vector_heading;
extended_calcium_analysis_output(cell_id).tuning_curve_HeadDirection = likelihood_heading;

  split_half_stability(cell_id,1) = heading_stability;
  tuning_vector(cell_id).likelihood = likelihood_heading;
 

end
save('extended_calcium_analysis_output.mat', 'extended_calcium_analysis_output');

end 


