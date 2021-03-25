function [extended_calcium_analysis_output] = extended_calcium_analysis_LT(ms,behav);

% Potentially run this by running select_folders_for_analysis.m
%%                             
 %clear all; 
 load 'behav.mat';
 load 'ms.mat';
 extended_calcium_analysis_output = [];
 
  poolobj = gcp('nocreate'); % If no pool, do not create new one.
if isempty(poolobj)
    poolsize = 0;
else
    poolsize = poolobj.NumWorkers
end

plotting = 0; 
save_figure = 0;
min_speed_threshold = 5; % 2 cm.s-1
params.smoothing = true;
numShuffles = 1000;
bin_size = 3; 
analyze_running = 1; 

if ~isfield(ms,'Binary');
    disp('Extracting binary information');
    ms = msExtractBinary(ms);
end

%% determine behav and ms timestamps and convert to seconds
behav_time = behav.time/1000; %convert to seconds 
behav_vec = behav.position(:,1);
[behav_time, idx]=unique(behav_time); %getting rid of duplicate timestamps
behav_vec = behav_vec(idx,:); %excluding the positions associated with these double timestamps
ca_time = ms.time/1000; % convert to seconds 
[ca_time, idx]=unique(ca_time); %getting rid of duplicate timestamps
sampling_frequency = 30; % This data set has been sampled at 30 images per second

%% Interpolate behavior
% In most cases, behavior data has to be interpolated to match neural temporal
% activity assuming it has a lower sampling frequency than neural activity
[interp_behav_vec] = interpolate_behavior(behav_vec, behav_time, ca_time);
interp_behav_vec(end) = interp_behav_vec(end-1); % Make sure to interpolate or remove every NaN so that every timepoint has a corresponding behavioral state

%% Compute velocities
% In this example, we will ignore moments when the mouse is immobile
[velocity] = extract_velocity(interp_behav_vec, ca_time);
acceleration = diff(velocity);
acceleration(end+1) = 0;

%% Find periods of immobility
% This will be usefull later to exclude immobility periods from analysis
running_ts = velocity > min_speed_threshold;

%% Compute occupancy and joint probabilities
% You can use 'min(behav_vec)' and 'max(behav_vec)'  to estimate the
% boundaries of the behavior vector (in our case, location in space)
min_X = 0;
max_X = 100; % Set these to the borders of the maze used. !!! Has to be the same for every mouse for consistent MI computation !!!

bin_vector = min_X:bin_size:max_X;
bin_centers_vector = bin_vector + bin_size/2;
bin_centers_vector(end) = [];

z_delta_x = diff(interp_behav_vec);
z_delta_x(end+1,1) = 0;
z_delta_x(isnan(z_delta_x)) = 0;
z_delta_x = zscore(z_delta_x);

% Right trajectories only
direction_indices = z_delta_x > 0.2;
inclusion_vector_right = direction_indices; % Only include right runs
inclusion_vector_right(running_ts == 0) = 0; % Only include periods of running

% Left trajectories only
direction_indices = z_delta_x < -0.2;
inclusion_vector_left = direction_indices; % Only include left runs
inclusion_vector_left(running_ts == 0) = 0; % Only include periods of running

% All trajectories
inclusion_vector = running_ts;

if analyze_running; 
included_velocity = velocity(inclusion_vector); 
included_acceleration = acceleration(inclusion_vector);

min_velocity = min(included_velocity);
min_acceleration  = min(included_acceleration)+.4; % I am exluding some of the extrimities as these bins are likely empty
max_velocity = max(included_velocity)-2; 
max_acceleration = max(included_acceleration)-.3;

velocity_binsize = (max_velocity-min_velocity)/33;
acceleration_binsize = (-(min_acceleration-max_acceleration))/33; 

velocity_bin_vector = min_velocity:velocity_binsize:max_velocity+velocity_binsize;
velocity_bin_centers_vector = velocity_bin_vector + velocity_binsize/2;
velocity_bin_vector(end) = [];

acceleration_bin_vector = min_acceleration:acceleration_binsize:max_acceleration+acceleration_binsize;
acceleration_bin_centers_vector = acceleration_bin_vector + acceleration_binsize/2;
acceleration_bin_vector(end) = [];
end

for cell_id = 1:ms.numNeurons
   
%% Binarize calcium trace
binarized_trace = ms.Binary(:, cell_id);
binarized_trace = binarized_trace(idx);

[KLD, MI, posterior, occupancy_vector, prob_being_active, likelihood ] = extract_1D_information(binarized_trace, interp_behav_vec, bin_vector, inclusion_vector);
tuning_curve_peak = max(likelihood);

[KLD_right, MI_right, posterior_right, occupancy_vector_right, ~, likelihood_right ] = extract_1D_information(binarized_trace, interp_behav_vec, bin_vector, inclusion_vector_right);
tuning_curve_peak_right = max(likelihood_right);

[KLD_left, MI_left, posterior_left, occupancy_vector_left, ~, likelihood_left ] = extract_1D_information(binarized_trace, interp_behav_vec, bin_vector, inclusion_vector_left);
tuning_curve_peak_left = max(likelihood_left);

if analyze_running;
[KLD_velocity, MI_velocity, posterior_velocity, occupancy_vector_velocity, ~, likelihood_velocity ] = extract_1D_information(binarized_trace, velocity, velocity_bin_vector, inclusion_vector);
[KLD_acceleration, MI_acceleration, posterior_acceleration, occupancy_vector_acceleration, ~, likelihood_acceleration] = extract_1D_information(binarized_trace, acceleration, acceleration_bin_vector, inclusion_vector);
end 

if plotting;
plotting_fig = figure;
subplot(3,1,1)
plot(bin_centers_vector,likelihood,'Color', [0 0.1 0.8])
title 'Tuning curve'
xlabel 'Location on the track (cm)'
ylabel 'Probability of firing in location'
subplot(3,1,2)
plot(bin_centers_vector,occupancy_vector,'Color', [0.1 0.8 0.1])
title 'Occupancy'
xlabel 'Location on the track (cm)'
ylabel 'Relative occupancy'
subplot(3,1,3)
plot(bin_centers_vector,posterior,'Color', [0.8 0.2 0])
title 'Posterior probability density function'
xlabel 'Location on the track (cm)'
ylabel 'Probability'
end

if save_figure;
    saveas(plotting_fig, sprintf('%d.jpg',cell_id));
end

%% Shuffle data
shuffled_tuning_curve = zeros(length(bin_centers_vector), numShuffles);

parfor k = 1:numShuffles
    random_ts = ceil(rand*length(ca_time));
    shuffled_binarized = zeros(length(binarized_trace),1);
    
    % Permute the trace
    shuffled_binarized(1:random_ts) = binarized_trace(end-random_ts+1:end);
    shuffled_binarized(random_ts+1:end) = binarized_trace(1:end-random_ts);
   
    % Compute tuning curve
    [~, shuffled_MI(:,k), shuffled_posterior(:,k), ~, ~, shuffled_likelihood(:,k)] = extract_1D_information(shuffled_binarized, interp_behav_vec, bin_vector, inclusion_vector);
    %shuffled_tuning_curve_peak(k) = max(shuffled_likelihood(:,k));    
    
    [~, shuffled_MI_left(:,k), shuffled_posterior_left(:,k), ~, ~, shuffled_likelihood_left(:,k)] = extract_1D_information(shuffled_binarized, interp_behav_vec, bin_vector, inclusion_vector_left);
    %shuffled_tuning_curve_peak_left(k) = max(shuffled_likelihood_left(:,k));    
    [~, shuffled_MI_right(:,k), shuffled_posterior_right(:,k), ~, ~, shuffled_likelihood_right(:,k)] = extract_1D_information(shuffled_binarized, interp_behav_vec, bin_vector, inclusion_vector_right);

    
    if analyze_running;
    [~, shuffled_MI_velocity(:,k), shuffled_posterior_velocity(:,k), ~ , ~, shuffled_likelihood_velocity(:,k)] = extract_1D_information(shuffled_binarized, velocity, velocity_bin_vector, inclusion_vector);
    %shuffled_peak_joint_probability_velocity(k) = max(shuffled_likelihood_velocity(:));  
    
    [~, shuffled_MI_acceleration(:,k), shuffled_posterior_acceleration(:,k), ~ , ~, shuffled_likelihood_acceleration(:,k)] = extract_1D_information(shuffled_binarized, acceleration, acceleration_bin_vector, inclusion_vector);
    %shuffled_peak_joint_probability_acceleration(k) = max(shuffled_likelihood_acceleration(:));  
    end 
end 

spatial_MI_pvalue = sum(shuffled_MI>MI)/numShuffles;
spatial_zMI = (MI-mean(shuffled_MI))/std(shuffled_MI);

spatial_MI_pvalue_left = sum(shuffled_MI_left>MI_left)/numShuffles;
spatial_zMI_left = (MI_left-mean(shuffled_MI_left))/std(shuffled_MI_left);

spatial_MI_pvalue_right = sum(shuffled_MI_right>MI_right)/numShuffles;
spatial_zMI_right = (MI_right-mean(shuffled_MI_right))/std(shuffled_MI_right);

%% Compute significance and determine if this cell is a place cell

pN = sum(shuffled_likelihood > likelihood,2)/numShuffles; %  pN, supra-threshold tests

pN_left = sum(shuffled_likelihood_left > likelihood_left,2)/numShuffles; %  pN, supra-threshold tests
pN_right = sum(shuffled_likelihood_right > likelihood_right,2)/numShuffles; %  pN, supra-threshold tests

significant_likelihood = likelihood;
significant_likelihood(pN > 0.01) = 0;

if analyze_running;
significant_likelihood_velocity = [];
significant_likelihood_acceleration = []; 
pN_velocity = [];
pN_acceleration = []; 

pN_velocity = sum(shuffled_likelihood_velocity > likelihood_velocity,2)/numShuffles; %  pN, supra-threshold tests
significant_likelihood_velocity = likelihood_velocity;
significant_likelihood_velocity(pN_velocity > 0.01) = NaN; %Determine whether this cell fires significantly different from chance in these particular bins 
significant_likelihood_velocity(end+1,1) = 0;

pN_acceleration = sum(shuffled_likelihood_acceleration > likelihood_acceleration,2)/numShuffles; %  pN, supra-threshold tests
significant_likelihood_acceleration = likelihood_acceleration;
significant_likelihood_acceleration(pN_acceleration > 0.01) = NaN; %Determine whether this cell fires significantly different from chance in these particular bins 
significant_likelihood_acceleration(end+1,1) = 0;
end


%% Compute split half stability index 
binarized_trace_early = [];
binarized_trace_late = [];
interp_behav_vec_early = [];
interp_behav_vec_late = [];

binarized_trace = ms.Binary(:, cell_id); 
binarized_trace(running_ts == 0) = 0; % only include periods of running 
nBinarizedEvents = sum(binarized_trace); % find total number of binarized events in inclusion vector 
half_recording_idx = nBinarizedEvents/2;

for frame = 1:length(binarized_trace);
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

inclusion_vector_early_right = inclusion_vector_right(1:length(binarized_trace_early),:);
inclusion_vector_late_right = inclusion_vector_right((length(binarized_trace_early)+1):length(interp_behav_vec),:);

inclusion_vector_early_left = inclusion_vector_left(1:length(binarized_trace_early),:);
inclusion_vector_late_left = inclusion_vector_left((length(binarized_trace_early)+1):length(interp_behav_vec),:);

[~, ~, ~, ~, ~, likelihood_early] = extract_1D_information(binarized_trace_early, interp_behav_vec_early, bin_vector, inclusion_vector_early);
[~, ~, ~, ~, ~, likelihood_early_right] = extract_1D_information(binarized_trace_early, interp_behav_vec_early, bin_vector, inclusion_vector_early_right);
[~, ~, ~, ~, ~, likelihood_early_left] = extract_1D_information(binarized_trace_early, interp_behav_vec_early, bin_vector, inclusion_vector_early_left);

[~, ~, ~, ~, ~, likelihood_late] = extract_1D_information(binarized_trace_late, interp_behav_vec_late, bin_vector, inclusion_vector_late);
[~, ~, ~, ~, ~, likelihood_late_right] = extract_1D_information(binarized_trace_late, interp_behav_vec_late, bin_vector, inclusion_vector_late_right);
[~, ~, ~, ~, ~, likelihood_late_left] = extract_1D_information(binarized_trace_late, interp_behav_vec_late, bin_vector, inclusion_vector_late_left);

%% left and right indices
if params.smoothing
    kernel_size = [bin_size bin_size];
    occupancy_std = 2;

[Xgrid,Ygrid]=meshgrid(-kernel_size(1)/2: kernel_size(1)/2, -kernel_size(2)/2:kernel_size(2)/2);
    Rgrid=sqrt((Xgrid.^2+Ygrid.^2));
    kernel = pdf('Normal', Rgrid, 0, occupancy_std);
    kernel = kernel./sum(sum(kernel));
    
    likelihood_early = conv2(likelihood_early, kernel, 'same'); % smoothing
    likelihood_late = conv2(likelihood_late, kernel, 'same'); % smoothing
    
    likelihood_early_right = conv2(likelihood_early_right, kernel, 'same'); % smoothing
    likelihood_late_right = conv2(likelihood_late_right, kernel, 'same'); % smoothing
    
    likelihood_early_left = conv2(likelihood_early_left, kernel, 'same'); % smoothing
    likelihood_late_left = conv2(likelihood_late_left, kernel, 'same'); % smoothing
    
    smoothed_likelihood = conv2(likelihood, kernel, 'same'); % smoothing
end

%correlate the tuning maps 
place_field_stability = corr2(likelihood_early,likelihood_late);
right_place_field_stability = corr2(likelihood_early_right, likelihood_late_right);
left_place_field_stability = corr2(likelihood_early_left, likelihood_late_left);

%% Output the data
% All trajectories 

extended_calcium_analysis_output(cell_id).KLD = KLD;
extended_calcium_analysis_output(cell_id).KLD_left = KLD_left;
extended_calcium_analysis_output(cell_id).KLD_right = KLD_right;

extended_calcium_analysis_output(cell_id).MI = MI;
extended_calcium_analysis_output(cell_id).MI_left = MI_left;
extended_calcium_analysis_output(cell_id).MI_right = MI_right;

extended_calcium_analysis_output(cell_id).MI_pvalue = spatial_MI_pvalue;
extended_calcium_analysis_output(cell_id).MI_pvalue_left = spatial_MI_pvalue_left;
extended_calcium_analysis_output(cell_id).MI_pvalue_right = spatial_MI_pvalue_right;

extended_calcium_analysis_output(cell_id).zMI = spatial_zMI;
extended_calcium_analysis_output(cell_id).zMI_left = spatial_zMI_left;
extended_calcium_analysis_output(cell_id).zMI_right = spatial_zMI_right;

extended_calcium_analysis_output(cell_id).prob_being_active = prob_being_active;

extended_calcium_analysis_output(cell_id).place_field_stability = place_field_stability;
extended_calcium_analysis_output(cell_id).right_place_field_stability = right_place_field_stability;
extended_calcium_analysis_output(cell_id).left_place_field_stability = left_place_field_stability;

extended_calcium_analysis_output(cell_id).likelihood = likelihood;
extended_calcium_analysis_output(cell_id).likelihood_left = likelihood_left;
extended_calcium_analysis_output(cell_id).likelihood_right = likelihood_right;

extended_calcium_analysis_output(cell_id).PDF = posterior;
extended_calcium_analysis_output(cell_id).PDF_left = posterior_left;
extended_calcium_analysis_output(cell_id).PDF_right = posterior_right;

extended_calcium_analysis_output(cell_id).pN = pN;
extended_calcium_analysis_output(cell_id).pN_left = pN_left;
extended_calcium_analysis_output(cell_id).pN_right = pN_right;

extended_calcium_analysis_output(cell_id).pN_velocity = pN_velocity;
extended_calcium_analysis_output(cell_id).pN_acceleration = pN_acceleration;

extended_calcium_analysis_output(cell_id).occupancy_vector = occupancy_vector;
extended_calcium_analysis_output(cell_id).occupancy_vector_left = occupancy_vector_left;
extended_calcium_analysis_output(cell_id).occupancy_vector_right = occupancy_vector_right;

 if analyze_running;
 extended_calcium_analysis_output(cell_id).MI_velocity = MI_velocity;
 extended_calcium_analysis_output(cell_id).MI_acceleration = MI_acceleration;
 
 extended_calcium_analysis_output(cell_id).KLD_velocity = KLD_velocity;
 extended_calcium_analysis_output(cell_id).KLD_acceleration = KLD_acceleration;
 
 extended_calcium_analysis_output(cell_id).PDF_velocity = posterior_velocity;
 extended_calcium_analysis_output(cell_id).PDF_acceleration = posterior_acceleration;

 extended_calcium_analysis_output(cell_id).pN_velocity = pN_velocity;
  extended_calcium_analysis_output(cell_id).pN_acceleration = pN_acceleration;

 extended_calcium_analysis_output(cell_id).occupancy_map_velocity = occupancy_vector_velocity;
  extended_calcium_analysis_output(cell_id).occupancy_map_acceleration = occupancy_vector_acceleration;

 extended_calcium_analysis_output(cell_id).tuning_curve_velocity = likelihood_velocity;
 extended_calcium_analysis_output(cell_id).tuning_curve_acceleration = likelihood_acceleration;
 end 

%pause;
end;

save('extended_calcium_analysis_output_pvalue.mat', 'extended_calcium_analysis_output');

end 


