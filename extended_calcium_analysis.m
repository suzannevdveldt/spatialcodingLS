function [extended_calcium_analysis_output] = extended_calcium_analysis;

% Potentially run this by running select_folders_for_analysis.m
%%                             
 clear all; 
  load 'behav.mat';
  load 'ms.mat';
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
min_speed_threshold = 5; % 2 cm.s-1
params.smoothing = true;
numShuffles = 1;
analyze_heading = 0; %this is analyzed in a separate head DIRECTION script (better) 
analyze_running = 0; %this is analyzed in the bootstrapping script
centroid_analysis = 1; 
place_field.centroid = zeros(1,2);

if ~isfield(ms,'Binary');
    disp('Extracting binary information');
    ms = msExtractBinary(ms);
end

%% determine behav and ms timestamps and convert to seconds
behav_time = behav.time/1000; %convert to seconds 
behav_vec = behav.position;
%behav_vec = behav.positionRescale;
[behav_time, idx]=unique(behav_time); %getting rid of duplicate timestamps
behav_vec = behav_vec(idx,:); %excluding the positions associated with these double timestamps
ca_time = ms.time/1000; % convert to seconds 
[ca_time, idx]=unique(ca_time); %getting rid of duplicate timestamps
% %behav_time = behav_time';


%% Binarize calcium trace
sampling_frequency = 30; % This data set has been sampled at 30 images per second

%  z_threshold = 2; % A 2 standard-deviation threshold is usually optimal to differentiate calcium ativity from background noise
% 
%  binarized_data = zeros(size(ca_data));
%   for cell_i = 1:size(ca_data,2)
%       binarized_data(:,cell_i) = extract_binary(ca_data(:,cell_i),sampling_frequency, z_threshold);
%  end

%% Interpolate behavior
% In most cases, behavior data has to be interpolated to match neural temporal
% activity assuming it has a lower sampling frequency than neural activity
interp_behav_vec(:,1) = interpolate_behavior(behav_vec(:,1), behav_time, ca_time); % in the X dimension
interp_behav_vec(:,2) = interpolate_behavior(behav_vec(:,2), behav_time, ca_time); % in the Y dimension
interp_behav_vec(end,:) = interp_behav_vec(end-1,:); % Make sure to interpolate or remove every NaN so that every timepoint has a corresponding behavioral state

%heading = mod(atan2(diff(interp_behav_vec(:,1)), diff(interp_behav_vec(:,2)))*180/pi, 360); % heading for each timestamp in degrees 
%heading(end+1)=NaN;
%% Compute velocities
% In this example, we will ignore moments when the mouse is immobile
[velocity] = extract_velocity(interp_behav_vec, ca_time);

%% Compute accelaration
acceleration = diff(velocity);
acceleration(end+1) = 0;
%% Find periods of immobility
% This will be usefull later to exclude immobility periods from analysis
running_ts = velocity > min_speed_threshold;

%% Create an inclusion vector to isolate only a specific behavior
inclusion_vector = running_ts; % Only include periods of running
included_velocity = velocity(inclusion_vector); 
included_acceleration = acceleration(inclusion_vector);

%% Compute occupancy and joint probabilities
% You can use 'min(behav_vec)' and 'max(behav_vec)'  to estimate the
% boundaries of the behavior vector (in our case, location in space)

min_X = 0;
min_Y = 0;
max_X = 105; % Set these to the borders of the maze used. !!! Has to be the same for every mouse for consistent MI computation !!!
max_Y = 55;

X_bin_vector = min_X:bin_size:max_X+bin_size;
X_bin_centers_vector = X_bin_vector + bin_size/2;
X_bin_centers_vector(end) = [];
Y_bin_vector = min_Y:bin_size:max_Y+bin_size;
Y_bin_centers_vector = Y_bin_vector + bin_size/2;
Y_bin_centers_vector(end) = [];

for cell_id = 1:ms.numNeurons %size(binarized_trace,2);
 
    cell_id
 ca_trace = ms.RawTraces(:,cell_id);
 
 ca_trace = ca_trace(idx); %excluding the positions associated with these double timestamps
 
 binarized_trace = ms.Binary(:, cell_id); 
 binarized_trace = binarized_trace(idx);

%binarized_trace = binarized_data(:,cell_id);

[MI, posterior, occupancy_map, prob_being_active, likelihood] = extract_2D_information(binarized_trace, interp_behav_vec, X_bin_vector, Y_bin_vector, inclusion_vector);
peak_joint_probability = max(likelihood(:));

if analyze_heading;
[MI_heading, posterior_heading, occupancy_vector_heading, prob_being_active, likelihood_heading ] = extract_1D_information(binarized_trace, heading, heading_bin_vector, inclusion_vector);
peak_joint_probability = max(likelihood_heading(:));
end

if analyze_running; 
  [MI_velocity, posterior_velocity, occupancy_vector_velocity, ~, likelihood_velocity ] = extract_1D_information(binarized_trace, velocity, velocity_bin_vector, inclusion_vector);
peak_joint_probability_velocity = max(likelihood_velocity(:));
[KLD_acceleration, MI_acceleration, posterior_acceleration, occupancy_vector_acceleration, ~, likelihood_acceleration ] = extract_1D_information(binarized_trace, acceleration, acceleration_bin_vector, inclusion_vector);
peak_joint_probability_acceleration = max(likelihood_acceleration(:));
end 

if plotting
plotting_fig = figure;
subplot(3,1,1)
imagesc(X_bin_centers_vector,Y_bin_centers_vector,likelihood)
daspect([1 1 1])
title 'Probability of firing in location'
xlabel 'Position (cm)'
ylabel 'Position (cm)'
colorbar
subplot(3,1,2)
imagesc(X_bin_centers_vector,Y_bin_centers_vector,occupancy_map)
daspect([1 1 1])
title 'Relative occupancy'
xlabel 'Position (cm)'
ylabel 'Position (cm)'
colorbar
subplot(3,1,3)
surf(X_bin_centers_vector,Y_bin_centers_vector,posterior)
daspect([1 1 0.01])
title 'Posterior probability density function'
xlabel 'Position (cm)'
ylabel 'Position (cm)'
colorbar
end

if save_figure;
    saveas(plotting_fig, sprintf('%d.jpg',cell_id));
end

%% Shuffle data
shuffled_data = zeros(length(Y_bin_centers_vector), length(X_bin_centers_vector), numShuffles);
shuffled_peak_joint_probability = [];

if analyze_heading;
shuffled_likelihood_heading = zeros(length(likelihood_heading), numShuffles);
shuffled_MI_heading = zeros(length(MI_heading), numShuffles);
shuffled_posterior_heading = zeros(length(posterior_heading), numShuffles);
end

if analyze_running;
    shuffled_likelihood_velocity = zeros(length(likelihood_velocity), numShuffles);
shuffled_MI_velocity = zeros(length(MI_velocity), numShuffles);
shuffled_posterior_velocity = zeros(length(posterior_velocity), numShuffles);

shuffled_likelihood_acceleration = zeros(length(likelihood_acceleration), numShuffles);
shuffled_MI_acceleration = zeros(length(MI_acceleration), numShuffles);
shuffled_posterior_acceleration = zeros(length(posterior_acceleration), numShuffles);
end

for k = 1:numShuffles;
    random_ts = ceil(rand*length(binarized_trace));
    shuffled_binarized = zeros(length(binarized_trace),1);
    
    % Permute the trace
    shuffled_binarized(1:random_ts) = binarized_trace(end-random_ts+1:end);
    shuffled_binarized(random_ts+1:end) = binarized_trace(1:end-random_ts);
    
    [shuffled_MI(k), shuffled_posterior(:,:,k), ~ , ~, shuffled_likelihood(:,:,k)] = extract_2D_information(shuffled_binarized, interp_behav_vec, X_bin_vector, Y_bin_vector, running_ts);
    %shuffled_peak_joint_probability(k) = max(shuffled_likelihood(:));
   if analyze_heading;
   [shuffled_MI_heading(:,k), shuffled_posterior_heading(:,k), ~ , ~, shuffled_likelihood_heading(:,k)] = extract_1D_information(shuffled_binarized, heading, heading_bin_vector, inclusion_vector);
    shuffled_peak_joint_probability(k) = max(shuffled_likelihood_heading(:)); 
    end 
    
    if analyze_running;
   [shuffled_MI_velocity(:,k), shuffled_posterior_velocity(:,k), ~ , ~, shuffled_likelihood_velocity(:,k)] = extract_1D_information(shuffled_binarized, velocity, velocity_bin_vector, inclusion_vector);
    shuffled_peak_joint_probability_velocity(k) = max(shuffled_likelihood_velocity(:));  
    
   [shuffled_MI_acceleration(:,k), shuffled_posterior_acceleration(:,k), ~ , ~, shuffled_likelihood_acceleration(:,k)] = extract_1D_information(shuffled_binarized, acceleration, acceleration_bin_vector, inclusion_vector);
    shuffled_peak_joint_probability_acceleration(k) = max(shuffled_likelihood_acceleration(:));  

    end 

end

%% Compute significance and determine if this cell is a place cell

pN = sum(shuffled_likelihood > likelihood,3)/numShuffles; %  pN, supra-threshold tests
significant_likelihood = likelihood;
significant_likelihood(pN > 0.01) = NaN; %Determine whether this cell fires significantly different from chance in these particular bins 

if analyze_heading;
pN_heading = sum(shuffled_likelihood_heading > likelihood_heading,2)/numShuffles; %  pN, supra-threshold tests
significant_likelihood_heading = likelihood_heading;
significant_likelihood_heading(pN_heading > 0.01) = 0;
end

if analyze_running;
significant_likelihood_velocity = [];
significant_likelihood_acceleration = []; 
pN_velocity = [];
pN_acceleration = []; 

pN_velocity = sum(shuffled_likelihood_velocity > likelihood_velocity,2)/numShuffles; %  pN, supra-threshold tests
significant_likelihood_velocity = likelihood_velocity;
significant_likelihood_velocity(pN_velocity > 0.01) = NaN; %Determine whether this cell fires significantly different from chance in these particular bins 
%pN_velocity(end+1) = 0;
significant_likelihood_velocity(end+1,1) = 0;

pN_acceleration = sum(shuffled_likelihood_acceleration > likelihood_acceleration,2)/numShuffles; %  pN, supra-threshold tests
significant_likelihood_acceleration = likelihood_acceleration;
significant_likelihood_acceleration(pN_acceleration > 0.01) = NaN; %Determine whether this cell fires significantly different from chance in these particular bins 
%pN_acceleration(end+1) = 0;
significant_likelihood_acceleration(end+1,1) = 0;
end

 if centroid_analysis;
 [mean_dispersion, centroid] = computeCentroidDispersion(binarized_trace, interp_behav_vec, inclusion_vector);
 end 

if plotting;
plotting_fig = figure;
subplot(3,1,1)
imagesc(X_bin_centers_vector,Y_bin_centers_vector, pN)
daspect([1 1 1])
title 'p-value map'
colorbar
subplot(3,1,2)
imagesc(X_bin_centers_vector,Y_bin_centers_vector,significant_likelihood);
daspect([1 1 1])
title 'significant tuning map'
colorbar
end

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

if analyze_heading;
heading_early = heading(1:length(binarized_trace_early),:);
heading_late = heading((length(binarized_trace_early)+1):length(heading),:);

[~, ~, ~, ~, ~, likelihood_heading_early] = extract_1D_information(binarized_trace_early, heading_early, heading_bin_vector, inclusion_vector_early);
[~, ~, ~, ~, ~, likelihood_heading_late] = extract_1D_information(binarized_trace_late, heading_late, heading_bin_vector, inclusion_vector_late);
end

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
    
    if analyze_heading;
        likelihood_heading_early = conv2(likelihood_heading_early, kernel, 'same'); % smoothing
        likelihood_heading_late = conv2(likelihood_heading_late, kernel, 'same'); % smoothing
    end

end

%correlate the tuning maps 
place_field_stability = corr2(likelihood_early,likelihood_late);

% if analyze_heading;
%     heading_stability = corr2(likelihood_heading_early, likelihood_heading_late);
% end

if plotting;
plotting_fig = figure;
subplot(3,1,1)
imagesc(X_bin_centers_vector,Y_bin_centers_vector, smoothed_likelihood)
caxis([0 .15])
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

%% Output the data
 extended_calcium_analysis_output(cell_id).MI = MI;
 extended_calcium_analysis_output(cell_id).prob_being_active = prob_being_active;
 extended_calcium_analysis_output(cell_id).place_field_stability = place_field_stability;
%  extended_calcium_analysis_output(cell_id).KLD = KLD;
 
 extended_calcium_analysis_output(cell_id).PDF = posterior;
extended_calcium_analysis_output(cell_id).pN = pN;
 extended_calcium_analysis_output(cell_id).occupancy_map = occupancy_map;
 extended_calcium_analysis_output(cell_id).tuning_map = likelihood;
 
  split_half_stability(cell_id,1) = place_field_stability;
  tuning_map(cell_id).likelihood = likelihood;
  
  extended_calcium_analysis_output(cell_id).mean_dispersion = mean_dispersion;
  extended_calcium_analysis_output(cell_id).centroid = centroid; 
  
  place_field.mean_dispersion(cell_id,1) = mean_dispersion;
  place_field.centroid(cell_id,:) = centroid; 
  
 if analyze_heading;
 extended_calcium_analysis_output(cell_id).MI_heading = MI_heading;
 extended_calcium_analysis_output(cell_id).KLD_heading = KLD_heading;
 extended_calcium_analysis_output(cell_id).heading_stability = heading_stability;
 extended_calcium_analysis_output(cell_id).PDF_heading = posterior_heading;
 extended_calcium_analysis_output(cell_id).pN_heading = pN_heading;
 extended_calcium_analysis_output(cell_id).occupancy_vector_heading = occupancy_vector_heading;
 extended_calcium_analysis_output(cell_id).tuning_curve_heading = likelihood_heading;
 end 

 if analyze_running;
 extended_calcium_analysis_output(cell_id).MI_velocity = MI_velocity;
 extended_calcium_analysis_output(cell_id).KLD_velocity = KLD_velocity;
 extended_calcium_analysis_output(cell_id).MI_acceleration = MI_acceleration;
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

end
save('extended_calcium_analysis_output_upd.mat', 'extended_calcium_analysis_output');

end 


