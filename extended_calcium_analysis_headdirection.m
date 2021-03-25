function [extended_calcium_analysis_output] = extended_calcium_analysis_headdirection_upd;

% Potentially run this by running select_folders_for_analysis.m
%%                             
 clear all; 
 load 'behav.mat';
 load 'ms.mat';
 load 'behavDLC.mat'; 
 
 load 'extended_calcium_analysis_output.mat'; 
 
 poolobj = gcp('nocreate'); % If no pool, do not create new one.
if isempty(poolobj)
    poolsize = 0;
else
    poolsize = poolobj.NumWorkers
end

plotting = 0; 
save_figure = 0;
min_speed_threshold = 5; % 2 cm.s-1
numShuffles = 1000;
params.smoothing = true;

if ~isfield(ms,'Binary');
    disp('Extracting bixsnary information');
    ms = msExtractBinary(ms);
end

%% determine behav and ms timestamps and convert to seconds

%% Re-aligning calcium and behavior times
behav_time = behav.time/1000; %convert to seconds 
behav_vec = behav.position;
heading_vec = behavDLC.headDirection;

[behav_time, idx]=unique(behav_time); %getting rid of duplicate timestamps
behav_vec = behav_vec(idx,:); %excluding the positions associated with these double timestamps
heading_vec = heading_vec(idx);
ca_time = ms.time/1000; % convert to seconds 
[ca_time, idx]=unique(ca_time); %getting rid of duplicate timestamps
heading_vec = heading_vec';

 %interpolate behav trajectory, this is needed to compute running speed
interp_behav_vec(:,1) = interpolate_behavior(behav_vec(:,1), behav_time, ca_time); % in the X dimension
interp_behav_vec(:,2) = interpolate_behavior(behav_vec(:,2), behav_time, ca_time); % in the Y dimension
interp_behav_vec(end,:) = interp_behav_vec(end-1,:); % Make sure to interpolate or remove every NaN so that every timepoint has a corresponding behavioral state

 %Interpolate head direction
[interp_heading_vec] = interpolate_HD(heading_vec, behav_time, ca_time);
interp_heading_vec(end) = interp_heading_vec(end-1); % Make sure to interpolate or remove every NaN so that every timepoint has a corresponding behavioral state

%% Compute velocities
% In this example, we will ignore moments when the mouse is immobile
[velocity] = extract_velocity(interp_behav_vec, ca_time);

%% Find periods of immobility
% This will be usefull later to exclude immobility periods from analysis
running_ts = velocity > min_speed_threshold;

%% Create an inclusion vector to isolate only a specific behavior
inclusion_vector = running_ts; 

for cell_id = 1:ms.numNeurons

%% Binarize calcium trace
sampling_frequency = 30; % This data set has been sampled at 30 images per second
binarized_trace = ms.Binary(:, cell_id); 
binarized_trace = binarized_trace(idx); %get rid of any duplicate timestamps

% Binning polar space
bin_size = 9; %open field with 3 cm bins in 45 cm environment has 256 bins total 
bin_vector = 0:bin_size:360;

% sorted_trace = zeros(length(binarized_trace),1); 
% 
%     [sorted_heading, sorted_heading_idx] = sort(deg2rad(heading));
%     sorted_trace = binarized_trace(sorted_heading_idx);
%     sorted_trace(sorted_trace < 0) = 0;

%% Compute occupancy and joint probabilities

[~, MI_heading, posterior_heading, occupancy_vector_heading, ~, likelihood_heading ] = extract_1D_information(binarized_trace, interp_heading_vec, bin_vector, inclusion_vector);
peak_joint_probability = max(likelihood_heading(:));

if plotting;
plotting_fig = figure;
subplot(3,2,1)
plot(likelihood_heading,'Color', [0 0.1 0.8])
xlim([0 256])
title 'Tuning curve'
xlabel 'direction (bin)'
ylabel 'Probability of firing in location'
subplot(3,2,3)
plot(occupancy_vector_heading,'Color', [0.1 0.8 0.1])
xlim([0 256])
title 'Occupancy'
xlabel 'Location on the track (cm)'
ylabel 'Relative occupancy'
subplot(3,2,5)
plot(posterior_heading,'Color', [0.8 0.2 0])
xlim([0 256])
title 'Posterior probability density function'
xlabel 'Location on the track (cm)'
ylabel 'Probability'
subplot(3,2,[2 4])
polarplot(likelihood_heading)
title 'Tuning curve'

end

if save_figure;
    saveas(plotting_fig, sprintf('%d.jpg',cell_id));
end

%% Shuffle data
shuffled_likelihood_heading = zeros(length(likelihood_heading), numShuffles);
shuffled_MI_heading = zeros(length(MI_heading), numShuffles);
shuffled_posterior_heading = zeros(length(posterior_heading), numShuffles);

parfor k = 1:numShuffles
    random_ts = ceil(rand*length(ca_time));
    shuffled_binarized = zeros(length(binarized_trace),1);
    
    % Permute the trace
    shuffled_binarized(1:random_ts) = binarized_trace(end-random_ts+1:end);
    shuffled_binarized(random_ts+1:end) = binarized_trace(1:end-random_ts);
    
    % Compute tuning curve
    % [~, ~, ~, ~, ~, shuffled_likelihood_heading(:,k)] = extract_1D_information(shuffled_binarized, interp_behav_vec, bin_vector, running_ts);

    [~, shuffled_heading_MI(k), ~, ~ , ~, shuffled_likelihood_heading(:,k)] = extract_1D_information(shuffled_binarized, interp_heading_vec, bin_vector, inclusion_vector);
    %shuffled_peak_joint_probability(k) = max(shuffled_likelihood_heading(:));   
end 


heading_MI_pvalue = sum(shuffled_heading_MI>MI_heading)/numShuffles;
heading_zMI = (MI_heading-mean(shuffled_heading_MI))/std(shuffled_heading_MI);


%% Compute significance and determine if this cell is a place cell

pN_heading = sum(shuffled_likelihood_heading > likelihood_heading,2)/numShuffles; %  pN, supra-threshold tests

if plotting;
plotting_fig = figure;
subplot(3,2,1)
plot(pN_heading)
ylim([0 .5]);
xlim([0 256])
title 'p-value map'
subplot(3,2,3)
plot(shuffled_likelihood_heading, 'b');
hold on
plot(likelihood_heading, 'r'); 
xlim([0 256])
ylim([0 .5]);
title 'significant tuning map'

subplot(3,2,[2 4])
polarplot(shuffled_likelihood_heading, 'b');
hold on
polarplot(likelihood_heading, 'r')
title 'Tuning curve'
end

significant_likelihood_heading = likelihood_heading;
significant_likelihood_heading(pN_heading > 0.01) = 0;

%% stability

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
% All trajectories 
extended_calcium_analysis_output(cell_id).MI_HeadDirection = MI_heading;
extended_calcium_analysis_output(cell_id).headdirection_stability = heading_stability;
 extended_calcium_analysis_output(cell_id).headdirection_MI_pvalue = heading_MI_pvalue;
 extended_calcium_analysis_output(cell_id).headdirection_zMI = heading_zMI;
 
extended_calcium_analysis_output(cell_id).PDF_HeadDirection = posterior_heading;
extended_calcium_analysis_output(cell_id).pN_HeadDirection = pN_heading;
extended_calcium_analysis_output(cell_id).occupancy_vector_HeadDirection = occupancy_vector_heading;
extended_calcium_analysis_output(cell_id).tuning_curve_HeadDirection = likelihood_heading;

  split_half_stability(cell_id,1) = heading_stability;
  tuning_vector(cell_id).likelihood = likelihood_heading;
  

end;


 % save('split_half_stability.mat', 'split_half_stability');
 % save('tuning_vector.mat', 'tuning_vector'); 

save('extended_calcium_analysis_output_HD.mat', 'extended_calcium_analysis_output');

end 


