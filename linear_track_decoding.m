%% This example is dedicate to spatial coding
% Used in figure XXXXX

function [] = linear_track_decoding;
   clear all;
   close all; 
  load 'behav.mat';
  load 'ms.mat';
  
  num_cell_to_use = 4; %% USED 40 cells for figure 

%% Load the data
ca_time = ms.time/1000;
behav_time = behav.time/1000;
behav_vec = behav.position(:,1);
ca_data = ms.Binary;

%% Control for unique timestamps (optional)
%In some rare cases, two or more timestamps with the same temporal value
%will prevent further interpolation. The following solves this.
[behav_time, IAbehav, ICbehav]=unique(behav_time);
[ca_time, IAms, ICms]=unique(ca_time);
ca_data = ca_data(IAms,:);
behav_vec = behav_vec(IAbehav,:);
%behav_time = behav_time; 
  
%  load 'ca_data.mat'
%  load 'ca_time.mat'
% load 'behav_vec.mat'
%  load 'behav_time.mat'
  
  shuffle = 30;

%% Binarize calcium trace
  sampling_frequency = 30; % This data set has been sampled at 30 images per second
  
binarized_data = ca_data; 

if shuffle
    shuffled_binarized = zeros(size(binarized_data));
    
  %  for cell_i = 1:ms.numNeurons;
        for cell_i = 1:size(ca_data,2)
random_ts = ceil(rand*length(ca_time));
     
shuffled_binarized(1:random_ts, cell_i) = binarized_data(end-random_ts+1:end, cell_i);
shuffled_binarized(random_ts+1:end,cell_i) = binarized_data(1:end-random_ts,cell_i);

    end
end

%% Interpolate behavior
% In most cases, behavior data has to be interpolated to match neural temporal
% activity assuming it has a lower sampling frequency than neural activity

interp_behav_vec(:,1) = interpolate_behavior(behav_vec(:,1), behav_time, ca_time); % in the X dimension

%% Compute velocities
% In this example, we will ignore moments when the mouse is immobile
[ velocity ] = extract_velocity(interp_behav_vec, ca_time);

%% Find periods of immobility
% This will be usefull later to exclude immobility periods from analysis
min_speed_threshold = 5; % 2 cm.s-1
running_ts = velocity > min_speed_threshold;

%% Create an inclusion vector to isolate only a specific behavior
inclusion_vector = running_ts; % Only include periods of running

%% Compute occupancy and joint probabilities
% You can use 'min(behav_vec)' and 'max(behav_vec)'  to estimate the
% boundaries of the behavior vector (in our case, location in space)

bin_size = 3;
bin_vector = 0:bin_size:behav.trackLength+bin_size; % start : bin_size : end
bin_size = bin_vector(2) - bin_vector(1);
bin_centers_vector = bin_vector + bin_size/2;
bin_centers_vector(end) = [];

%% Select the frames that are going to be used to train the decoder
training_set_creation_method = 'random'; % 'odd', odd timestamps; 'first_portion', first portion of the recording; 3, 'random' random frames
training_set_portion = 0.5; % Portion of the recording used to train the decoder for method 2 and 3


%% Decoding
for trial_i = 1:30
    disp(trial_i)
    
    training_ts = create_training_set(ca_time, training_set_creation_method, training_set_portion);
    training_ts(running_ts == 0) = 0; % Exclude periods of immobility from the traing set
    
    %% Create tuning curves for every cell
for cell_i = 1:size(binarized_data,2)
    [~, ~, ~, occupancy_vector, prob_being_active(cell_i), tuning_curve_data(:,cell_i) ] = extract_1D_information(binarized_data(:,cell_i), interp_behav_vec, bin_vector, training_ts);
end
    
    decoding_ts = ~training_ts; % Training timestamps are excluded
    decoding_ts(running_ts == 0) = 0; % Periods of immobility are excluded
    
    % Minimal a priori (use to remove experimental a priori)
    occupancy_vector = occupancy_vector./occupancy_vector*(1/length(occupancy_vector(:)));
    
    % Establish which cells are going to be used in the decoding process
     cell_used = logical(zeros(size(ca_data,2),1));
     while sum(cell_used) < num_cell_to_use
         randomidx = ceil(rand*(length(cell_used)));
         cell_used(randomidx) = 1;
     end
     
%cell_used = 0; 
    
    if shuffle
     [shuffled_decoded_probabilities] = bayesian_decode1D(shuffled_binarized, occupancy_vector, prob_being_active, tuning_curve_data, cell_used);
     [shuffled_decoded_probabilities] = bayesian_temporal_filter1D(shuffled_decoded_probabilities, ca_time, 0.5);
    end 
    
[decoded_probabilities] = bayesian_decode1D(binarized_data, occupancy_vector, prob_being_active, tuning_curve_data, cell_used);
[decoded_probabilities] = bayesian_temporal_filter1D(decoded_probabilities, ca_time, 0.5);

% figure;
% imagesc(decoded_probabilities);
% ax = gca;% ax.XLim = [20000 27181]; %to this for decoded probabilities 
% caxis([0 .1])
% 
% figure;
% imagesc(decoded_probabilities);
% ax = gca; ax.XLim = [5350 6750]%ax.XLim = [0 2000] % ; %to this for decoded probabilities 
% caxis([0 .1])
% colorbar

    %% Determine the decoded bin/location
[decoded_probabilities, decoded_bin] = max(decoded_probabilities,[],1);
decoded_position = bin_centers_vector(decoded_bin);

% figure;
% plot(decoded_position)
%  ax = gca;ax.XLim = [5350 6750] %ax.XLim = [0 2000] % 
%  set(gca, 'YDir','reverse')
%  hold on;
%  
if shuffle
    [shuffled_decoded_probabilities, shuffled_decoded_bin] = max(shuffled_decoded_probabilities,[],1);
    shuffled_decoded_position = bin_centers_vector(shuffled_decoded_bin);

    shuffled_decoded_bin(~decoding_ts) = nan;
    shuffled_decoded_position(~decoding_ts) = nan;
    shuffled_decoded_probabilities(:,:,~decoding_ts) = nan;
end 
    
    %% Remove timestamps used for training
    decoded_bin(~decoding_ts) = nan;
    decoded_position(~decoding_ts) = nan;
    decoded_probabilities(:,:,~decoding_ts) = nan;
    
    %% Decoding error
    % First we need to bin the actual data using the same bin vector used by
    % the decoder
    
   % Before looking at the error rate, we must first bin the actual data using the same bin vector used by
% the decoder
actual_bin = nan*interp_behav_vec;
actual_position = nan*interp_behav_vec;
for bin_i = 1:length(bin_vector)-1
    position_idx = find(interp_behav_vec>bin_vector(bin_i) & interp_behav_vec < bin_vector(bin_i+1));
    actual_bin(position_idx) = bin_i;
    actual_position(position_idx) = bin_centers_vector(bin_i);
end
%% Remove training timestamps to assess decoding error rate
 decoded_bin(~decoding_ts) = nan;
 decoded_position(~decoding_ts) = nan;
 decoded_probabilities(:,~decoding_ts) = nan;

%  plot(actual_position);
%   ax = gca;ax.XLim = [5350 6750];% ax.XLim = [0 2000] % 
%   
  
actual_bin(~decoding_ts) = nan;
actual_position(~decoding_ts) = nan;
actual_bin = actual_bin';
actual_position =  actual_position';


%% Compute decoding agreement
decoding_agreement_vector = double(decoded_bin == actual_bin);
decoding_agreement_vector(isnan(decoded_bin)) = nan;
decoding_agreement_vector(isnan(actual_bin)) = nan;
decoding_agreement_vector(isnan(decoding_agreement_vector)) = [];
decoding_agreement(trial_i) = sum(decoding_agreement_vector)./length(decoding_agreement_vector);
  
decoding_error = [];
decoding_error = actual_position - decoded_position;
mean_decoding_error(trial_i) = mean(abs(decoding_error), 'omitnan');

if shuffle
    shuffled_decoding_agreement_vector = double(shuffled_decoded_bin == actual_bin);
shuffled_decoding_agreement_vector(isnan(shuffled_decoded_bin)) = nan;
shuffled_decoding_agreement_vector(isnan(actual_bin)) = nan;
shuffled_decoding_agreement_vector(isnan(shuffled_decoding_agreement_vector)) = [];
shuffled_decoding_agreement(trial_i) = sum(shuffled_decoding_agreement_vector)./length(shuffled_decoding_agreement_vector);
  
shuffled_decoding_error = [];
shuffled_decoding_error = actual_position - shuffled_decoded_position;
shuffled_mean_decoding_error(trial_i) = mean(abs(shuffled_decoding_error), 'omitnan');
end 
    
 end
%     save('decoding_agreement.mat', 'decoding_agreement');
%     save('shuffled_decoding_agreement.mat', 'shuffled_decoding_agreement');
%     
%     save('decoding_error.mat', 'mean_decoding_error');
%     save('shuffled_decoding_error.mat', 'shuffled_mean_decoding_error'); 

 decoding.agreement = decoding_agreement;
decoding.error = mean_decoding_error;
decoding.shuffled_agreement = shuffled_decoding_agreement;
decoding.shuffled_error = shuffled_mean_decoding_error;
decoding.params.num_cell_to_use = num_cell_to_use;
decoding.params.shuffle = shuffle;

save('decoding.mat', 'decoding');
end