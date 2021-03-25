%% This example is dedicate to spatial coding
% Used in figure XXXXX

function [] = open_field_decoding;
  clear all;
 load 'behav.mat';
 load 'ms.mat';
 
 num_cell_to_use = 100;
 
 min_speed_threshold = 2; % 2 cm.s-1
 
 %% Load the data
ca_time = ms.time/1000;
behav_time = behav.time/1000;
 behav_vec = behav.position;
 %behav_vec = behav.positionRescale;
ca_data = ms.Binary;
 
ntrials = 30;
shuffle = 1;

% load('behav_time.mat');
% load('behav_vec.mat');
% load('behav_vec_rescaled.mat');
% %behav_vec = behav_vec_rescaled;
% load('ca_data.mat');
% load('ca_time.mat');
% load('ca_trace.mat');

%% Control for unique timestamps (optional)
% In some rare cases, two or more timestamps with the same temporal value
% will prevent further interpolation. The following solves this.
[behav_time, IAbehav, ICbehav]=unique(behav_time);
[ca_time, IAms, ICms]=unique(ca_time);
ca_data = ca_data(IAms,:);
behav_vec = behav_vec(IAbehav,:);
%behav_time = behav_time'; 

%% Binarize calcium trace
 sampling_frequency = 30; % This data set has been sampled at 30 images per second
z_threshold = 2; % A 2 standard-deviation threshold is usually optimal to differentiate calcium ativity from background noise

binarized_data = ca_data; 

if shuffle
    shuffled_binarized = zeros(size(binarized_data));
    
   % for cell_i = 1:ms.numNeurons;
   for cell_i = 1:size(ca_data,2); 
random_ts = ceil(rand*length(ca_time));
     
shuffled_binarized(1:random_ts, cell_i) = binarized_data(end-random_ts+1:end, cell_i);
shuffled_binarized(random_ts+1:end,cell_i) = binarized_data(1:end-random_ts,cell_i);

    end
end

%% Interpolate behavior
% In most cases, behavior data has to be interpolated to match neural temporal
% activity assuming it has a lower sampling frequency than neural activity

interp_behav_vec(:,1) = interpolate_behavior(behav_vec(:,1), behav_time, ca_time); % in the X dimension
interp_behav_vec(:,2) = interpolate_behavior(behav_vec(:,2), behav_time, ca_time); % in the Y dimension

%% Compute velocities
% In this example, we will ignore moments when the mouse is immobile
[ velocity ] = extract_velocity(interp_behav_vec, ca_time);

%% Find periods of immobility
% This will be usefull later to exclude immobility periods from analysis
running_ts = velocity > min_speed_threshold;

%% Create an inclusion vector to isolate only a specific behavior
inclusion_vector = running_ts; % Only include periods of running

%% Compute occupancy and joint probabilities
% You can use 'min(behav_vec)' and 'max(behav_vec)'  to estimate the
% boundaries of the behavior vector (in our case, location in space)

bin_size = 3;
X_bin_vector = 0:bin_size:behav.trackLength+bin_size; % start : bin_size : end
X_bin_size = X_bin_vector(2) - X_bin_vector(1);
X_bin_centers_vector = X_bin_vector + X_bin_size/2;
X_bin_centers_vector(end) = [];
Y_bin_size = X_bin_size;
Y_bin_vector = X_bin_vector;
Y_bin_centers_vector = X_bin_centers_vector; % this will bin space in a square

%% Select the frames that are going to be used to train the decoder
training_set_creation_method = 'random'; % 'odd', odd timestamps; 'first_portion', first portion of the recording; 3, 'random' random frames
training_set_portion = 0.5; % Portion of the recording used to train the decoder for method 2 and 3

%% Decoding
for trial_i = 1:ntrials 
    disp(trial_i)
    
    training_ts = create_training_set(ca_time, training_set_creation_method, training_set_portion);
    training_ts(running_ts == 0) = 0; % Exclude periods of immobility from the traing set
    
    %% Create tuning maps for every cell
    
    occupancy_vector = [];
    prob_being_active = [];
    tuning_map_data = []; 
    
    for cell_i = 1:size(binarized_data,2)
        [~, ~, occupancy_vector, prob_being_active(cell_i), tuning_map_data(:,:,cell_i)] = extract_2D_information(binarized_data(:,cell_i), interp_behav_vec, X_bin_vector, Y_bin_vector, training_ts);
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
    
    % use the shuffled traces to compute the decoded probabilities
    if shuffle
     
     [shuffled_decoded_probabilities] = bayesian_decode2D(shuffled_binarized, occupancy_vector, prob_being_active, tuning_map_data, cell_used);
     [shuffled_decoded_probabilities] = bayesian_temporal_filter2D(shuffled_decoded_probabilities, ca_time, 1.5);
    end 
    
    [decoded_probabilities] = bayesian_decode2D(binarized_data, occupancy_vector, prob_being_active, tuning_map_data, cell_used);
    [decoded_probabilities] = bayesian_temporal_filter2D(decoded_probabilities, ca_time, 1.5);

    
    %% Determine the decoded bin/location
    for step_i = 1:size(decoded_probabilities,3)
        [val(step_i,1),decoded_bin(step_i,1)] = max(max(decoded_probabilities(:,:,step_i),[],1,'omitnan'));
        [val(step_i,2),decoded_bin(step_i,2)] = max(max(decoded_probabilities(:,:,step_i),[],2,'omitnan'));
        
        decoded_position(step_i,1) = X_bin_centers_vector(decoded_bin(step_i,1));
        decoded_position(step_i,2) = Y_bin_centers_vector(decoded_bin(step_i,2));
    end
    
    if shuffle 
      for step_i = 1:size(shuffled_decoded_probabilities,3)
        [val(step_i,1),shuffled_decoded_bin(step_i,1)] = max(max(shuffled_decoded_probabilities(:,:,step_i),[],1,'omitnan'));
        [val(step_i,2),shuffled_decoded_bin(step_i,2)] = max(max(shuffled_decoded_probabilities(:,:,step_i),[],2,'omitnan'));
        
        shuffled_decoded_position(step_i,1) = X_bin_centers_vector(shuffled_decoded_bin(step_i,1));
        shuffled_decoded_position(step_i,2) = Y_bin_centers_vector(shuffled_decoded_bin(step_i,2));
      end
    end
    
    %% Remove timestamps used for training
    decoded_bin(~decoding_ts) = nan;
    decoded_position(~decoding_ts) = nan;
    decoded_probabilities(:,:,~decoding_ts) = nan;
    
    if shuffle
    shuffled_decoded_bin(~decoding_ts) = nan;
    shuffled_decoded_position(~decoding_ts) = nan;
    shuffled_decoded_probabilities(:,:,~decoding_ts) = nan;
    end
    
    %% Decoding error
    % First we need to bin the actual data using the same bin vector used by
    % the decoder
    actual_bin = nan*interp_behav_vec;
    actual_position = nan*interp_behav_vec;
    
    for y = 1:length(Y_bin_vector)-1
        for x = 1:length(X_bin_vector)-1
            position_idx = find(interp_behav_vec(:,1)>X_bin_vector(x) & interp_behav_vec(:,1) < X_bin_vector(x+1) & interp_behav_vec(:,2) > Y_bin_vector(y) & interp_behav_vec(:,2) < Y_bin_vector(y+1));
            actual_bin(position_idx,1) = x;
            actual_bin(position_idx,2) = y;
            
            actual_position(position_idx,1) = X_bin_centers_vector(x);
            actual_position(position_idx,2) = Y_bin_centers_vector(y);
        end
    end
    
    actual_bin(~decoding_ts) = nan;
    actual_position(~decoding_ts) = nan;
    
    %% Compute decoding agreement  
    decoding_agreement_vector = double(decoded_bin == actual_bin);
    decoding_agreement_vector(isnan(decoded_bin)) = nan;
    decoding_agreement_vector(isnan(actual_bin)) = nan;
    decoding_agreement_vector(isnan(decoding_agreement_vector)) = [];
    decoding_agreement(trial_i) = sum(decoding_agreement_vector)./length(decoding_agreement_vector);
   
    decoding_error = [];
    for step_i = 1:size(decoded_position,1)
        decoding_error(step_i) = sqrt((actual_position(step_i,1)-decoded_position(step_i,1)).^2 + (actual_position(step_i,2)-decoded_position(step_i,2)).^2);
    end 
    
    mean_decoding_error(trial_i) = mean(abs(decoding_error), 'omitnan');
    
     if shuffle
           shuffled_decoding_agreement_vector = double(shuffled_decoded_bin == actual_bin);
    shuffled_decoding_agreement_vector(isnan(shuffled_decoded_bin)) = nan;
    shuffled_decoding_agreement_vector(isnan(actual_bin)) = nan;
    shuffled_decoding_agreement_vector(isnan(shuffled_decoding_agreement_vector)) = [];
    shuffled_decoding_agreement(trial_i) = sum(shuffled_decoding_agreement_vector)./length(shuffled_decoding_agreement_vector);
    
     shuffled_decoding_error = [];
    for step_i = 1:size(shuffled_decoded_position,1)
        shuffled_decoding_error(step_i) = sqrt((actual_position(step_i,1)-shuffled_decoded_position(step_i,1)).^2 + (actual_position(step_i,2)-shuffled_decoded_position(step_i,2)).^2);
    end 
    
    shuffled_mean_decoding_error(trial_i) = mean(abs(shuffled_decoding_error), 'omitnan');
     end
    
    
end

    save('decoding_agreement.mat', 'decoding_agreement');    
    save('decoding_error.mat', 'mean_decoding_error');
    
    save('shuffled_decoding_error.mat', 'shuffled_mean_decoding_error');
    save('shuffled_decoding_agreement.mat', 'shuffled_decoding_agreement');

% decoding.agreement = decoding_agreement;
% decoding.error = mean_decoding_error;
% decoding.shuffled_agreement = shuffled_decoding_agreement;
% decoding.shuffled_error = shuffled_mean_decoding_error;
% decoding.params.num_cell_to_use = num_cell_to_use;
% decoding.params.shuffle = shuffle;

%save('decoding.mat', 'decoding');

end
