%% Control for unique timestamps (optional)
% In some rare cases, two or more timestamps with the same temporal value
% will prevent further interpolation. The following solves this.
%clear all
load('ms.mat');
load('behav.mat'); 
load('behavDLC.mat');

behav_time = behav.time/1000;
ca_time = ms.time/1000;
ca_data = ms.Binary;
HD_vec = behavDLC.headDirection;
position = behav.position;

[behav_time, IAbehav, ICbehav]=unique(behav_time);
[ca_time, IAms, ICms]=unique(ca_time);
ca_data = ca_data(IAms,:);
HD_vec = HD_vec(IAbehav);
position = position(IAbehav,:);


%% 

shuffle = 1;

num_cell_to_use = 20;

%% Binarize calcium trace
% sampling_frequency = 30; % This data set has been sampled at 30 images per second
% z_threshold = 2; % A 2 standard-deviation threshold is usually optimal to differentiate calcium ativity from background noise
% 
% binarized_data = zeros(size(ca_data));
% for cell_i = 1:size(ca_data,2)
%     binarized_data(:,cell_i) = extract_binary(ca_data(:,cell_i),sampling_frequency, z_threshold);
% end

binarized_data = ca_data; 

%% Interpolate behavior
% In most cases, behavior data has to be interpolated to match neural temporal
% activity assuming it has a lower sampling frequency than neural activity

%interp_HD_vec = interpolate_behavior(HD_vec', behav_time', ca_time);
interp_position = interpolate_behavior(position,behav_time, ca_time);

[interp_HD_vec] = interpolate_HD(HD_vec', behav_time, ca_time);
interp_HD_vec(end) = interp_HD_vec(end-1); % Make sure to interpolate or remove every NaN so that every timepoint has a corresponding behavioral state

% ERROR SCALE CORRECTION
%interp_position = interp_position*1.8;

if shuffle
    shuffled_binarized = zeros(size(binarized_data));
    
   % for cell_i = 1:ms.numNeurons;
   for cell_i = 1:size(ca_data,2); 
random_ts = ceil(rand*length(ca_time));
     
shuffled_binarized(1:random_ts, cell_i) = binarized_data(end-random_ts+1:end, cell_i);
shuffled_binarized(random_ts+1:end,cell_i) = binarized_data(1:end-random_ts,cell_i);

    end
end

%% Compute velocities
% In this example, we will ignore moments when the mouse is immobile
[ velocity ] = extract_velocity(interp_position, ca_time);

%% Find periods of immobility
% This will be usefull later to exclude immobility periods from analysis
min_speed_threshold = 5; % 2 cm.s-1
running_ts = velocity > min_speed_threshold;

%% Create an inclusion vector to isolate only a specific behavior
inclusion_vector = running_ts; % Only include periods of running
inclusion_vector(isnan(interp_HD_vec)) = 0;

%% Compute occupancy and joint probabilities
% You can use 'min(behav_vec)' and 'max(behav_vec)'  to estimate the
% boundaries of the behavior vector (in our case, location in space)
bin_size = 15;
deg_bin_vector = 0:bin_size:360; % start : bin_size : end
deg_bin_centers_vector = deg_bin_vector + bin_size/2;
deg_bin_centers_vector(end) = [];

%% Select the frames that are going to be used to train the decoder
training_set_creation_method = 'random'; % 'odd', odd timestamps; 'first_portion', first portion of the recording; 3, 'random' random frames
training_set_portion = 0.5; % Portion of the recording used to train the decoder for method 2 and 3

%% Decoding
% First let's binarize traces from all cells
for trial_i = 1:30
    disp(trial_i)
    
    training_ts = create_training_set(ca_time, training_set_creation_method, training_set_portion);
    training_ts(running_ts == 0) = 0; % Exclude periods of immobility from the traing set
    
    %% Create tuning maps for every cell
    for cell_i = 1:size(binarized_data,2)
        [~, ~, ~, occupancy_vector, prob_being_active(cell_i), tuning_curve_data(:,cell_i)] = extract_1D_information(binarized_data(:,cell_i), interp_HD_vec, deg_bin_vector, training_ts);
    end
    
    %% Plot the tunning curves
    tuning_plotting = false;
    if tuning_plotting
        [~,max_index] = max(tuning_curve_data,[],1);
        [~,sorted_index] = sort(max_index);
        sorted_tuning_curve_data = tuning_curve_data(:,sorted_index);
        
        figure
        imagesc(deg_bin_centers_vector,1:size(ca_data,2),sorted_tuning_curve_data')
        daspect([1 1 1])
        title 'Neuronal tuning curves'
        xlabel 'Head direction (�)'
        ylabel 'Cell ID'
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
    
  if shuffle
     [shuffled_decoded_probabilities] = bayesian_decode1D(shuffled_binarized, occupancy_vector, prob_being_active, tuning_curve_data, cell_used);

    end 
    
    [decoded_probabilities] = bayesian_decode1D(binarized_data, occupancy_vector, prob_being_active, tuning_curve_data, cell_used);
    
    %% Determine the decoded bin/head direction
    [max_decoded_prob, decoded_bin] = max(decoded_probabilities,[],1);
    decoded_HD = deg_bin_centers_vector(decoded_bin);
    
    if shuffle
    [shuffled_decoded_probabilities, shuffled_decoded_bin] = max(shuffled_decoded_probabilities,[],1);
    shuffled_decoded_HD = deg_bin_centers_vector(shuffled_decoded_bin);
     end 
    
    %% Remove timestamps used for training
    decoded_bin(~decoding_ts) = nan;
    decoded_HD(~decoding_ts) = nan;
    
    if shuffle; 
    shuffled_decoded_bin(~decoding_ts) = nan;
    shuffled_decoded_HD(~decoding_ts) = nan;
    end
    
    %% Decoding error
% Before looking at the error rate, we must first bin the actual data using the same bin vector used by
% the decoder
actual_bin = nan*interp_HD_vec;
actual_HD = nan*interp_HD_vec;
for bin_i = 1:length(deg_bin_vector)-1
    HD_idx = find(interp_HD_vec>deg_bin_vector(bin_i) & interp_HD_vec < deg_bin_vector(bin_i+1));
    actual_bin(HD_idx) = bin_i;
    actual_HD(HD_idx) = deg_bin_centers_vector(bin_i);
end
    
    %% Compute decoding agreement
    decoding_agreement_vector = double(decoded_bin' == actual_bin);
    decoding_agreement_vector(isnan(decoded_bin)) = nan;
    decoding_agreement_vector(isnan(actual_bin)) = nan;
    decoding_agreement_vector(isnan(decoding_agreement_vector)) = [];
    decoding_agreement(trial_i) = sum(decoding_agreement_vector)./length(decoding_agreement_vector);
    
    decoding_error = [];
    for step_i = 1:size(decoded_HD,2)
        decoding_error(step_i) = mod(decoded_HD(step_i)-actual_HD(step_i),360);
    end
    
    mean_decoding_error(trial_i) = mean(decoding_error, 'omitnan');
    
if shuffle
shuffled_decoding_agreement_vector = double(shuffled_decoded_bin' == actual_bin);
shuffled_decoding_agreement_vector(isnan(shuffled_decoded_bin)) = nan;
shuffled_decoding_agreement_vector(isnan(actual_bin)) = nan;
shuffled_decoding_agreement_vector(isnan(shuffled_decoding_agreement_vector)) = [];
shuffled_decoding_agreement(trial_i) = sum(shuffled_decoding_agreement_vector)./length(shuffled_decoding_agreement_vector);
  
shuffled_decoding_error = [];
    for step_i = 1:size(shuffled_decoded_HD,2)
        shuffled_decoding_error(step_i) = mod(shuffled_decoded_HD(step_i)-actual_HD(step_i),360);
    end
shuffled_mean_decoding_error(trial_i) = mean(abs(shuffled_decoding_error), 'omitnan');
end 
end


decoding.agreement = decoding_agreement;
decoding.error = mean_decoding_error;
decoding.shuffled_agreement = shuffled_decoding_agreement;
decoding.shuffled_error = shuffled_mean_decoding_error;
decoding.params.num_cell_to_use = num_cell_to_use;
decoding.params.shuffle = shuffle;

%decoding_test = decoding; 

save('decoding.mat', 'decoding');



