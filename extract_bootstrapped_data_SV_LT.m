function [cell_ID, mean_MI, SEM_MI, mean_MI_shuffled, SEM_MI_shuffled] = extract_bootstrapped_data;

clear all;
load('ms.mat');
load('behav.mat');

ca_data = ms.RawTraces;
ca_time = ms.time/1000;
behav_vec = behav.position(:,1);
behav_time = behav.time/1000;

%% Parameters
sampling_frequency = 30;
training_set_portion = 0.5;
numShuffles = 30;
bin_size = 3; 
analyze_running = 1; 

% Creating binnedvectors
min_X = 0;
max_X = 100; % Set these to the borders of the maze used. !!! Has to be the same for every mouse for consistent MI computation !!!

X_bin_vector = min_X:bin_size:max_X+bin_size;
X_bin_centers_vector = X_bin_vector + bin_size/2;
X_bin_centers_vector(end) = [];

%% Clean data: keep only unique points
[behav_time, IAbehav, ICbehav]=unique(behav_time);
[ca_time, IAms, ICms]=unique(ca_time);
ca_data = ca_data(IAms,:);
behav_vec = behav_vec(IAbehav,:);

%% Binarize data
binarized_data = ms.Binary;
binarized_data = binarized_data(IAms,:); %keep only unique points 

%% Interpolate behavior
% In most cases, behavior data has to be interpolated to match neural temporal
% activity assuming it has a lower sampling frequency than neural activity
[interp_behav_vec] = interpolate_behavior(behav_vec, behav_time, ca_time);
interp_behav_vec(end) = interp_behav_vec(end-1); % Make sure to interpolate or remove every NaN so that every timepoint has a corresponding behavioral state

z_delta_x = diff(interp_behav_vec);
z_delta_x(end+1,1) = 0;
z_delta_x(isnan(z_delta_x)) = 0;
z_delta_x = zscore(z_delta_x);

%% Compute velocities
% In this example, we will ignore moments when the mouse is immobile
[velocity] = extract_velocity(interp_behav_vec, ca_time);

%% Compute accelaration
acceleration = diff(velocity);
acceleration(end+1) = 0;

%% Find periods of immobility
% This will be usefull later to exclude immobility periods from analysis
min_speed_threshold = 5; % cm.s-1
running_ts = velocity > min_speed_threshold;

%% Create an inclusion vector to isolate only a specific behavior
inclusion_vector = running_ts; % Only include periods of running
included_velocity = velocity(inclusion_vector); 
included_acceleration = acceleration(inclusion_vector);

%% right trajectories only 
direction_indices = z_delta_x > 0.2;
inclusion_vector_right = direction_indices; % Only include right runs
inclusion_vector_right(running_ts == 0) = 0; % Only include periods of running

%% left trajectories only
direction_indices = z_delta_x < -0.2;
inclusion_vector_left = direction_indices; % Only include left runs
inclusion_vector_left(running_ts == 0) = 0; % Only include periods of running

%%
if analyze_running 
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

%% Perform bootstrapping here

for shuffle_i = 1:numShuffles
    bootstrap_ts = zeros(1,length(ca_time));
    numFrames = length(ca_time);
    half_ts = ceil(numFrames*training_set_portion);
    bootstrap_ts(1:half_ts) = 1;
    bootstrap_ts = logical(bootstrap_ts(randperm(numFrames)));
    
    bootstrap_ts = bootstrap_ts == 1 & inclusion_vector' == 1; % Exclude periods of immobility from the traing set
    
    bootstrap_ts_right = bootstrap_ts == 1 & inclusion_vector_right' == 1;
    bootstrap_ts_left = bootstrap_ts == 1 & inclusion_vector_left' == 1;
    
    random_ts = ceil(rand*length(ca_time));
    shuffled_binarized = zeros(size(binarized_data));
    
    % Permute the trace
    shuffled_binarized(1:random_ts,:) = binarized_data(end-random_ts+1:end,:);
    shuffled_binarized(random_ts+1:end,:) = binarized_data(1:end-random_ts,:);
    
    for cell_i = 1:size(binarized_data,2)
        [~, MI(shuffle_i,cell_i), ~, ~, ~, ~] = extract_1D_information(binarized_data(:,cell_i), interp_behav_vec, X_bin_vector, bootstrap_ts);
        [~, MI_right(shuffle_i,cell_i), ~, ~, ~, ~] = extract_1D_information(binarized_data(:,cell_i), interp_behav_vec, X_bin_vector, bootstrap_ts_right);
        [~, MI_left(shuffle_i,cell_i), ~, ~, ~, ~] = extract_1D_information(binarized_data(:,cell_i), interp_behav_vec, X_bin_vector, bootstrap_ts_left);

        [~, shuffled_MI(shuffle_i,cell_i), ~, ~, ~, ~] = extract_1D_information(shuffled_binarized(:,cell_i), interp_behav_vec, X_bin_vector, bootstrap_ts);
        [~,shuffled_MI_right(shuffle_i,cell_i), ~, ~, ~, ~] = extract_1D_information(shuffled_binarized(:,cell_i), interp_behav_vec, X_bin_vector, bootstrap_ts_right);
        [~,shuffled_MI_left(shuffle_i,cell_i), ~, ~, ~, ~] = extract_1D_information(shuffled_binarized(:,cell_i), interp_behav_vec, X_bin_vector, bootstrap_ts_left);
        
        if analyze_running;
        [~,MI_velocity(shuffle_i,cell_i), ~, ~, ~, ~] = extract_1D_information(binarized_data(:,cell_i), velocity, velocity_bin_vector, bootstrap_ts);
        [~, shuffled_MI_velocity(shuffle_i,cell_i), ~, ~, ~, ~] = extract_1D_information(shuffled_binarized(:,cell_i), velocity, velocity_bin_vector, bootstrap_ts);
    
        [~, MI_acceleration(shuffle_i,cell_i), ~, ~, ~, ~] = extract_1D_information(binarized_data(:,cell_i), acceleration, acceleration_bin_vector, bootstrap_ts);
        [~, shuffled_MI_acceleration(shuffle_i,cell_i), ~, ~, ~, ~] = extract_1D_information(shuffled_binarized(:,cell_i), acceleration, acceleration_bin_vector, bootstrap_ts);
        end 
    end
end

%%

cell_ID = 1:size(ca_data,2);
mean_MI = mean(MI);
SEM_MI = std(MI)/sqrt(numShuffles);
mean_MI_shuffled = mean(shuffled_MI);
SEM_MI_shuffled = std(shuffled_MI)/sqrt(numShuffles);

mean_MI_right = mean(MI_right);
SEM_MI_right = std(MI_right)/sqrt(numShuffles);
mean_MI_shuffled_right = mean(shuffled_MI_right);
SEM_MI_shuffled_right = std(shuffled_MI_right)/sqrt(numShuffles);

mean_MI_left = mean(MI_left);
SEM_MI_left = std(MI_left)/sqrt(numShuffles);
mean_MI_shuffled_left = mean(shuffled_MI_left);
SEM_MI_shuffled_left = std(shuffled_MI_left)/sqrt(numShuffles);

bootstrapped_data.cell_ID = cell_ID';
bootstrapped_data.mean_MI = mean_MI';
bootstrapped_data.SEM_MI = SEM_MI';
bootstrapped_data.mean_MI_shuffled = mean_MI_shuffled';
bootstrapped_data.SEM_MI_shuffled = SEM_MI_shuffled';

bootstrapped_data.mean_MI_right = mean_MI_right';
bootstrapped_data.SEM_MI_right = SEM_MI_right';
bootstrapped_data.mean_MI_shuffled_right = mean_MI_shuffled_right';
bootstrapped_data.SEM_MI_shuffled_right = SEM_MI_shuffled_right';

bootstrapped_data.mean_MI_left = mean_MI_left';
bootstrapped_data.SEM_MI_left = SEM_MI_left';
bootstrapped_data.mean_MI_shuffled_left = mean_MI_shuffled_left';
bootstrapped_data.SEM_MI_shuffled_left = SEM_MI_shuffled_left';

if analyze_running;
    
mean_MI_velocity = mean(MI_velocity);
SEM_MI_velocity = std(MI_velocity)/sqrt(numShuffles);
mean_MI_velocity_shuffled = mean(shuffled_MI_velocity);
SEM_MI_velocity_shuffled = std(shuffled_MI_velocity)/sqrt(numShuffles);

mean_MI_acceleration = mean(MI_acceleration);
SEM_MI_acceleration = std(MI_acceleration)/sqrt(numShuffles);
mean_MI_acceleration_shuffled = mean(shuffled_MI_acceleration);
SEM_MI_acceleration_shuffled = std(shuffled_MI_acceleration)/sqrt(numShuffles);

bootstrapped_data.mean_MI_velocity = mean_MI_velocity';
bootstrapped_data.SEM_MI_velocity = SEM_MI_velocity';
bootstrapped_data.mean_MI_velocity_shuffled = mean_MI_velocity_shuffled';
bootstrapped_data.SEM_MI_velocity_shuffled = SEM_MI_velocity_shuffled';

bootstrapped_data.mean_MI_acceleration = mean_MI_acceleration';
bootstrapped_data.SEM_MI_acceleration = SEM_MI_acceleration';
bootstrapped_data.mean_MI_acceleration_shuffled = mean_MI_acceleration_shuffled';
bootstrapped_data.SEM_MI_acceleration_shuffled = SEM_MI_acceleration_shuffled';
end 


save('bootstrapped_data_upd.mat', 'bootstrapped_data');

end
