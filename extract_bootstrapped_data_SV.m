function [cell_ID, mean_MI, SEM_MI, mean_MI_shuffled, SEM_MI_shuffled] = extract_bootstrapped_data;
% % 
  clear all;
   load('ms.mat');
   load('behav.mat');
 
%% Parameters
sampling_frequency = 30;
training_set_portion = 0.5;
numShuffles = 30;
bin_size = 3;

%% Load the open field data
  ca_data = ms.RawTraces;
  ca_time = ms.time/1000;
  behav_vec = behav.position;
  %behav_vec = behav.positionRescale;
  behav_time = behav.time/1000;

%% Creating binnedvectors
min_X = 0;
min_Y = 0;

X_bin_vector = min_X:bin_size:49+bin_size; 
X_bin_centers_vector = X_bin_vector + bin_size/2;
X_bin_centers_vector(end) = [];
Y_bin_vector = min_Y:bin_size:49+bin_size; % same here
Y_bin_centers_vector = Y_bin_vector + bin_size/2;
Y_bin_centers_vector(end) = [];

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
interp_behav_vec(:,1) = interp1(behav_time,behav_vec(:,1),ca_time);
interp_behav_vec(:,2) = interp1(behav_time,behav_vec(:,2),ca_time);
interp_behav_vec(end,:) = interp_behav_vec(end-1,:); % Make sure to interpolate or remove every NaN so that every timepoint has a corresponding behavioral state

%% Compute velocities
% In this example, we will ignore moments when the mouse is immobile
velocity(1,:) = 0; % Velocity at t1 (coordinates x1,y1) to be zero. Delta_v1 is estimated for a displacement to coordinates x2,y2.
    for i=2:length(ca_time)
        velocity(i,:) = sqrt((interp_behav_vec(i,1)-interp_behav_vec(i-1,1)).^2 + (interp_behav_vec(i,2)-interp_behav_vec(i-1,2)).^2)/(ca_time(i)-ca_time(i-1));
    end
velocity = smooth(velocity,round(1/mode(diff(ca_time))));

%% Compute acceleration
acceleration = diff(velocity);
acceleration(end+1) = 0;

%% Find periods of immobility
% This will be usefull later to exclude immobility periods from analysis
min_speed_threshold = 5; % cm.s-1
running_ts = velocity > min_speed_threshold;

%% Create an inclusion vector to isolate only a specific behavior
inclusion_vector = running_ts; % Only include periods of running

% inclusion vector for velocity and acceleration
min_speed_threshold = 0; % cm.s-1
running_ts = velocity > min_speed_threshold;
inclusion_vector_velocity = running_ts;

%% Compute bins for velocity and acceleration
%here compute the min and max values after kicking out the non-running
%timestamps 
min_velocity = 0;
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

%% Perform bootstrapping here

for shuffle_i = 1:numShuffles
    bootstrap_ts = zeros(1,length(ca_time));
    numFrames = length(ca_time);
    half_ts = ceil(numFrames*training_set_portion);
    bootstrap_ts(1:half_ts) = 1;
    bootstrap_ts = logical(bootstrap_ts(randperm(numFrames)));
    
    bootstrap_ts = bootstrap_ts == 1 & inclusion_vector' == 1; % Exclude periods of immobility from the traing set
    bootstrap_ts_velocity = bootstrap_ts == 1 & inclusion_vector_velocity' == 1; % Exclude periods of immobility from the traing set
    
    random_ts = ceil(rand*length(ca_time));
    shuffled_binarized = zeros(size(binarized_data));
    
    % Permute the trace
    shuffled_binarized(1:random_ts,:) = binarized_data(end-random_ts+1:end,:);
    shuffled_binarized(random_ts+1:end,:) = binarized_data(1:end-random_ts,:);
    
    for cell_i = 1:size(binarized_data,2)
        %% spatial 
        [MI(shuffle_i,cell_i), ~, ~, ~, ~] = extract_2D_information(binarized_data(:,cell_i), interp_behav_vec, X_bin_vector, Y_bin_vector, bootstrap_ts);
        [shuffled_MI(shuffle_i,cell_i), ~, ~, ~, ~] = extract_2D_information(shuffled_binarized(:,cell_i), interp_behav_vec, X_bin_vector, Y_bin_vector, bootstrap_ts);
    
        %% velocity
        [~, MI_velocity(shuffle_i,cell_i), ~, ~, ~, ~] = extract_1D_information(binarized_data(:,cell_i), velocity, velocity_bin_vector,   bootstrap_ts_velocity);
        [~, shuffled_MI_velocity(shuffle_i,cell_i), ~, ~, ~, ~] = extract_1D_information(shuffled_binarized(:,cell_i), velocity, velocity_bin_vector,   bootstrap_ts_velocity);
    
        % acceleration
        [~, MI_acceleration(shuffle_i,cell_i), ~, ~, ~, ~] = extract_1D_information(binarized_data(:,cell_i), acceleration, acceleration_bin_vector,   bootstrap_ts_velocity);
        [~, shuffled_MI_acceleration(shuffle_i,cell_i), ~, ~, ~, ~] = extract_1D_information(shuffled_binarized(:,cell_i), acceleration, acceleration_bin_vector,   bootstrap_ts_velocity);

    end
end

cell_ID = 1:size(ca_data,2);
mean_MI = mean(MI);
SEM_MI = std(MI)/sqrt(numShuffles);
mean_MI_shuffled = mean(shuffled_MI);
SEM_MI_shuffled = std(shuffled_MI)/sqrt(numShuffles);

bootstrapped_data.cell_ID = cell_ID';
bootstrapped_data.mean_MI = mean_MI';
bootstrapped_data.SEM_MI = SEM_MI';
bootstrapped_data.mean_MI_shuffled = mean_MI_shuffled';
bootstrapped_data.SEM_MI_shuffled = SEM_MI_shuffled';

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

save('bootstrapped_data_upd_2021.mat', 'bootstrapped_data');


end
