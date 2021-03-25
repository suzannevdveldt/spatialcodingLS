function [cell_ID, mean_MI, SEM_MI, mean_MI_shuffled, SEM_MI_shuffled] = extract_bootstrapped_data_SV_heading;

clear all;
load('ms.mat');
load('behav.mat');
load('behavDLC.mat');

%% Load the open field data
ca_data = ms.RawTraces;
ca_time = ms.time/1000;
behav_vec = behav.position;
behav_time = behav.time/1000;
heading_vec = behavDLC.headDirection';

%% Parameters
sampling_frequency = 30;
training_set_portion = 0.5;
numShuffles = 30;
min_speed_threshold = 5; % cm.s-1

%% Clean data: keep only unique points

[behav_time, idx]=unique(behav_time); %getting rid of duplicate timestamps
behav_vec = behav_vec(idx,:); %excluding the positions associated with these double timestamps
heading_vec = heading_vec(idx);

[ca_time, IAms, ICms]=unique(ca_time);
ca_data = ca_data(IAms,:);

behav_time = behav_time';
heading_vec = heading_vec';

%% Binarize data
binarized_data = ms.Binary;
binarized_data = binarized_data(IAms,:); %keep only unique points 

%% Interpolate behavior
% In most cases, behavior data has to be interpolated to match neural temporal
% activity assuming it has a lower sampling frequency than neural activity
interp_behav_vec(:,1) = interp1(behav_time,behav_vec(:,1),ca_time);
interp_behav_vec(:,2) = interp1(behav_time,behav_vec(:,2),ca_time);
interp_behav_vec(end,:) = interp_behav_vec(end-1,:); % Make sure to interpolate or remove every NaN so that every timepoint has a corresponding behavioral state

 %Interpolate head direction 
[interp_heading_vec] = interpolate_HD(heading_vec, behav_time, ca_time);
interp_heading_vec(end) = interp_heading_vec(end-1); % Make sure to interpolate or remove every NaN so that every timepoint has a corresponding behavioral state

%% Compute velocities
% In this example, we will ignore moments when the mouse is immobile
velocity(1,:) = 0; % Velocity at t1 (coordinates x1,y1) to be zero. Delta_v1 is estimated for a displacement to coordinates x2,y2.
    for i=2:length(ca_time)
        velocity(i,:) = sqrt((interp_behav_vec(i,1)-interp_behav_vec(i-1,1)).^2 + (interp_behav_vec(i,2)-interp_behav_vec(i-1,2)).^2)/(ca_time(i)-ca_time(i-1));
    end
velocity = smooth(velocity,round(1/mode(diff(ca_time))));

%% Find periods of immobility
% This will be usefull later to exclude immobility periods from analysis
running_ts = velocity > min_speed_threshold;

%% Create an inclusion vector to isolate only a specific behavior
inclusion_vector = running_ts; % Only include periods of running

%% Compute bins for velocity and acceleration
%here compute the min and max values after kicking out the non-running
%timestamps 
bin_vector = 0:9:360;

%% Perform bootstrapping here

for shuffle_i = 1:numShuffles
    bootstrap_ts = zeros(1,length(ca_time));
    numFrames = length(ca_time);
    half_ts = ceil(numFrames*training_set_portion);
    bootstrap_ts(1:half_ts) = 1;
    bootstrap_ts = logical(bootstrap_ts(randperm(numFrames)));
    
    bootstrap_ts = bootstrap_ts == 1 & inclusion_vector' == 1; % Exclude periods of immobility from the traing set
    
    random_ts = ceil(rand*length(ca_time));
    shuffled_binarized = zeros(size(binarized_data));
    
    % Permute the trace
    shuffled_binarized(1:random_ts,:) = binarized_data(end-random_ts+1:end,:);
    shuffled_binarized(random_ts+1:end,:) = binarized_data(1:end-random_ts,:);
    
    for cell_i = 1:size(binarized_data,2)
        [~, MI_heading(shuffle_i,cell_i), ~, ~, ~, ~] = extract_1D_information(binarized_data(:,cell_i), interp_heading_vec, bin_vector, bootstrap_ts);
        [~, shuffled_MI_heading(shuffle_i,cell_i), ~, ~, ~, ~] = extract_1D_information(shuffled_binarized(:,cell_i), interp_heading_vec, bin_vector, bootstrap_ts);
    end
end

mean_MI_heading = mean(MI_heading);
SEM_MI_heading = std(MI_heading)/sqrt(numShuffles);
mean_MI_heading_shuffled = mean(shuffled_MI_heading);
SEM_MI_heading_shuffled = std(shuffled_MI_heading)/sqrt(numShuffles);

bootstrapped_data.mean_MI_HeadDirection = mean_MI_heading';
bootstrapped_data.SEM_MI_HeadDirection = SEM_MI_heading';
bootstrapped_data.mean_MI_HeadDirection_shuffled = mean_MI_heading_shuffled';
bootstrapped_data.SEM_MI_HeadDirection_shuffled = SEM_MI_heading_shuffled';

save('bootstrapped_data_HD.mat', 'bootstrapped_data');

end
