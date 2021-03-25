dt = median(diff(ms.time))/1000; % Conversion from ms to s
Fs = 1/dt;
z_threshold = 2;

[bFilt,aFilt] = butter(2,  2/(30/2), 'low');

raw_trace = ms.RawTraces(:,2);
filt_trace = zscore(filtfilt(bFilt,aFilt,raw_trace))';

%% Load the open field data
ca_data = ms.RawTraces;
ca_time = ms.time/1000;
behav_vec = behav.position;
behav_time = behav.time/1000;
%heading_vec = behavDLC.headDirection';

%% Parameters
sampling_frequency = 30;
training_set_portion = 0.5;
numShuffles = 30;
min_speed_threshold = 5; % cm.s-1

%% Clean data: keep only unique points

[behav_time, idx]=unique(behav_time); %getting rid of duplicate timestamps
behav_vec = behav_vec(idx,:); %excluding the positions associated with these double timestamps

[ca_time, IAms, ICms]=unique(ca_time);
ca_data = ca_data(IAms,:);

%behav_time = behav_time';


%% Heading
heading_vec = heading_vec';
heading_vec = heading_vec(idx);
[interp_heading_vec] = interpolate_HD(heading_vec, behav_time, ca_time);
interp_heading_vec = interp_heading_vec';

filt_HD = filtfilt(bFilt,aFilt,interp_heading_vec);

close all;
s = 3800
e = 6000
x1 = filt_trace(1,s:e);
y1 = filt_HD(1,s:e);
z1 = ms.Binary(s:e,79);

figure;
plot(x1) 
hold on;
plot(y1./360)
hold on 
plot(z1, 'LineWidth',2) 


%% Speed

interp_behav_vec(:,1) = interpolate_behavior(behav_vec(:,1), behav_time, ca_time); % in the X dimension
interp_behav_vec(:,2) = interpolate_behavior(behav_vec(:,2), behav_time, ca_time); % in the Y dimension
interp_behav_vec(end,:) = interp_behav_vec(end-1,:); % Make sure to interpolate or remove every NaN so that every timepoint has a corresponding behavioral state

close all;
s = 18000
e = 20000
x1 = filt_trace(1,s:e);

z1 = ms.Binary(s:e,29);

[velocity] = extract_velocity(interp_behav_vec, ca_time);

y1 = velocity(s:e);

figure;
plot(x1) 
hold on
plot(y1)
hold on 
plot(z1.*2, 'LineWidth',2) 

%% acceleration

interp_behav_vec(:,1) = interpolate_behavior(behav_vec(:,1), behav_time, ca_time); % in the X dimension
interp_behav_vec(:,2) = interpolate_behavior(behav_vec(:,2), behav_time, ca_time); % in the Y dimension
interp_behav_vec(end,:) = interp_behav_vec(end-1,:); % Make sure to interpolate or remove every NaN so that every timepoint has a corresponding behavioral state

[velocity] = extract_velocity(interp_behav_vec, ca_time);

acceleration = diff(velocity);
acceleration(end+1) = 0;


close all;
s = 5800
e = 7800
x1 = filt_trace(1,s:e);

z1 = ms.Binary(s:e,2);

y1 = acceleration(s:e);

figure;
plot(x1) 
hold on
plot(y1)
hold on 
plot(z1, 'LineWidth',2) 

acceleration_smooth = smooth(acceleration,round(1/mode(diff(ca_time))));


close all;
s = 5700
e = 6500
x1 = filt_trace(1,s:e);
z1 = ms.Binary(s:e,2);
y1 = acceleration_smooth(s:e);

figure;
plot(x1) 
hold on
plot(y1)
hold on 
plot(z1, 'LineWidth',2) 

