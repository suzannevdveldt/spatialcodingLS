%% Extracting firing probabilities (BETA)
%% Guillaume Etter 2017-2018

load('ms.mat');

%  load 'ca_data.mat'
%  load 'ca_time.mat'

if ~isfield(ms,'Binary')
   ms = msExtractBinary(ms); 
end

%% Here's the probability any cell is active
prob_being_active = sum(ms.Binary)./ms.numFrames;

%% Then you can look at the probability that a cell will be active on the next frame provided it was active on the current frame:

prob_transitioning_from_active_to_active = [];
prob_transitioning_from_active_to_inactive = [];
prob_transitioning_from_inactive_to_active = [];
prob_transitioning_from_inactive_to_inactive = []; 

for cell_i = 1:ms.numNeurons;
    following_state = [];
    spikes_idx = find(ms.Binary(:,cell_i) == 1);
    for spike_i = 1:length(spikes_idx);        
        if spikes_idx(spike_i)+1<ms.numFrames;
        following_state(spike_i) = ms.Binary((spikes_idx(spike_i)+1),cell_i);
        end
    following_inactive_state = []; 
    silent_idx = find(ms.Binary(:,cell_i) == 0);
    for silent_i = 1:length(silent_idx);
        if silent_idx(silent_i)+1<ms.numFrames;
        following_inactive_state(silent_i) = ms.Binary((silent_idx(silent_i)+1),cell_i);
        end
    end
end
    
prob_transitioning_from_active_to_active(cell_i,1) = sum(following_state)./length(spikes_idx); %dividing the number of active following states by total number of cell being active 
prob_transitioning_from_active_to_inactive(cell_i,1) = 1-prob_transitioning_from_active_to_active(cell_i);
prob_transitioning_from_inactive_to_active(cell_i,1) = sum(following_inactive_state)./length(silent_idx);
prob_transitioning_from_inactive_to_inactive(cell_i,1) = 1-prob_transitioning_from_inactive_to_active(cell_i); 
probability_being_active = prob_being_active';
end

firing_probabilities.prob_being_active = probability_being_active;
firing_probabilities.prob_transitioning_from_active_to_active = prob_transitioning_from_active_to_active;
firing_probabilities.prob_transitioning_from_active_to_silent = prob_transitioning_from_active_to_inactive;
firing_probabilities.prob_transitioning_from_inactive_to_active = prob_transitioning_from_inactive_to_active;
firing_probabilities.prob_transitioning_from_inactive_to_inactive = prob_transitioning_from_inactive_to_inactive;

save('firing_probabilities.mat', 'firing_probabilities'); 

