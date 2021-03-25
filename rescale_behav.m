

a = min(behav.position(:,1))
b = max(behav.position(:,1))

c = b-a

scale = 45./ c

behav.position(:,1) = behav.position(:,1) .* scale 


 a = min(behav.position(:,2))
 b = max(behav.position(:,2))
 
 c = b-a
 
 scale = 45./ c

%behav_vec_rescaled(:,2) = behav_vec(:,2) .* scale 
behav.position(:,2) = behav.position(:,2) .* scale 

% a = min(behav.positionRescale(:,1))
% b = max(behav.positionRescale(:,1))
% a = min(behav.positionRescale(:,2))
% b = max(behav.positionRescale(:,2))

%%
% 

 a = min(behav_vec(:,1))
 b = max(behav_vec(:,1))

c = b-a

scale = 45./ c
  behav_vec_rescaled(:,1) = behav_vec(:,1) .* scale 

%  a = min(behav_vec(:,2))
%  b = max(behav_vec(:,2))
 behav_vec_rescaled(:,2) = behav_vec(:,2) .* scale 


a = min(behav_vec_rescaled(:,1))
b = max(behav_vec_rescaled(:,1))
a = min(behav_vec_rescaled(:,2))
b = max(behav_vec_rescaled(:,2))
