%suzanne script

% bins behavioral data and displays how long the mouse was in each bin. 
%You can also see the total distance the mouse travel in one session. 

%--------------USER INPUTS--------------

binSize = 3; %in cm 

%---------------------------------------
load('ms.mat'); 
load('behav'); 
%finds the trajectory of the mouse: 
[unique_behavtime, IAbehav, ICbehav]=unique(behav.time);
[unique_mstime, IAms, ICms]=unique(ms.time);

% Interpolate the X position of mouse. 
interpolated_X=interp1(behav.time(IAbehav),behav.position(IAbehav,1),ms.time(IAms),'pchip');
interpolated_Y=interp1(behav.time(IAbehav),behav.position(IAbehav,2),ms.time(IAms),'pchip');

%the length of the open field
length_of_field_x = round(max(interpolated_X)/binSize)*binSize;
length_of_field_y = round(max(interpolated_Y)/binSize)*binSize;
bins_x = length_of_field_x/binSize;  
bins_y = length_of_field_y/binSize; 

%now that we know the size of the bins etc, we can calculate how much time
%the mouse spends in each: 

%if we know that each frame is 30F/s ~ 1F = 0.0333sec
%we count each frame that the mouse spends in each bin

bin = zeros(bins_y, bins_x); 
binLocation = strings(10,10); 
for i=1: bins_x
    for j = 1: bins_y        
        x = find((i-1)*binSize < interpolated_X & interpolated_X < i*binSize );
        y = find((j-1)*binSize < interpolated_Y  & interpolated_Y < j*binSize );
        frames = intersect(x,y);
        bin(j,i) = length(frames);
        string = sprintf("%d,%d",i, j); 
        binLocations(j,i) = string; 
    end 
end 

j = 1; 
i = 1; 
for binNum =1: bins_y*bins_x      
    LinearBin(binNum) = bin(j,i);
    LinearBinLocations(binNum) = binLocations(j,i);
    i = i+1;
    if binNum == bins_x*j
        i = 1;
        j = j+1; 
    end  
end 


%converting time in bin to probability of mouse in bin 
NormalizedTimeInBin = (LinearBin(1,:) ./ length(ms.time)) .* 100;

%now to convert frames to seconds:
LinearBin = LinearBin./30; 

variance = var(NormalizedTimeInBin);

% %and bar graphing that: 
% figure
% subplot(3,2,[1,2]);
% bar(LinearBin)
% set(gca, 'xticklabel', LinearBinLocations)
% title('The amount of time the mouse spends in each bin')
% xlabel('Bin Number')
% ylabel('Time in Seconds')
% 
% subplot(3,2,[3,4]);
% bar(NormalizedTimeInBin)
% set(gca, 'xticklabel', LinearBinLocations)
% title('Percentage of time animal spends in bin')
% xlabel('Bin Number')
% ylabel('Percentage of Time')
% 
% subplot(3,2,5);
% plot(interpolated_X,interpolated_Y)
% title('Lil mouse running around lil open field'); 
% xlabel('length (cm)');
% ylabel('width (cm)');

% Calculating the total distance the mouse traveled: 
distance =0;
for i=2: length(interpolated_X)    
    distance = sqrt((interpolated_X(i) - interpolated_X(i-1))^2 + (interpolated_Y(i) - interpolated_Y(i-1))^2) + distance ;    
end 

% calculate mean speed
meanspeed = mean(behav.speed);

behavior.NormalizedTimeInBin = NormalizedTimeInBin;
behavior.distance = distance;
behavior.variance = variance;
behavior.meanspeed = meanspeed; 

save('behavior.mat', 'behavior');

 