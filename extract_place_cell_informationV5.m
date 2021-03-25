function [place_cell_data] = Basic_Calcium_Analysis(ms,behav);

%% You can run this function by running the function select_folders_for_analysis to automatically 
% select the folders you want to look at.
% This function call functions:
    % msExtractTrans
    % msSpatialFiringV2  (uses stability threshod .5, threshold .8, infield activity threshold 3) 
    % msSpatialFiringV3  (uses place field threshold z = 2, 
                             %place field variance .01, 200 shuffles and 1000 iterations )
    % msCalciumRunningSpeed 
                             
%%                             
clear all; 
 
 load('ms.mat');
 load('behav.mat');

%% Cell characteristics analysis 
 cell_firing_characteristics = msExtractTrans(ms);
 save('cell_firing_characteristics.mat', 'cell_firing_characteristics');
 
 %% Place Cell Analysis 
for cell_i = 1:ms.numNeurons;
     % Determine possible place cells using stability, in field activity
     % ratio and firing using 
    [place_cell_data] = msSpatialFiringV2( ms, behav, cell_i,  0, 1);
    AllSpatialFiringData.V2(cell_i) = place_cell_data;
    Cell_ID(cell_i) = place_cell_data.Cell_ID;
    Place_cellV2(cell_i) = place_cell_data.IsPlaceCell;
    FiringProb(cell_i) = place_cell_data.CellFiringProb;
    SplitHalfStability(cell_i) = place_cell_data.PlaceFieldStability;
    PlaceFieldArea(cell_i) = place_cell_data.PlaceFieldArea;
    InFieldActivityRatio(cell_i) = place_cell_data.InFieldActivityRatio;
    PlaceField{cell_i} = place_cell_data.PlaceFieldMask;
    PlaceFieldCentroid = place_cell_data.PlaceFieldCentroid;
    PlaceFieldCentroidX(cell_i) = PlaceFieldCentroid(1);
    if isnan(PlaceFieldCentroidX(cell_i));
         PlaceFieldCentroidY(cell_i) = NaN;
    else
           PlaceFieldCentroidY(cell_i) = PlaceFieldCentroid(2);
    end

    % Confirm place cell activity is 2sd above chance, and the average percent variance >.05. 
    [place_field_data] = msSpatialFiringV4(ms,behav,cell_i, 1);
    AllSpatialFiringData.V4(cell_i) = place_field_data;
    Place_cellV4(cell_i) = place_field_data.IsPlaceCell;
    PlaceFieldMaskV4{cell_i} = place_field_data.PlaceFieldMask;  
    Zmap{cell_i} = place_field_data.Z_map;
    MeanPercentVariance{cell_i} = place_field_data.MeanPercentVariance;
   %pause
   
    if Place_cellV2(cell_i) == 1 & Place_cellV4(cell_i) == 1; 
            Place_cell(cell_i) = 1;
        else 
            Place_cell(cell_i) = 0; 
    end
     
    mean_rise_time(cell_i) = cell_firing_characteristics{1,cell_i}.mean_rise_time; 
    mean_decay_time(cell_i) = cell_firing_characteristics{1,cell_i}.mean_decay_time;
    mean_prominence(cell_i) = cell_firing_characteristics{1,cell_i}.mean_prominence;
    mean_width(cell_i) =  cell_firing_characteristics{1, cell_i}.mean_width;
    mean_rise_area(cell_i) = cell_firing_characteristics{1, cell_i}.mean_rise_area;
    mean_decay_area(cell_i) = cell_firing_characteristics{1,cell_i}.mean_decay_area;
   
    close all 
end

save('AllSpatialFiringData.mat', 'AllSpatialFiringData');

%% Analyze the behavior of the session   
[behavioral_data] = msBehavioralAnalysis(ms,behav);
save('behavioral_data.mat', 'behavioral_data');

%% Speed Correlation and Acceleration / Decelleration 

% I think I would use msCalciumRunningSpeed here but double check 


%% Save results in a file that can be easily accessed 

summary(:,1) = Cell_ID;
summary(:,2) = Place_cell;
summary(:,3) = SplitHalfStability;
summary(:,4) = InFieldActivityRatio;
summary(:,5) = FiringProb;
summary(:,6) = PlaceFieldArea;
summary(:,7) = PlaceFieldCentroidX;
summary(:,8) = PlaceFieldCentroidY;
summary(:,9) = mean_rise_time;
summary(:,10) = mean_decay_time;
summary(:,11) = mean_prominence;
summary(:,12) = mean_width;
summary(:,13) = mean_rise_area;
summary(:,14) = mean_decay_area;

save('summary.mat', 'summary');
 
clear all
close all
end
