function [centroidShift, meanCentroidShift] = extract_centroidShift_8days(mouse_session_folder, only_place_cells, only_stable_cells)
%% mouse_session_folder: folder to analyze (contains folders for each session days)
% sig_threshold: correlation threshold above which correlation is rejected
% and transformed into NaN. Recommended threshold: 0.05

% Identify the mouse number
split_path = strsplit(mouse_session_folder,filesep);
mouse_number = split_path{end};

%% Make a list of folders in chronological order
session_folders = dir([mouse_session_folder filesep strcat('*',mouse_number,'*')]);
for file_i = 1:length(session_folders)
    split_session = strsplit(session_folders(file_i).name,'_');
    session_dates(file_i) = str2num(split_session{end});
end

[session_dates] = sort(session_dates,'ascend');

%% Consider only place cells
load([mouse_session_folder filesep strcat(mouse_number,'_',num2str(session_dates(1))) filesep 'SignCells.mat'])
place_cells{1} = SignCells;
load([mouse_session_folder filesep strcat(mouse_number,'_',num2str(session_dates(2))) filesep 'SignCells.mat'])
place_cells{2} = SignCells;
load([mouse_session_folder filesep strcat(mouse_number,'_',num2str(session_dates(3))) filesep 'SignCells.mat'])
place_cells{3} = SignCells;
load([mouse_session_folder filesep strcat(mouse_number,'_',num2str(session_dates(4))) filesep 'SignCells.mat'])
place_cells{4} = SignCells;

if ~only_place_cells % Include all cells
    place_cells{1} = ones(length(place_cells{1}),1);
    place_cells{2} = ones(length(place_cells{2}),1);
    place_cells{3} = ones(length(place_cells{3}),1);
    place_cells{4} = ones(length(place_cells{4}),1);
end

%% Consider only stable cells 
% 
% load([mouse_session_folder filesep strcat(mouse_number,'_',num2str(session_dates(1))) filesep 'split_half_stability.mat'])
% place_cells{1} = split_half_stability >= .5;
% load([mouse_session_folder filesep strcat(mouse_number,'_',num2str(session_dates(2))) filesep 'split_half_stability.mat'])
% place_cells{2} = split_half_stability >= .5;
% load([mouse_session_folder filesep strcat(mouse_number,'_',num2str(session_dates(3))) filesep 'split_half_stability.mat'])
% place_cells{3} = split_half_stability >= .5;
% load([mouse_session_folder filesep strcat(mouse_number,'_',num2str(session_dates(4))) filesep 'split_half_stability.mat'])
% place_cells{4} = split_half_stability >= .5;
% 
% if ~only_stable_cells % Include all cells
%     place_cells{1} = ones(length(place_cells{1}),1);
%     place_cells{2} = ones(length(place_cells{2}),1);
%     place_cells{3} = ones(length(place_cells{3}),1);
%     place_cells{4} = ones(length(place_cells{4}),1);
% end

%% Open the cell_index file and register cell pairs present on different days
load([mouse_session_folder filesep 'CellReg' filesep 'cell_to_index_map.mat'])

cellIndex = [];
for day = 1:4
    ct=1;
    for cell_i = 1:size(cell_to_index_map)
       if cell_to_index_map(cell_i,1) > 0 && cell_to_index_map(cell_i,day) > 0
            if place_cells{1}(cell_to_index_map(cell_i,1)) == 1 && place_cells{day}(cell_to_index_map(cell_i,day)) == 1
                cellIndex{day}(ct) = cell_i;
                ct=ct+1;
            end
       end
    end
end

%% Compute centroid distances here
meanCentroidShift=[];
centroidShift={};
for day = 2:4
    
    % Reference day
    load([mouse_session_folder filesep strcat(mouse_number,'_',num2str(session_dates(1))) filesep 'tuning_map.mat'])
    load([mouse_session_folder filesep strcat(mouse_number,'_',num2str(session_dates(1))) filesep 'place_field.mat'])
    for cell_i = 1:length(cellIndex{day})
            cell2pick = cell_to_index_map(cellIndex{day}(cell_i),1);
            current_map = tuning_map(cell2pick).likelihood;
            current_map = smoothPlaceField(current_map);
            centroidRef(cell_i,:) = place_field.centroid(cell_i,:);
%             [~,peakX] = max(max(current_map,[],1,'omitnan'));
%             [~,peakY] = max(max(current_map,[],2,'omitnan'));
%             peakX=peakX*3; % Convert to cm
%             peakY=peakY*3; % Convert to cm
%             
%             centroidRef(cell_i,1) = peakX;
%             centroidRef(cell_i,2) = peakY;
    end

    load([mouse_session_folder filesep strcat(mouse_number,'_',num2str(session_dates(day))) filesep 'tuning_map.mat'])
     load([mouse_session_folder filesep strcat(mouse_number,'_',num2str(session_dates(day))) filesep 'place_field.mat'])
    for cell_i = 1:length(cellIndex{day})
            cell2pick = cell_to_index_map(cellIndex{day}(cell_i),day);
            current_map = tuning_map(cell2pick).likelihood;
            current_map = smoothPlaceField(current_map);
            centroidComp(cell_i,:) = place_field.centroid(cell_i,:);
%             [~,peakX] = max(max(current_map,[],1,'omitnan'));
%             [~,peakY] = max(max(current_map,[],2,'omitnan'));
%             peakX=peakX*3; % Convert to cm
%             peakY=peakY*3; % Convert to cm
            
%             centroidComp(cell_i,1) = peakX;
%             centroidComp(cell_i,2) = peakY;
    end
    
    %% Compute centroid distance
    for cell_i = 1:length(cellIndex{day})
        centroidShift{day-1}(cell_i) = sqrt((centroidComp(cell_i,1)-centroidRef(cell_i,1)).^2 + (centroidComp(cell_i,2)-centroidRef(cell_i,2)).^2);
    end
    meanCentroidShift(day-1)=mean(centroidShift{day-1},'omitnan');
    
end

save('centroidShift.mat', 'centroidShift', 'meanCentroidShift'); 

end
