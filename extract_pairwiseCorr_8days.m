function [PWCorr1_8, PWCorr2_8, PWCorr3_8] = extract_pairwiseCorr_8days(mouse_session_folder, sig_threshold, only_place_cells, only_stable_cells)
%% mouse_session_folder: folder to analyze (contains folders for each session days)
% sig_threshold: correlation threshold above which correlation is rejected
% and transformed into NaN. Recommended threshold: 0.05

nshuffles = 30;

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

cellIndex1_2 = [];
cellIndex2_3 = []; 
cellIndex1_8 = [];
cellIndex2_8 = [];
cellIndex3_8 = [];

for cell_i = 1:size(cell_to_index_map)
    % cellIndex1_2
   if cell_to_index_map(cell_i,1) > 0 && cell_to_index_map(cell_i,2) > 0
        if place_cells{1}(cell_to_index_map(cell_i,1)) == 1 || place_cells{2}(cell_to_index_map(cell_i,2)) == 1
            cellIndex1_2(end+1,:) = cell_i;
        end
   end
   %cellIndex2_3
       if cell_to_index_map(cell_i,2) > 0 && cell_to_index_map(cell_i,3) > 0
        if place_cells{2}(cell_to_index_map(cell_i,2)) == 1 || place_cells{3}(cell_to_index_map(cell_i,3)) == 1
            cellIndex2_3(end+1,:) = cell_i;
        end
    end
    % cellIndex1_8
    if cell_to_index_map(cell_i,1) > 0 && cell_to_index_map(cell_i,4) > 0
        if place_cells{1}(cell_to_index_map(cell_i,1)) == 1 || place_cells{4}(cell_to_index_map(cell_i,4)) == 1
            cellIndex1_8(end+1,:) = cell_i;
        end
    end
     % cellIndex2_8
    if cell_to_index_map(cell_i,2) > 0 && cell_to_index_map(cell_i,4) > 0
        if place_cells{2}(cell_to_index_map(cell_i,2)) == 1 || place_cells{4}(cell_to_index_map(cell_i,4)) == 1
            cellIndex2_8(end+1,:) = cell_i;
        end
    end
      % cellIndex3_8
    if cell_to_index_map(cell_i,3) > 0 && cell_to_index_map(cell_i,4) > 0
        if place_cells{3}(cell_to_index_map(cell_i,3)) == 1 || place_cells{4}(cell_to_index_map(cell_i,4)) == 1
            cellIndex3_8(end+1,:) = cell_i;
        end
    end
end

%% Re-arrange all tuning maps into linear vector

%% Day 1-2
if ~isempty(cellIndex1_2)
    %Day 1
    day = 1;
    
    %Open tuning map here
    load([mouse_session_folder filesep strcat(mouse_number,'_',num2str(session_dates(day))) filesep 'tuning_map.mat'])
    
    % Compute current day population vector
    PV1 = [];
    for cell_i = 1:length(cellIndex1_2)
        cell2pick = cell_to_index_map(cellIndex1_2(cell_i),day);
        current_map = tuning_map(cell2pick).likelihood;
        current_map = smoothPlaceField(current_map);
        PV1(:,cell_i) = current_map(:);
    end
    
    % Day 2
    day = 2;
    
    % Open tuning map here
    load([mouse_session_folder filesep strcat(mouse_number,'_',num2str(session_dates(day))) filesep 'tuning_map.mat'])
    
    % Compute current day population vector
    PV2 = [];
    for cell_i = 1:length(cellIndex1_2)
        cell2pick = cell_to_index_map(cellIndex1_2(cell_i),day);
        current_map = tuning_map(cell2pick).likelihood;
        current_map = smoothPlaceField(current_map);
        PV2(:,cell_i) = current_map(:);
    end
    
    %% PVcorr for day 1-2
    for cell_i = 1:length(cellIndex1_2)
        [PWCorr1_2(cell_i), p_value] = corr(PV1(:,cell_i),PV2(:,cell_i));
        if ~isempty(sig_threshold)
            if p_value > sig_threshold
                PWCorr1_2(cell_i) = NaN;
            end
        end
    end
    
     %% Compute shuffled
    for shuffle = 1:nshuffles;
       Shuffled_CellIndex1_2 = cellIndex1_2(randperm(length(cellIndex1_2)));
        
    PV2 = [];
     for cell_i = 1:length(Shuffled_CellIndex1_2)
         cell2pick = cell_to_index_map(Shuffled_CellIndex1_2(cell_i),day);
         current_map = tuning_map(cell2pick).likelihood;
        current_map = smoothPlaceField(current_map);
         PV2(:,cell_i) = current_map(:);
     end
     
     %% PVcorr for day 1-2
     for cell_i = 1:length(cellIndex1_2)
         [Shuffled_PWCorr1_2(cell_i,shuffle), p_value] = corr(PV1(:,cell_i),PV2(:,cell_i));
         if ~isempty(sig_threshold)
              if p_value > sig_threshold
                  Shuffled_PWCorr1_2(cell_i,shuffle) = NaN;
              end
         end
     end
    end
    
    Shuffled_PWCorr1_2 = mean(Shuffled_PWCorr1_2,2);
else
    PWCorr1_2 = NaN;
end



%% Day 2-3
if ~isempty(cellIndex2_3)
    %Day 1
    day = 2;
    
    %Open tuning map here
    load([mouse_session_folder filesep strcat(mouse_number,'_',num2str(session_dates(day))) filesep 'tuning_map.mat'])
    
    % Compute current day population vector
    PV2 = [];
    for cell_i = 1:length(cellIndex2_3)
        cell2pick = cell_to_index_map(cellIndex2_3(cell_i),day);
        current_map = tuning_map(cell2pick).likelihood;
        current_map = smoothPlaceField(current_map);
        PV2(:,cell_i) = current_map(:);
    end
    
    % Day 3
    day = 3;
    
    % Open tuning map here
    load([mouse_session_folder filesep strcat(mouse_number,'_',num2str(session_dates(day))) filesep 'tuning_map.mat'])
    
    % Compute current day population vector
    PV3 = [];
    for cell_i = 1:length(cellIndex2_3)
        cell2pick = cell_to_index_map(cellIndex2_3(cell_i),day);
        current_map = tuning_map(cell2pick).likelihood;
        current_map = smoothPlaceField(current_map);
        PV3(:,cell_i) = current_map(:);
    end
    
    %% PVcorr for day 2-3
    for cell_i = 1:length(cellIndex2_3)
        [PWCorr2_3(cell_i), p_value] = corr(PV2(:,cell_i),PV3(:,cell_i));
        if ~isempty(sig_threshold)
            if p_value > sig_threshold
                PWCorr2_3(cell_i) = NaN;
            end
        end
    end
    
    
     %% Compute shuffled
    for shuffle = 1:nshuffles;
       Shuffled_CellIndex2_3 = cellIndex2_3(randperm(length(cellIndex2_3)));
        
    PV3 = [];
     for cell_i = 1:length(Shuffled_CellIndex2_3)
         cell2pick = cell_to_index_map(Shuffled_CellIndex2_3(cell_i),day);
         current_map = tuning_map(cell2pick).likelihood;
        current_map = smoothPlaceField(current_map);
         PV3(:,cell_i) = current_map(:);
     end
     
     %% PVcorr for day 2-3
     for cell_i = 1:length(cellIndex2_3)
         [Shuffled_PWCorr2_3(cell_i,shuffle), p_value] = corr(PV2(:,cell_i),PV3(:,cell_i));
         if ~isempty(sig_threshold)
              if p_value > sig_threshold
                  Shuffled_PWCorr2_3(cell_i,shuffle) = NaN;
              end
         end
     end
    end
    
    Shuffled_PWCorr2_3 = mean(Shuffled_PWCorr2_3,2);
    
else
    PWCorr2_3 = NaN;
end

%% Day 1-8
if ~isempty(cellIndex1_8)
    %Day 1
    day = 1;
    
    %Open tuning map here
    load([mouse_session_folder filesep strcat(mouse_number,'_',num2str(session_dates(day))) filesep 'tuning_map.mat'])
    
    % Compute current day population vector
    PV1 = [];
    for cell_i = 1:length(cellIndex1_8)
        cell2pick = cell_to_index_map(cellIndex1_8(cell_i),day);
        current_map = tuning_map(cell2pick).likelihood;
        current_map = smoothPlaceField(current_map);
        PV1(:,cell_i) = current_map(:);
    end
    
    % Day 8
    day = 4;
    
    % Open tuning map here
    load([mouse_session_folder filesep strcat(mouse_number,'_',num2str(session_dates(day))) filesep 'tuning_map.mat'])
    
    % Compute current day population vector
    PV8 = [];
    for cell_i = 1:length(cellIndex1_8)
        cell2pick = cell_to_index_map(cellIndex1_8(cell_i),day);
        current_map = tuning_map(cell2pick).likelihood;
        current_map = smoothPlaceField(current_map);
        PV8(:,cell_i) = current_map(:);
    end
    
    %% PVcorr for day 1-8
    for cell_i = 1:length(cellIndex1_8)
        [PWCorr1_8(cell_i), p_value] = corr(PV1(:,cell_i),PV8(:,cell_i));
        if ~isempty(sig_threshold)
            if p_value > sig_threshold
                PWCorr1_8(cell_i) = NaN;
            end
        end
    end
    
    
    % Compute shuffled
    for shuffle = 1:nshuffles;
       Shuffled_CellIndex1_8 = cellIndex1_8(randperm(length(cellIndex1_8)));
        
    PV8 = [];
     for cell_i = 1:length(Shuffled_CellIndex1_8)
         cell2pick = cell_to_index_map(Shuffled_CellIndex1_8(cell_i),day);
         current_map = tuning_map(cell2pick).likelihood;
        current_map = smoothPlaceField(current_map);
         PV8(:,cell_i) = current_map(:);
     end
     
     %% PVcorr for day 1-2
     for cell_i = 1:length(cellIndex1_8)
         [Shuffled_PWCorr1_8(cell_i,shuffle), p_value] = corr(PV1(:,cell_i),PV8(:,cell_i));
         if ~isempty(sig_threshold)
              if p_value > sig_threshold
                  Shuffled_PWCorr1_8(cell_i,shuffle) = NaN;
              end
         end
     end
    end
    
    Shuffled_PWCorr1_8 = mean(Shuffled_PWCorr1_8,2);
else
    PWCorr1_8 = NaN;
end

%% Day 2-8
if ~isempty(cellIndex2_8)
    %Day 2
    day = 2;
    %Open tuning map here
    load([mouse_session_folder filesep strcat(mouse_number,'_',num2str(session_dates(day))) filesep 'tuning_map.mat'])
    
    % Compute current day population vector
    PV2 = [];
    for cell_i = 1:length(cellIndex2_8)
        cell2pick = cell_to_index_map(cellIndex2_8(cell_i),day);
        current_map = tuning_map(cell2pick).likelihood;
        current_map = smoothPlaceField(current_map);
        PV2(:,cell_i) = current_map(:);
    end
    
    % Day 8
    day = 4;
    
    % Open tuning map here
    load([mouse_session_folder filesep strcat(mouse_number,'_',num2str(session_dates(day))) filesep 'tuning_map.mat'])
    
    % Compute current day population vector
    PV8 = [];
    for cell_i = 1:length(cellIndex2_8)
        cell2pick = cell_to_index_map(cellIndex2_8(cell_i),day);
        current_map = tuning_map(cell2pick).likelihood;
        current_map = smoothPlaceField(current_map);
        PV8(:,cell_i) = current_map(:);
    end
    
    %% PVcorr for day 2-8
    for cell_i = 1:length(cellIndex2_8)
        [PWCorr2_8(cell_i), p_value] = corr(PV2(:,cell_i),PV8(:,cell_i));
        if ~isempty(sig_threshold)
            if p_value > sig_threshold
                PWCorr2_8(cell_i) = NaN;
            end
        end
    end
    
    
    % Compute shuffled
    for shuffle = 1:nshuffles;
       Shuffled_CellIndex2_8 = cellIndex2_8(randperm(length(cellIndex2_8)));
        
    PV8 = [];
     for cell_i = 1:length(Shuffled_CellIndex2_8)
         cell2pick = cell_to_index_map(Shuffled_CellIndex2_8(cell_i),day);
         current_map = tuning_map(cell2pick).likelihood;
        current_map = smoothPlaceField(current_map);
         PV8(:,cell_i) = current_map(:);
     end
     
     %% PVcorr
     for cell_i = 1:length(cellIndex2_8)
         [Shuffled_PWCorr2_8(cell_i,shuffle), p_value] = corr(PV2(:,cell_i),PV8(:,cell_i));
         if ~isempty(sig_threshold)
              if p_value > sig_threshold
                  Shuffled_PWCorr2_8(cell_i,shuffle) = NaN;
              end
         end
     end
    end
    
    Shuffled_PWCorr2_8 = mean(Shuffled_PWCorr2_8,2);
    
else
    PWCorr2_8 = NaN;
end

%% Day 3-8
if ~isempty(cellIndex3_8)
    %Day 1
    day = 3;
    
    %Open tuning map here
    load([mouse_session_folder filesep strcat(mouse_number,'_',num2str(session_dates(day))) filesep 'tuning_map.mat'])
    
    % Compute current day population vector
    PV3 = [];
    for cell_i = 1:length(cellIndex3_8)
        cell2pick = cell_to_index_map(cellIndex3_8(cell_i),day);
        current_map = tuning_map(cell2pick).likelihood;
        current_map = smoothPlaceField(current_map);
        PV3(:,cell_i) = current_map(:);
    end
    
    % Day 8
    day = 4;
    
    % Open tuning map here
    load([mouse_session_folder filesep strcat(mouse_number,'_',num2str(session_dates(day))) filesep 'tuning_map.mat'])
    
    % Compute current day population vector
    PV8 = [];
    for cell_i = 1:length(cellIndex3_8)
        cell2pick = cell_to_index_map(cellIndex3_8(cell_i),day);
        current_map = tuning_map(cell2pick).likelihood;
        current_map = smoothPlaceField(current_map);
        PV8(:,cell_i) = current_map(:);
    end
    
    %% PVcorr for day 3-8
    for cell_i = 1:length(cellIndex3_8)
        [PWCorr3_8(cell_i), p_value] = corr(PV3(:,cell_i),PV8(:,cell_i));
        if ~isempty(sig_threshold)
            if p_value > sig_threshold
                PWCorr3_8(cell_i) = NaN;
            end
        end
    end
    % Compute shuffled
    for shuffle = 1:nshuffles;
       Shuffled_CellIndex3_8 = cellIndex3_8(randperm(length(cellIndex3_8)));
        
    PV8 = [];
     for cell_i = 1:length(Shuffled_CellIndex3_8)
         cell2pick = cell_to_index_map(Shuffled_CellIndex3_8(cell_i),day);
         current_map = tuning_map(cell2pick).likelihood;
        current_map = smoothPlaceField(current_map);
         PV8(:,cell_i) = current_map(:);
     end
     
     %% PVcorr for day 3-8
     for cell_i = 1:length(cellIndex3_8)
         [Shuffled_PWCorr3_8(cell_i,shuffle), p_value] = corr(PV3(:,cell_i),PV8(:,cell_i));
         if ~isempty(sig_threshold)
              if p_value > sig_threshold
                  Shuffled_PWCorr3_8(cell_i,shuffle) = NaN;
              end
         end
     end
    end
    
    Shuffled_PWCorr3_8 = mean(Shuffled_PWCorr3_8,2);
else
    PWCorr3_8 = NaN;
end

PWCorr.PWCorr1_2 = PWCorr1_2';
PWCorr.PWCorr2_3 = PWCorr2_3';
PWCorr.PWCorr1_8 = PWCorr1_8';
PWCorr.PWCorr2_8 = PWCorr2_8';
PWCorr.PWCorr3_8 = PWCorr3_8';

PWCorr.Shuffled_PWCorr1_2 = Shuffled_PWCorr1_2;
PWCorr.Shuffled_PWCorr2_3 = Shuffled_PWCorr2_3;
PWCorr.Shuffled_PWCorr1_8 = Shuffled_PWCorr1_8;
PWCorr.Shuffled_PWCorr2_8 = Shuffled_PWCorr2_8;
PWCorr.Shuffled_PWCorr3_8 = Shuffled_PWCorr3_8;


save('PWCorr.mat', 'PWCorr'); 

end
