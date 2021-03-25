function [PW_matrix] = extract_pairwiseCorr_fullMatrix(mouse_session_folder, sig_threshold, only_place_cells, only_stable_cells)
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


cellIndex = {};

for day_i = 1:4
    for day_j = 1:4
            ct=1;
            for cell_i = 1:size(cell_to_index_map)
                if cell_to_index_map(cell_i,day_i)>0 && cell_to_index_map(cell_i,day_j)>0
                    if place_cells{day_i}(cell_to_index_map(cell_i,day_i)) == 1 && place_cells{day_j}(cell_to_index_map(cell_i,day_j)) == 1
                        cellIndex{day_i,day_j}(ct) = cell_i;
                        ct=ct+1;
                    end
                end
            end
    end
end

%% Re-arrange all tuning maps into linear vector

for day_i=1:4
    for day_j=1:4
        current_dayPairPWCorr = [];
%        if day_j<=day_i
        %Open tuning map here
        load([mouse_session_folder filesep strcat(mouse_number,'_',num2str(session_dates(day_i))) filesep 'tuning_map.mat'])
        tuning_map_A = tuning_map;
        
        load([mouse_session_folder filesep strcat(mouse_number,'_',num2str(session_dates(day_j))) filesep 'tuning_map.mat'])
        tuning_map_B = tuning_map;
        
        for cell_i = 1:length(cellIndex{day_i,day_j})
            cell2pick = cell_to_index_map(cellIndex{day_i,day_j}(cell_i),day_i);
            current_map = tuning_map_A(cell2pick).likelihood;
            current_map = smoothPlaceField(current_map);
            PV_A = current_map(:);

            cell2pick = cell_to_index_map(cellIndex{day_i,day_j}(cell_i),day_j);
            current_map = tuning_map_B(cell2pick).likelihood;
            current_map = smoothPlaceField(current_map);
            PV_B = current_map(:);

            %% PVcorr for day 1-2
            [current_dayPairPWCorr(cell_i), p_value] = corr(PV_A,PV_B);
            if ~isempty(sig_threshold)
                if p_value > sig_threshold
                    current_dayPairPWCorr(cell_i) = NaN;
                end
            end
        end
        
        PW_matrix(day_i,day_j) = mean(current_dayPairPWCorr,'omitnan');
        PW_matrix_SD(day_i,day_j) = std(current_dayPairPWCorr, 'omitnan');
        PW_matrix_N(day_i,day_j) = length(current_dayPairPWCorr); 
%       end
    end
end

save('PWCorr_fullMatrix.mat', 'PW_matrix');
save('PWCorr_fullMatrix_SD.mat', 'PW_matrix_SD'); 
save('PWCorr_fullMatrix_N.mat', 'PW_matrix_N');

%% Plot the results
labels = {'Day 1', 'Day 2','Day 3', 'Day 8'}
figure
imagesc(PW_matrix)
colormap Parula
daspect([1 1 1])
colorbar
set(gca,'xtick',[1:4],'ytick',[1:4],'xticklabel',labels,'yticklabels',labels,'CLim',[0 1])


end
