function [place_cell_data] = extract_place_cell_informationoverdays;


topLevelFolder = uigetdir;
cd (topLevelFolder);
matrix = zeros(15);

load('StabilityPlaceCellsZmap.mat');
load('alligned_cellsfirst2days.mat');

nsessions = size(alligned_cellsfirst2days,2);


for iteration = 1:nsessions;
    disp('Select directory of session to be loaded');
    directory_name = uigetdir;
    directories{iteration} = directory_name
end


for trial = 1:length(alligned_cellsfirst2days);
    
    iteration = 1;
    for iteration = 1:nsessions;
    directory_name = directories{iteration};
    cd (directory_name);
    load('ms.mat')
    load('behav.mat')
    try
        cell_j = alligned_cellsfirst2days(trial,iteration);
        [place_cell_data] = msSpatialFiringV2( ms, behav, cell_j,  1, true);
        %msSpatialFiringV4( ms, behav, cell_j, true);
    catch 
        fprintf('no cell recorded');
    end 
    
    end
    pause
    trial = trial +1;
    close all
end
   
cd (topLevelFolder);

%% safe the important information in a variable

%alligned_cell_analysis.isPlaceCell 
%alligned_cell_analysis.PlaceFieldCorr
%alligned_cell_analysis.Stability

%% this originally plotted the spatial footprints for each cell 
% folder = dir('*.mat'); 
% 
% for i = 1: length(folder)               %Going through every object in the folder to extract index maping 
%     if(contains(folder(i).name,'cellRegistered'))
%         cellIndex = load(folder(i).name,'cell_registered_struct');          %extracting the cell index matrix across sessions
%         cellIndex = cellIndex.cell_registered_struct.cell_to_index_map;     %Saving the matrix within the function
%         load(folder(i).name,'cell_registered_struct');  
%     end
% end
% 
% 
% cell_row = find(cell_registered_struct.cell_to_index_map(:,1) == cell{1,1})
% 
% 
% figure; 
% set(gcf,'Position',[100 100 2000 500])
%     subplot(1,nsessions+1,1)
%     imagesc(ms.CorrProj);
%     daspect([1 1 1]);
%     ax2=gca;
%     colormap(ax2, gray);
%     title('Correlation');
%     
%     for session = 1:nsessions;
%     if session == 1;
%     subplot(1,nsessions+1,session+1);
%     cell_id = cell_registered_struct.cell_to_index_map(cell_row,session);
%     imagesc(permute(cell_registered_struct.spatial_footprints_corrected{session,1}(cell_id,:,:),[2 3 1]))
%     hold on
%     else 
%     subplot(1,nsessions+1,session+1);
%     cell_id_registered = cell_registered_struct.cell_to_index_map(cell_row,session);
%         if cell_id_registered ~= 0;
%         subplot(1,nsessions+1,session+1);
%         imagesc(permute(cell_registered_struct.spatial_footprints_corrected{session,1}(cell_id_registered,:,:),[2 3 1]))
%         else 
%             continue 
%         end
%         
%     end
%     end 
%     
%     pause
%     close all
% end
% 
%     
%     cd (topLevelFolder);
end

    




%% Plot results
%  figure
%  subplot(4,1,1);
%  bar(Cell_ID,Stability);
%  title 'Place Field Stability'
%  
%  subplot(4,1,2)
%  bar(Cell_ID,InFieldActivityRatio);
%  title 'In Field Activity Ratio'
%  
%  subplot(4,1,3);
%  bar(Cell_ID,Place_cell);
%  title 'Place cell'
%  
%  subplot(4,1,4);
%  bar(Cell_ID,MI);
%  title 'MI value'


%%

% 
% summary(:,1) = Cell_ID;
% summary(:,2) = Place_cell;
% summary(:,3) = Stability;
% summary(:,4) = InFieldActivityRatio;
% summary(:,5) = FiringProb;
% summary(:,6) = PlaceFieldArea;
% summary(:,7) = MI;
% 
% save('AllSpatialFiringData.mat', 'AllSpatialFiringData');
% save('PlaceField.mat', 'PlaceField');
% save('summary.mat', 'summary');
