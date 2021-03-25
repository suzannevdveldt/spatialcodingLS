function [StabilityCellsOverDays] = msStabilityCellsOverDays;


%% smoothing parameters 
bin_size = 3; %cm
kernel_size = [bin_size bin_size];
occupancy_std = 2;

[Xgrid,Ygrid]=meshgrid(-kernel_size(1)/2: kernel_size(1)/2, -kernel_size(2)/2:kernel_size(2)/2);
Rgrid=sqrt((Xgrid.^2+Ygrid.^2));
kernel = pdf('Normal', Rgrid, 0, occupancy_std);
kernel = kernel./sum(sum(kernel));    

%% inputting the folders 
disp('Input folder with your CellAllignment.mat file:   ');
topLevelFolder = uigetdir;
cd (topLevelFolder)
load('CellAllignment.mat'); 

%% load the directories of all files of interest; 

prompt = 'Total number of open field recordings analyzed =       ';
nsessions = input(prompt);

prompt = 'Input session code: 3 if only A1-A3, %4 if until B1, 8 if A1-B3, 13 if until C3  =       ';
sessioncode = input(prompt);

for loadedSession = 1:nsessions;
    loadedSession
    disp('Select directory of session to be loaded');
    directory_name = uigetdir;
    directories{loadedSession} = directory_name;
end 

%% Load the cellIndex of the cells that you alligned 
for i = 1:sessioncode; %3 if only A1-A3, %4 if until B1, 8 if A1-B3, 13 if until C3
    if i == 1;   
        cellIndex = CellAllignment.A1A2.CellIDs;
        firstsession = 1;
        lastsession = 2;
    elseif i == 2;   %CellAllignment.A2A3 = session_alligned_cells;
        cellIndex = CellAllignment.A2A3.CellIDs;
        firstsession = 2;
        lastsession = 3; 
    elseif i == 3;   %CellAllignment.A1A3 = session_alligned_cells;
        cellIndex = CellAllignment.A1A3.CellIDs;
        firstsession = 1;
        lastsession = 3;
%     elseif i == 4;   %CellAllignment.A3B1 = session_alligned_cells;
%         cellIndex = CellAllignment.A3B1.CellIDs;
%         firstsession = 3;
%         lastsession = 4;
%      elseif i == 5;   %CellAllignment.A1B1 = session_alligned_cells;
%          cellIndex = CellAllignment.A1B1.CellIDs;
%          firstsession = 1;
%          lastsession = 4;
%        
%     elseif i == 6;   %CellAllignment.B2B3 = session_alligned_cells;
%         cellIndex = CellAllignment.B2B3.CellIDs;
%         firstsession = 5;
%         lastsession = 6;
%     elseif i == 7;   %CellAllignment.B1B3 = session_alligned_cells;
%         cellIndex = CellAllignment.B1B3.CellIDs;
%         firstsession = 4;
%         lastsession = 6; 
%     elseif i == 8;   %CellAllignment.A1B3 = session_alligned_cells;
%         cellIndex = CellAllignment.A1B3.CellIDs;
%         firstsession = 1;
%         lastsession = 6;
%  elseif i == 9;   %CellAllignment.B1B2 = session_alligned_cells;
%          cellIndex = CellAllignment.B1B2.CellIDs;
%          firstsession = 4;
%          lastsession = 5;
%          
%  %   elseif i == 9;   %CellAllignment.B3C1 = session_alligned_cells;
%  %       cellIndex = CellAllignment.B3C1.CellIDs;
%  %       firstsession = 6;
%  %       lastsession = 7;
%     elseif i == 10;  %CellAllignment.C1C2 = session_alligned_cells;
%         cellIndex = CellAllignment.C1C2.CellIDs;
%         firstsession = 7;
%         lastsession = 8; 
%     elseif i == 11;  %CellAllignment.C2C3 = session_alligned_cells;
%         cellIndex = CellAllignment.C2C3.CellIDs;
%         firstsession = 8;
%         lastsession = 9;
%     elseif i == 12;  %CellAllignment.C1C3 = session_alligned_cells;
%         cellIndex = CellAllignment.C1C3.CellIDs;
%         firstsession = 7;
%         lastsession = 9;
%     elseif i == 13;  %CellAllignment.A1C3 = session_alligned_cells;
%         cellIndex = CellAllignment.A1C3.CellIDs;
%         firstsession = 1;
%         lastsession = 9;
 else
 end
 
nsessions = length(firstsession:lastsession);
sessionorder = 1;

%This loads all the Zmaps for the alligned cells and puts them in a separate variable 
for SessionLoading = firstsession:lastsession;
    directory_name = directories{SessionLoading};
    cd (directory_name);
    load('extended_calcium_analysis_output.mat') %this loads the file that has our tuning maps 
    for cell_j = 1:length(extended_calcium_analysis_output);       
        TuningMap{cell_j,sessionorder} = conv2(extended_calcium_analysis_output(cell_j).tuning_map,kernel, 'same'); %load and smooth tuning map
        %PDF{cell_j,sessionorder} = extended_calcium_analysis_output(cell_j).PDF; 
        %pN{cell_j,sessionorder} = extended_calcium_analysis_output(cell_j).pN;
    end 
    sessionorder = sessionorder +1; 
    cd (topLevelFolder);
end

%% here we are making a matrix that only has the cell_IDs of all the alligned cells 
cell_ID = cell(1,  nsessions);
tuningmap_stability = zeros(length(cellIndex),length(nsessions));
%PDF_stability = zeros(length(cellIndex),length(nsessions));
%pN_stability = zeros(length(cellIndex),length(nsessions));
TuningMapCorrelationRandom = zeros(length(cellIndex),length(nsessions));
 
for cell_i = 1 : length(cellIndex(:,1));

    % here we load all the cell IDs into the temporary variable cell_ID 
    for session = 1: nsessions;
        cell_ID{1,session} = cellIndex(cell_i,session);
    end 

%% here i combine all the place fields of the sessions in one mat file 

% now correlate the place fields for this cell with the next session
    for session = 1:(nsessions-1);
        cellA = cell2mat(cell_ID(1,session));
        cellB = cell2mat(cell_ID(1,session+1));

        if cellA~= 0 && cellB ~=0;
      
       %%now doing the same for the tuningmaps by extended_place_cell_analysis 
        tuning_map_dayA = TuningMap(cellA,session);
        tuning_map_dayB = TuningMap(cellB,session+1);
        tuning_map_dayB{1,1} = imresize(tuning_map_dayB{1,1},size(tuning_map_dayA{1,1}));
        tuning_map_dayA{1,1}(isnan(tuning_map_dayA{1,1}))=0;
        tuning_map_dayB{1,1}(isnan(tuning_map_dayB{1,1}))=0;
        tuningmap_stability(cell_i, session) = corr2(tuning_map_dayA{1,1},tuning_map_dayB{1,1});     
        
        %% now doing the same for the PDFs
%         PDF_dayA = PDF(cellA,session);
%         PDF_dayB = PDF(cellB,session+1);
%         PDF_dayB{1,1} = imresize(PDF_dayB{1,1},size(PDF_dayA{1,1}));
%         PDF_dayA{1,1}(isnan(PDF_dayA{1,1}))=0;
%         PDF_dayB{1,1}(isnan(PDF_dayB{1,1}))=0;
%         PDF_stability(cell_i, session) = corr2(PDF_dayA{1,1},PDF_dayB{1,1});
        
        %% now doing the same for the pNs
%         pN_dayA = PDF(cellA,session);
%         pN_dayB = PDF(cellB,session+1);
%         pN_dayB{1,1} = imresize(pN_dayB{1,1},size(pN_dayA{1,1}));
%         pN_dayA{1,1}(isnan(pN_dayA{1,1}))=0;
%         pN_dayB{1,1}(isnan(pN_dayB{1,1}))=0;
%         pN_stability(cell_i, session) = corr2(pN_dayA{1,1},pN_dayB{1,1});        
        
        % now I want to compute the tuningmap stability to a random other
        % cell 
        ix  = randi([1,length(extended_calcium_analysis_output)]); % instead of picking a random cell from the alignment index I want a cell from the ms file 
        %cellC = cellIndex(ix,session+1);
        cellC = ix;
        if cellC ~= cellB
        tuning_map_dayC = TuningMap(cellC,session+1);
        tuning_map_dayC{1,1} = imresize(tuning_map_dayC{1,1},size(tuning_map_dayA{1,1}));
        tuning_map_dayC{1,1}(isnan(tuning_map_dayC{1,1}))=0;
        TuningMapCorrelationRandom(cell_i, session) = corr2(tuning_map_dayA{1,1},tuning_map_dayC{1,1});     
        else
         ix  = randi([1,length(cellIndex)]);
         cellC = cellIndex(ix,session+1);
         tuning_map_dayC = TuningMap(cellC,session+1);
        tuning_map_dayC{1,1} = imresize(tuning_map_dayC{1,1},size(tuning_map_dayA{1,1}));
        tuning_map_dayC{1,1}(isnan(tuning_map_dayC{1,1}))=0;
        TuningMapCorrelationRandom(cell_i, session) = corr2(tuning_map_dayA{1,1},tuning_map_dayC{1,1});     
        end
    end 
end

%% save the place field stability in a temporary variable
cellIndex_size = length(cellIndex);

TuningMapCorrelation = tuningmap_stability;
%PDFCorrelation = PDF_stability; 
%pNCorrelation = pN_stability; 
TuningMapCorrelationRandomPair = TuningMapCorrelationRandom;
 
%% Read out first colum and check if it is a place cell, if so keep it in another 
count = 1;
iteration = 1;
round = 1; 

MeanRiseTime = [];
MeanDecayTime = [];
MeanWidth = [];
SplitHalfStability = [];

PlaceCells = [];
InFieldActivityRatio = [];
AlignedPlaceCells = [];
InFieldActivityRatioPlaceCells = [];

AlignedPlaceCellsPrep = [];
SplitHalfStabilityPlaceCells = [];
TuningMapCorrelationPlaceCells = [];

MI = [];
Prob_being_active = [];
Prob_being_activePlaceCells = [];
Prob_active_to_active = [];
Prob_inactive_to_active = []; 

for loadsummary = firstsession:lastsession;
 cd (directories{loadsummary});
 load 'summary.mat'
 load 'extended_calcium_analysis_output.mat'
 load 'firing_probabilities.mat'
      for cell_j = 1:length(summary(:,1));
          
          row = find(cellIndex(:,iteration) == cell_j);
           if ~isnan(row); % if this cell is part of the alligned cells; 
               
          % General firing properties of alligned cells
          MeanRiseTime(row,iteration) = summary(cell_j,9);
          MeanDecayTime(row,iteration) = summary(cell_j,10);
          MeanWidth(row,iteration) = summary(cell_j,12);
          
          MI(row,iteration) = extended_calcium_analysis_output(cell_j).MI;
          Prob_being_active(row,iteration) = extended_calcium_analysis_output(cell_j).prob_being_active;
          Prob_active_to_active(row,iteration) = firing_probabilities.prob_transitioning_from_active_to_active(cell_j);
          Prob_inactive_to_active(row,iteration) = firing_probabilities.prob_transitioning_from_inactive_to_active(cell_j);

          InFieldActivityRatio(row,iteration) = summary(cell_j,4);    
          SplitHalfStability(row,iteration) = extended_calcium_analysis_output(cell_j).place_field_stability;
          
          round = round + 1; 
           else
           end
      end
      iteration = iteration + 1; 
end

iteration = 1; 
PlaceCells_Binary = zeros(size(cellIndex));

for loadsummary = firstsession:lastsession;
 cd (directories{loadsummary});
 load 'extended_calcium_analysis_output.mat'
      for cell_j = 1:length(summary(:,1));
        
                if extended_calcium_analysis_output(cell_j).SignSpatialMI == 1; % if this cell is a place cell % now i kind of want to put a 1 in my 
                    row = find(cellIndex(:,iteration) == cell_j);
                    if ~isnan(row); %if this place cell is part of the alligned cells
                        PlaceCells_Binary(row,iteration) = 1;
                    
                     % Firing properties of alligned place cells 
                    PlaceCells(end+1,:) = cell_j ; %keeps track of the place cells
                    AlignedPlaceCells(end+1,:) = cellIndex(row,:); %keeps track of place cells and aligned partners
                    TuningMapCorrelationPlaceCells(end+1,:) = TuningMapCorrelation(row,:);
                    
                    InFieldActivityRatioPlaceCells(end+1,:) = InFieldActivityRatio(row,:);    
                    SplitHalfStabilityPlaceCells(end+1,:) = SplitHalfStability(row,:);
                    Prob_being_activePlaceCells(end+1,:) = Prob_being_active(row,:);
                    count = count + 1;
                    
                    else
                    end
                end
                cell_j = cell_j + 1;
        end
       iteration = iteration + 1; 
       count = 1;
      % PlaceCells = PlaceCells PlaceCellsSession; 
end

InFieldActivityRatioPlaceCells = unique(InFieldActivityRatioPlaceCells, 'rows', 'stable');
SplitHalfStabilityPlaceCells = unique(SplitHalfStabilityPlaceCells, 'rows', 'stable');
Prob_being_activePlaceCells = unique(Prob_being_activePlaceCells, 'rows', 'stable');
TuningMapCorrelationPlaceCells = unique(TuningMapCorrelationPlaceCells, 'rows', 'stable');

AlignedPlaceCells = unique(AlignedPlaceCells, 'rows','stable'); 

 %ZmapStabilityPlaceCells = ZmapStabilityPlaceCells';
 %PlaceFieldStabilitySummaryPlaceCells = PlaceFieldStabilitySummaryPlaceCells';

        if i == 1   CellAllignmentAnalyzed.A1A2.CellIDs = CellAllignment.A1A2.CellIDs;
                    CellAllignmentAnalyzed.A1A2.TuningMapCorrelation = TuningMapCorrelation;
                    CellAllignmentAnalyzed.A1A2.TuningMapCorrelationRandomPair = TuningMapCorrelationRandomPair;

                   
                   % CellAllignmentAnalyzed.A1A2.PDFCorrelation = PDFCorrelation;
                   % CellAllignmentAnalyzed.A1A2.pNCorrelation = pNCorrelation;
 
                    CellAllignmentAnalyzed.A1A2.MeanRiseTime = MeanRiseTime;
                    CellAllignmentAnalyzed.A1A2.MeanDecayTime = MeanDecayTime;
                    CellAllignmentAnalyzed.A1A2.MeanWidth = MeanWidth;
                    
                    CellAllignmentAnalyzed.A1A2.MI = MI;
                    CellAllignmentAnalyzed.A1A2.Prob_being_active = Prob_being_active;
                    CellAllignmentAnalyzed.A1A2.Prob_active_to_active = Prob_active_to_active;
                    CellAllignmentAnalyzed.A1A2.Prob_inactive_to_active = Prob_inactive_to_active;
                    CellAllignmentAnalyzed.A1A2.InFieldActivityRatio = InFieldActivityRatio;
                    CellAllignmentAnalyzed.A1A2.SplitHalfStability = SplitHalfStability;
                    
                    CellAllignmentAnalyzed.A1A2.PlaceCells = PlaceCells; 
                    CellAllignmentAnalyzed.A1A2.AlignedPlaceCells = AlignedPlaceCells;
                    CellAllignmentAnalyzed.A1A2.PlaceCells_Binary = PlaceCells_Binary;
                    
                    CellAllignmentAnalyzed.A1A2.TuningMapCorrelationPlaceCells = TuningMapCorrelationPlaceCells;
                    CellAllignmentAnalyzed.A1A2.Prob_being_activePlaceCells = Prob_being_activePlaceCells;
                    CellAllignmentAnalyzed.A1A2.InFieldActivityRatioPlaceCells = InFieldActivityRatioPlaceCells;
                    CellAllignmentAnalyzed.A1A2.SplitHalfStabilityPlaceCells = SplitHalfStabilityPlaceCells;       
                 
                   
    elseif i == 2;  CellAllignmentAnalyzed.A2A3.CellIDs = CellAllignment.A2A3.CellIDs;
                    CellAllignmentAnalyzed.A2A3.TuningMapCorrelation = TuningMapCorrelation;
                    CellAllignmentAnalyzed.A2A3.TuningMapCorrelationRandomPair = TuningMapCorrelationRandomPair;

                  %  CellAllignmentAnalyzed.A2A3.PDFCorrelation = PDFCorrelation;
                 %   CellAllignmentAnalyzed.A2A3.pNCorrelation = pNCorrelation;
                    
                    CellAllignmentAnalyzed.A2A3.MeanRiseTime = MeanRiseTime;
                    CellAllignmentAnalyzed.A2A3.MeanDecayTime = MeanDecayTime;
                    CellAllignmentAnalyzed.A2A3.MeanWidth = MeanWidth;
                    
                    CellAllignmentAnalyzed.A2A3.MI = MI;
                    CellAllignmentAnalyzed.A2A3.Prob_being_active = Prob_being_active;
                    CellAllignmentAnalyzed.A2A3.Prob_active_to_active = Prob_active_to_active;
                    CellAllignmentAnalyzed.A2A3.Prob_inactive_to_active = Prob_inactive_to_active;
                     CellAllignmentAnalyzed.A2A3.InFieldActivityRatio = InFieldActivityRatio
                     CellAllignmentAnalyzed.A2A3.SplitHalfStability = SplitHalfStability;
                    
                    CellAllignmentAnalyzed.A2A3.PlaceCells = PlaceCells;
                    CellAllignmentAnalyzed.A2A3.AlignedPlaceCells = AlignedPlaceCells;
                    CellAllignmentAnalyzed.A2A3.PlaceCells_Binary = PlaceCells_Binary;
                    CellAllignmentAnalyzed.A2A3.TuningMapCorrelationPlaceCells = TuningMapCorrelationPlaceCells;
                    CellAllignmentAnalyzed.A2A3.Prob_being_activePlaceCells = Prob_being_activePlaceCells;
                    CellAllignmentAnalyzed.A2A3.InFieldActivityRatioPlaceCells = InFieldActivityRatioPlaceCells;
                    CellAllignmentAnalyzed.A2A3.SplitHalfStabilityPlaceCells = SplitHalfStabilityPlaceCells;       
                    
    elseif i == 3;  CellAllignmentAnalyzed.A1A3.CellIDs = CellAllignment.A1A3.CellIDs;
                    CellAllignmentAnalyzed.A1A3.TuningMapCorrelation = TuningMapCorrelation;
                   CellAllignmentAnalyzed.A1A3.TuningMapCorrelationRandomPair = TuningMapCorrelationRandomPair;

                 %   CellAllignmentAnalyzed.A1A3.PDFCorrelation = PDFCorrelation;
                %    CellAllignmentAnalyzed.A1A3.pNCorrelation = pNCorrelation;
                    
                    CellAllignmentAnalyzed.A1A3.MeanRiseTime = MeanRiseTime;
                    CellAllignmentAnalyzed.A1A3.MeanDecayTime = MeanDecayTime;
                    CellAllignmentAnalyzed.A1A3.MeanWidth = MeanWidth;
                    
                    CellAllignmentAnalyzed.A1A3.MI = MI;
                    CellAllignmentAnalyzed.A1A3.Prob_being_active = Prob_being_active;
                    CellAllignmentAnalyzed.A1A3.Prob_active_to_active = Prob_active_to_active;
                    CellAllignmentAnalyzed.A1A3.Prob_inactive_to_active = Prob_inactive_to_active;
                    CellAllignmentAnalyzed.A1A3.InFieldActivityRatio = InFieldActivityRatio;
                    CellAllignmentAnalyzed.A1A3.SplitHalfStability = SplitHalfStability;
                    
                    CellAllignmentAnalyzed.A1A3.PlaceCells = PlaceCells; 
                    CellAllignmentAnalyzed.A1A3.AlignedPlaceCells = AlignedPlaceCells;
                    CellAllignmentAnalyzed.A1A3.PlaceCells_Binary = PlaceCells_Binary;
                    CellAllignmentAnalyzed.A1A3.TuningMapCorrelationPlaceCells = TuningMapCorrelationPlaceCells;
                    CellAllignmentAnalyzed.A1A3.Prob_being_activePlaceCells = Prob_being_activePlaceCells;
                    CellAllignmentAnalyzed.A1A3.InFieldActivityRatioPlaceCells = InFieldActivityRatioPlaceCells;
                    CellAllignmentAnalyzed.A1A3.SplitHalfStabilityPlaceCells = SplitHalfStabilityPlaceCells;       

                    
%     elseif i == 4;   CellAllignment.A3B1.ZmapCorr = ZmapStabilitySummary;
%                     CellAllignment.A3B1.PlaceFieldStabilitySummary = PlaceFieldStabilitySummary;
%                     CellAllignment.A3B1.TuningMapStabilitySummary = tuningmapStabilitySummary;
%                     CellAllignment.A3B1.PDFStabilitySummary = PDFStabilitySummary;
%                      CellAllignment.A3B1.pNStabilitySummary = pNStabilitySummary;
%                     
%                     
%                     CellAllignment.A3B1.FiringProbability = FiringProbability;
%                     CellAllignment.A3B1.MeanRiseTime = MeanRiseTime;
%                     CellAllignment.A3B1.MeanDecayTime = MeanDecayTime;
%                     CellAllignment.A3B1.MeanWidth = MeanWidth;
%                     
%                     CellAllignment.A3B1.KLD = KLD;
%                     CellAllignment.A3B1.MI = MI;
%                     CellAllignment.A3B1.Prob_being_active = Prob_being_active;
%             
%                     
%                     CellAllignment.A3B1.PlaceCells = PlaceCells; 
%                     CellAllignment.A3B1.AlignedPlaceCells = AlignedPlaceCells;
%                     CellAllignment.A3B1.InFieldActivityRatio = InFieldActivityRatio;
%                     CellAllignment.A3B1.FiringProbabilityPlaceCell = FiringProbabilityPlaceCell; 
%                     CellAllignment.A3B1.PlaceFieldArea = PlaceFieldArea; 
%                     CellAllignment.A3B1.SplitHalfStability = SplitHalfStability;
%                     
%                     CellAllignment.A3B1.ZmapCorrPlaceCells = ZmapStabilityPlaceCells;
%                     CellAllignment.A3B1.PlaceFieldStabilitySummaryPlaceCells = PlaceFieldStabilitySummaryPlaceCells;
%                     
%     
%          elseif i == 5; 
%                      CellAllignment.A1B1.ZmapCorr = ZmapStabilitySummary;
%                      CellAllignment.A1B1.PlaceFieldStabilitySummary = PlaceFieldStabilitySummary;
%                      CellAllignment.A1B1.TuningMapStabilitySummary = tuningmapStabilitySummary;
%                      CellAllignment.A1B1.PDFStabilitySummary = PDFStabilitySummary;
%                       CellAllignment.A1B1.pNStabilitySummary = pNStabilitySummary;
%                      
%                      
%                      CellAllignment.A1B1.FiringProbability = FiringProbability;
%                      CellAllignment.A1B1.MeanRiseTime = MeanRiseTime;
%                      CellAllignment.A1B1.MeanDecayTime = MeanDecayTime;
%                      CellAllignment.A1B1.MeanWidth = MeanWidth;
%                      
%                      CellAllignment.A1B1.KLD = KLD;
%                      CellAllignment.A1B1.MI = MI;
%                      CellAllignment.A1B1.Prob_being_active = Prob_being_active;
%                      
%                      CellAllignment.A1B1.PlaceCells = PlaceCells;
%                      CellAllignment.A1B1.AlignedPlaceCells = AlignedPlaceCells;
%                      CellAllignment.A1B1.InFieldActivityRatio = InFieldActivityRatio;
%                      CellAllignment.A1B1.FiringProbabilityPlaceCell = FiringProbabilityPlaceCell; 
%                      CellAllignment.A1B1.PlaceFieldArea = PlaceFieldArea; 
%                      CellAllignment.A1B1.SplitHalfStability = SplitHalfStability;
%                      
%                      CellAllignment.A1B1.ZmapCorrPlaceCells = ZmapStabilityPlaceCells;
%                      CellAllignment.A1B1.PlaceFieldStabilitySummaryPlaceCells = PlaceFieldStabilitySummaryPlaceCells;
%         
%     elseif i == 6; CellAllignment.B2B3.ZmapCorr = ZmapStabilitySummary;
%                     CellAllignment.B2B3.PlaceFieldStabilitySummary = PlaceFieldStabilitySummary;
%                     CellAllignment.B2B3.TuningMapStabilitySummary = tuningmapStabilitySummary;
%                     CellAllignment.B2B3.PDFStabilitySummary = PDFStabilitySummary;
%                     CellAllignment.B2B3.pNStabilitySummary = pNStabilitySummary;
%                     
%                     
%                     
%                     CellAllignment.B2B3.FiringProbability = FiringProbability;
%                     CellAllignment.B2B3.MeanRiseTime = MeanRiseTime;
%                     CellAllignment.B2B3.MeanDecayTime = MeanDecayTime;
%                     CellAllignment.B2B3.MeanWidth = MeanWidth;
%                     
%                      CellAllignment.B2B3.KLD = KLD;
%                      CellAllignment.B2B3.MI = MI;
%                      CellAllignment.B2B3.Prob_being_active = Prob_being_active;
%                     
%                     CellAllignment.B2B3.PlaceCells = PlaceCells; 
%                     CellAllignment.B2B3.AlignedPlaceCells = AlignedPlaceCells;
%                     CellAllignment.B2B3.InFieldActivityRatio = InFieldActivityRatio;
%                     CellAllignment.B2B3.FiringProbabilityPlaceCell = FiringProbabilityPlaceCell; 
%                     CellAllignment.B2B3.PlaceFieldArea = PlaceFieldArea; 
%                     CellAllignment.B2B3.SplitHalfStability = SplitHalfStability;
%                     
%                     CellAllignment.B2B3.ZmapCorrPlaceCells = ZmapStabilityPlaceCells;
%                     CellAllignment.B2B3.PlaceFieldStabilitySummaryPlaceCells = PlaceFieldStabilitySummaryPlaceCells;
%                     
%     elseif i == 7; CellAllignment.B1B3.ZmapCorr = ZmapStabilitySummary;
%                     CellAllignment.B1B3.PlaceFieldStabilitySummary = PlaceFieldStabilitySummary;
%                     CellAllignment.B1B3.TuningMapStabilitySummary = tuningmapStabilitySummary;
%                     CellAllignment.B1B3.PDFStabilitySummary = PDFStabilitySummary;
%                      CellAllignment.B1B3.pNStabilitySummary = pNStabilitySummary;
%                     
%                     CellAllignment.B1B3.FiringProbability = FiringProbability;
%                     CellAllignment.B1B3.MeanRiseTime = MeanRiseTime;
%                     CellAllignment.B1B3.MeanDecayTime = MeanDecayTime;
%                     CellAllignment.B1B3.MeanWidth = MeanWidth;
%                     
%                     CellAllignment.B1B3.KLD = KLD;
%                     CellAllignment.B1B3.MI = MI;
%                     CellAllignment.B1B3.Prob_being_active = Prob_being_active;
%                     
%                     CellAllignment.B1B3.PlaceCells = PlaceCells; 
%                     CellAllignment.B1B3.AlignedPlaceCells = AlignedPlaceCells;
%                     CellAllignment.B1B3.InFieldActivityRatio = InFieldActivityRatio;
%                     CellAllignment.B1B3.FiringProbabilityPlaceCell = FiringProbabilityPlaceCell; 
%                     CellAllignment.B1B3.PlaceFieldArea = PlaceFieldArea; 
%                     CellAllignment.B1B3.SplitHalfStability = SplitHalfStability;
%                     
%                     CellAllignment.B1B3.ZmapCorrPlaceCells = ZmapStabilityPlaceCells;
%                     CellAllignment.B1B3.PlaceFieldStabilitySummaryPlaceCells = PlaceFieldStabilitySummaryPlaceCells;
%                     
%                     
%     elseif i == 8;  
%                     CellAllignment.A1B3.ZmapCorr = ZmapStabilitySummary;
%                     CellAllignment.A1B3.PlaceFieldStabilitySummary = PlaceFieldStabilitySummary;
%                     CellAllignment.A1B3.TuningMapStabilitySummary = tuningmapStabilitySummary;
%                     CellAllignment.A1B3.PDFStabilitySummary = PDFStabilitySummary;
%  CellAllignment.A1B3.pNStabilitySummary = pNStabilitySummary;
%                     
%                     CellAllignment.A1B3.FiringProbability = FiringProbability;
%                     CellAllignment.A1B3.MeanRiseTime = MeanRiseTime;
%                     CellAllignment.A1B3.MeanDecayTime = MeanDecayTime;
%                     CellAllignment.A1B3.MeanWidth = MeanWidth;
%                     
%                     CellAllignment.A1B3.KLD = KLD;
%                     CellAllignment.A1B3.MI = MI;
%                     CellAllignment.A1B3.Prob_being_active = Prob_being_active;
%                     
%                     CellAllignment.A1B3.PlaceCells = PlaceCells; 
%                     CellAllignment.A1B3.AlignedPlaceCells = AlignedPlaceCells;
%                     CellAllignment.A1B3.InFieldActivityRatio = InFieldActivityRatio;
%                     CellAllignment.A1B3.FiringProbabilityPlaceCell = FiringProbabilityPlaceCell; 
%                     CellAllignment.A1B3.PlaceFieldArea = PlaceFieldArea; 
%                     CellAllignment.A1B3.SplitHalfStability = SplitHalfStability;
%                     
%                     CellAllignment.A1B3.ZmapCorrPlaceCells = ZmapStabilityPlaceCells;
%                     CellAllignment.A1B3.PlaceFieldStabilitySummaryPlaceCells = PlaceFieldStabilitySummaryPlaceCells;
%                     
%      
%         elseif i == 9;  
%                     CellAllignment.B1B2.ZmapCorr = ZmapStabilitySummary;
%                     CellAllignment.B1B2.PlaceFieldStabilitySummary = PlaceFieldStabilitySummary;
%                     CellAllignment.B1B2.TuningMapStabilitySummary = tuningmapStabilitySummary;
%                     CellAllignment.B1B2.PDFStabilitySummary = PDFStabilitySummary;
%  CellAllignment.B1B2.pNStabilitySummary = pNStabilitySummary;
%                     
%                     CellAllignment.B1B2.FiringProbability = FiringProbability;
%                     CellAllignment.B1B2.MeanRiseTime = MeanRiseTime;
%                     CellAllignment.B1B2.MeanDecayTime = MeanDecayTime;
%                     CellAllignment.B1B2.MeanWidth = MeanWidth;
%                     
%                     CellAllignment.B1B2.KLD = KLD;
%                     CellAllignment.B1B2.MI = MI;
%                     CellAllignment.B1B2.Prob_being_active = Prob_being_active;
%                     
%                     CellAllignment.B1B2.PlaceCells = PlaceCells;
%                     CellAllignment.B1B2.AlignedPlaceCells = AlignedPlaceCells;
%                     CellAllignment.B1B2.InFieldActivityRatio = InFieldActivityRatio;
%                     CellAllignment.B1B2.FiringProbabilityPlaceCell = FiringProbabilityPlaceCell; 
%                     CellAllignment.B1B2.PlaceFieldArea = PlaceFieldArea; 
%                     CellAllignment.B1B2.SplitHalfStability = SplitHalfStability;
%                     
%                     CellAllignment.B1B2.ZmapCorrPlaceCells = ZmapStabilityPlaceCells;
%                     CellAllignment.B1B2.PlaceFieldStabilitySummaryPlaceCells = PlaceFieldStabilitySummaryPlaceCells;
%                      
%                     
%     elseif i == 9;  CellAllignment.B3C1.ZmapCorr = ZmapStabilitySummary;
%                     CellAllignment.B3C1.PlaceFieldStabilitySummary = PlaceFieldStabilitySummary;
%                     CellAllignment.B3C1.TuningMapStabilitySummary = tuningmapStabilitySummary;
%                     CellAllignment.B3C1.PDFStabilitySummary = PDFStabilitySummary;
%  CellAllignment.B3C1.pNStabilitySummary = pNStabilitySummary;
%                     
%                     CellAllignment.B3C1.FiringProbability = FiringProbability;
%                     CellAllignment.B3C1.MeanRiseTime = MeanRiseTime;
%                     CellAllignment.B3C1.MeanDecayTime = MeanDecayTime;
%                     CellAllignment.B3C1.MeanWidth = MeanWidth;
%                     
%                     CellAllignment.B3C1.KLD = KLD;
%                     CellAllignment.B3C1.MI = MI;
%                     CellAllignment.B3C1.Prob_being_active = Prob_being_active;
%                     
%                     
%                     CellAllignment.B3C1.PlaceCells = PlaceCells; 
%                     CellAllignment.B3C1.AlignedPlaceCells = AlignedPlaceCells;
%                     CellAllignment.B3C1.InFieldActivityRatio = InFieldActivityRatio;
%                     CellAllignment.B3C1.FiringProbabilityPlaceCell = FiringProbabilityPlaceCell; 
%                     CellAllignment.B3C1.PlaceFieldArea = PlaceFieldArea; 
%                     CellAllignment.B3C1.SplitHalfStability = SplitHalfStability;
%                     
%                     CellAllignment.B3C1.ZmapCorrPlaceCells = ZmapStabilityPlaceCells;
%                     CellAllignment.B3C1.PlaceFieldStabilitySummaryPlaceCells = PlaceFieldStabilitySummaryPlaceCells;
% 
%                     
%     elseif i == 10;
%         
%                     CellAllignment.C1C2.ZmapCorr = ZmapStabilitySummary;
%                     CellAllignment.C1C2.PlaceFieldStabilitySummary = PlaceFieldStabilitySummary;
%                     CellAllignment.C1C2.TuningMapStabilitySummary = tuningmapStabilitySummary;
%                     CellAllignment.C1C2.PDFStabilitySummary = PDFStabilitySummary;
%                     CellAllignment.C1C2.pNStabilitySummary = pNStabilitySummary;
%                     
%                     CellAllignment.C1C2.FiringProbability = FiringProbability;
%                     CellAllignment.C1C2.MeanRiseTime = MeanRiseTime;
%                     CellAllignment.C1C2.MeanDecayTime = MeanDecayTime;
%                     CellAllignment.C1C2.MeanWidth = MeanWidth;
%                     
%                     CellAllignment.C1C2.KLD = KLD;
%                     CellAllignment.C1C2.MI = MI;
%                     CellAllignment.C1C2.Prob_being_active = Prob_being_active;
%                     
%                     
%                     CellAllignment.C1C2.PlaceCells = PlaceCells;
%                     CellAllignment.C1C2.AlignedPlaceCells = AlignedPlaceCells;
%                     CellAllignment.C1C2.InFieldActivityRatio = InFieldActivityRatio;
%                     CellAllignment.C1C2.FiringProbabilityPlaceCell = FiringProbabilityPlaceCell; 
%                     CellAllignment.C1C2.PlaceFieldArea = PlaceFieldArea; 
%                     CellAllignment.C1C2.SplitHalfStability = SplitHalfStability;
%                     
%                     CellAllignment.C1C2.ZmapCorrPlaceCells = ZmapStabilityPlaceCells;
%                     CellAllignment.C1C2.PlaceFieldStabilitySummaryPlaceCells = PlaceFieldStabilitySummaryPlaceCells;
%                     
%     elseif i == 11; 
%         
%                     CellAllignment.C2C3.ZmapCorr = ZmapStabilitySummary;
%                     CellAllignment.C2C3.PlaceFieldStabilitySummary = PlaceFieldStabilitySummary;
%                     CellAllignment.C2C3.TuningMapStabilitySummary = tuningmapStabilitySummary;
%                     CellAllignment.C2C3.PDFStabilitySummary = PDFStabilitySummary;
%  CellAllignment.C2C3.pNStabilitySummary = pNStabilitySummary;
%                     
%                     CellAllignment.C2C3.FiringProbability = FiringProbability;
%                     CellAllignment.C2C3.MeanRiseTime = MeanRiseTime;
%                     CellAllignment.C2C3.MeanDecayTime = MeanDecayTime;
%                     CellAllignment.C2C3.MeanWidth = MeanWidth;
%                     
%                     CellAllignment.C2C3.KLD = KLD;
%                     CellAllignment.C2C3.MI = MI;
%                     CellAllignment.C2C3.Prob_being_active = Prob_being_active;
%                     
%                     
%                     CellAllignment.C2C3.PlaceCells = PlaceCells; 
%                     CellAllignment.C2C3.AlignedPlaceCells = AlignedPlaceCells;
%                     CellAllignment.C2C3.InFieldActivityRatio = InFieldActivityRatio;
%                     CellAllignment.C2C3.FiringProbabilityPlaceCell = FiringProbabilityPlaceCell; 
%                     CellAllignment.C2C3.PlaceFieldArea = PlaceFieldArea; 
%                     CellAllignment.C2C3.SplitHalfStability = SplitHalfStability;
%                     
%                     CellAllignment.C2C3.ZmapCorrPlaceCells = ZmapStabilityPlaceCells;
%                     CellAllignment.C2C3.PlaceFieldStabilitySummaryPlaceCells = PlaceFieldStabilitySummaryPlaceCells;
%                    
%     elseif i == 12; 
%                     CellAllignment.C1C3.ZmapCorr = ZmapStabilitySummary;
%                     CellAllignment.C1C3.PlaceFieldStabilitySummary = PlaceFieldStabilitySummary;
%                     CellAllignment.C1C3.TuningMapStabilitySummary = tuningmapStabilitySummary;
%                     CellAllignment.C1C3.PDFStabilitySummary = PDFStabilitySummary;
%                     CellAllignment.C1C3.pNStabilitySummary = pNStabilitySummary;
%                     
%                     CellAllignment.C1C3.FiringProbability = FiringProbability;
%                     CellAllignment.C1C3.MeanRiseTime = MeanRiseTime;
%                     CellAllignment.C1C3.MeanDecayTime = MeanDecayTime;
%                     CellAllignment.C1C3.MeanWidth = MeanWidth;
%                     
%                     CellAllignment.C1C3.KLD = KLD;
%                     CellAllignment.C1C3.MI = MI;
%                     CellAllignment.C1C3.Prob_being_active = Prob_being_active;
%                     
%                     CellAllignment.C1C3.PlaceCells = PlaceCells; 
%                     CellAllignment.C1C3.AlignedPlaceCells = AlignedPlaceCells;
%                     CellAllignment.C1C3.InFieldActivityRatio = InFieldActivityRatio;
%                     CellAllignment.C1C3.FiringProbabilityPlaceCell = FiringProbabilityPlaceCell; 
%                     CellAllignment.C1C3.PlaceFieldArea = PlaceFieldArea; 
%                     CellAllignment.C1C3.SplitHalfStability = SplitHalfStability;
%                     
%                     CellAllignment.C1C3.ZmapCorrPlaceCells = ZmapStabilityPlaceCells;
%                     CellAllignment.C1C3.PlaceFieldStabilitySummaryPlaceCells = PlaceFieldStabilitySummaryPlaceCells;
%                                        
%     elseif i == 13; 
%                     CellAllignment.A1C3.ZmapCorr = ZmapStabilitySummary;
%                     CellAllignment.A1C3.PlaceFieldStabilitySummary = PlaceFieldStabilitySummary;
%                     CellAllignment.A1C3.TuningMapStabilitySummary = tuningmapStabilitySummary;
%                     CellAllignment.A1C3.PDFStabilitySummary = PDFStabilitySummary;
%                     CellAllignment.A1C3.pNStabilitySummary = pNStabilitySummary;
%                     
%                     CellAllignment.A1C3.FiringProbability = FiringProbability;
%                     CellAllignment.A1C3.MeanRiseTime = MeanRiseTime;
%                     CellAllignment.A1C3.MeanDecayTime = MeanDecayTime;
%                     CellAllignment.A1C3.MeanWidth = MeanWidth;
%                     
%                     CellAllignment.A1C3.KLD = KLD;
%                     CellAllignment.A1C3.MI = MI;
%                     CellAllignment.A1C3.Prob_being_active = Prob_being_active;
%                     
%                     
%                     CellAllignment.A1C3.PlaceCells = PlaceCells; 
%                     CellAllignment.A1C3.AlignedPlaceCells = AlignedPlaceCells;
%                     CellAllignment.A1C3.InFieldActivityRatio = InFieldActivityRatio;
%                     CellAllignment.A1C3.FiringProbabilityPlaceCell = FiringProbabilityPlaceCell; 
%                     CellAllignment.A1C3.PlaceFieldArea = PlaceFieldArea; 
%                     CellAllignment.A1C3.SplitHalfStability = SplitHalfStability;
%                     
%                     CellAllignment.A1C3.ZmapCorrPlaceCells = ZmapStabilityPlaceCells;
%                     CellAllignment.A1C3.PlaceFieldStabilitySummaryPlaceCells = PlaceFieldStabilitySummaryPlaceCells;

    else
        end
    
end
cd (topLevelFolder)
save('CellAllignmentAnalyzed.mat', 'CellAllignmentAnalyzed');
end 


