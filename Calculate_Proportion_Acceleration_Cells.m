%% This function calculates the overall proportion of place cells, as well as the proportion of place cells per bin

function [Proportion] = calculate_proportion_acceleration_cells;

disp('Select toplevel folder');
topLevelFolder = uigetdir;
cd (topLevelFolder);
matrix = zeros(15);

%% input number of sessions to look at 
prompt = 'Input number of recordings to combine:   ';
nsessions = input(prompt);
%% this loads all the directories of the sessions to look at

for iteration = 1:nsessions;
    disp('Select directory of session to be loaded');
    directory_name = uigetdir;
    directories{iteration} = directory_name;
    directory_name;
end

%% now we go through all directories to extract place cell information 

session = 1;
for session = 1:nsessions;
    directory_name = directories{session};
    cd (directory_name);
    
    load('extended_calcium_analysis_output_2021.mat')

    summary = [];
    
  %% calculate the overall proportion of place cells 
  NCells = size(extended_calcium_analysis_output,2);

  NPlaceCells = 0;
    
 for cell_i = 1:NCells;
 if extended_calcium_analysis_output(cell_i).SignMI_acceleration == 1
 NPlaceCells = NPlaceCells + 1; 
 end
 end

  Proportion.PercentagePlaceCells = (NPlaceCells/NCells) * 100;
  ProportionPlaceCells(session,:) = Proportion.PercentagePlaceCells;
        
   save('ProportionAcceleration.mat', 'Proportion');
      
    session
end 

  cd (topLevelFolder);
 
  ProportionsSummary.Session = ProportionPlaceCells;
  save('ProportionsAcceleration.mat', 'ProportionsSummary');
  
    
end

    



