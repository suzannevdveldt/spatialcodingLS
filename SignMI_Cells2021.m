% Check significance MI value

function [] = SignMI_Cells(SignMI_LS);

disp('Select directory of open field session');
topLevelFolder = uigetdir;
cd(topLevelFolder);

prompt = 'Input number of recordings to combine:   ';
nsessions = input(prompt);

for iteration = 1:nsessions;
    cd(topLevelFolder);
    disp('Select directory of session to be loaded');
    directory_name = uigetdir;
    cd(directory_name);
        
prompt = 'Input animal ID:   '
ID = input(prompt);

SignCells = []; 
Cells = [];

load('extended_calcium_analysis_output_2021.mat');

 TF = isfield(extended_calcium_analysis_output, 'SignBootstrappedMI');
 if TF == 1; 
 extended_calcium_analysis_output = rmfield(extended_calcium_analysis_output, 'SignBootstrappedMI');
 end

SignMI_LS(:,1) == ID; %%Animal ID 
Cells = ans;
for row = 1:length(Cells);
    if Cells(row) == 1, %if this neuron is recorded in that particular animal
        cell_id = SignMI_LS(row,2);
        SignCells(cell_id,1) = 1;
        extended_calcium_analysis_output(cell_id).SignMI_HD = 1; 

    end
end

for row = 1:length(extended_calcium_analysis_output);
    if isempty(extended_calcium_analysis_output(row).SignMI_HD)
        extended_calcium_analysis_output(row).SignMI_HD  = 0;
       
    end
end 
save('extended_calcium_analysis_output_2021.mat', 'extended_calcium_analysis_output'); 

end
end

 