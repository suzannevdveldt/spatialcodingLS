function msConjunctiveCoding;

%clear all;

disp('Select directory of open field session');
topLevelFolder = uigetdir;
cd(topLevelFolder);

ExtendedSummary = 0;

prompt = 'Input number of recordings to combine:   ';
nsessions = input(prompt);

%this will ask you to open the folder 
for iteration = 1:nsessions;
    cd(topLevelFolder);
    disp('Select directory of session to be loaded');
    directory_name = uigetdir;
    cd(directory_name);
    %fprintf(directory_name); % so you can keep track of last one opened 
    
    %% here load the variables that I want to combine 
   % load('summary.mat');
    load('ms.mat');
    load('extended_calcium_analysis_output.mat');
 %   load('firing_probabilities.mat');
    
prompt = 'Input animal number:   ';
animal_ID = input(prompt);

%animal_ID = ID_list(iteration,1);

extended_summary = [];

    for cell_i = 1:ms.numNeurons; 
         extended_summary(cell_i,1) = string(animal_ID); %animal id 
        extended_summary(cell_i,2) = cell_i;  %cell_id
        
if extended_calcium_analysis_output(cell_i).spatial_MI_pvalue <= .01
    extended_summary(cell_i,3) = 1; 
else 
 extended_summary(cell_i,3) = 0;
end 

if extended_calcium_analysis_output(cell_i).velocity_MI_pvalue <= .01
    extended_summary(cell_i,4) = 1; 
else 
   extended_summary(cell_i,4) = 0;
end 

if extended_calcium_analysis_ou5tput(cell_i).acceleration_MI_pvalue <= .01
    extended_summary(cell_i,5) = 1; 
else 
 extended_summary(cell_i,5) = 0;
end 

if extended_calcium_analysis_output(cell_i).headdirection_MI_pvalue <= .01
    extended_summary(cell_i,6) = 1; 
else 
 extended_summary(cell_i,6) = 0;
end 
    end
    
    if iteration == 1;
        ExtendedSummary = extended_summary;
    else
        ExtendedSummary = [ExtendedSummary; extended_summary];
    end
  
end

    %% now I want to combine all the summaries together;
cd(topLevelFolder);

ConjunctiveCoding.Coding = ExtendedSummary;
ConjunctiveCoding.Info = ["spatial"; "velocity"; "acceleration"; "headdirection"]

save('ConjunctiveCoding.mat', 'ConjunctiveCoding');

end

