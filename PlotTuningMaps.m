figure;
for i = 1:30;
    subplot(5,6,i);
    smoothed_map = smoothPlaceField(extended_calcium_analysis_output(i).tuning_map);
    imagesc(smoothed_map);
   % imagesc(extended_calcium_analysis_output(i).tuning_map);
    
    
end

smoothed_map = smoothPlaceField(tuning_map(18).likelihood); 
imagesc(smoothed_map); 