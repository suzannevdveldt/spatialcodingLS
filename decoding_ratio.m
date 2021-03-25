function decoding_ratio

 %% here load the variables that I want to combine 
    load('decoding.mat');
 
decoding_error = decoding.error;
shuffled_decoding_error = decoding.shuffled_error;

decoding_agreement = decoding.agreement;
shuffled_decoding_agreement = decoding.shuffled_agreement;

% compute means 

mean_decoding_error = mean(decoding_error);
mean_shuffled_decoding_error= mean(shuffled_decoding_error);

mean_decoding_agreement = mean(decoding_agreement); 
mean_shuffled_decoding_agreement = mean(shuffled_decoding_agreement);

% compute scores 
decoding_error_score = mean_decoding_error/mean_shuffled_decoding_error;
decoding_agreement_score = mean_decoding_agreement/mean_shuffled_decoding_agreement;

decoding.mean_decoding_error = mean_decoding_error;
decoding.mean_decoding_agreement = mean_decoding_agreement; 
decoding.error_score = decoding_error_score;
decoding.agreement_score = decoding_agreement_score;

save('decoding.mat', 'decoding')
end

