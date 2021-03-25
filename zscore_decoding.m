clear all;

load('decoding_agreement.mat');
load('decoding_error.mat');
load('shuffled_decoding_agreement.mat');
load('shuffled_decoding_error.mat');

nshuffles = 30;

mean_shuffled_error = mean(shuffled_mean_decoding_error);
SD_shuffled_error = std(shuffled_mean_decoding_error);

mean_shuffled_agreement = mean(shuffled_decoding_agreement);
SD_shuffled_agreement = std(shuffled_decoding_agreement);

for n = 1:nshuffles;

Z_ActualvsShuffleDecodingError(n,1) = (mean_decoding_error(1,n) - mean_shuffled_error)./ SD_shuffled_error;
Z_ActualvsShuffleDecodingAgreement(n,1) = (decoding_agreement(1,n) - mean_shuffled_agreement)./ SD_shuffled_agreement;

end 

MeanZ_ActualvsShuffleDecodingError = mean(Z_ActualvsShuffleDecodingError);
MeanZ_ActualvsShuffleDecodingAgreement = mean(Z_ActualvsShuffleDecodingAgreement);