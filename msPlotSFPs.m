function [] = msPlotSFPs(SFP)
load('SFP.mat');
figure;imagesc(max(permute(SFP,[2 3 1]),[],3));
end
