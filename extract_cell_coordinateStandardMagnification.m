function [cell_coordinate] = extract_cell_coordinateStandardMagnification;
clear all

diameterGRIN = .5; 

load('ms.mat')
load('SFP.mat')
%load('summary.mat')

NCells = ms.numNeurons; 

prompt = 'Input GRIN lens AP coordinate:   ';
GRINcoordinate.AP = input(prompt);

prompt = 'Input GRIN lens ML coordinate:   ';
GRINcoordinate.ML = input(prompt);

prompt = 'Input GRIN lens DV coordinate:   ';
GRINcoordinate.DV = input(prompt);

Anatomy.GRIN = GRINcoordinate;

%% Plot the SFPs to allow for manual identification of zoomfactor
figure;imagesc(max(permute(SFP,[2 3 1]),[],3));
% 
 fprintf('Select the upper boundary of the grin     ');
 [XTopGrin,YTopGrin] = ginput(1); 
 fprintf('Select the lower boundary of the grin      ');
 [XBottomGrin,YBottomGrin] = ginput(1);

%the difference between these numbers should be equal to .5 mm 

GrinPlotSize = 160; %YBottomGrin - YTopGrin; %% This value will change depending on what downsampling you use 
PixelsperMM = GrinPlotSize / diameterGRIN;

% manually select where the middle of the GRIN would be 
fprintf('Select the approximate middle of the GRIN lens     ');
[XMiddleGrin,YMiddleGrin] = ginput(1);

for cell_id = 1:NCells;
%imagesc(permute(cell_registered_struct.spatial_footprints_corrected{session,1}(cell_id,:,:),[2 3 1]))
figure;imagesc(permute(SFP(cell_id,:,:),[2 3 1]))
[Xcell,Ycell] = ginput(1); 

%% Calculate how far the cell is from the middle of the GRIn 

XDistance = XMiddleGrin - Xcell;
YDistance = YMiddleGrin - Ycell; 

XDistanceMM = XDistance / PixelsperMM; 
YDistanceMM = YDistance / PixelsperMM; 

Anatomy(cell_id).APcoordinate = GRINcoordinate.AP + XDistanceMM;
Anatomy(cell_id).MLcoordinate = GRINcoordinate.ML + YDistanceMM; 
Anatomy(cell_id).DVcoordinate = GRINcoordinate.DV;

%% Save results in a file that can be easily accessed 
%summary(cell_id,17) = Anatomy(cell_id).APcoordinate;
%summary(cell_id,18) = Anatomy(cell_id).MLcoordinate;
%summary(cell_id,19) = Anatomy(cell_id).DVcoordinate; 

close all;
end 

save('Anatomy.mat', 'Anatomy');
%save('summary.mat', 'summary');
end
 