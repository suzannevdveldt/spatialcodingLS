function [Proportion] = BinCellCoordinates;
%% this function calculates the number of cells per anatomical bin, the number of place cells and the proportion of place cells 

APbinsize = .2;
MLbinsize = .2;
DVbinsize = .3;

load('Anatomy.mat');
load('summary.mat');
load('extended_calcium_analysis_output.mat');

%Creating bins
APbins = -0.2:APbinsize:1.3;
MLbins = -0.2:MLbinsize:.8;
DVbins = 2.4:DVbinsize:3.6;

% Binning the AP,ML and DV coordinates 
for i = 1:size(Anatomy,2);
Anatomy(i).APBin = discretize(Anatomy(i).APcoordinate, APbins);
tempAPbin(i) = Anatomy(i).APBin;
Anatomy(i).MLBin = discretize(Anatomy(i).MLcoordinate, MLbins);
tempMLbin(i) = Anatomy(i).MLBin;
Anatomy(i).DVBin = discretize(Anatomy(i).DVcoordinate, DVbins);
tempDVbin(i) = Anatomy(i).DVBin;
SignSpatialMI(i) =  extended_calcium_analysis_output(i).SignSpatialMI; 
end

% % Now for each bin, calculate the total number of cells in that bin and
% % calculate the proportion of place cells 
    %go through AP
for nbin =  1:(length(APbins)-1); 
    nplacecellsAP(1,nbin) = 0;
     ncellsAP(1,nbin) = sum(tempAPbin==nbin);
        for cell_i = 1:size(tempAPbin,2); 
            if tempAPbin(cell_i) == nbin && SignSpatialMI(cell_i) == 1; 
                nplacecellsAP(1,nbin) =  nplacecellsAP(1,nbin) + 1; 
            else
            end
        end
        proportionplacecellsAP(1,nbin) = ((nplacecellsAP(1,nbin)/ncellsAP(1,nbin)) .* 100);
end

% got through ML 
for nbin =  1:(length(MLbins)-1); 
    nplacecellsML(1,nbin) = 0;
     ncellsML(1,nbin) = sum(tempMLbin==nbin);
        for cell_i = 1:size(tempMLbin,2); 
            if tempMLbin(cell_i) == nbin && SignSpatialMI(cell_i) == 1; 
                nplacecellsML(1,nbin) =  nplacecellsML(1,nbin) + 1; 
            else
            end
        end
        proportionplacecellsML(1,nbin) = ((nplacecellsML(1,nbin)/ncellsML(1,nbin)) .* 100);
end
    
    %go through DV
for nbin =  1:(length(DVbins)-1); 
    nplacecellsDV(1,nbin) = 0;
     ncellsDV(1,nbin) = sum(tempDVbin==nbin);
        for cell_i = 1:size(tempDVbin,2); 
            if tempDVbin(cell_i) == nbin && SignSpatialMI(cell_i) == 1; 
                nplacecellsDV(1,nbin) =  nplacecellsDV(1,nbin) + 1; 
            else
            end
        end
        proportionplacecellsDV(1,nbin) = ((nplacecellsDV(1,nbin)/ncellsDV(1,nbin)) .* 100);
end
    
  Proportion.NCells.AP = ncellsAP;
  Proportion.NPlaceCells.AP = nplacecellsAP;
  Proportion.NCells.ML = ncellsML;
  Proportion.NPlaceCells.ML = nplacecellsML;
  Proportion.NCells.DV = ncellsDV;
  Proportion.NPlaceCells.DV = nplacecellsDV;
  
  Proportion.PlaceCells.AP = proportionplacecellsAP';
  Proportion.PlaceCells.ML = proportionplacecellsML';
  Proportion.PlaceCells.DV = proportionplacecellsDV';
  
 save('Proportion.mat', 'Proportion');
 save('Anatomy.mat', 'Anatomy');
 save('summary.mat', 'summary');
end
