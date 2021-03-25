function [interp_HD_vec] = interpolate_HD(HD_vec, behav_time, ca_time)
%INTERP_BEHAV Interpolates head direction angles
% This function uses angles containing missing values as inputs, and
% outputs interpolated angles using a windowing approach and unwrapping
% angles

% Copyright (C) 2020 by Guillaume Etter
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or any
% later version.
% Contact: etterguillaume@gmail.com

% behav_vec: m x 1  vector or a m x 2 matrix representing head direction signals m is corresponds to the total number of frames and has to be the
% same size as calcium_time
%
% behav_time: vector containing the timestamps (is s) corresponding to each
% behav frame.

% ca_time: vector containing the timestamps (is s) corresponding to each
% calcium imaging frame. ca_time should be longer or the same length as behav_time

   %% Fill missing values
   HD_vec = unwrap(deg2rad(HD_vec));
   HD_vec = fillmissing(HD_vec,'pchip');
   HD_vec = rad2deg(wrapTo2Pi(HD_vec));
   
   %% Perform linear interpolation
   interp_HD_vec = interp1(behav_time,HD_vec,ca_time,'nearest');

end


