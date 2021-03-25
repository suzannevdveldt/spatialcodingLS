function smoothed_map = smooth_map(input_map)
    % Used to smooth place fields
   kernel_size = [3 3];
   occupancy_std = 2;

[Xgrid,Ygrid]=meshgrid(-kernel_size(1)/2: kernel_size(1)/2, -kernel_size(2)/2:kernel_size(2)/2);
Rgrid=sqrt((Xgrid.^2+Ygrid.^2));
kernel = pdf('Normal', Rgrid, 0, occupancy_std);
kernel = kernel./sum(sum(kernel));
smoothed_map = conv2(input_map, kernel, 'same'); % smoothing 

end

