function msPlotTrajectory

load('behav.mat');
load('ms.mat');

%behav.position = behav.positionRescale

[unique_behavtime, IAbehav, ICbehav]=unique(behav.time);
[unique_mstime, IAms, ICms]=unique(ms.time);

interpolated_X=interp1(behav.time(IAbehav),behav.position(IAbehav,1),ms.time(IAms));
interpolated_Y=interp1(behav.time(IAbehav),behav.position(IAbehav,2),ms.time(IAms));

figure;
    plot(interpolated_X,interpolated_Y,'color',[0 0.2 0.5]); hold on;
    %scatter(interpolated_X(logical(binarized_trace)),interpolated_Y(logical(binarized_trace)), [], [0.8 0 0],'.');
    daspect([1 1 1])
    ax = gca;
    ax.YDir = 'Reverse';
    set(ax, 'Visible','off');
    
end
