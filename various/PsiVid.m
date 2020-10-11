%SPHEREVID creates a video of the values for usage
% Author: Samuel Dudley - Queensland Institute of Medical Research 2018
clc, clear, close all

% Set how many frames per second.
fps = 10;

% SET THE OUTFILE NAME HERE
outfilename = fullfile('Videos', 'test_vid.mp4');

writerObj = VideoWriter(outfilename, 'MPEG-4');
writerObj.FrameRate = fps;
writerObj.Quality = 100;
open(writerObj)

outwidth = 1024;
outheight = 576;

hf = figure(888); clf
figpos = [50 50 outwidth outheight];
set(hf, 'position', figpos, 'color', [1 1 1], 'paperPositionMode', 'auto');

ax = axes('visible', 'off', ...
    'units', 'normalized', ...
    'position', [0.05 0.05 0.9 0.9], ...
    'color', 'w');

% MAKE SURE FOR RELATIVE SATURATION: 0 1
% head (h) -15 0
% H dunno.
caxis([0 100])
cbar = colorbar;
set(cbar,...
    'units', 'normalized', ...
    'FontSize', 14, ...
    'TickLabelInterpreter', 'LaTeX')

% SET THE NAME OF COLORBAR HERE:
ylabel(cbar, 'Name me bitch', ...
    'Interpreter', 'Latex', 'Fontsize', 20)

% LOAD DATA
load(fullfile('Analysis','55_Test4.mat'));

for jj = 1:round(0.7*length(t_total))
    % Calculate data
    [S, ~] = calcSk(h(:,jj), cfg);
    psi = calcPsi(h(:,jj), S, cfg);
    
    % Plot the stuff here. It should automatically go to figure 888
    surf(cfg.X, cfg.Z, vec2mat(psi./cfg.PsatNode*100, cfg.r)')
    caxis([0 100])
    cbar = colorbar;
    set(cbar,...
        'units', 'normalized', ...
        'FontSize', 14, ...
        'TickLabelInterpreter', 'LaTeX')
    ylabel(cbar, 'Saturation (\%)', ...
        'Interpreter', 'LaTeX', 'Fontsize', 20)
    view(2) % top down view
    title(['Percentage Saturation, Years: ', num2str(t_total(jj)/365)],...
        'Interpreter', 'LaTeX', 'FontSize', 20)
    drawnow
    
    % Now we write the video
    writeVideo(writerObj, getframe(hf));
end

close(writerObj);



    
    