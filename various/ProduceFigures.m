% Produce Figures
% matlab script that we will use to produce all the figures.

clear; clc; 
addpath(genpath('export_fig'))
addpath(genpath('Analysis'))
addpath(genpath('Data'))
addpath(genpath('Functions'))

% Say for example we wanted to create a figure of the pressure head at a
% given point in time.

loadString = 'MitchTesting';
folder = 'Analysis';

% LOAD -> 7_MitchTesting.mat

nameString = loadName(loadString, folder, 14);
load(nameString); clear loadString folder nameString

% Calculate S and ignore k.

[S, ~] = calcSk(h(:,7), cfg);

% Calculate psi and convert to psi Relative
psi = calcPsi(h, S, cfg);
psiRel = psi./cfg.PsatNode;

% reconstruct matrix size. 

psiRel = vec2mat(psiRel, cfg.r)';

% create beautiful plot:

surf(cfg.X, cfg.Z, psiRel)
view(2)

cbar = colorbar;
caxis([0 1])

cbar.TickLabelInterpreter = 'Latex';

ylabel(cbar, '$\psi/\psi_{\mbox{sat}}$', 'Interpreter', 'Latex', 'FontSize', 22)

ax = gca;
ax.TickLabelInterpreter = 'Latex';
ax.FontSize = 18;

ylabel('$z$ [m]', 'Interpreter', 'latex', 'FontSize', 22);
xlabel('$x$ [m]', 'Interpreter', 'latex', 'FontSize', 22);

set(gcf, 'Units', 'points', ...
    'Position', [10 10 1100 600], ... % THIS 1100 and 600 are what determine the size.
    'PaperUnits', 'centimeters', ...
    'PaperSize', [14 7], ...
    'Color', 'w')

%           Folder      / Name.pdf         quality (setting to 101 is best)
export_fig FiguresEXPORT/ExampleFigure.pdf -q101

function nameString = loadName(loadString, folder, ID)

    filesFull = dir(fullfile(folder, ['*', loadString, '*.mat']));

    fileCell = {filesFull.name};

    if nargin == 3
        warning('Specified ID chosen. This is not reliable, please check loaded file.')
    end

    if length(fileCell) > 1 && nargin < 3
        fprintf('\nMultiple files found.\n')
        fprintf('Please select from the following\n\n')
        for j = 1:length(fileCell)
            fprintf('\tID %g: %s\n', j, fileCell{j})
        end
        prompt = '\nPlease select a file ID Number:  ';
        fileID = input(prompt);

        while fileID > length(fileCell)
            prompt = 'Please select a valid file ID Number:  ';
            fileID = input(prompt);
        end

        fileToLoad = fileCell{fileID};
    elseif nargin == 3
        fileToLoad = fileCell{ID};
    else
        fileToLoad = fileCell{1};
    end

    fprintf('\n\nYou chose: %s\n\n', fileToLoad)

    nameString = fileToLoad;

end
