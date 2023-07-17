%% calculate msd 
% will write for a specific file, then refactor into helper functions?
close all; clf; clear
set(0,'DefaultFigureWindowStyle','docked')
filename = "test8.msd";
data = load(filename);
time = data(:,1);
% cellPos is (# timesteps, d*N), where d=2 for my 2D simulations
cellPos = data(:,2:end);
numCells = size(cellPos, 2) / 2;
numSteps = size(cellPos, 1);

maxLag = round(numSteps/4);  % Maximum time lag to consider
msd_per_particle_per_lag = zeros(numCells, maxLag);

for lag = 1:maxLag
    displacements = cellPos(lag+1:end, :) - cellPos(1:end-lag, :);
    msd_per_particle = sum(displacements.^2, 2);
    msd_per_particle_per_lag(:, lag) = mean(msd_per_particle, 1);
end
MSD = mean(msd_per_particle_per_lag, 1);
figure()
plot(time(2:maxLag), MSD(2:end),'k')

%% plot cell tracks
% divide cellPos (# timesteps x 2N) into xCoords and yCoords, which are (#
% timesteps x N)
xCoords = cellPos(:,1:2:end);
yCoords = cellPos(:,2:2:end);

% Plot cell tracks
figure;
hold on;
for i = 1:size(xCoords, 2)
    plot(xCoords(:, i), yCoords(:, i), '-','linewidth', 2);
end
hold off;
axis equal off;

%% plot neighbor exchange count vs time
filename = "test8.cellContact";
delimiter = ','; % Delimiter used within each matrix

% Read the entire file as a single string
fileContent = fileread(filename);

% Split the string into individual matrices using newline separations
matricesStr = strsplit(fileContent, '\n\n');

% Remove empty cell from last newline character
nonemptyIndices = ~cellfun(@isempty, matricesStr);
matricesStr = matricesStr(nonemptyIndices);

NE = zeros(1, size(matricesStr,2));

% Process each matrix
for i = 1:numel(matricesStr)
    % Split the matrix string into rows
    rowsStr = strsplit(matricesStr{i}, '\n');
    rowsStr(cellfun('isempty', rowsStr)) = [];
    
    % Split rows into values
    values = cellfun(@(x) strsplit(x, delimiter), rowsStr, 'UniformOutput', false);
    
    % Convert values to a numeric matrix
    % matrix(i,j) = # vertex contacts between cell i and cell j at a
    % particular time 
    matrix = str2double(vertcat(values{:}));
    matrix(matrix >= 1) = 1;

    % move to next iteration if we're on first, since we want differences
    % in time
    if i == 1
        continue;
    end

    % Process the matrix
    newmat = matrix - oldmatrix;

    % store matrix to compare next timestep
    oldmatrix = matrix;
end
%%
