%% calculate msd and other dynamic quantities
% will write for a specific file, then refactor into helper functions?
close all; clf; clear
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