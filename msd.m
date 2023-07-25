%% calculate msd 
% will write for a specific file, then refactor into helper functions?
close all; clf; clear
set(0,'DefaultFigureWindowStyle','docked')
set(groot, 'defaultAxesFontSize', 12)
%folder = "test8";
v0_arr=["0.01"];% "0.02" "0.04" "0.08" "0.16"];
k_ecm_arr=["0.005"];% "0.05" "0.5" "5"];
k_off_arr=["1.0"];
att_arr=["0.001"];% "0.01" "0.1"];
for ii=1:length(v0_arr)
    v0 = v0_arr(ii);
    for jj=1:length(k_ecm_arr)
        k_ecm = k_ecm_arr(jj);
        for kk=1:length(k_off_arr)
            k_off = k_off_arr(kk);
            for ll=1:length(att_arr)
                att = att_arr(ll);
                figure(1); clf; hold on; tiledlayout(2, 2);
                folder = "pipeline\cells\psm" + "\psm_calA01.0_phi0.74_tm10.0_v0"+v0+"_t_abp1.0k_ecm"+k_ecm+...
                    "k_off"+k_off+"/_N40_dur1000_att"+att+"_start1_end1_sd1";
                output_folder = "output\"+extractAfter(folder, 9);
                
                filename = folder + ".msd";
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
                %figure(1)
                nexttile;
                plot(time(2:maxLag), MSD(2:end),'k')
                xlabel("Time (min)")
                ylabel("MSD(t)/(25\pi \mum^2)")
                
                %% plot cell tracks
                % divide cellPos (# timesteps x 2N) into xCoords and yCoords, which are (#
                % timesteps x N)
                xCoords = cellPos(:,1:2:end);
                yCoords = cellPos(:,2:2:end);
                
                % Plot cell tracks
                %figure(2);clf;hold on;
                nexttile; hold on;
                for i = 1:size(xCoords, 2)
                    plot(xCoords(:, i), yCoords(:, i), '-','linewidth', 2);
                end
                axis equal off;

                %% for each timestep, plot the voronoi diagram to make sure it works
                % might need to remove last cell which is boundary, not
                % sure right now.
                nexttile; hold on;
                for i = 1:size(xCoords,1)
                    n = size(xCoords,2);
                    [v, c] = voronoin([xCoords(i,:)' yCoords(i,:)']);
                    vn = zeros(n, n);
                    for j=1:n
                        for k=j+1:n
                            s = size(intersect(c{j}, c{k}));
                            if (1 < s(2))
                                vn(j,k) = 1;
                                vn(k,j) = 1;
                            end
                        end
                    end
                end
                
                %% plot neighbor exchange count vs time
                filename = folder + ".cellContact";
                delimiter = ','; % Delimiter used within each matrix
                
                % Read the entire file as a single string
                fileContent = fileread(filename);
                
                % Split the string into individual matrices using newline separations
                matricesStr = strsplit(fileContent, '\n\n');
                
                % Remove empty cell from last newline character
                nonemptyIndices = ~cellfun(@isempty, matricesStr);
                matricesStr = matricesStr(nonemptyIndices);
                
                NE = zeros(1, numel(matricesStr));
                
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
                
                    % binarize matrix so that any vertex contact is counted as a cell
                    % contact, and there can be no more than 1 cell contact between cells 
                    matrix(matrix >= 1) = 1;
                
                    % skip to next iteration if on first, since we want differences in time
                    if i == 1
                        continue;
                    else
                        % store matrix to compare next timestep
                        oldmatrix = matrix;
                    end
                
                    % Process the matrix
                    newmat = matrix - oldmatrix;
                    NE(i) = sum(newmat, 'all');
                end
                
                %figure()
                nexttile;
                plot(time, NE,'k')
                xlabel("Time (min)")
                ylabel("NE count")
                
                print(output_folder+"msd_tracks.eps",'-depsc')
            end
        end
    end
end