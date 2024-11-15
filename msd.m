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
                % timesteps x N-1), -1 because removing boundary
                xCoords = cellPos(:,1:2:end-2);
                yCoords = cellPos(:,2:2:end-2);
                
                % Plot cell tracks
                %figure(2);clf;hold on;
                nexttile; hold on;
                for i = 1:size(xCoords, 2)
                    plot(xCoords(:, i), yCoords(:, i), '-','linewidth', 2);
                end
                axis equal off;

                %% for each timestep, plot the voronoi diagram to make sure it works
                nexttile; hold on;
                changeInContacts = zeros(1, size(xCoords,1)-1);
                neighborDistances = [];
                videoWriter = VideoWriter(output_folder+"video.mp4", 'MPEG-4');
                videoWriter.FrameRate = 5;
                open(videoWriter);
                n = size(xCoords,2);
                voronoiCenterPositions = zeros(size(xCoords,1), 2, n);
                for timeii = 1:size(xCoords,1)
                    [vx, vy] = voronoi(xCoords(timeii,:), yCoords(timeii,:));
                    figure(2);
                    plot(vx,vy,'k')
                    axis square equal
                    xlim([-2, 20])
                    ylim([-2, 20])
                    writeVideo(videoWriter, getframe(gcf));
                    [v, c] = voronoin([xCoords(timeii,:)' yCoords(timeii,:)']);
                    v(isinf(v)) = nan;
                    vn = zeros(n, n);
                    % calculate c.o.m. of each polygon j, 
                    % make a dictionary for polygon j -> cellID[j]
                    % Calculate the center of mass (COM) for each Voronoi cell
                    cellComs = cell2mat(cellfun(@(indices) mean(v(indices, :),"omitnan"),...
                        c, 'UniformOutput', false))';

                    voronoiCenterPositions(timeii,:,:) = cellComs;

                    % Duplicate point coordinates to match with cellComs for vectorized distance calculation
                    pointCoords = [xCoords(timeii,:); yCoords(timeii,:)];
                    
                    % Calculate the Euclidean distance between each point and each cell's COM
                    %distances = sqrt(sum((pointCoords - cellComs).^2, 1));
                    distances = pdist2(pointCoords', cellComs');
                    
                    % Find the index of the closest cell (minimum distance) for each point
                    [~, pointCellDict] = min(distances, [], 1);
                    for j=1:n
                        for k=j+1:n
                            s = size(intersect(c{j}, c{k}));
                            if (1 < s(2))
                                %instead of j,k here, should be cellID[j]
                                %and cellID[k]
                                %vn(pointCellDict(j),pointCellDict(k)) = 1;
                                %vn(pointCellDict(k),pointCellDict(j)) = 1;
                                vn(j,k) = 1;
                                vn(k,j) = 1;
                                distance = sqrt((xCoords(timeii, j) - xCoords(timeii,k))^2 + ...
                                    (yCoords(timeii,j) - yCoords(timeii,k))^2);
                                neighborDistances = [neighborDistances distance];
                            end
                        end
                    end
                    if (timeii > 1)
                        %changeInContacts(i-1) = sum(vn - oldContactMatrix,"all");
                        changeInContacts(i-1) = sum(vn ~= oldVn,"all");
                    end
                    oldC = c;
                    oldV = v;
                    oldVn = vn;
                end
                % close(videoWriter);
                figure(1);
                plot(time(2:end), cumsum(abs(changeInContacts)))
                nexttile();
                histogram(neighborDistances);
                %print(output_folder+"msd_tracks.eps",'-depsc')
            end
        end
    end
end

figure(); clf; hold on;
for ii=1:39
    plot(voronoiCenterPositions(:,1,ii),voronoiCenterPositions(:,2,ii));
end