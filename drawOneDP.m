%%Draw sims for laserAblation.cpp
% output is a movie made from stitching the position file frames together
% different from drawLoadingSims.m because it plots psi information
%pwd should give ~/Documents/YalePhd/projects/dpm
clear all
close all
clc
% function drawWoundSims(N, NV, ndelete, calA, kl, att, v0, B,...
%     Dr0,duration, boolCIL, showPeriodicImages, showverts, isTestData)
isTestData = true;
addpath('/Users/AndrewTon/Documents/YalePhD/projects/dpm/bash')
addpath('/Users/AndrewTon/Documents/YalePhD/projects/dpm/matlab_funcs')

%CHANGE THESE PARAMETERS                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   NEEDED

%PC directory
pc_dir = "/Users/AndrewTon/Documents/YalePhD/projects/dpm/";
showverts = 0;
showSprings = 0;
makeAMovie = 1;
set(0,'DefaultFigureWindowStyle','docked')
showPeriodicImages = 0;
fnum=1;
showGlobalIndex = 0;
showcirculoline = 0;
showQuiver = 0;
walls=0;

showCustomView = 0;
viewLeft = 0.8;
viewRight = 1.6;
viewTop = 1.8;
viewBottom = 1.0;

showVoid=0;
showCornersOrEdges = 0;
showPurseString = 0;
showShapeHistogram = 0;

%txt = 'N = '+N+', NV = '+NV+', calA_o='+calA+', att='+att+', B='+B;
txt='test';

figure(1), clf, hold on, box on;  

nvestr = pc_dir + "oneDP.pos";
%%
[trajectoryData, cell_count] = readDPMClassPosOutput(nvestr);

% get number of frames
NFRAMES = trajectoryData.NFRAMES;
NCELLS = trajectoryData.NCELLS;
nv = trajectoryData.nv;
flag = trajectoryData.flag;
flagX = trajectoryData.flagPosX;
flagY = trajectoryData.flagPosY;
time = trajectoryData.time;

% get frames to plot
FSTART = 1;
FSTEP = 1;
FEND = NFRAMES;


if makeAMovie == 1
    moviestr = pc_dir+'oneDP.mp4';
    if exist(moviestr, 'file')==2
      delete(moviestr);
    end
    vobj = VideoWriter(moviestr, 'MPEG-4');
    vobj.Quality = 100;
        
    vobj.FrameRate = 5;
    open(vobj);
end

if showPeriodicImages == 1
    itLow = -1;
    itHigh = 1;
    boxAxLow = -0.25;
    boxAxHigh = 1.25;
else
    itLow = 0;
    itHigh = 0;
    boxAxLow = 0;
    boxAxHigh = 1;
end

for ff = FSTART:FSTEP:FEND
    %nv can change, so recompute color map each frame
    [nvUQ, ~, IC] = unique(nonzeros(nv(ff,:).*(flag(ff,:)+1)));
    %IC = IC * 0 + 1; % <- use for single colored configurations
    NUQ = length(nvUQ);
    %NUQ = 8; % 1<- use for single colored configurations
    %IC = [1;2]; % <- use for 2 cell testing simulations
    %cellCLR = jet(2); % <- use for 2 cell testing simulations
    cellCLR = jet(NUQ);
    NCELLS = cell_count(ff);
    
    area = trajectoryData.area(ff,:);
    perimeter = trajectoryData.perimeter(ff,:);
    shape = 1/(4*pi)*perimeter.^2./area;

    if ~makeAMovie
        fnum = fnum+1;
    end
    figure(fnum), clf, hold on, box on;

    fprintf('printing frame ff = %d/%d\n',ff,FEND);

    % get cell positions
    xpos = trajectoryData.xpos(ff,:);
    ypos = trajectoryData.ypos(ff,:);
    gi = trajectoryData.gi(ff,:);
    l0 = trajectoryData.l0(ff,:);
    vrad = trajectoryData.vrad(ff,:);
    psi = trajectoryData.psi(ff,:);
    

    %if L is not constant, use the next 3 lines
    L = trajectoryData.L(ff,:);
    Lx = L(1);
    Ly = L(2);

    for nn = 1:NCELLS
        xtmp = xpos{nn};
        ytmp = ypos{nn};
        gitmp = gi{nn};
        l0tmp = l0{nn};
        vradtmp = vrad{nn};
        psitmp = psi(nn);
        costmp = cos(psitmp);
        sintmp = sin(psitmp);

        clr = cellCLR(IC(nn),:);

        cx = mean(xtmp);
        cy = mean(ytmp);
        if showverts == 1
            for vv = 1:nv(ff,nn)
                xplot = xtmp(vv) - vradtmp(vv);
                yplot = ytmp(vv) - vradtmp(vv);
                if showSprings
                    [xs,ys] = spring(xtmp(vv), ytmp(vv),...
                        xtmp(mod(vv, nv(ff,nn))+1), ytmp(mod(vv, nv(ff,nn))+1)...
                        , 3, l0tmp(vv), vradtmp(vv)/2);
                    plot(xs,ys,'LineWidth', 1,'Color' ,'black');
                    plot(xtmp(vv),ytmp(vv), '.r', 'color', 'blue', 'markersize', 50)
                else 
                    for xx = itLow:itHigh
                        for yy = itLow:itHigh
                            rectangle('Position',[xplot+xx*Lx, yplot + yy*Ly, 2*vradtmp(vv), 2*vradtmp(vv)],'Curvature',[1 1],'EdgeColor','k','FaceColor',clr);
                            if showGlobalIndex
                                %text(xtmp(vv), ytmp(vv), num2str(gitmp(vv)), 'FontSize', 6);
                                text(xtmp(vv)-0.02, ytmp(vv)+0.005, 'x', 'FontSize', 20,'color','red');
                            end
                        end
                    end
                end
            end

            for vv = 1:nv(ff,nn)
                xplot = xtmp(vv) - vradtmp(vv);
                yplot = ytmp(vv) - vradtmp(vv);
                if showcirculoline == 1% calculate coordinates of a rectangle representing the line segment between successive vertices in a DP
                    vnext = mod(vv, nv(ff,nn))+1;
                    xtmpnext = xtmp(vnext);
                    ytmpnext = ytmp(vnext);
                    rx = xtmpnext - xtmp(vv);
                    ry = ytmpnext - ytmp(vv);
                    if (rx == 0) % if line is vertical, perpendicular is <1,0>
                       perp_x = 1;
                       perp_y = 0;
                    else
                         % dot product of r and perp = 0, so perp is perpendicular to r
                        perp_x = -ry/rx;
                        perp_y = 1;
                    end
                    norm = sqrt(perp_x^2 + perp_y^2);
                    perp_x = perp_x / norm;
                    perp_y = perp_y / norm;
                    % calculate 4 coordinates of a rectangle
                    % for the segment
                    offsetx = vradtmp(vv)*perp_x;
                    offsety = vradtmp(vv)*perp_y;
                    cornerx = [xtmp(vv)-offsetx, xtmp(vv)+offsetx, ...
                        xtmpnext+offsetx, xtmpnext-offsetx];
                    cornery = [ytmp(vv)-offsety, ytmp(vv)+offsety,...
                        ytmpnext+offsety, ytmpnext-offsety];
                    %patch(cornerx, cornery, cornerx./cornerx, 'black','EdgeColor','blue', 'LineWidth',2)
                    patch(cornerx, cornery, cornerx./cornerx, clr,'EdgeColor',' black', 'LineWidth',2)
                end
            end
        else
            rx = xtmp - cx;
            ry = ytmp - cy;
            rads = sqrt(rx.^2 + ry.^2);
            xtmp = xtmp + 0.4*l0tmp(1)*(rx./rads);
            ytmp = ytmp + 0.4*l0tmp(1)*(ry./rads);
            for xx = itLow:itHigh
                for yy = itLow:itHigh
                    vpos = [xtmp + xx*Lx, ytmp + yy*Ly];
                    finfo = [1:nv(ff,nn) 1];
                    %disp("finfo is "+ finfo)
                    patch('Faces',finfo,'vertices',vpos,'FaceColor',clr,'EdgeColor','k','linewidth',1);
                end
            end
        end
        


        %plot arrows representing directors
        if (showQuiver == 1)
            for xx = itLow:itHigh
                for yy = itLow:itHigh
                    quiver(cx + xx*Lx, cy + yy*Ly, costmp, sintmp,...
                        'r', 'LineWidth', 2, 'MaxHeadSize', 2);
                end
            end
        end
    end

    axis equal;
    ax = gca;
    %ax.XTick = [];
    %ax.YTick = [];
    if showPeriodicImages == 1
        ax.XLim = [boxAxLow boxAxHigh]*Lx;
        ax.YLim = [boxAxLow boxAxHigh]*Ly;
        % plot box
        plot([0 Lx Lx 0 0], [0 0 Ly Ly 0], 'k-', 'linewidth', 1.5);
    elseif walls == 1
        ax.XLim = [0 1]*Lx;
        ax.YLim = [0 1]*Ly;
        % plot box
        plot([0 Lx Lx 0 0], [0 0 Ly Ly 0], 'k-', 'linewidth', 1.5);
    elseif showCustomView == 1
        ax.XLim = [viewLeft viewRight];
        ax.YLim = [viewBottom viewTop];
    else
        if (ff == 1)
            xlimit = xlim;
            ylimit = ylim;
            diff = 0.1* (xlimit(2) - xlimit(1));
            xlimit(1) = xlimit(1) - diff;
            xlimit(2) = xlimit(2) + diff;
            ylimit(1) = ylimit(1) - diff;
            ylimit(2) = ylimit(2) + diff;
        end
        ax.XLim = xlimit;
        ax.YLim = ylimit;
        axis off;
%         viewScale = 1.2;
%         viewLx = viewScale*Lx;
%         viewLxLow = -(viewScale-1)*Lx;
%         viewLy = viewScale*Ly;
%         viewLyLow = -(viewScale-1)*Ly;
%         ax.XLim = [-(viewScale-1) viewScale]*Lx;
%         ax.YLim = [-(viewScale-1) viewScale]*Ly;
        % plot box
        %plot([viewLxLow viewLx viewLx viewLxLow viewLxLow], [viewLyLow viewLyLow viewLy viewLy viewLyLow], 'k-', 'linewidth', 1.5);
    end
    
    annotationStr = "$$t/\tau$$ = "+time(ff);
    %annotationStr = "frame = "+ff;
    %annotation('textbox',[0.48, 0.5, 0, 0],...
    %    'interpreter', 'latex', 'String', annotationStr, 'Edgecolor','none', 'FitBoxToText','on');
    if showVoid
        if showVoidBlack 
            scatter(voidLocations{ff}(:,1), voidLocations{ff}(:,2),...
            20, 'black', 's','filled');
        else
            scatter(voidLocations{ff}(:,1), voidLocations{ff}(:,2),...
            10, 'red', '_');
        end
        offset = 0;
        nff = ff- offset;
        if (nff > length(voidLocations))
            continue
        end
        if (nff > 0)
            %plot(voidLocations{nff}(:,1)...
            %    , voidLocations{nff}(:,2)...
            %    ,'k-', 'linewidth', 3);
            %voidArea(ff) = polyarea(voidLocations{nff}(:,1), voidLocations{nff}(:,2));
        end
    end
    if showCornersOrEdges
        % scatter(edgeLocations{ff}(:,1), edgeLocations{ff}(:,2), 30, 'black', 'x');
        for k = 1:length(edgeLocations{ff}(:,1))
            % plot text centered on an (x,y) point with label equal to the 3rd column value 
            text(edgeLocations{ff}(k,1), edgeLocations{ff}(k,2),...
                num2str(edgeLocations{ff}(k,3)),...
                'HorizontalAlignment', 'Center',...
                'VerticalAlignment', 'Middle', 'FontSize', 8);
        end
    end
    if showPurseString
        scatter(purseLocations{ff}(1:2:end,1),...
            purseLocations{ff}(1:2:end,2), 20, 'blue', 'o','MarkerFaceColor','blue');
        scatter(purseLocations{ff}(2:2:end,1),...
            purseLocations{ff}(2:2:end,2), 10, 'red', 'x');
        purseLocs = purseLocations{ff};
        purseLength = length(purseLocs(1:2:end,1));
        for i=1:purseLength
            lowOdd = 2*(i-1)+1;
            %plot([purseLocs(lowOdd,1) purseLocs(lowOdd+1,1)]...
            %    ,[purseLocs(lowOdd,2) purseLocs(lowOdd+1,2)], 'LineWidth', 2, 'Color' ,'black');
            
            %[xs,ys] = spring(purseLocs(lowOdd,1),purseLocs(lowOdd,2)...
            %    ,purseLocs(lowOdd+1,1), purseLocs(lowOdd+1,2)...
            %    ,2, 0, 0.1);
            %plot(xs,ys,'LineWidth', 1,'Color' ,'black');
        end
    end

    if (showShapeHistogram)
        histEdges = 1.0:0.04:3.0;
        histAx = axes('Position', [.7 .7   .25 .25]);
        box on
        histogram(histAx, shape, histEdges);
        xlabel('$\mathcal{A}$','Interpreter','LaTeX','Fontsize',16)
        ylabel('P($\mathcal{A}$)','Interpreter','LaTeX','Fontsize',16)
        ylim([0 10])
        xticks([1.0:0.4:3.0])
        annotation('textbox',[.75 .8 .1 .1], 'edgecolor', 'none', 'string', "mean="+num2str(mean(shape)))
        annotation('textbox',[.75 .77 .1 .1], 'edgecolor', 'none', 'string', "std="+num2str(std(shape)))
    end
    % if making a movie, save frame
    if makeAMovie == 1
        set(gca,'visible','off')
        currframe = getframe(gcf);
        writeVideo(vobj,currframe);
    end
end
    
% close video object
if makeAMovie == 1
    close(vobj);
end
