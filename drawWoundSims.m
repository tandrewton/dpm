%%Draw sims for laserAblation.cpp
% output is a movie made from stitching the position file frames together
% different from drawLoadingSims.m because it plots psi information
%pwd should give ~/Documents/YalePhd/projects/dpm

%function drawWoundSims(N, strainRate_ps, calA0, smooth, deltaSq, d_flag, att, boundaryType) %uncomment if using function call to pipeline data
%isTestData = false; %uncomment if using function call to pipeline data

isTestData = true; %uncomment if using test data
addpath('/Users/AndrewTon/Documents/YalePhD/projects/dpm/bash')
addpath('C:\Users\atata\projects\dpm\bash')
addpath('/Users/AndrewTon/Documents/YalePhD/projects/dpm/matlab_funcs')
addpath('C:\Users\atata\projects\dpm\matlab_funcs')

%CHANGE THESE PARAMETERS                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   NEEDED

runType = "ablate";
%N="40";
ndelete="6";
%calA0="1.10";
%strainRate_ps="0.001";
%deltaSq = "2.0";
%k_a = "1.0";
k_l = "1.0";
k_ps = "1.0"; %purse-string spring constant
k_lp = "4.0"; %lamellipodia spring constant
%smooth = "1";
tau_lp = "1.0"; %lamellipodia lifetime
%d_flag = "0.0"; %lamellipodia max length
%boundaryType = "0"; 
%att="0.2";
B="1.0";
Dr0="0.5";
boolCIL="0";
Duration="800";
FSKIP = 1;

etaStr = " ";
startSeed = 1;
max_seed = 1;
no_plots = 1;
makeAMovie = 0; %if makeAMovie is 0, then plot every frame separately and dont save a movie object
%plotCells = makeAMovie; % if plotCells is 0, then skip plotting altogether
plotCells = 1;
set(0,'DefaultFigureWindowStyle','docked')
showPeriodicImages = 0;
showWoundAndShapeProperties = 1; 


showverts = 1;
showBoundaries = 0;
showcirculoline = 0; % show line segments of circulo-lines
isReadAndPlotTrajectoryQualities = 1; % read nvestr and plot associated quantities
att_range = 0.0;
showArea = 1;
showQuiver = 0;
walls = 0;
showCustomView = 0; % specific choice of coordinates to zoom in on for movie
viewLeft = 0.8;
viewRight = 1.6;
viewTop = 1.8;
viewBottom = 1.0;

%disable showVoid if using printConfig on its own, outside of
%dampedNVE/dampedNP0 routines
showGlobalIndex = 1;
showVoid = 0;
showVoidBlack = 0; % print void in larger black circles to see easier
showVoidLite = 1; % print void, but in a way that works with printConfiguration on its own
showCornersOrEdges = 0;
showPurseString = 0;
showProtrusion = 1;
showShapeHistogram = 0;
 
%PC directory
%pc_dir = "/Users/AndrewTon/Documents/YalePhD/projects/dpm/";
pc_dir="C:\Users\atata\projects\dpm\";
%pipeline is the location of data generated during simulations
subdir_pipeline = pc_dir + "pipeline/cells/"+runType+"/";

%output is location of results of this postprocessing
subdir_output = pc_dir + "output/cells/"+runType+"/";
mkdir(subdir_pipeline);
mkdir(subdir_output);

load("polyBoundary.txt"); % load boundaries of polygon walls
%txt = 'N = '+N+', NV = '+NV+', calA_o='+calA+', att='+att+', B='+B;
txt='test';

fnum = 1;
figure(13), clf, hold on, box on;
for seed = startSeed:max_seed
    if (isTestData)
        run_name = runType+txt;     
        pipeline_dir =  subdir_pipeline + run_name + "/";
        output_dir = subdir_output + run_name + "/";
        mkdir(pipeline_dir)
        mkdir(output_dir)
        fileheader=run_name + "_seed" + seed;
        fileheader_short = fileheader;
        nvestr = pc_dir+'test.pos';
        energystr = pc_dir+'test.energy';
        stressstr = pc_dir+'test.stress';
        boundaryStr = pc_dir+'test.void';
        edgeStr = pc_dir+'test.edge';
        purseStr = pc_dir+'test.purseString';
        voidAreaStr = pc_dir+'test.voidArea';
        innerStr = pc_dir+ 'test.innerCellShape';
        bulkStr = pc_dir+ 'test.bulkCellShape';
        woundPropertiesStr = pc_dir+ 'test.woundProperties';
        innerAndBulkCellIDStr = pc_dir+'test.cellID';
    else
        run_name =runType+"_A0"+calA0+"_k_l"+k_l+"_w_ps"+strainRate_ps+ ...
            "_dsq"+deltaSq+"_k_ps"+k_ps+"_k_lp"+k_lp+...
            "_t_lp"+tau_lp+"_d_flag"+d_flag+"_bd"+boundaryType+"_sm"+smooth;
        pipeline_dir =  subdir_pipeline + run_name + "/";
        output_dir = subdir_output + run_name + "/";
        mkdir(pipeline_dir)
        mkdir(output_dir)
        fileheader=run_name +"_N"+N+"_Dur"+Duration+"_att"+att+"_sd"+ ...
            startSeed+"_sd"+max_seed+"_sd"+seed;
        fileheader_short = "_N"+N+"_Dur"+Duration+"_att"+att+"_sd"+seed;
        nvestr = pipeline_dir+fileheader+'.pos';
        energystr = pipeline_dir+fileheader+'.energy';
        stressstr = pipeline_dir+fileheader+'.stress';
        boundaryStr = pipeline_dir+fileheader+ ".void";
        edgeStr = pipeline_dir+fileheader+ '.edge';
        purseStr = pipeline_dir+fileheader+ '.purseString';
        voidAreaStr = pipeline_dir+fileheader+ '.voidArea';
        innerStr = pipeline_dir+fileheader+ '.innerCellShape';
        bulkStr = pipeline_dir+fileheader+ '.bulkCellShape';
        woundPropertiesStr = pipeline_dir+fileheader+ '.woundProperties';
        innerAndBulkCellIDStr = pipeline_dir+fileheader+ '.cellID';
    end
    
    cd(output_dir)
    
    if (~no_plots)
        figure(11); clf; hold on;
        stress = load(stressstr);
        plot(stress(:,1), stress(:,3), 'r-', 'linewidth',2, 'DisplayName',...
            '$S_{xx}$');
        plot(stress(:,1), stress(:,4), 'b-', 'linewidth',2, 'DisplayName',...
            '$S_{yy}$');
        plot(stress(:,1), stress(:,5),'k-','linewidth',2, 'DisplayName',...
            '$S_{xy}$');
        xlabel('$\tau$','Interpreter','latex');
        ylabel('Stress','Interpreter','latex');
        legend('Location', 'southeast', 'Interpreter', 'latex');
        ax = gca;
        ax.FontSize = 24;
        if seed == max_seed 
         %annotation('textbox',[.2 .5 .5 .3],'String', txt, 'Edgecolor','none')
         saveas(gcf, 'Stress'+runType+fileheader_short+'_'+max_seed+ ...
             '.eps', 'epsc')
        end
    
        figure(12); clf;hold on;
        plot(stress(:,1), stress(:,6), 'r--', 'linewidth',2, 'DisplayName',...
            '$Sh_{xx}$');
        plot(stress(:,1), stress(:,7), 'b--', 'linewidth',2, 'DisplayName',...
            '$Sh_{yy}$');
        plot(stress(:,1), stress(:,8),'k--','linewidth',2, 'DisplayName',...
            '$Sh_{xy}$');
        xlabel('$\tau$','Interpreter','latex');
        ylabel('Stress','Interpreter','latex');
        legend('Location', 'southeast', 'Interpreter', 'latex');
        ax = gca;
        ax.FontSize = 24;
        if seed == max_seed 
         %annotation('textbox',[.2 .5 .5 .3],'String', txt, 'Edgecolor','none')
         saveas(gcf, 'ShapeStress'+runType+fileheader_short+'_'+max_seed+ ...
             '.eps', 'epsc')
        end
    
        figure(14); clf; hold on 
        energy = load(energystr);
        U = energy(:,3);
        K = energy(:,4);
        %U_ps = energy(:,5);
        %U_crawling = energy(:,6);
        plot(energy(:,1), K, 'r-', 'linewidth',2, 'DisplayName',...
            '$K$');
        %plot(energy(:,1), U_ps, 'b-', 'linewidth',2, 'DisplayName',...
        %    '$U_{ps}$');
        %plot(energy(:,1), U_crawling,'k-','linewidth',2, 'DisplayName',...
        %    '$U_{crawling}$');
        plot(energy(:,1), U,'r--','linewidth',2, 'DisplayName',...
           '$U$');
        plot(energy(:,1), K+U,'r.','linewidth',2, 'DisplayName',...
           '$K+U$');
        xlabel('$\tau$','Interpreter','latex');
        ylabel('Energy','Interpreter','latex');
        legend('Location', 'southeast', 'Interpreter', 'latex');
        ax = gca;
        ax.FontSize = 24;
        if seed == max_seed 
         saveas(gcf, 'Energy'+runType+fileheader_short+'_'+max_seed+ ...
             '.eps', 'epsc')
        end
    
        if showArea
            figure(15)
            voidArea = load(voidAreaStr);
            plot(voidArea(:,1), voidArea(:,2), 'linewidth', 4)
            %ylim([0 voidArea(1,2)])
            xlabel('$t/\tau$','Interpreter','latex','fontsize', 24);
            ylabel('Area','Interpreter','latex','fontsize', 24);
            %set(gca,'Yscale','log')
            if seed == max_seed 
                saveas(gcf, 'VoidArea'+runType+fileheader_short+'_'+max_seed+'.eps', 'epsc')
            end
        end

        if showWoundAndShapeProperties
            innerCellShape = load(innerStr);
            bulkCellShape = load(bulkStr);
            woundProperties = load(woundPropertiesStr)
            cellID = load(innerAndBulkCellIDStr);
            figure(16); clf; hold on;
            %plot(innerCellShape(:,1), mean(innerCellShape(:,2:end),2),  ...
            %    'linewidth', 4, 'DisplayName', "inner shapes")
            %plot(bulkCellShape(:,1),mean(bulkCellShape(:,2:end),2), ...
            %    'linewidth', 4, 'DisplayName', "bulk shapes")
            % cellID row = [ci inInitialWoundNeighbors inFinalWoundNeighbors]
            % we want to access inFinalWoundNeighbors of bulkCellShape
            % bulkCellShape row = [time shape(0) shape(1) ... shape(NCELLS)]
            timeAndInnerShapes = bulkCellShape(:,[1; cellID(:,3)]==1);
            timeAndOuterShapes = bulkCellShape(:,[1; ~cellID(:,3)]==1);
            plot(timeAndInnerShapes(:,1), mean(timeAndInnerShapes(:,2:end),2),  ...
                'linewidth', 4, 'DisplayName', "inner shapes")
            plot(timeAndOuterShapes(:,1),mean(timeAndOuterShapes(:,2:end),2), ...
                'linewidth', 4, 'DisplayName', "bulk shapes")
            xlabel('$t/\tau$','Interpreter','latex','fontsize', 24);
            ylabel('Shape','Interpreter','latex','fontsize', 24);
            legend('location','northwest','fontsize', 14)
            if seed == max_seed 
                saveas(gcf, 'cellShapes'+runType+fileheader_short+'_'+max_seed+'.eps', 'epsc')
            end
        end
    end

    if (isReadAndPlotTrajectoryQualities)
        nvestr
        % read in position data
        [trajectoryData, cell_count] = readDPMClassPosOutput(nvestr);
    
        % get number of frames
        NFRAMES = trajectoryData.NFRAMES;
        NCELLS = trajectoryData.NCELLS;
        nv = trajectoryData.nv;
        flag = trajectoryData.flag;
        flagX = trajectoryData.flagPosX;
        flagY = trajectoryData.flagPosY;
        time = trajectoryData.time;
        %NVTOT = sum(nv);
        
        %if L is constant, use the next 3 lines
        %L = trajectoryData.L(1,:);
        %Lx = L(1);
        %Ly = L(2);
    
        % get cell colors
        %[nvUQ, ~, IC] = unique(nv);
       % NUQ = length(nvUQ);
        %cellCLR = jet(NUQ);
    
        % get frames to plot
        FSTART = 1;
        FSTEP = FSKIP;
        FEND = NFRAMES;
    
        if makeAMovie == 1
            moviestr = runType+fileheader_short+'.mp4';
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
    
        if (~no_plots)
            figure(fnum), clf, hold on, box on;
            
            if showVoid
                voidLocations = readDataBlocks(boundaryStr, 2);
                voidArea = zeros(NFRAMES,1);
            end
            if showCornersOrEdges
                edgeLocations = readDataBlocks(edgeStr, 3);
            end
            if showPurseString
                purseLocations = readDataBlocks(purseStr, 2);
            end
        
    %         zc = trajectoryData.zc; % # cell contacts
    %         zv = trajectoryData.zv; % # vertex contacts
    %         figure(2)
    %         yyaxis left
    %         plot(time-time(1), mean(zc,2), 'r', 'linewidth', 4,'DisplayName', 'zc')
    %         xlabel('$t/\tau$','Interpreter','latex','fontsize', 36);
    %         ylabel('$\langle zc \rangle$','Interpreter','LaTeX','Fontsize',36)
    %     
    %     
    %         yyaxis right
    %         plot(time-time(1), mean(zv,2), 'b', 'linewidth', 4,'DisplayName', 'zv')
    %         %xlabel('$t/\tau$','Interpreter','latex','fontsize', 36);
    %         ylabel('$\langle zv \rangle$','Interpreter','LaTeX','Fontsize',36)
    %         legend('location','Southeast')
    %     
    %         saveas(gcf, 'contacts'+runType+fileheader_short+'_'+max_seed+ ...
    %              '.eps', 'epsc')
        
            figure(3)
            area = mean(trajectoryData.area,2);
            perimeter = mean(trajectoryData.perimeter,2);
            shape = 1/(4*pi)*perimeter.^2./area;
            size(shape)
            plot(time-time(1), shape, 'k', 'linewidth', 4)
            xlabel('$t/\tau$','Interpreter','latex','fontsize', 36);
            ylabel('$\mathcal{A}$','Interpreter','LaTeX','Fontsize',36)
            saveas(gcf, 'shapes'+runType+fileheader_short+'_'+max_seed+...
                '.eps','epsc')
        end
        if (plotCells == 1)
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
                    vradtmp = vrad{nn}*(1 + att_range);
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
                            for xx = itLow:itHigh
                                for yy = itLow:itHigh
                                    rectangle('Position',[xplot+xx*Lx, yplot + yy*Ly, 2*vradtmp(vv), 2*vradtmp(vv)],'Curvature',[1 1],'EdgeColor','k','FaceColor',clr);
                                    if showGlobalIndex
                                        text(xtmp(vv), ytmp(vv), num2str(gitmp(vv)), 'FontSize', 6);
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
                                patch(cornerx, cornery, cornerx./cornerx, 'black','EdgeColor','blue', 'LineWidth',2)
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
                                patch('Faces',finfo,'vertices',vpos,'FaceColor',clr,'EdgeColor','k','linewidth',2);
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
                if (showProtrusion)
                    plot(flagX(ff,:)./flag(ff,:), flagY(ff,:)./flag(ff,:), 'ro', 'linewidth', 2);
                end
                %plot([polyBoundary(1,1:2:end) polyBoundary(1,1)],...
                %    [polyBoundary(1,2:2:end) polyBoundary(1,2)], 'k','linewidth', 4)
                ann = annotation('textbox', [.42 .05 .6 .05],'interpreter', 'latex',...
                    'String', etaStr, 'Edgecolor','none');
                ann.FontSize = 30;
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
                    viewScale = 1.2;
                    viewLx = viewScale*Lx;
                    viewLxLow = -(viewScale-1)*Lx;
                    viewLy = viewScale*Ly;
                    viewLyLow = -(viewScale-1)*Ly;
                    ax.XLim = [-(viewScale-1) viewScale]*Lx;
                    ax.YLim = [-(viewScale-1) viewScale]*Ly;
                    % plot box
                    plot([viewLxLow viewLx viewLx viewLxLow viewLxLow], [viewLyLow viewLyLow viewLy viewLy viewLyLow], 'k-', 'linewidth', 1.5);
                end
                
                annotationStr = "$$t/\tau$$ = "+time(ff);
                %annotationStr = "frame = "+ff;
                annotation('textbox',[0.48, 0.5, 0, 0],...
                    'interpreter', 'latex', 'String', annotationStr, 'Edgecolor','none', 'FitBoxToText','on');
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
                        purseLocations{ff}(1:2:end,2), 50, 'blue', 'o','MarkerFaceColor','blue');
                    scatter(purseLocations{ff}(2:2:end,1),...
                        purseLocations{ff}(2:2:end,2), 25, 'red', 'x');
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
                    currframe = getframe(gcf);
                    writeVideo(vobj,currframe);
                end
            end
        end
    
    
        % close video object
        if makeAMovie == 1
            close(vobj);
        end
        cd ../../../../
        if showVoid
            figure(15); clf, hold on, box on
            time= time - time(1);
            plot(time,voidArea, 'DisplayName', 'A', 'linewidth', 3);
            plot(time, gradient(voidArea(:))./gradient(time(:)), 'DisplayName', 'dA/dt','linewidth', 3);
            xlabel('$\tau$','Interpreter','latex');
            ylabel('A, dA/dt');
            saveas(gcf, 'Area'+runType+fileheader_short+'_'+max_seed+'.eps', 'epsc')
            legend
        end
    end
end
