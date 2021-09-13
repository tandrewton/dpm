%%Draw sims for laserAblation.cpp
% output is a movie made from stitching the position file frames together
% different from drawLoadingSims.m because it plots psi information
%pwd should give ~/Documents/YalePhd/projects/dpm

% function drawWoundSims(N, NV, ndelete, calA, kl, att, v0, B,...
%     Dr0,duration, boolCIL, showPeriodicImages, showverts, isTestData)
isTestData = true;
addpath('/Users/AndrewTon/Documents/YalePhD/projects/dpm/bash')
addpath('/Users/AndrewTon/Documents/YalePhD/projects/dpm/matlab_funcs')

%CHANGE THESE PARAMETERS                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   NEEDED

runType = "ablate";
N="48";
ndelete="10";
calA="1.08";
att="0.3";
v0="0.5";
B="1.0";
Dr0="0.5";
boolCIL="0";
duration="5000";
FSKIP = 1;
eta = str2num(v0)^2 / str2num(Dr0) / str2num(att) / 2;
eta = num2str(eta);
%etaStr = "$$\eta = $$" + eta;
%if (isTestData)
etaStr = " ";
%end

startSeed = 1;
max_seed = 1;
makeAMovie = 1;
showPeriodicImages = 0;
showverts = 0;
showBoundaries = 0;
showQuiver = 0;
walls = 0;
showVoid = 1;
 
%PC directory
pc_dir = "/Users/AndrewTon/Documents/YalePhD/projects/dpm/";
%pipeline is the location of data generated during simulations
subdir_pipeline = pc_dir + "pipeline/cells/"+runType+"/";

%output is location of results of this postprocessing
subdir_output = pc_dir + "output/cells/"+runType+"/";
mkdir(subdir_pipeline);
mkdir(subdir_output);


%txt = 'N = '+N+', NV = '+NV+', calA_o='+calA+', att='+att+', B='+B;
txt='test';

figure(13), clf, hold on, box on;
for seed = startSeed:max_seed
    if (isTestData)
        run_name = runType+txt;     
        pipeline_dir =  subdir_pipeline + run_name + "/";
        output_dir = subdir_output + run_name + "/";
        mkdir(pipeline_dir)
        mkdir(output_dir)
        fileheader=run_name + "_seed" + seed;
        nvestr = pc_dir+'/pos.test';
        energystr = pc_dir+'/energy.test';
        stressstr = pc_dir+'/stress.test';
        boundaryStr = pc_dir+'void.test';
    else
        run_name =runType+"_N"+N+"_ndel"+ndelete+"_calA0"+calA+...
            "_att"+att+"_v0"+v0+"_B"+B+"_Dr0"+Dr0+"_CIL"+boolCIL+"_duration"+duration;     
        pipeline_dir =  subdir_pipeline + run_name + "/";
        output_dir = subdir_output + run_name + "/";
        mkdir(pipeline_dir)
        mkdir(output_dir)
        fileheader=run_name + "_seed" + seed;
        nvestr = pipeline_dir+fileheader+'.pos';
        energystr = pipeline_dir+fileheader+'.energy';
        stressstr = pipeline_dir+fileheader+'.stress';
        boundaryStr = pipeline_dir+fileheader+ ".void";
    end
    

    figure(13);
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
     saveas(gcf, output_dir + 'Stress'+startSeed+'_'+max_seed, 'epsc')
    end

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
        moviestr = output_dir + runType+fileheader+'seed_'+seed+'.mp4';
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

    fnum = 1;
    figure(fnum), clf, hold on, box on;
    if showVoid
        voidLocations = readDataBlocks(boundaryStr, 2);
        voidArea = zeros(NFRAMES,1);
    end
    
    for ff = FSTART:FSTEP:FEND
        %nv can change, so recompute color map each frame
        [nvUQ, ~, IC] = unique(nonzeros(nv(ff,:).*(flag(ff,:)+1)));
        NUQ = length(nvUQ);
        cellCLR = jet(NUQ);

        NCELLS = cell_count(ff);
        figure(fnum), clf, hold on, box on;
        fprintf('printing frame ff = %d/%d\n',ff,FEND);

        % get cell positions
        xpos = trajectoryData.xpos(ff,:);
        ypos = trajectoryData.ypos(ff,:);
        l0 = trajectoryData.l0(ff,:);
        psi = trajectoryData.psi(ff,:);

        %if L is not constant, use the next 3 lines
        L = trajectoryData.L(ff,:);
        Lx = L(1);
        Ly = L(2);

        for nn = 1:NCELLS
            xtmp = xpos{nn};
            ytmp = ypos{nn};
            l0tmp = l0(nn);
            psitmp = psi(nn);
            costmp = cos(psitmp);
            sintmp = sin(psitmp);

            clr = cellCLR(IC(nn),:);

            cx = mean(xtmp);
            cy = mean(ytmp);
            if showverts == 1
                for vv = 1:nv(ff,nn)
                    xplot = xtmp(vv) - 0.5*l0tmp;
                    yplot = ytmp(vv) - 0.5*l0tmp;
                    for xx = itLow:itHigh
                        for yy = itLow:itHigh
                            rectangle('Position',[xplot + xx*Lx, yplot + yy*Ly, l0tmp, l0tmp],'Curvature',[1 1],'EdgeColor','k','FaceColor',clr);
                        end
                    end
                end
            else
                rx = xtmp - cx;
                ry = ytmp - cy;
                rads = sqrt(rx.^2 + ry.^2);
                xtmp = xtmp + 0.4*l0tmp*(rx./rads);
                ytmp = ytmp + 0.4*l0tmp*(ry./rads);
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

        plot(flagX(ff,:)./flag(ff,:), flagY(ff,:)./flag(ff,:), 'ro', 'linewidth', 5);
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
        
        annotationStr = "$$t/\tau$$ = "+round(time(ff));
        annotation('textbox',[0.48, 0.5, 0, 0],...
            'interpreter', 'latex', 'String', annotationStr, 'Edgecolor','none', 'FitBoxToText','on');
        if showVoid
            %scatter(voidLocations{ff}(:,1), voidLocations{ff}(:,2),30, 'black', 's','filled');
            offset = 0;
            nff = ff- offset;
            if (nff > 0)
                plot([voidLocations{nff}(:,1)]...
                    , [voidLocations{nff}(:,2)]...
                    ,'k-', 'linewidth', 3);
                voidArea(ff) = polyarea(voidLocations{nff}(:,1), voidLocations{nff}(:,2));
            end
        end

        % if making a movie, save frame
        if makeAMovie == 1
            currframe = getframe(gcf);
            writeVideo(vobj,currframe);
        end
    end


    % close video object
    if makeAMovie == 1
        close(vobj);
    end
    if showVoid
        figure(14); clf, hold on, box on
        time= time - time(1);
        plot(time,voidArea, 'DisplayName', 'A', 'linewidth', 3);
        plot(time, gradient(voidArea(:))./gradient(time(:)), 'DisplayName', 'dA/dt','linewidth', 3);
        xlabel('$\tau$','Interpreter','latex');
        ylabel('A, dA/dt');
        legend
    end
end
