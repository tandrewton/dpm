%%Draw sims for laserAblation.cpp
% output is a movie made from stitching the position file frames together
% different from drawLoadingSims.m because it plots psi information
%pwd should give ~/Documents/YalePhd/projects/dpm

%function drawCellSim(N, att, initialPressure, prate, adhrate, Duration)
%isTestData = false; %uncomment if using function call to pipeline data

isTestData = true; %uncomment if using test data
addpath('/Users/AndrewTon/Documents/YalePhD/projects/dpm/bash')
addpath('C:\Users\atata\projects\dpm\bash')
addpath('/Users/AndrewTon/Documents/YalePhD/projects/dpm/matlab_funcs')
%CHANGE THESE PARAMETERS                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   NEEDED

%NT_calA01.0_initPressure-0.01_prate0.002_adhrate0.0
% _NCELLS10_Duration100_att0.01_startseed1_endseed1_seed1
runType = "NT";
%N="40";
calA0="1.0";
%att="0.2";
%Duration="800";
FSKIP = 1;

startSeed = 1;
max_seed = 1;
makeAMovie = 1; %if makeAMovie is 0, then plot every frame separately
set(0,'DefaultFigureWindowStyle','docked')
showPeriodicImages = 0;

showverts = 0;
showcirculoline = 0;
walls = 0;
att_range = 0.3;
 
%PC directory
%pc_dir = "/Users/AndrewTon/Documents/YalePhD/projects/dpm/";
pc_dir="C:\Users\atata\projects\dpm\";
%pipeline is the location of data generated during simulations
subdir_pipeline = pc_dir + "pipeline/cells/"+runType+"/";

%output is location of results of this postprocessing
subdir_output = pc_dir + "output/cells/"+runType+"/";
mkdir(subdir_pipeline);
mkdir(subdir_output);


%txt = 'N = '+N+', NV = '+NV+', calA_o='+calA+', att='+att+', B='+B;
txt='test';
load("polyBoundary.txt"); % load boundaries of polygon walls
load("initPosSP.txt");
load("initPosSP2.txt");

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
        nvestr = pc_dir+'test.pos';
        energystr = pc_dir+'test.energy';
        stressstr = pc_dir+'test.stress';
    else
        %NT_calA01.0_initPressure-0.01_prate0.002_adhrate0.0
        % _NCELLS10_Duration100_att0.01_startseed1_endseed1_seed1
        run_name =runType+"_calA0"+calA0+'_initPressure'+initialPressure...
            +'_prate'+prate+'_adhrate'+adhrate;
        pipeline_dir =  subdir_pipeline + run_name + "/";
        output_dir = subdir_output + run_name + "/";
        mkdir(pipeline_dir)
        mkdir(output_dir)
        fileheader=run_name +"_NCELLS"+N+"_Duration"+Duration+"_att"+att+"_startseed"+ ...
            startSeed+"_endseed"+max_seed+"_seed"+seed;
        nvestr = pipeline_dir+fileheader+'.pos';
        energystr = pipeline_dir+fileheader+'.energy';
        stressstr = pipeline_dir+fileheader+'.stress';
    end

    % read in position data
    nvestr
    [trajectoryData, cell_count] = readCellClassPosOutput(nvestr);

    % get number of frames
    NFRAMES = trajectoryData.NFRAMES;
    NCELLS = trajectoryData.NCELLS;
    nv = trajectoryData.nv;
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
        movieName = output_dir + runType+fileheader+'seed_'+seed+'.mp4';
        if exist(movieName, 'file')==2
          delete(movieName);
        end
        moviestr = movieName;
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

    figure(fnum), clf, hold on, box on;

    for ff = FSTART:FSTEP:FEND
        %nv can change, so recompute color map each frame
        [nvUQ, ~, IC] = unique(nonzeros(nv(ff,:)));
        %IC = IC * 0 + 1; % <- use for single colored configurations
        NUQ = length(nvUQ);
        %NUQ = 8; % 1<- use for single colored configurations
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
        cellID = trajectoryData.cellID(ff,:);
        xpos = trajectoryData.xpos(ff,:);
        ypos = trajectoryData.ypos(ff,:);
        gi = trajectoryData.gi(ff,:);
        l0 = trajectoryData.l0(ff,:);
        vrad = trajectoryData.vrad(ff,:);
        psi = trajectoryData.psi(ff,:);
        

        %if L is not constant, use the next 3 lines
        L = trajectoryData.L(ff,:);
        L_left = L(1);
        L_bottom = L(2);
        Lx = L(3);
        Ly = L(4);

        for nn = 1:NCELLS
            cellIDtmp = cellID(nn);
            xtmp = xpos{nn};
            ytmp = ypos{nn};
            gitmp = gi{nn};
            l0tmp = l0{nn};
            vradtmp = vrad{nn}*1.0;
            psitmp = psi(nn);
            costmp = cos(psitmp);
            sintmp = sin(psitmp);

            %clr = cellCLR(IC(nn),:);
            colors = ['r','g','b','c','m','y','k'];
            clr = colors(cellIDtmp+1);

            cx = mean(xtmp);
            cy = mean(ytmp);
            if showverts == 1
                for vv = 1:nv(ff,nn)
                    xplot = xtmp(vv) - vradtmp(vv);
                    yplot = ytmp(vv) - vradtmp(vv);
                    for xx = itLow:itHigh
                        for yy = itLow:itHigh
                            rectangle('Position',[xplot+xx*Lx, yplot + yy*Ly, 2*vradtmp(vv), 2*vradtmp(vv)],'Curvature',[1 1],'EdgeColor','k','FaceColor',clr);
                            %text(xplot-0.25,yplot-0.25,num2str(vv-1))
                        end
                    end
                end
                for vv = 1:nv(ff,nn)
                    if nn == 1
                        continue
                    end
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
                text(cx,cy,num2str(nn))
                for xx = itLow:itHigh
                    for yy = itLow:itHigh
                        vpos = [xtmp + xx*Lx, ytmp + yy*Ly];
                        finfo = [1:nv(ff,nn) 1];
                        %disp("finfo is "+ finfo)
                        patch('Faces',finfo,'vertices',vpos,'FaceColor',clr,'EdgeColor','k','linewidth',2);
                    end
                end
            end

        end

        %for kk=1:4
        %    plot([polyBoundary(kk,1:2:end) polyBoundary(kk,1)],...
        %        [polyBoundary(kk,2:2:end) polyBoundary(kk,2)], 'k','linewidth', 1)
        %end
        %scatter(initPosSP(:,1), initPosSP(:,2),'ro')
        %scatter(initPosSP2(:,1),initPosSP2(:,2),'ko')

        axis equal;
        ax = gca;
        % since boundaries are last in the file, flip children order to get
        % boundary patch objects below cell patch objects
        set(gca,'children',flipud(get(gca,'children'))) 
        %ax.XTick = [];
        %ax.YTick = [];
        if showPeriodicImages == 1
            ax.XLim = [boxAxLow boxAxHigh]*Lx;
            ax.YLim = [boxAxLow boxAxHigh]*Ly;
            % plot box
            plot([0 Lx Lx 0 0], [0 0 Ly Ly 0], 'k-', 'linewidth', 1.5);
        elseif walls == 1
            ax.XLim = [L_left Lx];
            ax.YLim = [L_bottom Ly];
            % plot box
            plot([L_left Lx Lx L_left L_left], [L_bottom L_bottom Ly Ly L_bottom], 'k-', 'linewidth', 1.5);
        else
            viewScale = 1.5;
            viewLx = viewScale*Lx;
            viewLxLow = -(viewScale-1)*Lx;
            viewLy = viewScale*Ly;
            viewLyLow = -(viewScale-1)*Ly;
            ax.XLim = [-(viewScale-1) viewScale]*Lx;
            ax.YLim = [-(viewScale-1) viewScale]*Ly;
            % plot box
            %plot([viewLxLow viewLx viewLx viewLxLow viewLxLow], [viewLyLow viewLyLow viewLy viewLy viewLyLow], 'k-', 'linewidth', 1.5);
        end
        
        %annotationStr = "$$t/\tau$$ = "+time(ff);
        framenum = ff-1;
        annotationStr = "frame="+framenum;
        annotation('textbox',[0.4, 0.4, 0, 0],...
            'interpreter', 'latex', 'String', annotationStr, 'Edgecolor','none', 'FitBoxToText','on');

        % if making a movie, save frame
        if makeAMovie == 1
            currframe = getframe(gcf);
            writeVideo(vobj,currframe);
        end
        %while (ff == 77)
        %    disp("hi")
        %end

    end


    % close video object
    if makeAMovie == 1
        close(vobj);
    end

end
