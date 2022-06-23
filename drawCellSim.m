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
walls = 1;
 
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
    
%     figure(13); clf; hold on;
%     stress = load(stressstr);
%     plot(stress(:,1), stress(:,3), 'r-', 'linewidth',2, 'DisplayName',...
%         '$S_{xx}$');
%     plot(stress(:,1), stress(:,4), 'b-', 'linewidth',2, 'DisplayName',...
%         '$S_{yy}$');
%     plot(stress(:,1), stress(:,5),'k-','linewidth',2, 'DisplayName',...
%         '$S_{xy}$');
%     xlabel('$\tau$','Interpreter','latex');
%     ylabel('Stress','Interpreter','latex');
%     legend('Location', 'southeast', 'Interpreter', 'latex');
%     ax = gca;
%     ax.FontSize = 24;
%     if seed == max_seed 
%      %annotation('textbox',[.2 .5 .5 .3],'String', txt, 'Edgecolor','none')
%      saveas(gcf, output_dir + 'Stress'+runType+fileheader+'_'+max_seed+ ...
%          '.eps', 'epsc')
%     end

%     figure(14); clf;hold on;
%     plot(stress(:,1), stress(:,6), 'r--', 'linewidth',2, 'DisplayName',...
%         '$Sh_{xx}$');
%     plot(stress(:,1), stress(:,7), 'b--', 'linewidth',2, 'DisplayName',...
%         '$Sh_{yy}$');
%     plot(stress(:,1), stress(:,8),'k--','linewidth',2, 'DisplayName',...
%         '$Sh_{xy}$');
%     xlabel('$\tau$','Interpreter','latex');
%     ylabel('Stress','Interpreter','latex');
%     legend('Location', 'southeast', 'Interpreter', 'latex');
%     ax = gca;
%     ax.FontSize = 24;
%     if seed == max_seed 
%      %annotation('textbox',[.2 .5 .5 .3],'String', txt, 'Edgecolor','none')
%      saveas(gcf, output_dir + 'ShapeStress'+runType+fileheader+'_'+max_seed+ ...
%          '.eps', 'epsc')
%     end

%     figure(11); clf; hold on 
%     energy = load(energystr);
%     U = energy(:,3);
%     K = energy(:,4);
%     plot(energy(:,1), K, 'r-', 'linewidth',2, 'DisplayName',...
%         '$K$');
%      plot(energy(:,1), U,'--','linewidth',2, 'DisplayName',...
%         '$U$');
%     xlabel('$\tau$','Interpreter','latex');
%     ylabel('Energy','Interpreter','latex');
%     legend('Location', 'southeast', 'Interpreter', 'latex');
%     ax = gca;
%     ax.FontSize = 24;
%     if seed == max_seed 
%      saveas(gcf, output_dir + 'Energy'+runType+fileheader+'_'+max_seed+ ...
%          '.eps', 'epsc')
%     end

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
                    for xx = itLow:itHigh
                        for yy = itLow:itHigh
                            rectangle('Position',[xplot+xx*Lx, yplot + yy*Ly, 2*vradtmp(vv), 2*vradtmp(vv)],'Curvature',[1 1],'EdgeColor','k','FaceColor',clr);
                        end
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
            ax.XLim = [L_left Lx];
            ax.YLim = [L_bottom Ly];
            % plot box
            plot([L_left Lx Lx L_left L_left], [L_bottom L_bottom Ly Ly L_bottom], 'k-', 'linewidth', 1.5);
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
        
        %annotationStr = "$$t/\tau$$ = "+time(ff);
        %annotation('textbox',[0.48, 0.5, 0, 0],...
        %    'interpreter', 'latex', 'String', annotationStr, 'Edgecolor','none', 'FitBoxToText','on');


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

end
