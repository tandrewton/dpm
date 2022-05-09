isTestData = true; %uncomment if using test data
addpath('/Users/AndrewTon/Documents/YalePhD/projects/dpm/bash')
addpath('/Users/AndrewTon/Documents/YalePhD/projects/dpm/matlab_funcs')

% this code is good for running dampedNVE2D(), which doesn't have
% purse-string or crawling yet. other code will expect those functions

%CHANGE THESE PARAMETERS                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   NEEDED

runType = "ablate";
%N="40";
ndelete="10";
calA0="1.10";
%strainRate_ps="0.01";
%deltaSq = "2.0";
k_ps = "1.0"; %purse-string spring constant
k_lp = "2.0"; %lamellipodia spring constant
%tau_lp = "1.0"; %lamellipodia lifetime
%d_flag = "0.0"; %lamellipodia max length
prate = "0.00"; %perimeter relaxation rate
%att="0.2";
B="1.0";
Dr0="0.5";
boolCIL="0";
Duration="800";
FSKIP = 1;

etaStr = " ";
startSeed = 1;
max_seed = 1;
makeAMovie = 1; %if makeAMovie is 0, then plot every frame separately
set(0,'DefaultFigureWindowStyle','docked')
showPeriodicImages = 0;

showverts = 1;
showBoundaries = 0;
showArea = 1;
showQuiver = 0;
walls = 0;
%disable showVoid if using printConfig on its own, outside of
%dampedNVE/dampedNP0 routines
showGlobalIndex = 0;
showVoid = 0;
showVoidBlack = 0; % print void in larger black circles to see easier
showCornersOrEdges = 0;
showPurseString = 1;
showShapeHistogram = 0;
 
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
        boundaryStr = pc_dir+'test.void';
        edgeStr = pc_dir+'test.edge';
        purseStr = pc_dir+'test.purseString';
        voidAreaStr = pc_dir+'test.voidArea';
    else
        run_name =runType+"_calA0"+calA0+"_strainRate_ps"+strainRate_ps+ ...
            "_deltaSq"+deltaSq+"_k_ps"+k_ps+"_k_lp"+k_lp+...
            "_tau_lp"+tau_lp+"_d_flag"+d_flag+"_prate"+prate;
        pipeline_dir =  subdir_pipeline + run_name + "/";
        output_dir = subdir_output + run_name + "/";
        mkdir(pipeline_dir)
        mkdir(output_dir)
        fileheader=run_name +"_NCELLS"+N+"_Duration"+Duration+"_att"+att+"_startseed"+ ...
            startSeed+"_endseed"+max_seed+"_seed"+seed;
        nvestr = pipeline_dir+fileheader+'.pos';
        energystr = pipeline_dir+fileheader+'.energy';
        stressstr = pipeline_dir+fileheader+'.stress';
        boundaryStr = pipeline_dir+fileheader+ ".void";
        edgeStr = pipeline_dir+fileheader+ '.edge';
        purseStr = pipeline_dir+fileheader+ '.purseString';
        voidAreaStr = pipeline_dir+fileheader+ '.voidArea';
    end

    figure(11); clf; hold on 
    stress = load(stressstr);
    energy = load(energystr);
    U = energy(:,3);
    K = energy(:,4);
    U_ps = energy(:,5);
    U_crawling = energy(:,6);
    plot(energy(:,1), K, 'r-', 'linewidth',2, 'DisplayName',...
        '$K$');
    plot(energy(:,1), U_ps, 'b-', 'linewidth',2, 'DisplayName',...
        '$U_{ps}$');
    plot(energy(:,1), U_crawling,'k-','linewidth',2, 'DisplayName',...
        '$U_{crawling}$');
    plot(energy(:,1), U,'--','linewidth',2, 'DisplayName',...
        '$U$');
    xlabel('$\tau$','Interpreter','latex');
    ylabel('Energy','Interpreter','latex');
    legend('Location', 'southeast', 'Interpreter', 'latex');
    ax = gca;
    ax.FontSize = 24;

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

    for ff = FSTART:FSTEP:FEND
        %nv can change, so recompute color map each frame
        [nvUQ, ~, IC] = unique(nonzeros(nv(ff,:).*(flag(ff,:)+1)));
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
                    for xx = itLow:itHigh
                        for yy = itLow:itHigh
                            rectangle('Position',[xplot+xx*Lx, yplot + yy*Ly, 2*vradtmp(vv), 2*vradtmp(vv)],'Curvature',[1 1],'EdgeColor','k','FaceColor',clr);
                            if showGlobalIndex
                                text(xtmp(vv), ytmp(vv), num2str(gitmp(vv)));
                            end
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

        plot(flagX(ff,:)./flag(ff,:), flagY(ff,:)./flag(ff,:), 'ro', 'linewidth', 2);
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
    end
end