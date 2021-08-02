%%Draw sims for notchTest.cpp
% output is a movie made from stitching the position file frames together
%pwd should give ~/Documents/YalePhd/projects/dpm
%works on cluster to write an avi file, but avi is a terrible format...
isTestData = true;
addpath("/Users/AndrewTon/Documents/YalePhD/projects/Jamming/CellSim/cells/bash/seq/")
%addpath("/home/at965/cells/bash/seq")

%CHANGE THESE PARAMETERS AS NEEDED
runType = "notch";
N="96";
NV="24";
calA="1.08";
kl="1.0";
kb="0";
att="0.5";
B="0.1";
strain="0.5";
strainRate="0.01";
duration="100";
loadingType="uniaxial";
FSKIP = 1; %# frames skipped to lower sampling freq

%PC directory
pc_dir = "/Users/AndrewTon/Documents/YalePhD/projects/dpm/";
%pipeline is the location of data generated during simulations
subdir_pipeline = pc_dir + "pipeline/cells/"+runType+"/";

%output is location of results of this postprocessing
subdir_output = pc_dir + "output/cells/"+runType+"/";
mkdir(subdir_pipeline);
mkdir(subdir_output);

startSeed = 1;
max_seed = 1;
makeAMovie = 1;
showPeriodicImages = 1;
plotVoronoi = 0;

% show vertices or not
showverts = 0;

%txt = 'N = '+N+', NV = '+NV+', calA_o='+calA+', att='+att+', B='+B+', load='+loadingType;
txt = 'test';
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
    else
        run_name =runType+"_N"+N+"_NV"+NV+"_calA0"+calA+"_kl"+kl+...
            "_att"+att+"_B"+B+"_strain"+strain+"_strainRate"+strainRate+...
            "_duration"+duration+"_loading_"+loadingType;     
        pipeline_dir =  subdir_pipeline + run_name + "/";
        output_dir = subdir_output + run_name + "/";
        mkdir(pipeline_dir)
        mkdir(output_dir)
        fileheader=run_name + "_seed" + seed;
        nvestr = pipeline_dir+fileheader+'.pos';
        energystr = pipeline_dir+fileheader+'.energy';
        stressstr = pipeline_dir+fileheader+'.stress';
    end
    

    figure(13);
    stress = load(stressstr);
    plot(stress(:,2), stress(:,3), 'r-', 'linewidth',2, 'DisplayName',...
         '$S_{xx}$');
    plot(stress(:,2), stress(:,4), 'b-', 'linewidth',2, 'DisplayName',...
         '$S_{yy}$');
    plot(stress(:,2), stress(:,5),'k-','linewidth',2, 'DisplayName',...
         '$S_{xy}$');
    xlabel('$\epsilon_x$ (tensile strain)','Interpreter','latex');
    ylabel('Stress','Interpreter','latex');
    %plot(stress(:,1), stress(:,2), 'k-', 'linewidth', 2);
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
    FEND = NFRAMES/2;

    if makeAMovie == 1
       
        moviestr = output_dir + runType+fileheader+'seed_'+seed+'.mp4';
        voronoi_movie_str = output_dir + "voronoi"+runType+fileheader+'seed_'+seed+'.mp4';
        vobj = VideoWriter(moviestr,'MPEG-4');
        vor_vobj = VideoWriter(voronoi_movie_str, 'MPEG-4');

        vobj.Quality = 100;
        vobj.FrameRate = 15;
        open(vobj);
        vor_vobj.Quality = 100;
        vor_vobj.FrameRate = 15;
        open(vor_vobj);
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
    for ff = FSTART:FSTEP:FEND
        vor_cx = [];
        vor_cy = [];
        %nv can change, so recompute color map each frame
        [nvUQ, ~, IC] = unique(nonzeros(nv(ff,:)));
        NUQ = length(nvUQ);
        cellCLR = jet(NUQ);

        NCELLS = cell_count(ff);
        figure(fnum), clf, hold on, box on;
        fprintf('printing frame ff = %d/%d\n',ff,FEND);

        % get cell positions
        xpos = trajectoryData.xpos(ff,:);
        ypos = trajectoryData.ypos(ff,:);
        l0 = trajectoryData.l0(ff,:);

        %if L is not constant, use the next 3 lines
        L = trajectoryData.L(ff,:);
        Lx = L(1);
        Ly = L(2);

        %idea: to keep the box size constant visually,
        % just store Lx0, Ly0. Every frame, compute delta = Lx(new) - Lx(old),
        % Ly(new) - Ly(old). Scale particle sizes by (1-delta). 

        for nn = 1:NCELLS
            xtmp = xpos{nn};
            ytmp = ypos{nn};
            l0tmp = l0(nn);
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
                        patch('Faces',finfo,'vertices',vpos,'FaceColor',clr,'EdgeColor','k','linewidth',2);
                    end
                end
            end
        end

        % plot box
        plot([0 Lx Lx 0 0], [0 0 Ly Ly 0], 'k-', 'linewidth', 1.5);
        axis equal;
        ax = gca;
        ax.XTick = [];
        ax.YTick = [];
        if showPeriodicImages == 1
            ax.XLim = [boxAxLow boxAxHigh]*Lx;
            ax.YLim = [boxAxLow boxAxHigh]*Ly;
        else
            ax.XLim = [0 1]*Lx;
            ax.YLim = [0 1]*Ly;
        end

        % if making a movie, save frame
        if makeAMovie == 1
            currframe = getframe(gcf);
            writeVideo(vobj,currframe);
        end
        if plotVoronoi == 1 
            mat_cx = mean(cell2mat(xpos));
            mat_cy = mean(cell2mat(ypos));
            for xx=itLow:itHigh
                for yy=itLow:itHigh
                    vor_cx = [vor_cx, mat_cx + xx*Lx];
                    vor_cy = [vor_cy, mat_cy + yy*Ly];
                end
            end
            figure(fnum+20), clf, hold on, box on;
            plot([0 Lx Lx 0 0], [0 0 Ly Ly 0], 'k-', 'linewidth', 1.5);
            axis equal;
            ax = gca;
            ax.XTick = [];
            ax.YTick = [];
            if showPeriodicImages == 1
                ax.XLim = [boxAxLow boxAxHigh]*Lx;
                ax.YLim = [boxAxLow boxAxHigh]*Ly;
            else
                ax.XLim = [0 1]*Lx;
                ax.YLim = [0 1]*Ly;
            end
            voronoi(ax, vor_cx,vor_cy);    
            if makeAMovie == 1
                currframe = getframe(gcf);
                writeVideo(vor_vobj, currframe);
            end
        end
    end


    % close video object
    if makeAMovie == 1
        close(vobj);
        if plotVoronoi == 1
            close(vor_vobj);
        end
    end
end