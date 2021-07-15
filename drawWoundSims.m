%%Draw sims for laserAblation.cpp
% output is a movie made from stitching the position file frames together
% different from drawLoadingSims.m because it plots psi information
%pwd should give ~/Documents/YalePhd/projects/dpm
%works on cluster to write an avi file, but avi is a terrible format...
isTestData = true;
addpath("/Users/AndrewTon/Documents/YalePhD/projects/Jamming/CellSim/cells/bash/seq/")
%addpath("/home/at965/cells/bash/seq")

%CHANGE THESE PARAMETERS AS NEEDED

runType = "ablate";
N="96";
NV="24";
calA="1.08";
kl="1.0";
kb="0";
att="0.5";
B="0.5";
Dr0="1.0";

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


txt = 'N = '+N+', NV = '+NV+', calA_o='+calA+', att='+att+', B='+B;

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
            "_att"+att+"_B"+B+"_Dr0"+Dr0;     
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
    %NVTOT = sum(nv);
    L = trajectoryData.L(1,:);
    Lx = L(1);
    Ly = L(2);


    % show vertices or not
    showverts = 1;

    % get cell colors
    %[nvUQ, ~, IC] = unique(nv);
   % NUQ = length(nvUQ);
    %cellCLR = jet(NUQ);

    % get frames to plot
    if showverts == 0
        FSTART = 1;
        FSTEP = 1;
        FEND = NFRAMES;
    else
        FSTART = 1;
        FSTEP = 1;
        FEND = NFRAMES;
    end

    if makeAMovie == 1
       
        moviestr = output_dir + runType+fileheader+'seed_'+seed+'.mp4';
        vobj = VideoWriter(moviestr,'MPEG-4');

        vobj.FrameRate = 5;
        open(vobj);
    %end

        fnum = 1;
        figure(fnum), clf, hold on, box on;
        for ff = FSTART:FSTEP:FEND
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
            psi = trajectoryData.psi(ff,:);
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
                %plot arrows representing directors
                for xx = -1:1
                    for yy = -1:1
                        quiver(cx + xx*Lx, cy + yy*Ly, costmp, sintmp,...
                            'r', 'LineWidth', 2, 'MaxHeadSize', 2);
                    end
                end
                if showverts == 1
                    for vv = 1:nv(ff,nn)
                        xplot = xtmp(vv) - 0.5*l0tmp;
                        yplot = ytmp(vv) - 0.5*l0tmp;
                        for xx = -1:1
                            for yy = -1:1
                                %if (ff == 5 || ff == 6) 
                                %    disp(ff+","+l0tmp)
                                %end
                                rectangle('Position',[xplot + xx*Lx, yplot + yy*Ly, l0tmp, l0tmp],'Curvature',[1 1],'EdgeColor','k','FaceColor',clr);
                            end
                        end
                    end
                else
                    %cx = mean(xtmp);
                    %cy = mean(ytmp);
                    rx = xtmp - cx;
                    ry = ytmp - cy;
                    rads = sqrt(rx.^2 + ry.^2);
                    xtmp = xtmp + 0.4*l0tmp*(rx./rads);
                    ytmp = ytmp + 0.4*l0tmp*(ry./rads);
                    for xx = -1:1
                        for yy = -1:1
                            vpos = [xtmp + xx*Lx, ytmp + yy*Ly];
                            finfo = [1:nv(ff,nn) 1];
                            disp("finfo is ", finfo)
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
            ax.XLim = [-0.25 1.25]*Lx;
            ax.YLim = [-0.25 1.25]*Ly;

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
end
