%%Draw sims for laserAblation.cpp
% output is a movie made from stitching the position file frames together
% different from drawLoadingSims.m because it plots psi information
%pwd should give ~/Documents/YalePhd/projects/dpm

%function drawCellSim(N, att, initialPressure, prate, adhrate, Duration)
%isTestData = false; %uncomment if using function call to pipeline data

isTestData = true; %uncomment if using test data
testDataID = num2str(testDataii);
addpath('/Users/AndrewTon/Documents/YalePhD/projects/dpm/bash')
addpath('C:\Users\atata\projects\dpm\bash')
addpath('/Users/AndrewTon/Documents/YalePhD/projects/dpm/matlab_funcs')
%CHANGE THESE PARAMETERS                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   NEEDED

%psm/psm_calA01.05_t_maxwell25.0_v00.05_t_abp1.0_sm1
% /_NCELLS10_dur100_att0.1_startsd1_endsd1_sd1.tissue

runType = "psm";
N="100";
calA0="1.05";
t_maxwell = "25.0";
v0 = "0.05";
t_abp = "1.0";
sm = "1";
att="0.1";
Duration="400";
FSKIP = 1;

startSeed = 1;
max_seed = 1;
makeAMovie = 1; %if makeAMovie is 0, then plot every frame separately
set(0,'DefaultFigureWindowStyle','docked')
showPeriodicImages = 0;

showverts = 0;
showGlobalIndex = 0;
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

phi_array = [];
fnum = 1;
figure(13), clf, hold on, box on;
for seed = startSeed:max_seed
    if (isTestData)
        run_name = runType+txt;     
        pipeline_dir =  subdir_pipeline + run_name + "/";
        output_dir = subdir_output + run_name + "/";
        mkdir(pipeline_dir)
        mkdir(output_dir)
        fileheader=run_name+testDataID+"_seed" + seed;
        nvestr = pc_dir+'test'+testDataID+'.pos';
        energystr = pc_dir+'test'+testDataID+'.energy';
        stressstr = pc_dir+'test'+testDataID+'.stress';
        tissuestr = pc_dir+'test'+testDataID+'.tissue';
    else
        %psm/psm_calA01.05_t_maxwell25.0_v00.05_t_abp1.0_sm1
        % /_NCELLS10_dur100_att0.1_startsd1_endsd1_sd1.tissue
        run_name =runType+"_calA0"+calA0+'_t_maxwell'+t_maxwell...
            +'_v0'+v0+'_t_abp'+t_abp+'_sm'+sm;
        pipeline_dir =  subdir_pipeline + run_name + "/";
        output_dir = subdir_output + run_name + "/";
        mkdir(pipeline_dir)
        mkdir(output_dir)
        fileheader="_NCELLS"+N+"_dur"+Duration+"_att"+att+"_startsd"+ ...
            startSeed+"_endsd"+max_seed+"_sd"+seed;
        nvestr = pipeline_dir+fileheader+'.pos';
        energystr = pipeline_dir+fileheader+'.energy';
        stressstr = pipeline_dir+fileheader+'.stress';
        tissuestr = pipeline_dir+fileheader+'.tissue';
    end

%         A = load(tissuestr);
%         %plot(A(:,1), A(:,2), 'displayname', 'boundary area')
%         %plot(A(:,1), A(:,3), 'displayname', 'total cell area')
%         plot(A(:,1), (A(:,2) - A(:,3))./A(:,2), 'displayname',...
%             "$\phi$, $<A>$ = " + mean(A(:,2)), 'linewidth', 2)
%         legend()
%         xlabel('$\tau$','Interpreter','latex');
%         ylabel('Extracellular Fraction','Interpreter','latex');
%         legend('Location', 'northeast', 'Interpreter', 'latex');
%         ylim([0 1])
%         ax = gca;
%         ax.FontSize = 20;

    % read in position data
    [trajectoryData, cell_count] = readCellClassPosOutput(nvestr);

    % get number of frames
    NFRAMES = trajectoryData.NFRAMES;
    NCELLS = trajectoryData.NCELLS;
    nv = trajectoryData.nv;
    time = trajectoryData.time;

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

    %FSTART=FEND; % hack to skip to the last frame to save the boundary info

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
        cellarea = trajectoryData.area(ff,:);
        

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
            areatmp = cellarea(nn);

            %clr = cellCLR(IC(nn),:);
            colors = ['r','g','b','c','m','y','k'];
            %clr = colors(cellIDtmp+1);
            clr='k';

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
                            if showGlobalIndex
                                text(xtmp(vv), ytmp(vv), num2str(gitmp(vv)), 'FontSize', 6);
                            end
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
                %xtmp = xtmp + 0.4*l0tmp(1)*(rx./rads);
                %ytmp = ytmp + 0.4*l0tmp(1)*(ry./rads);
                %text(cx,cy,num2str(nn)) % plot cell # on each cell
                for xx = itLow:itHigh
                    for yy = itLow:itHigh
                        vpos = [xtmp + xx*Lx, ytmp + yy*Ly];
                        finfo = [1:nv(ff,nn) 1];
                        %disp("finfo is "+ finfo)
                        % if cellID is boundary, have it be black exterior
                        % with white interior
                        if (nn == NCELLS)
                            % if cellID is a boundary, have it be black
                            % exterior with white interior
                            patch('Faces',finfo,'vertices',vpos,'FaceColor','w','EdgeColor','k','linewidth',0.001);
                        else
                            % if cellID is a real cell, have it be black
                            % exterior with black interior
                            %patch('Faces',finfo,'vertices',vpos,'FaceColor','k','EdgeColor','k','linewidth',0.001);
                            patch('Faces',finfo,'vertices',vpos,'FaceColor',clr, ...
                                'FaceAlpha', 'flat','FaceVertexAlphaData', 0, 'EdgeColor','k',...
                                'linewidth',0.001);
                        end
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
            viewScale = 1.3;
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
        %annotation('textbox',[0.4, 0.4, 0, 0],...
        %    'interpreter', 'latex', 'String', annotationStr, 'Edgecolor','none', 'FitBoxToText','on');

        % if making a movie, save frame
        if makeAMovie == 1
            axis off;
            currframe = getframe(gcf);
            writeVideo(vobj,currframe);
        end

        if (ff == FEND)
            axis off;
            %exportgraphics(gcf, 'last_frame_PSM_sim_att'+att+'_sd'+seed+'.tif', 'Resolution', 1000);
            %writematrix(vpos, "last_frame_PSM_images/" + ...
            %    "last_frame_PSM_sim_att"+att+"_sd"+seed+".txt");
            %exportgraphics(gcf, 'last_frame_PSM_sim_att'+att+'_sd'+seed+'_bd.tif', 'Resolution', 1000);
            phi_array = [phi_array sum(cellarea(1:end-1))/cellarea(end)];
        end
    end


    % close video object
    if makeAMovie == 1
        close(vobj);
    end

end
