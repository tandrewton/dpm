%pwd should give ~/Documents/YalePhd/projects/dpm
function drawCellSim(N, calA0, phi, ka, kb, att, att2, v0, gamma)
%close all; clear
isTestData = false; %uncomment if using function call to pipeline data

%isTestData = true; %uncomment if using test data
%testDataIDs = ["a_0.0006_a2_0.0012_p_0.75_t_1"
%"a_0.0006_a2_0.012_p_0.75_t_1"
%"a_0.006_a2_0.0012_p_0.75_t_1"
%"a_0.006_a2_0.012_p_0.75_t_1"
%"a_0.06_a2_0.0012_p_0.75_t_1"
%"a_0.06_a2_0.012_p_0.75_t_1"];

%for i=1:length(testDataIDs)
%    testDataID = testDataIDs(i);

%testDataID = "a_0.05_a2_0.05_p_0.8_t_1.0_gamma_0.01";
%testDataID = "9";

addpath('/Users/AndrewTon/Documents/YalePhD/projects/dpm/bash')
addpath('C:\Users\atata\projects\dpm\bash')
addpath('/Users/AndrewTon/Documents/YalePhD/projects/dpm/matlab_funcs')
addpath('C:\Users\atata\projects\dpm\matlab_funcs')
set(0,'DefaultFigureWindowStyle','docked')

%psm/psm_calA01.05_t_maxwell25.0_v00.05_t_abp1.0_sm1
% /_NCELLS10_dur100_att0.1_startsd1_endsd1_sd1.tissue

runType = "psm";
%N="100";
%calA0="1.05";
t_maxwell = "0";
%v0 = "0.05";
t_abp = "1.0";
kl = "1.0";
%att="0.1";
Duration="200";
FSKIP = 1;
startSeed = 1;
max_seed = 1;
%att_range = 0.3;

%if makeAMovie is 0, then plot every frame separately
%forImageAnalysis = ~isTestData;
%forImageAnalysis = true;
forImageAnalysis=false;
skipPlottingCells = false;
if (forImageAnalysis)
    showCatchBonds = 0;
    showverts = 1;
    showcirculoline = 1;
    makeAMovie = 1;
else
    showCatchBonds = 0;
    showverts = 1;
    showcirculoline = 1;
    makeAMovie = 1;
end

lengthscale = sqrt(25 * pi); % micron
timescale = 3; % minute
forceScale = 1; % nanonewton
speedConversionFactor = lengthscale / timescale; % micron/min

%PC directory
%pipeline is the location of data generated during simulations
subdir_pipeline = "pipeline/cells/"+runType+"/";

%output is location of results of this postprocessing
subdir_output = "output/cells/"+runType+"/";
mkdir(subdir_pipeline);
mkdir(subdir_output);
txt='test';

theta = linspace(0, 2*pi, 100); % Adjust the number of points as needed
fnum = 1;
fnum_boundary = 1000;
for seed = startSeed:max_seed
    if (isTestData)
        run_name = runType+txt;     
        pipeline_dir =  subdir_pipeline + run_name + "/";
        output_dir = subdir_output + run_name + "/";
        fileheader=run_name+testDataID+"_seed" + seed;
        fileheader_short = fileheader;
        nvestr = "test"+testDataID+'.pos';
        energystr = "test"+testDataID+'.energy';
        stressstr = "test"+testDataID+'.stress';
        tissuestr = "test"+testDataID+'.tissue';
        catchBondStr = "test"+testDataID+'.catchBond';
    else
        %psm/psm_calA01.05_tm0.0_v00.1_t_abp50.0
        %k_off1000.0/_N40_dur1000_att0_start1_end1_sd1.tissue
        run_name =runType+"_calA0"+calA0+'_phi'+phi+'_tm'+t_maxwell...
            +'_v0'+v0+'_t_abp'+t_abp+'_gamma'+gamma+'_kl'+kl+'_ka'+ka+'_kb'+kb;
        pipeline_dir =  subdir_pipeline + run_name + "/";
        output_dir = subdir_output + run_name + "/";
        fileheader="_N"+N+"_dur"+Duration+"_att"+att+"_att2"+att2+"_start"+ ...
            startSeed+"_end"+max_seed+"_sd"+seed;
        fileheader_short = "_N"+N+"_dur"+Duration+"_att"+att+"_att2"+att2+"_sd"+seed;
        nvestr = pipeline_dir+fileheader+'.pos';
        energystr = pipeline_dir+fileheader+'.energy';
        stressstr = pipeline_dir+fileheader+'.stress';
        tissuestr = pipeline_dir+fileheader+'.tissue';
        catchBondStr = pipeline_dir+fileheader+'.catchBond';
    end
    mkdir(pipeline_dir)
    mkdir(output_dir)
    % read in position data
    disp(nvestr)

    [trajectoryData, cell_count] = readCellClassPosOutput(nvestr);

    % set up, clear files for writing
    shapeFile = output_dir + fileheader + 'shape.csv';
    fopen(shapeFile,'w');
    speedFile = output_dir + fileheader + 'speed.csv';
    fopen(speedFile,'w');
    shapes = [];
    speeds = [];

    if showCatchBonds
        catchBondLocations = readDataBlocks(catchBondStr, 4);
    end
    cd(output_dir)

    % get number of frames
    NFRAMES = trajectoryData.NFRAMES;
    NCELLS = trajectoryData.NCELLS;
    nv = trajectoryData.nv;
    time = trajectoryData.time;

    % get frames to plot
    FSTART = 1;
    FSTEP = FSKIP;
    FEND = NFRAMES;

    if makeAMovie && ~skipPlottingCells
        movieName = runType+fileheader_short+'.mp4';
        exist(movieName, 'file')
        runType+fileheader+movieName
        if exist(movieName, 'file')==2
          delete(movieName);
        end
        moviestr = movieName;
        vobj = VideoWriter(moviestr, 'MPEG-4');
        vobj.Quality = 100;
        vobj.FrameRate = 5;
        open(vobj);
    end
    itLow = 0;
    itHigh = 0;
    boxAxLow = 0;
    boxAxHigh = 1;


    for ff = FSTART:FSTEP:FEND
        %nv can change, so recompute color map each frame
        [nvUQ, ~, IC] = unique(nonzeros(nv(ff,:)));
        %IC = IC * 0 + 1; % <- use for single colored configurations
        NUQ = length(nvUQ);
        %NUQ = 8; % 1<- use for single colored configurations
        cellCLR = jet(NUQ);
        NCELLS = cell_count(ff);


        % calculate the cell speeds
        if (ff~=FEND && ff~=1) 
            xpos_current = trajectoryData.xpos(ff,:);
            xpos_next = trajectoryData.xpos(ff+1,:);
            ypos_current = trajectoryData.ypos(ff,:);
            ypos_next = trajectoryData.ypos(ff+1,:);
            speed_ff = [];
            for nn=1:NCELLS
                cx_current = mean(xpos_current{nn});
                cy_current = mean(ypos_current{nn});
                cx_next = mean(xpos_next{nn});
                cy_next = mean(ypos_current{nn});
                speed = sqrt((cx_next - cx_current)^2 + (cy_next - cy_current)^2);
                speed = speed / (time(ff+1) - time(ff));
                speed_ff = [speed_ff; speed];
            end
            speeds = [speeds; mean(speed_ff)];
            % make a histogram of speeds, see how it looks here.
            %figure(2)
            %histogram(speeds*speedConversionFactor)
            %xlabel('v (μm/min)')
            %ylabel('# cells')
        end
        
        area = trajectoryData.area(ff,:);
        perimeter = trajectoryData.perimeter(ff,:);
        shape = 1/(4*pi)*perimeter.^2./area;
        shapes = [shapes; shape];
        %writematrix(shape(1:end-1), shapeFile, 'WriteMode','append');

        if skipPlottingCells
            continue
        end

        if ~makeAMovie
            fnum = fnum+1;
        end
        figure(fnum), clf, hold on, axis off;
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

        if (showCatchBonds)
            % calculate line between catch bond anchor points
            catchBond = catchBondLocations{ff};
            rx = catchBond(:,1) - catchBond(:,3);
            ry = catchBond(:,2) - catchBond(:,4);
            [offsetx, offsety] = patchRectangleOffsets(rx, ry, mean(vrad{1}/2));
            cornerx = [catchBond(:,1) - offsetx, catchBond(:,1) + offsetx, catchBond(:,3) + offsetx, catchBond(:,3) - offsetx];
            cornery = [catchBond(:,2) - offsety, catchBond(:,2) + offsety, catchBond(:,4) + offsety, catchBond(:,4) - offsety];
            patch(cornerx', cornery', 'red', 'linestyle', 'none')
        end

        for nn = 1:NCELLS
            cellIDtmp = cellID(nn);
            xtmp = xpos{nn};
            ytmp = ypos{nn};
            gitmp = gi{nn};
            l0tmp = l0{nn};
            vradtmp = vrad{nn};
            psitmp = psi(nn);
            costmp = cos(psitmp);
            sintmp = sin(psitmp);
            areatmp = area(nn);

            %clr = cellCLR(IC(nn),:);
            colors = ['r','g','b','c','m','y','k'];
            %clr = colors(cellIDtmp+1);
            clr='b';

            cx = mean(xtmp);
            cy = mean(ytmp);
            %text(cx, cy, num2str(nn));
            if (~showverts || (showverts && showcirculoline))
                rx = xtmp - cx;
                ry = ytmp - cy;
                rads = sqrt(rx.^2 + ry.^2);
                vpos = [xtmp, ytmp];
                finfo = [1:nv(ff,nn) 1];
                % if cellID is boundary, have it be black exterior
                % with white interior
                if (cellID(nn) == 1)
                    % if cellID is a boundary, have it be blue
                    % exterior with white interior
                    patch('Faces',finfo,'vertices',vpos,'FaceColor','w','EdgeColor','b','linewidth',0.001);
                    if (forImageAnalysis)
                        % switch to bd figure, plot bd, switch back
                        figure(fnum_boundary); clf; axis off;
                        patch('Faces',finfo,'vertices',vpos,'FaceColor','k','EdgeColor','k','linewidth',0.001)
                        figure(fnum);
                    end
                else
                    % if cellID is a real cell, give it an interior color
                    patch('Faces',finfo,'vertices',vpos,'FaceColor','cyan','linestyle','none');
                end
            end
            if showverts == 1
                if (cellID(nn) == 0)
                    boundaryX = xtmp + vradtmp * cos(theta);
                    boundaryY = ytmp + vradtmp * sin(theta);
                    patch(boundaryX', boundaryY', clr, 'linestyle', 'none')                 
                end
                % calculate coordinates of a rectangle representing 
                % the line segment between successive vertices in a DP
                if (cellID(nn) == 0 && showcirculoline == 1)
                    [cornerx, cornery] = patchConnectedRectanglesCorners(xtmp, ytmp, vradtmp);
                    patch(cornerx', cornery', clr,'LineStyle', 'none')
                end
            end
        end
        % last patch is the boundary, flip order so it's plotted underneath
        set(gca,'children',flipud(get(gca,'children')))
        axis equal;
        ax = gca;
        viewScale = 1.3;
        viewLx = viewScale*Lx;
        viewLxLow = -(viewScale-1)*Lx;
        viewLy = viewScale*Ly;
        viewLyLow = -(viewScale-1)*Ly;
        ax.XLim = [-(viewScale-1) viewScale]*Lx;
        ax.YLim = [-(viewScale-1) viewScale]*Ly;
        % plot box
        %plot([viewLxLow viewLx viewLx viewLxLow viewLxLow], [viewLyLow viewLyLow viewLy viewLy viewLyLow], 'k-', 'linewidth', 1.5);

        % if making a movie, save frame
        if makeAMovie && ~skipPlottingCells
            currframe = getframe(gcf);
            writeVideo(vobj,currframe);
        end

        if (forImageAnalysis)
            % fix fnum_boundary to have same axes as fnum
            figure(fnum_boundary)
            ax_boundary = gca;
            ax_boundary.XLim = ax.XLim;
            ax_boundary.YLim = ax.YLim;
            ax_boundary.DataAspectRatio = ax.DataAspectRatio;
            ax_boundary.Visible = ax.Visible;

            figure(fnum)
            exportgraphics(gcf, runType+fileheader_short+'fr'+ff+'.tif', 'Resolution', 200);
            figure(fnum_boundary);
            exportgraphics(gcf, runType+fileheader_short+'fr'+ff+'_bd.tif', 'Resolution', 200);
        end
    end

    % close video object
    if makeAMovie && ~skipPlottingCells
        close(vobj);
    end
    cd ../../../../
    writematrix(shapes(:,1:end-1), shapeFile, 'WriteMode','append');
    writematrix(speeds, speedFile, 'WriteMode', 'append');
    fclose('all');
end
end
%end

function [cornerx, cornery] = patchConnectedRectanglesCorners(midptx, midpty, width)
%INPUT: midptx, midpty are N x 1 vectors representing N coordinates
%           where midpt, circshift(midpt, -1) are the midpoints of two opposite
%           sides of a rectangle
%       width is N x 1 representing width of each rectangle

%OUTPUT: [cornerx, cornery] are the locations of the corners of the
%       rectangle for the patch object. The rectangles will all be
%       connected, which is useful for drawing cells.
    
    shift_midptx = circshift(midptx, -1);
    shift_midpty = circshift(midpty, -1);
    
    rx = midptx - shift_midptx;
    ry = midpty - shift_midpty;
    [offsetx, offsety] = patchRectangleOffsets(rx, ry, width);

    cornerx = [midptx-offsetx, midptx+offsetx, shift_midptx+offsetx, shift_midptx-offsetx];
    cornery = [midpty-offsety, midpty+offsety, shift_midpty+offsety, shift_midpty-offsety];
end

function [offsetx, offsety] = patchRectangleOffsets(rx, ry, width)
% INPUT: Nx1 vectors rx, ry which point from one midpoint to an opposite
%           midpoint on the same rectangle, for N rectangles
%        Nx1 vector width of the rectangle width
% OUTPUT: offsetx, offsety which represent displacements to the original
%           midpoints of the rectangle, used for a patch object

    % dot product of r and perp = 0, so perp is perpendicular to r
    perp_x = -ry ./ rx;
    perp_y = ones(size(rx));

    % if line is vertical, perpendicular is <1,0>
    isVertical = (rx == 0);
    perp_x(isVertical) = 1;
    perp_y(isVertical) = 0;

    norm = sqrt(perp_x.^2 + perp_y.^2);
    perp_x = perp_x ./ norm;
    perp_y = perp_y ./ norm;
    
    % calculate 4 coordinates of a rectangle
    % for the segment
    offsetx = width .* perp_x;
    offsety = width .* perp_y;
end