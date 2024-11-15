%pwd should give ~/Documents/YalePhd/projects/dpm
function drawCellSim(N, calA0, phi, kl, kb, att, att2, v0, t_maxwell, gamma, k_on, k_off, k_ecm, numSeeds)
%close all; clear
isTestData = false; %uncomment if using function call to pipeline data

%isTestData = true; %uncomment if using test data
%testDataIDs = ["testa_0.05_a2_0.05_tm_10000.0_p_0.8_t_1.0_gamma_0_k_on_1.0_k_off_0.5_k_ecm_0.05"];

%for i=1:length(testDataIDs)
%    testDataID = testDataIDs(i);

%testDataID = "a_0.05_a2_0.05_tm_10000.0_p_0.8_t_1.0_gamma_0_k_on_1.0_k_off_0.5_k_ecm_0.05";
%testDataID = "9";

addpath('/Users/AndrewTon/Documents/YalePhD/projects/dpm/bash')
addpath('C:\Users\atata\projects\dpm\bash')
addpath('/home/at965/dpm/bash')
addpath('/Users/AndrewTon/Documents/YalePhD/projects/dpm/matlab_funcs')
addpath('C:\Users\atata\projects\dpm\matlab_funcs')
addpath('/home/at965/dpm/matlab_funcs')
if ~isunix
    set(0,'DefaultFigureWindowStyle','docked')
end

%psm/psm_calA01.05_t_maxwell25.0_v00.05_t_abp1.0_sm1
% /_NCELLS10_dur100_att0.1_startsd1_endsd1_sd1.tissue

runType = "psm";
%N="100";
%calA0="1.05";
%t_maxwell = "0";
%v0 = "0.05";
t_abp = "1.0";
ka = "2.5";
%kl = "1.0";
%att="0.1";
Duration="300";
FSKIP = 1;
startSeed = 1;
max_seed = numSeeds; % gets overwritten if OS is unix (I use slurm in unix)
%att_range = 0.3;

%if makeAMovie is 0, then plot every frame separately
%forImageAnalysis = ~isTestData;
forImageAnalysis = true;
%forImageAnalysis=false;
skipPlottingCells = false;
if (forImageAnalysis)
    showBonds = 1;
    showverts = 1;
    showcirculoline = 1;
    makeAMovie = 1;
else
    showBonds = 1;
    showverts = 1;
    showcirculoline = 0;
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

if isunix
    % for cluster runs
    subdir_pipeline = "/gpfs/gibbs/pi/ohern/at965/dpm/"+runType+"/";
    subdir_output = "/gpfs/gibbs/pi/ohern/at965/dpm/"+runType+"/output/";
    max_seed = numSeeds;
end

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
        bondStr = "test"+testDataID+'.bond';
    else
        %psm/psm_calA01.05_tm0.0_v00.1_t_abp50.0
        %k_off1000.0/_N40_dur1000_att0_start1_end1_sd1.tissue
        run_name =runType+"_calA0"+calA0+'_phi'+phi+'_tm'+t_maxwell...
            +'_v0'+v0+'_t_abp'+t_abp+'_gamma'+gamma+'_k_on_'+k_on+'_k_off_'+k_off+'_k_ecm_'+k_ecm+'_kl'+kl+'_ka'+ka+'_kb'+kb;
        pipeline_dir =  subdir_pipeline + run_name + "/";
        output_dir = subdir_output + run_name + "/";
        fileheader="_N"+N+"_dur"+Duration+"_att"+att+"_att2"+att2+"_start"+ ...
            startSeed+"_end"+max_seed+"_sd"+seed;
        fileheader_short = "_N"+N+"_dur"+Duration+"_att"+att+"_att2"+att2+"_sd"+seed;
        nvestr = pipeline_dir+fileheader+'.pos';
        energystr = pipeline_dir+fileheader+'.energy';
        stressstr = pipeline_dir+fileheader+'.stress';
        tissuestr = pipeline_dir+fileheader+'.tissue';
        bondStr = pipeline_dir+fileheader+'.bond';
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

    if showBonds
        bondLocations = readDataBlocks(bondStr, 4);
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
        if isunix
            movieName = runType+fileheader_short+'.avi';
        end
        runType+fileheader+movieName
        if exist(movieName, 'file')==2
          delete(movieName);
        end
        moviestr = movieName;
        if ~isunix
            vobj = VideoWriter(moviestr, 'MPEG-4');
            vobj.Quality = 100;
        else
            vobj = VideoWriter(moviestr, 'Uncompressed AVI');
        end
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
        NUQ = length(nvUQ);
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

        % uncomment to have frames print each frame to different figure numbers
        %if ~makeAMovie
        %    fnum = fnum+1;
        %end
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

            
            clr = [0.5 0 0.5];

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
                    %patch('Faces',finfo,'vertices',vpos,'FaceColor','w','EdgeColor',[0 1 0],'linewidth',1);
                    patch('Faces',finfo,'vertices',vpos,'FaceColor',[0 0.4078 0.2157],'FaceAlpha', 0.86, 'EdgeColor',[0.4784 0.7882 0.2627],'linewidth',5);
                    if (forImageAnalysis)
                        % switch to bd figure, plot bd, switch back
                        figure(fnum_boundary); clf; axis off;
                        patch('Faces',finfo,'vertices',vpos,'FaceColor','k','EdgeColor','k','linewidth',0.001)
                        figure(fnum);
                    end
                else
                    % if cellID is a real cell, give it an interior color
                    patch('Faces',finfo,'vertices',vpos,'FaceColor',[1 0 1],'linestyle','none');
                end
            end
            if showverts == 1
                if (cellID(nn) == 0 || (cellID(nn) == 1 && showcirculoline == 0))
                    boundaryX = xtmp + vradtmp * cos(theta);
                    boundaryY = ytmp + vradtmp * sin(theta);
                    patch(boundaryX', boundaryY', 'k', 'linestyle', 'none')
                elseif (cellID(nn) == 1)
                    boundaryX = xtmp + vradtmp * cos(theta);
                    boundaryY = ytmp + vradtmp * sin(theta);
                    patch(boundaryX', boundaryY', [0 1 0], 'linestyle', 'none')
                end
                % calculate coordinates of a rectangle representing 
                % the line segment between successive vertices in a DP
                if showcirculoline == 1
                    [cornerx, cornery] = patchConnectedRectanglesCorners(xtmp, ytmp, vradtmp);
                    if cellID(nn) == 0
                        patch(cornerx', cornery', 'k','LineStyle', 'none')
                    elseif cellID(nn) == 1
                        patch(cornerx', cornery', [0 1 0],'LineStyle', 'none')
                    end
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

        if (forImageAnalysis)
            % fix fnum_boundary to have same axes as fnum
            figure(fnum_boundary)
            ax_boundary = gca;
            ax_boundary.XLim = ax.XLim;
            ax_boundary.YLim = ax.YLim;
            ax_boundary.DataAspectRatio = ax.DataAspectRatio;
            ax_boundary.Visible = ax.Visible;

            figure(fnum)
            if ff == 60
                exportgraphics(gcf, runType+fileheader_short+'fr'+ff+'.tif', 'Resolution', 600);
            else
                exportgraphics(gcf, runType+fileheader_short+'fr'+ff+'.tif', 'Resolution', 200);
            end
            figure(fnum_boundary);
            exportgraphics(gcf, runType+fileheader_short+'fr'+ff+'_bd.tif', 'Resolution', 200);
            
            %switch views back to cells to plot bonds
            figure(fnum)
        end

        % want bonds in movie but not in tif, so write tif before movie
        if (showBonds)
            % calculate line between catch bond anchor points
            bond = bondLocations{ff};
            rx = bond(:,1) - bond(:,3);
            ry = bond(:,2) - bond(:,4);
            [offsetx, offsety] = patchRectangleOffsets(rx, ry, mean(vrad{1}/2));
            cornerx = [bond(:,1) - offsetx, bond(:,1) + offsetx, bond(:,3) + offsetx, bond(:,3) - offsetx];
            cornery = [bond(:,2) - offsety, bond(:,2) + offsety, bond(:,4) + offsety, bond(:,4) - offsety];
            patch(cornerx', cornery', 'red', 'linestyle', 'none')
        end

        % if making a movie, save frame
        if makeAMovie && ~skipPlottingCells
            currframe = getframe(gcf);
            writeVideo(vobj,currframe);
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