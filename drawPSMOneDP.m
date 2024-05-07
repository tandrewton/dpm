close all; clear; clc
filename = "_N40_dur300_att0.01_att20.01_start1_end10_sd5.xStream";
figure(); hold on;
cellpos = load(filename);
N = cellpos(1,1);
NCELLS = N-1;
NVTOT = cellpos(1,2);
% drop header row
cellpos = cellpos(2:end,:);
NV = 30;
plotVertices = 1;
plotCirculoLines = 1;
%clrs = jet(5);
%clr = clrs(1,:);
clr = [0 0.9 0.9];

% frameNum should be a scalar
%frameNum = 1;
%cellNums = [1,30,40,13,35,38];
frameNum = 1;
%cellNums = [35];
cellNums = [1:40];

for cellNum=cellNums
    numVerts = NV;
    cellRowInd = 1+cellNum*NV + frameNum * NVTOT;
    cellRows = cellRowInd:cellRowInd+NV-1;

    % if cellNum is the boundary, need to increase numVerts
    isBoundary = mod(cellNum, NCELLS) == 0;
    if isBoundary
        numVerts = NVTOT - NCELLS*NV;
        cellRowInd = 1+cellNum*NV + frameNum * NVTOT;
        cellRows = cellRowInd:cellRowInd+numVerts-1;
    end

    
    xtmp = cellpos(cellRows,1);
    ytmp = cellpos(cellRows, 2);
    vradtmp = cellpos(cellRows, 3);
    
    for vv=1:numVerts
        xplot = xtmp(vv) - vradtmp(vv);
        yplot = ytmp(vv) - vradtmp(vv);
        vnext = mod(vv, numVerts)+1;
        xtmpnext = xtmp(vnext);
        ytmpnext = ytmp(vnext);
        rx = xtmpnext - xtmp(vv);
        ry = ytmpnext - ytmp(vv);

        if plotVertices
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

            if plotCirculoLines
                % plots the circulo line part of all cells
                patch(cornerx, cornery, cornerx./cornerx, [0 0 0],'linestyle', 'none', 'LineWidth',1)
                
                % plots the vertices part of all cells
                circle2(xplot, yplot, vradtmp(vv), [0 0 0], 'none');
            else
                patch(cornerx, cornery, cornerx./cornerx, clr, 'LineWidth',1)
                circle2(xplot, yplot, vradtmp(vv), clr, [0 0 0]);
            end
            
            
            %text(xtmp(vv)-0.02, ytmp(vv)+0.005, num2str(vv), 'FontSize', 5,'color','red');
            %text(xtmp(vv)-0.02, ytmp(vv)+0.005, num2str(cellNum), 'FontSize', 5,'color','red');
        end
        if plotCirculoLines
            if isBoundary
                % plots the vertex part of the boundary cells
                circle2(xplot, yplot, vradtmp(vv), [0 0 0], 'none');
            else
                vpos = [xtmp, ytmp];
                finfo = [1:numVerts 1];

                % plots the interior of the non-boundary cells
                patch('Faces',finfo,'vertices',vpos,'FaceColor','cyan','linestyle','none');
                
                
                %plots the circulo-line part of the outer boundary of
                %non-boundary cells
                [cornerx, cornery] = patchConnectedRectanglesCorners(xtmp, ytmp, vradtmp);
                patch(cornerx', cornery', clr,'LineStyle', 'none')
            end
        end
    end
end
axis off
%print("cells",'-depsc2');

function h = circle2(x,y,r, color, EdgeColor)
    % plot a circle at x,y with radius r
    d = r*2;
    px = x+r;
    py = y+r;
    t = linspace(0, 2*pi);
    %h = rectangle('Position',[px py d d],'Curvature',[1,1], 'FaceColor',...
    %    color, 'LineWidth', 1);
    patch(r*cos(t)+px, r*sin(t)+py, ones(size(t))*2, color, 'EdgeColor', EdgeColor, 'LineWidth', 1)
    daspect([1,1,1])
end

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
