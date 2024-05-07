filename = "_N40_dur300_att0.01_att20.01_start1_end10_sd5.xStream";
cellpos = load(filename);
% drop header row
N = cellpos(1,1);
NVTOT = cellpos(1,2);
NV = 30;

frameNum = 59;
cellNum = 19;
cellRowInd = 2+cellNum*NV + frameNum * NVTOT;
cellRows = cellRowInd:cellRowInd+NV;

xtmp = cellpos(cellRows,1);
ytmp = cellpos(cellRows, 2);
vradtmp = cellpos(cellRows, 3);

cx = mean(xtmp);
cy = mean(ytmp);

figure(); hold on;

for vv=1:NV
    xplot = xtmp(vv) - vradtmp(vv);
    yplot = ytmp(vv) - vradtmp(vv);
    if showcirculoline == 1% calculate coordinates of a rectangle representing the line segment between successive vertices in a DP
        vnext = mod(vv, NV)+1;
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
        patch(cornerx, cornery, cornerx./cornerx, clr,'EdgeColor','black', 'LineWidth',1)
    end
    circle2(xplot+xx*Lx, yplot + yy*Ly, vradtmp(vv), [0 0.9 0.9]);
end

function h = circle2(x,y,r, color)
    % plot a circle at x,y with radius r
    d = r*2;
    px = x+r;
    py = y+r;
    t = linspace(0, 2*pi);
    %h = rectangle('Position',[px py d d],'Curvature',[1,1], 'FaceColor',...
    %    color, 'LineWidth', 1);
    patch(r*cos(t)+px, r*sin(t)+py, ones(size(t))*2, color, 'EdgeColor', 'black', 'LineWidth', 1)
    daspect([1,1,1])
end