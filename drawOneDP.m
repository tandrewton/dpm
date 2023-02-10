%%Draw sims for laserAblation.cpp
% output is a movie made from stitching the position file frames together
% different from drawLoadingSims.m because it plots psi information
%pwd should give ~/Documents/YalePhd/projects/dpm
clear all
close all
clc
% function drawWoundSims(N, NV, ndelete, calA, kl, att, v0, B,...
%     Dr0,duration, boolCIL, showPeriodicImages, showverts, isTestData)
isTestData = true;
addpath('/Users/AndrewTon/Documents/YalePhD/projects/dpm/bash')
addpath('/Users/AndrewTon/Documents/YalePhD/projects/dpm/matlab_funcs')

%CHANGE THESE PARAMETERS                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   NEEDED

%PC directory
pc_dir = "/Users/AndrewTon/Documents/YalePhD/projects/dpm/";
showverts = 1;


%txt = 'N = '+N+', NV = '+NV+', calA_o='+calA+', att='+att+', B='+B;
txt='test';

figure(1), clf, hold on, box on;  

nvestr = pc_dir + "oneDP.pos";

%%

positions = load(nvestr);
nv = max(positions(:,1))+1;
xpos = positions(:,2);
ypos = positions(:,3);
vrad = positions(:,4);
l0 = positions(:,5);

[nvUQ, ~, IC] = unique(nv);
NUQ = length(nvUQ);
cellCLR = jet(NUQ);

nn=1; % only 1 cell
xtmp = xpos;
ytmp = ypos;
l0tmp = l0;
vradtmp = vrad;
clr = cellCLR(IC(nn),:);

cx = mean(xtmp);
cy = mean(ytmp);
for vv = 1:nv(nn)
    xplot = xtmp(vv) - vradtmp(vv);
    yplot = ytmp(vv) - vradtmp(vv);
    rectangle('Position',[xplot, yplot, 2*vradtmp(vv), 2*vradtmp(vv)],'Curvature',[1 1],'EdgeColor','k','FaceColor',clr);
end

axis equal
