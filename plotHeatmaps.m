clear;close all;
set(0,'DefaultFigureWindowStyle','docked')
isPlottingAreaVelocity = false;
isPlottingShapes = ~isPlottingAreaVelocity;
assert((isPlottingAreaVelocity && isPlottingShapes) == false);
%%
folderPath = "/Users/AndrewTon/Documents/YalePhD/projects/dpm/output/cells/ablate/array_output_figures/heatmaps/";
if (isPlottingAreaVelocity)
    heatmapFile = "arealVelocityPStrue-calA01.20_tau-ka.mat";
else
    heatmapFile = "rosetteCellsFinalShapePStrue-calA01.20_tau-ka.mat";
end
%shapeFile = "rosetteCellsFinalShapePStrue-activity-ka.mat";
%heatmapFile = "rosetteCellsFinalShapePStrue-calA01.20_tau-ka.mat";

% physical units here
cellArea = 25; % microns^2
constrictionRate = 0.3; % micron/sec
adhesionForce = 1e-9; % Newton

% simulation to real unit conversion
bulkModulusConvert = adhesionForce/cellArea*1e12/1000; %kPa
timeConvert = sqrt(cellArea)/constrictionRate/60; % minutes
areaVelocityConvert = constrictionRate*sqrt(cellArea)*60; %micron^2/min

hfile = load(folderPath + heatmapFile);
hfile = hfile.h;
hfile.ColorData(isnan(hfile.ColorData)) = 0;
xpoints = cellfun(@str2num,convertCharsToStrings(hfile.XData));
%xinds = 1:length(xpoints);

ypoints = cellfun(@str2num,convertCharsToStrings(hfile.YData(1:end)));
%ypoints = [ypoints; inf];
%yinds = 1:length(ypoints);

im = zeros(length(xpoints), length(ypoints));
im(im==0) = nan;
for ii=1:length(xpoints)
    for jj=1:length(ypoints)
        im(ii,jj) = hfile.ColorData(jj, ii);
    end
end

figure(); clf;
% if using areal velocity, we need to convert to physical units here
if (isPlottingAreaVelocity)
    im = im*areaVelocityConvert;
end
heatmap(im)

figure(); clf;

%h = heatmap(xlabs,ylabs,interpIm');
%h = heatmap(xpoints, 1.0./(ypoints), im');
h = heatmap(xpoints, flip(ypoints), flip(im'));
h.FontSize = 20;
h.GridVisible = 'off';
h.CellLabelColor = 'none';
cdl = h.XDisplayLabels; % current display labels
skipNum = round(length(h.XDisplayLabels)/10); % remove all but 10 labels
for ii=1:skipNum-1
    cdl(ii+1:skipNum:end) = {' '};
end

for ii=1:length(cdl)
    cdl(ii) = {num2str(bulkModulusConvert*str2num(cdl{ii}))};
end
h.XDisplayLabels = cdl;

for ii=1:length(h.YDisplayLabels)
    if (h.YDisplayLabels{ii} == "Inf")
        h.YDisplayLabels{ii} = '\infty';
    end
end

cdl = h.YDisplayLabels; 
skipNum = round(length(h.YDisplayLabels)/10);
for ii=1:skipNum-1
    cdl(ii+1:skipNum:end) = {' '};
end

for ii=1:length(cdl)
    cdl(ii) = {num2str(timeConvert*str2num(cdl{ii}))};
end

h.YDisplayLabels = cdl;
%h.YDisplayLabels = compose('%.4g',str2double(h.YDisplayLabels));
h.GridVisible = 'off';

%% Show interpolated heatmap and scatter data samples on linear axes
figure(3); clf; hold on;

% interpolate between points to determine an underlying color
interpX = interpn(xpoints*bulkModulusConvert, 3);
interpY = interpn(ypoints*timeConvert, 3);
[interpMeshX, interpMeshY] = meshgrid(interpX, interpY);
interpData = interp2(im', 3);

% after setting up interpData, use it as a heatmap underlying the real data
imagesc(interpX, interpY, interpData); colormap jet(16); axis xy;
colorbar;

inputCombos = combvec((xpoints*bulkModulusConvert)', (ypoints*timeConvert)');

scatter(inputCombos(1,:), inputCombos(2,:),30, 'k', 'o', 'filled', 'MarkerFaceColor', 'white', 'MarkerEdgeColor', 'k', 'linewidth', 1);

xlabel('B (kPa)', 'FontSize', 20, 'interpreter','latex')
ylabel('$\tau$ (minutes)', 'FontSize', 20, 'Interpreter', 'latex')
xlim([min(xpoints) max(xpoints)]*bulkModulusConvert);
ylim([min(ypoints) max(ypoints)]*timeConvert);
set(gca, 'YScale', 'log', 'XScale', 'log')
set(gca,'FontSize', 20)
% 
% %integer format for x and y ticks
% yt = get(gca,'ytick');
% xt = get(gca,'xtick');
% for j=1:length(yt)
%     YTL{1,j} = num2str(yt(j),'%d');
% end
% yticklabels(YTL);
% for j=1:length(xt)
%     if (xt(j) >= 1)
%         XTL{1,j} = num2str(xt(j),'%d');
%     elseif (xt(j) >= 0.1)
%         XTL{1,j} = num2str(xt(j), '%.1f');
%     elseif (xt(j) >= 0.01)
%         XTL{1,j} = num2str(xt(j), '%.2f');
%     end
% end
% xticklabels(XTL);
