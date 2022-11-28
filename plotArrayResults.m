close all; clear
addpath('/Users/AndrewTon/Documents/YalePhD/projects/dpm/bash')
addpath('C:\Users\atata\projects\dpm\bash')
addpath('/Users/AndrewTon/Documents/YalePhD/projects/dpm/matlab_funcs')
addpath('C:\Users\atata\projects\dpm\matlab_funcs')

% first identify all the strings I'll need to run over

%    "C:\Users\atata\projects\dpm\pipeline/cells/ablate/
% ablate_A01.10_k_l1.0_w_ps0.005_dsq0.0_k_ps4.0_k_lp4.0_t_lp1.0_d_flag3.0_bd0_sm1/
% ablate_A01.10_k_l1.0_w_ps0.005_dsq0.0_k_ps4.0_k_lp4.0_t_lp1.0_d_flag3.0_bd0_sm1_N40_Dur1000_att0.20_sd1_sd1_sd1.pos"

% simulation parameters go here
runType = "ablate";
N="30";
%ndelete="6";
%calA0="1.10";
strainRate_ps="0.005";
%deltaSq = "2.0";
%k_a = "1.0";
k_l = "1.0";
k_ps = "4.0"; %purse-string spring constant
k_lp = "4.0"; %lamellipodia spring constant
%smooth = "0";
tau_lp = "1.0"; %lamellipodia lifetime
%d_flag = "0.0"; %lamellipodia max length
boundaryType = "0"; 
%att="0.2";
B="1.0";
Dr0="0.5";
boolCIL="0";
Duration="1000";

numSeeds = 1;
startSeed = 1;
max_seed = numSeeds;
set(0,'DefaultFigureWindowStyle','docked')

%PC directory
%pc_dir = "/Users/AndrewTon/Documents/YalePhD/projects/dpm/";
pc_dir="C:\Users\atata\projects\dpm\";
%pipeline is the location of data generated during simulations
subdir_pipeline = pc_dir + "pipeline/cells/"+runType+"/";
%output is location of results of this postprocessing
subdir_output = pc_dir + "output/cells/"+runType+"/";

calA0_arr = ["1.0" "1.05" "1.10"];
att_arr = ["0.04" "0.08" "0.12" "0.16" "0.20"];
sm_arr = ["0" "1"];
k_a_arr = ["0.5" "1.0" "2.0" "4.0"];

isCrawling = true;

if (isCrawling)
    dsq = "0.0"; % for C
    d_flag = "3.0"; % for C
    numPlots = length(calA0_arr)*length(sm_arr)*length(k_a_arr); %for C
else
    dsq = "4.0"; % for PS
    d_flag = "0.0"; % for PS
    numPlots = length(calA0_arr)*length(sm_arr)*length(k_a_arr); %for PS
end

% set up plotting windows


% 3 values of calA0, 2 values of smoothness, 4 values of yield length
numPlotTypes = 2; % voidArea and shape vs time
for i=1:numPlots
    for j=1:numPlotTypes
        figure(j+numPlotTypes*(i-1)); clf; hold on;    
    end
end
heatmap_fig_num = numPlots*numPlotTypes+1;
figure(heatmap_fig_num); clf; %shape_max/shape_bulk
figure(heatmap_fig_num+1); clf;%healing time
figure(heatmap_fig_num+2); clf; % rosette number

heatmap1 = zeros(length(calA0_arr), length(calA0_arr)); 
heatmap2 = zeros(length(calA0_arr), length(calA0_arr));
heatmap3 = zeros(length(calA0_arr), length(calA0_arr)); 

for shapeii=1:length(calA0_arr)
    calA0=calA0_arr(shapeii);
    for i=1:length(att_arr)
        att = att_arr(i);
        for j=1:length(k_a_arr)
            k_a = k_a_arr(j);
            for k=1:length(sm_arr)
                sm = sm_arr(k);
                % might also be looping over m=numSeeds to accumulate some results
                voidArea = zeros(0,2);
                meanInnerShapes = NaN(0,1);
                timeInnerShapes = zeros(0,1);
                timestep = 0; % determine on the fly to pad with zeros
                innerShapeArr = []; % fill with meanInnerShape and dynamically pad rows with nans
                for m=1:numSeeds
                    % construct filenames to find the right simulation
                    bd = "0";
                    seed = m;
                    run_name =runType+"_A0"+calA0+"_k_l"+k_l+"k_a"+k_a+"_w_ps"+ ...
                        strainRate_ps+ "_dsq"+deltaSq+"_k_ps"+k_ps+"_k_lp"+ ...
                        k_lp+"_d_flag"+d_flag+"_bd"+boundaryType+"_sm"+sm;
                    pipeline_dir =  subdir_pipeline + run_name + "/";
                    output_dir = subdir_output + run_name + "/";

                    if ~exist(pipeline_dir, 'dir')
                        mkdir(pipeline_dir)
                    end
                    if ~exist(output_dir, 'dir')
                        mkdir(output_dir)
                    end

                    fileheader=run_name +"_N"+N+"_Dur"+Duration+"_att"+att+"_sd"+ ...
                        startSeed+"_sd"+max_seed+"_sd"+seed;
                    fileheader_short = "_N"+N+"_Dur"+Duration+"_att"+att+"_sd"+seed;
                    nvestr = pipeline_dir+fileheader+'.pos';
                    energystr = pipeline_dir+fileheader+'.energy';
                    stressstr = pipeline_dir+fileheader+'.stress';
                    boundaryStr = pipeline_dir+fileheader+ ".void";
                    edgeStr = pipeline_dir+fileheader+ '.edge';
                    purseStr = pipeline_dir+fileheader+ '.purseString';
                    voidAreaStr = pipeline_dir+fileheader+ '.voidArea';
                    innerStr = pipeline_dir+fileheader+ '.innerCellShape';
                    bulkStr = pipeline_dir+fileheader+ '.bulkCellShape';
                    woundPropertiesStr = pipeline_dir+fileheader+ '.woundProperties';
                    innerAndBulkCellIDStr = pipeline_dir+fileheader+ '.cellID';
    
                    voidArea_sd = load(voidAreaStr);
                    voidArea_sd(voidArea_sd == 1e10) = NaN;
                    bulkCellShape_sd = load(bulkStr);
                    woundProperties_sd = load(woundPropertiesStr);
                    cellID = load(innerAndBulkCellIDStr);
                    innerShapes_sd = bulkCellShape_sd(:,[1; cellID(:,3)]==1);
                    meanInnerShapes_sd = nanmean(innerShapes_sd(:,2:end),2);
                    timeInnerShapes_sd = innerShapes_sd(:,1);

                    % pad shorter voidArea with zeros to add them together
                    % pad meanInnerShapes end+1:length with NaNs to nanmean them later 
                    if (length(voidArea) < length(voidArea_sd))
                        voidArea(length(voidArea_sd),:) = 0;
                        meanInnerShapes(end+1:length(meanInnerShapes_sd),:) = nan;
                        voidArea(:,1) = voidArea_sd(:,1); %extend time column
                        timeInnerShapes = timeInnerShapes_sd; %extend time column
                    elseif (length(voidArea_sd) < length(voidArea))
                        voidArea_sd(length(voidArea),:) = 0;
                        meanInnerShapes_sd(end+1:length(meanInnerShapes),:) = nan;
                        voidArea_sd(:,1) = voidArea(:,1); %extend time column
                    end
                    % otherwise they have exactly the same length, so don't
                    % adjust lengths

                    voidArea(:,2) = voidArea(:,2) + voidArea_sd(:,2);
                    %meanInnerShapes = meanInnerShapes + meanInnerShapes_sd;
                    sizeInnerShapeArr = size(innerShapeArr);
                    differenceSize = sizeInnerShapeArr(2) - length(meanInnerShapes_sd);
                    if (differenceSize < 0)
                        innerShapeArr = padarray(innerShapeArr, [0 -differenceSize], nan, 'post');
                    elseif (differenceSize > 0)
                        meanInnerShapes_sd = padarray(meanInnerShapes_sd, [0 differenecSize], nan, 'post');
                    end
                    innerShapeArr(end+1, :) = meanInnerShapes_sd;
                end

                voidArea = voidArea / numSeeds;
                %meanInnerShapes = meanInnerShapes / numSeeds;
                meanInnerShapes = nanmean(innerShapeArr, 1);

                
                if (isCrawling)
                    % plot area vs time for C
                    figure((shapeii-1)*length(sm_arr)*length(k_a_arr) ...
                        + (k-1)*length(k_a_arr) ...
                        + l)
                else 
                    %plot area vs time for PS
                    figure((shapeii-1)*length(sm_arr)*length(k_a_arr) ...
                        + (k-1)*length(k_a_arr) ...
                        + j)
                end
                if (isCrawling)
                    displayStr = "C: A0="+calA0+",att="+att+",sm="+sm+",ka="+ka;
                else
                    displayStr = "P: A0="+calA0+",att="+att+",sm="+sm+",ka="+ka;
                end

                plot(voidArea(:,1), voidArea(:,2), 'linewidth', 4, 'DisplayName', displayStr)
                xlabel('$t/\tau$','Interpreter','latex','fontsize', 24);
                ylabel('Area','Interpreter','latex','fontsize', 24);
                %set(gca,'Yscale','log')
                ylim([0 inf])
                legend('location','northeast','fontsize', 8)

%                     % plot shape vs time
%                     figure(length(calA0_arr)+shapeii + 100*(k-1) + 1000*(j-1))
%                    
%                     % cellID row = [ci inInitialWoundNeighbors inFinalWoundNeighbors]
%                     % we want to access inFinalWoundNeighbors of bulkCellShape
%                     % bulkCellShape row = [time shape(0) shape(1) ... shape(NCELLS)]
%                     %innerShapes = bulkCellShape(:,[1; cellID(:,3)]==1);
%                     %outerShapes = bulkCellShape(:,[1; ~cellID(:,3)]==1);
%                     %plot(innerShapes(:,1), nanmean(innerShapes(:,2:end),2),  ...
%                     %    'linewidth', 4, 'DisplayName', "inner shapes",...
%                     %     'Color', color_array(l),'LineStyle', style_array(j))
%                     %plot(timeAndOuterShapes(:,1),nanmean(timeAndOuterShapes(:,2:end),2), ...
%                     %    'linewidth', 4, 'DisplayName', "bulk shapes")
%                     plot(timeInnerShapes, meanInnerShapes, 'linewidth', 4, 'DisplayName', "A0="+calA0+",att="+att+",sm="+sm+",dsq="+deltaSq)
%                     xlabel('$t/\tau$','Interpreter','latex','fontsize', 24);
%                     ylabel('Shape','Interpreter','latex','fontsize', 24);
%                     legend('location','southeast','fontsize', 6)

                heatmap1(shapeii,j) = max(meanInnerShapes)/min(meanInnerShapes);
                heatmap2(shapeii,j) = max(innerShapes_sd(:,1));
                heatmap3(shapeii,j) = woundProperties_sd(2);
            end
        end
    end
end

if (isCrawling)
    % % how to save figures:
    % need to get them saved in order
    FolderName = "output/cells/ablate/array_output_figures/activity_sweep";   % Your destination folder
    FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
    for iFig = 1:length(FigList)
      FigHandle = FigList(iFig);
      FigName   = num2str(get(FigHandle, 'Number'));
      set(0, 'CurrentFigure', FigHandle);
      saveas(FigHandle, fullfile(FolderName, "C"+FigName+".png"))
    end
else
    % % how to save figures:
    % need to get them saved in order
    FolderName = "output/cells/ablate/array_output_figures/activity_sweep";   % Your destination folder
    FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
    for iFig = 1:length(FigList)
      FigHandle = FigList(iFig);
      FigName   = num2str(get(FigHandle, 'Number'));
      set(0, 'CurrentFigure', FigHandle);
      saveas(FigHandle, fullfile(FolderName, "PS"+FigName+".png"))
    end
end

%plot array results now : heatmaps
heatmap_filepath = "output/cells/ablate/array_output_figures";
mkdir(heatmap_filepath);


figure(heatmap_fig_num)
h = heatmap(att_arr, calA0_arr, heatmap1);
h.XLabel = 'att';
h.YLabel = '$\mathcal{A}_0$';
h.Title = '$\langle\mathcal{A}_{max}\rangle / \langle\mathcal{A}_{0}\rangle$';
h.NodeChildren(3).XAxis.Label.Interpreter = 'latex';
h.NodeChildren(3).YAxis.Label.Interpreter = 'latex';
h.NodeChildren(3).Title.Interpreter = 'latex';
saveas(gcf, heatmap_filepath+"/relativeMaxShapeDeformationSmooth"+ sm_arr(1) + ...
             '.eps', 'epsc')

figure(heatmap_fig_num+1)
h = heatmap(att_arr, calA0_arr, heatmap2);
h.XLabel = 'att';
h.YLabel = '$\mathcal{A}_0$';
h.Title = 'Healing time';
h.NodeChildren(3).XAxis.Label.Interpreter = 'latex';
h.NodeChildren(3).YAxis.Label.Interpreter = 'latex';
h.NodeChildren(3).Title.Interpreter = 'latex';
saveas(gcf, heatmap_filepath+"/healingTimesSmooth"+ sm_arr(1) + ...
             '.eps', 'epsc')

figure(heatmap_fig_num+2)
heatmap3(heatmap3==0) = nan;
h = heatmap(att_arr, calA0_arr, heatmap3);
h.XLabel = 'att';
h.YLabel = '$\mathcal{A}_0$';
h.Title = 'Rosette number';
h.NodeChildren(3).XAxis.Label.Interpreter = 'latex';
h.NodeChildren(3).YAxis.Label.Interpreter = 'latex';
h.NodeChildren(3).Title.Interpreter = 'latex';
saveas(gcf, heatmap_filepath+"/rosetteNumberSmooth"+ sm_arr(1) + ...
             '.eps', 'epsc')


