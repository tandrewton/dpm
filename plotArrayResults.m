close all; clear
addpath('/Users/AndrewTon/Documents/YalePhD/projects/dpm/bash')
addpath('C:\Users\atata\projects\dpm\bash')
addpath('/Users/AndrewTon/Documents/YalePhD/projects/dpm/matlab_funcs')
addpath('C:\Users\atata\projects\dpm\matlab_funcs')

% first identify all the strings I'll need to run over

%    "C:\Users\atata\projects\dpm\pipeline/cells/ablate/

% ablate_A01.05_t_stress125.0k_l1.0_k_a1.0_kb_0.1_w_ps0.005_dsq4.0_k_ps4.0_k_lp4.0_d_flag0.0_bd0_sm1
% / _N50_Dur500_att0.1_sd1_sd10_sd1.pos

% simulation parameters go here
runType = "ablate";
N="50";
%ndelete="6";
%calA0="1.10";
strainRate_ps="0.005";
%deltaSq = "4.0";
k_a = "1.0";
k_l = "1.0";
k_b = "0.01";
k_ps = "4.0"; %purse-string spring constant
k_lp = "4.0"; %lamellipodia spring constant
sm = "1";
tau_lp = "1.0"; %lamellipodia lifetime
%d_flag = "0.0"; %lamellipodia max length
boundaryType = "0"; 
%att="0.2";
B="1.0";
%bd = "1";
boolCIL="0";
Duration="1000";

numSeeds = 10;
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
array_output_dir = subdir_output + "array_output_figures/";

% all possible parameter variations go here
% for numPlots, just multiply the lengths of all the parameter arrays

% select 2 arrays to have length > 1. script will generate plots for
% variations in these two parameters

N_arr = ["50"];                 %i
calA0_arr = ["1.05"];           %ii
t_stress_arr = ["1.0" "5.0" "25.0" "125.0" "625.0"]; %iii
%t_stress_arr = ["125.0"]; %iii
%att_arr = ["0.05" "0.1" "0.15" "0.2" "0.25" "0.29"]; %j
att_arr = ["0.1"]; % j
om_arr = ["0.005"];             %jj
%om_arr = ["0.001" "0.005" "0.01" "0.05"];             %jj
kl_arr = ["1.0"]; %jjj
%kl_arr = ["0.1" "0.5" "1.0" "5.0" "10.0"]; %jjj
%ka_arr = ["1.0"];               %k
ka_arr = ["0.1" "0.5" "1.0" "5.0" "10.0"];    %k
kb_arr = ["0.01"]; %kk
%kb_arr = ["0" "0.001" "0.01" "0.1"]; %kk
deltaSq_arr = ["4.0"];          %kkk
d_flag_arr = ["0.0"];           %l


% loop logic selects the 2 arrays above with length > 1.
%  it then identifies the correct parameter name, and the iterator in the
%  big nested loop that runs all of the plotting routines.
pm1 = []; pm2 = [];
param_cell_list = {t_stress_arr att_arr om_arr kl_arr ka_arr kb_arr};
param_id_list = ["tau" "adhesion" "activity" "kl" "ka" "kb"];
param_ind_list = [1 2 3 4 5 6];
for i=1:length(param_cell_list)
    if (length(param_cell_list{i}) > 1)
        if (isempty(pm1))
            pm1 = param_cell_list{i};
            pm1_str = param_id_list(i);
            pm1_ind_num = param_ind_list(i);
        elseif (isempty(pm2))
            pm2 = param_cell_list{i};
            pm2_str = param_id_list(i);
            pm2_ind_num = param_ind_list(i);
        else
            disp("error, more than 2 elements of cell list are nonempty");
            assert(false);
        end
    end
end

assert(~isempty(pm1) && ~isempty(pm2));

pm1pm2_folder = "cfg_"+pm1_str+"_"+pm2_str+"/";
mkdir(array_output_dir+pm1pm2_folder);

bigproduct = length(N_arr)*length(calA0_arr)*length(t_stress_arr)*...
    length(att_arr)*length(om_arr)*length(kl_arr)*...
    length(ka_arr)*length(kb_arr)*length(deltaSq_arr)*length(d_flag_arr);
numPlots = bigproduct;

showLastFrameOfSimulations = true;
isCrawling = false;
% set up plotting windows

numPlotTypes = 2; % parameter 1 (t_stress), parameter 2 (kb)
for i=1:numPlots
    for j=1:numPlotTypes
        figure(j+numPlotTypes*(i-1)); clf; hold on;    
    end
end

heatmap_fig_num = numPlots*numPlotTypes+1;

% 4 model outputs
figure(heatmap_fig_num); clf; %shape_max(any time)/shape(starting) of rosette cells
figure(heatmap_fig_num+1); clf;%healing time  
figure(heatmap_fig_num+2); clf; % rosette number
figure(heatmap_fig_num+3); clf; % final shape parameters of rosette cells

numHeatmaps = 4;
heatmap1 = zeros(length(pm1), length(pm2)); 
heatmap2 = zeros(length(pm1), length(pm2));
heatmap3 = zeros(length(pm1), length(pm2));
heatmap4 = zeros(length(pm1), length(pm2)); 
heatmap4_std = zeros(length(pm1), length(pm2));

numItsTotal = 1;

for i=1:length(N_arr)
    N = N_arr(i);
    for ii=1:length(calA0_arr)
        calA0 = calA0_arr(ii);
        for iii=1:length(t_stress_arr)
            t_stress = t_stress_arr(iii);
            for j=1:length(att_arr)
                att = att_arr(j);
                for jj=1:length(om_arr)
                    strainRate_ps = om_arr(jj);
                    for jjj=1:length(kl_arr)
                        k_l = kl_arr(jjj);
                        for k=1:length(ka_arr)
                            k_a = ka_arr(k);
                            for kk=1:length(kb_arr)
                                k_b = kb_arr(kk);
                                for kkk=1:length(deltaSq_arr)
                                    deltaSq = deltaSq_arr(kkk);
                                    for l=1:length(d_flag_arr)
                                        d_flag = d_flag_arr(l);

                                        iterator_arr = [iii j jj jjj k kk];

                                        pm1_ind = iterator_arr(pm1_ind_num);
                                        pm2_ind = iterator_arr(pm2_ind_num);

                                        voidArea = zeros(0,2);
                                        meanInnerShapes = NaN(0,1);
                                        timeInnerShapes = zeros(0,1);
                                        timestep = 0; % determine on the fly to pad with zeros
                                        innerShapeArr = []; % fill with meanInnerShape and dynamically pad rows with nans
                                        stdMeanShapeArr = []; % fill with mean, std of shapes of rosette cells in last frame of each seed
                                        healingTime = 0;
                                        rosetteNumber = 0;

                                        if (isCrawling)
                                            displayStr = "C__A0="+calA0+",att="+att+",sm="+sm+",t_{stress}="+t_stress;
                                        else
                                            displayStr = "P"+"N="+N+",tau="+t_stress+",att="+att+",om="+...
                                                strainRate_ps+",kl="+k_l+",ka="+k_a+",kb="+k_b+",dsq="+deltaSq;
                                        end

                                        numGoodSeeds = numSeeds;
                                        for m=1:numSeeds
                                            %for m=1:1
                                            % construct filenames to find the right simulation
                                            seed = m;
                                            run_name=runType+"_A0"+calA0+"_t_stress"+t_stress+"k_l"+...
                                                k_l+"_k_a"+k_a+"_k_b"+k_b+"_w_ps"+strainRate_ps+ ...
                                                "_dsq"+deltaSq+"_k_ps"+k_ps+"_k_lp"+k_lp+...
                                                "_d_flag"+d_flag+"_bd"+boundaryType+"_sm"+sm;
                                            pipeline_dir =  subdir_pipeline + run_name + "/";
                                            output_dir = subdir_output + run_name + "/";

                                            if ~exist(pipeline_dir, 'dir')
                                                mkdir(pipeline_dir)
                                            end
                                            if ~exist(output_dir, 'dir')
                                                mkdir(output_dir)
                                            end

                                            fileheader="_N"+N+"_Dur"+Duration+"_att"+att+"_sd"+ ...
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

                                            % if load fails, go to next
                                            % loop iteration
                                            try 
                                                voidArea_sd = load(voidAreaStr);
                                                bulkCellShape_sd = load(bulkStr);
                                                woundProperties_sd = load(woundPropertiesStr);
                                                cellID = load(innerAndBulkCellIDStr);
                                            catch
                                                disp('Did not load file: '+pipeline_dir+fileheader);
                                                numGoodSeeds = numGoodSeeds - 1;
                                                continue;
                                            end
                                            voidArea_sd(voidArea_sd == 1e10) = NaN;
                                            innerShapes_sd = bulkCellShape_sd(:,[1; cellID(:,3)]==1);
                                            meanInnerShapes_sd = nanmean(innerShapes_sd(:,2:end),2);
                                            timeInnerShapes_sd = innerShapes_sd(:,1);

                                            if (showLastFrameOfSimulations && seed == 1)
                                                figure(1000 + numItsTotal); hold on;
                                                [trajectoryData, cell_count] = readDPMClassPosOutput(nvestr);
                                                ff = trajectoryData.NFRAMES;

                                                % get cell positions
                                                xpos = trajectoryData.xpos(ff,:);
                                                ypos = trajectoryData.ypos(ff,:);
                                                gi = trajectoryData.gi(ff,:);
                                                l0 = trajectoryData.l0(ff,:);
                                                vrad = trajectoryData.vrad(ff,:);
                                                psi = trajectoryData.psi(ff,:);
                                                NCELLS = cell_count(ff);
                                                nv = trajectoryData.nv;

                                                [nvUQ, ~, IC] = unique(nonzeros(nv(ff,:)));
                                                NUQ = length(nvUQ);
                                                cellCLR = jet(NUQ);

                                                for nn = 1:NCELLS
                                                    xtmp = xpos{nn};
                                                    ytmp = ypos{nn};
                                                    gitmp = gi{nn};
                                                    l0tmp = l0{nn};
                                                    vradtmp = vrad{nn};
                                                    psitmp = psi(nn);
                                                    costmp = cos(psitmp);
                                                    sintmp = sin(psitmp);

                                                    clr = cellCLR(IC(nn),:);

                                                    cx = mean(xtmp);
                                                    cy = mean(ytmp);

                                                    rx = xtmp - cx;
                                                    ry = ytmp - cy;
                                                    rads = sqrt(rx.^2 + ry.^2);
                                                    xtmp = xtmp + 0.4*l0tmp(1)*(rx./rads);
                                                    ytmp = ytmp + 0.4*l0tmp(1)*(ry./rads);
                                                    vpos = [xtmp, ytmp];
                                                    finfo = [1:nv(ff,nn) 1];
                                                    %disp("finfo is "+ finfo)
                                                    patch('Faces',finfo,'vertices',vpos,'FaceColor',clr,'EdgeColor','k','linewidth',2);
                                                end

                                                % show simulation
                                                % parameters on the figure
                                                ann = annotation('textbox', [.1 .9 0.1 0.1],'interpreter', 'latex',...
                                                'String', displayStr,'FitBoxToText','on','Edgecolor','none',...
                                                'FaceAlpha', 0.5, 'backgroundcolor','white','interpreter', 'latex');
                                                ann.FontSize = 10;

                                                axis equal; axis off;
                                                saveas(gcf, array_output_dir+'/'+pm1pm2_folder+displayStr+'.eps', 'epsc')
                                            end

                                            % pad shorter voidArea with zeros to add them together
                                            % pad shorted voidArea with its last known value to add
                                            % them together?
                                            % pad meanInnerShapes end+1:length with NaNs to nanmean them later
                                            if (length(voidArea) < length(voidArea_sd))
                                                if (voidArea > 0)
                                                    voidArea(length(voidArea_sd),:) = voidArea(length(voidArea), :);
                                                else
                                                    voidArea(length(voidArea_sd),:) = 0;
                                                end
                                                meanInnerShapes(end+1:length(meanInnerShapes_sd),:) = nan;
                                                voidArea(:,1) = voidArea_sd(:,1); %extend time column
                                                timeInnerShapes = timeInnerShapes_sd; %extend time column
                                            elseif (length(voidArea_sd) < length(voidArea))
                                                if (voidArea_sd > 0)
                                                    voidArea_sd(length(voidArea),:) = voidArea_sd(length(voidArea_sd),:);
                                                else
                                                    voidArea_sd(length(voidArea),:) = 0;
                                                end
                                                meanInnerShapes_sd(end+1:length(meanInnerShapes),:) = nan;
                                                voidArea_sd(:,1) = voidArea(:,1); %extend time column
                                            end
                                            % otherwise they have exactly the same length, so don't
                                            % adjust lengths

                                            voidArea(:,2) = voidArea(:,2) + voidArea_sd(:,2);
                                            healingTime = healingTime + woundProperties_sd(1);
                                            rosetteNumber = rosetteNumber + woundProperties_sd(2);
                                            % get the distribution quantities of the rosette cells in their last frame
                                            [stdShapes_sd, meanShapes_sd] = std(innerShapes_sd(end, 2:end));

                                            %voidArea = voidArea + voidArea_sd;
                                            %meanInnerShapes = meanInnerShapes + meanInnerShapes_sd;
                                            sizeInnerShapeArr = size(innerShapeArr);
                                            differenceSize = sizeInnerShapeArr(2) - length(meanInnerShapes_sd);
                                            if (differenceSize < 0)
                                                innerShapeArr = padarray(innerShapeArr, [0 -differenceSize], nan, 'post');
                                            elseif (differenceSize > 0)
                                                meanInnerShapes_sd = padarray(meanInnerShapes_sd, [0 differenceSize], nan, 'post');
                                            end
                                            % add meanInnerShapes_sd as a
                                            % new row to innerShapeArr
                                            innerShapeArr(end+1, :) = meanInnerShapes_sd;
                                            stdMeanShapeArr(end+1,:) = [stdShapes_sd meanShapes_sd];
                                        end

                                        voidArea(:,2) = voidArea(:,2) / numGoodSeeds;
                                        healingTime = healingTime / numGoodSeeds;
                                        rosetteNumber = rosetteNumber / numGoodSeeds;
                                        minInnerShapes = mean(min(innerShapeArr, [], 2),'omitnan');
                                        maxInnerShapes = mean(max(innerShapeArr, [], 2),'omitnan');
                                        meanRosetteShape = mean(stdMeanShapeArr(:,2));
                                        stdRosetteShape = mean(stdMeanShapeArr(:,1).^2);

                                        % plot area vs time 
                                        figure(pm1_ind)

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

                                        %heatmap1(shapeii,j) = max(meanInnerShapes)/min(meanInnerShapes);
                                        %heatmap2(shapeii,j) = max(innerShapes_sd(:,1));
                                        %heatmap3(shapeii,j) = woundProperties_sd(2);

                                        heatmap1(pm1_ind,pm2_ind) = maxInnerShapes/minInnerShapes;
                                        %heatmap2(i,j) = max(innerShapes_sd(:,1));
                                        heatmap2(pm1_ind,pm2_ind) = healingTime;
                                        heatmap3(pm1_ind,pm2_ind) = rosetteNumber;
                                        heatmap4(pm1_ind, pm2_ind) = meanRosetteShape;
                                        heatmap4_std(pm1_ind, pm2_ind) = stdRosetteShape;

                                        numItsTotal = numItsTotal + 1;
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

if (isCrawling)
    % % how to save figures:
    % need to get them saved in order
    FolderName = "output/cells/ablate/array_output_figures/activity_sweep";   % Your destination folder
    mkdir(FolderName);
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
heatmap_filepath = "output/cells/ablate/array_output_figures/heatmaps";
mkdir(heatmap_filepath);

figure(heatmap_fig_num + numHeatmaps)
h = heatmap(pm2, pm1, heatmap1);
h.YLabel = pm1_str;
h.XLabel = pm2_str;

%h.YLabel = '$\mathcal{A}_0$';
h.Title = '$\langle\mathcal{A}_{max}\rangle / \langle\mathcal{A}_{0}\rangle$';
h.NodeChildren(3).XAxis.Label.Interpreter = 'latex';
h.NodeChildren(3).YAxis.Label.Interpreter = 'latex';
h.NodeChildren(3).Title.Interpreter = 'latex';
fontsize(gca, scale=1.5)
colormap parula
saveas(gcf, heatmap_filepath+"/relativeMaxShapeDeformation"+ ...
             "PS"+~isCrawling+"-"+pm1_str+"-"+pm2_str+".eps", 'epsc')

figure(heatmap_fig_num + numHeatmaps + 1)
h = heatmap(pm2, pm1, heatmap2);
h.YLabel = pm1_str;
h.XLabel = pm2_str;

h.Title = 'Healing time';
h.NodeChildren(3).XAxis.Label.Interpreter = 'latex';
h.NodeChildren(3).YAxis.Label.Interpreter = 'latex';
h.NodeChildren(3).Title.Interpreter = 'latex';
fontsize(gca, scale=1.5)
colormap parula
saveas(gcf, heatmap_filepath+"/healingTimesSmooth"+...
             "PS"+~isCrawling+"-"+pm1_str+"-"+pm2_str+".eps", 'epsc')

figure(heatmap_fig_num + numHeatmaps + 2)
heatmap3(heatmap3==0) = nan;
h = heatmap(pm2, pm1, heatmap3);
h.YLabel = pm1_str;
h.XLabel = pm2_str;

h.Title = 'Rosette number';
h.NodeChildren(3).XAxis.Label.Interpreter = 'latex';
h.NodeChildren(3).YAxis.Label.Interpreter = 'latex';
h.NodeChildren(3).Title.Interpreter = 'latex';
fontsize(gca, scale=1.5)
colormap parula
saveas(gcf, heatmap_filepath+"/rosetteNumberSmooth" + ...
             "PS"+~isCrawling+"-"+pm1_str+"-"+pm2_str+".eps", 'epsc')

for i=1:2
    % hack - not sure why I need to reset figure to get right fontsize, but
    % not going to try to fix it now. loop iterator i runs the same code
    % twice, which gets around this fontsize problem.

    figure(heatmap_fig_num + numHeatmaps + 3)
    %convert matrices heatmap4, stdevs
    clabel = arrayfun(@(x,y){sprintf('%0.2f +/- %0.2f',x,y)}, heatmap4, heatmap4_std);
    
    heatmap_custom(heatmap4, cellstr(pm2), cellstr(pm1), clabel,'Colorbar',true,'FontSize', 6);
    ylabel(pm1_str);
    xlabel(pm2_str);
    title('$\langle\mathcal{A}_{rosette}\rangle$', 'Interpreter','latex');
    fontsize(gca, scale=1.5)
    colormap parula
    saveas(gcf, heatmap_filepath+"/rosetteCellsFinalShape"+ ...
                 "PS"+~isCrawling+"-"+pm1_str+"-"+pm2_str+".eps", 'epsc')
end
