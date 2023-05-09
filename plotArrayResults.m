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
strainRate_ps="1.0  ";
%deltaSq = "4.0";
k_a = "1.0";
k_l = "1.0";
k_b = "0.001";
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
Duration="2000";

numSeeds = 5;
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
calA0_arr = ["1.20"];           %ii
%t_stress_arr = ["9830.4" "39321.6"]; %iii
%t_stress_arr = ["19.2" "4915.2" "9830.4"];
t_stress_arr = ["40.0" "100.0" "1000.0" "1500.0" "5500.0"];
%t_stress_arr=["2.4" "4.8" "9.6" "19.2" "76.8" "307.2" "1228.8" "4915.2" "9830.4" "39321.6"];
%t_stress_arr=["76.8" "9830.4"];
att_arr = ["0.1"]; % j
om_arr = ["1.0"]; %jj
kl_arr = ["1.0"]; %jjj
%ka_arr = ["0.25" "1.0" "5.0" "10.0" "50.0"];               %k
%ka_arr=["0.25" "0.5" "1.0" "2.0" "4.0" "8.0" "16.0" "32.0" "64.0" "128.0" "256.0"]; %k
%ka_arr=["0.5" "1.0" "2.5" "5.0" "12.5" "25.0" "50.0"];
%ka_arr=["16.0" "64.0"];
ka_arr=["16.0" "32.0"];
kb_arr = ["0.01"]; %kk
deltaSq_arr = ["4.0"];          %kkk
%d_flag_arr = ["0.0"];           %l
%tau_s_arr=["19.2" "76.8" "307.2" "1228.8" "4915.2" "9830.4"]; %l
%tau_s_arr=["153.6" "307.2" "614.4" "1228.8"];
%tau_s_arr=["500.0"];
tau_s_arr=["0.25"];

d_flag = "0.0"; 

% loop logic selects the 2 arrays above with length > 1.
%  it then identifies the correct parameter name, and the iterator in the
%  big nested loop that runs all of the plotting routines.
pm1 = []; pm2 = [];
param_cell_list = {t_stress_arr att_arr om_arr kl_arr ka_arr kb_arr tau_s_arr};
param_id_list = ["tau" "adhesion" "activity" "kl" "ka" "kb" "tau_s"];
param_ind_list = [1 2 3 4 5 6 7];
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

if (length(pm1) == length(pm2) && length(pm1) == 2)
    disp("comparing parameter pair! for publication")
    isProductionRun = true;
    isUsePhysicalQuantites = true;
    
    % physical units here
    cellArea = [25,16]; % microns^2
    constrictionRate = 0.3; % micron/sec 
    adhesionForce = 1e-9; % Newton
    
    % simulation to real unit conversion
    bulkModulusConvert = adhesionForce./cellArea*1e12/1000; %kPa
    timeConvert = sqrt(cellArea)./constrictionRate/60; % minutes
    areaVelocityConvert = constrictionRate.*sqrt(cellArea)*60; %micron^2/min
else
    isUsePhysicalQuantities = false;
    isProductionRun=false;
end

pm1pm2_folder = "cfg_"+pm1_str+"_"+pm2_str+"/";
mkdir(array_output_dir+pm1pm2_folder);

bigproduct = length(N_arr)*length(calA0_arr)*length(t_stress_arr)*...
    length(att_arr)*length(om_arr)*length(kl_arr)*...
    length(ka_arr)*length(kb_arr)*length(deltaSq_arr)*length(tau_s_arr);
numPlots = bigproduct;

showLastFrameOfSimulations = false;
showPlots = true; % if false, don't show area vs time
showPhysicalUnits = 1;
isCrawling = false;

% set up plotting windows

numPlotTypes = 2; % parameter 1 (t_stress), parameter 2 (kb)
for i=1:numPlots
    for j=1:numPlotTypes
        figure(j+numPlotTypes*(i-1)); clf; hold on;    
    end
end

heatmap_fig_num = numPlots*numPlotTypes+1;

% 6 model outputs
figure(heatmap_fig_num); clf; %shape_max(any time)/shape(starting) of rosette cells
figure(heatmap_fig_num+1); clf;%healing time  
figure(heatmap_fig_num+2); clf; % rosette number
figure(heatmap_fig_num+3); clf; % final shape parameters of rosette cells
figure(heatmap_fig_num+4); clf; % initial shape parameters of rosette cells
figure(heatmap_fig_num+5); clf; % areal velocity from A_max-A_0.05 / (t_0.05 - t_max)

numHeatmaps = 6;
heatmap1 = zeros(length(pm1), length(pm2)); 
heatmap2 = zeros(length(pm1), length(pm2));
heatmap3 = zeros(length(pm1), length(pm2));
heatmap4 = zeros(length(pm1), length(pm2)); 
heatmap4_std = zeros(length(pm1), length(pm2));
heatmap5 = zeros(length(pm1), length(pm2)); 
heatmap5_std = zeros(length(pm1), length(pm2));
heatmap6 = zeros(length(pm1), length(pm2));

numItsTotal = 1;

for i=1:length(N_arr)
    N = N_arr(i);
    for ii=1:length(calA0_arr)
        calA0 = calA0_arr(ii);
        for iii=1:length(t_stress_arr)
            t_stress = t_stress_arr(iii);
            if (str2num(t_stress) >= 1e5 && length(t_stress_arr) > 1)
                pm1(end) = "\infty";
            end
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
                                    for l=1:length(tau_s_arr)
                                        tau_s = tau_s_arr(l);
                                        if (isProductionRun)
                                            %hardcoding some things
                                            %depending on if its for
                                            %production
                                            if (t_stress == "19.2")
                                                tau_s = "0";
                                            else
                                                %Duration = "2000";
                                            end
                                        end

                                        iterator_arr = [iii j jj jjj k kk l];

                                        pm1_ind = iterator_arr(pm1_ind_num);
                                        pm2_ind = iterator_arr(pm2_ind_num);

                                        voidArea = zeros(0,4);
                                        meanInnerShapes = NaN(0,1);
                                        timeInnerShapes = zeros(0,1);
                                        timestep = 0; % determine on the fly to pad with zeros
                                        innerShapeArr = []; % fill with meanInnerShape and dynamically pad rows with nans
                                        stdMeanShapeArr = []; % fill with mean, std of shapes of rosette cells in last frame of each seed
                                        stdMeanShapeInitialArr = [];
                                        healingTime = 0;
                                        rosetteNumber = 0;

                                        if (isCrawling)
                                            displayStr = "C__A0="+calA0+",att="+att+",sm="+sm+",t_{stress}="+t_stress;
                                        else
                                            displayStr = "P"+"N="+N+",tau="+t_stress+",att="+att+...
                                                ",kl="+k_l+",ka="+k_a+",kb="+k_b+",taus="+tau_s;
                                        end

                                        numGoodSeeds = numSeeds;
                                        for m=1:numSeeds
                                            %for m=1:1
                                            % construct filenames to find the right simulation
                                            seed = m;
                                            run_name=runType+"_A0"+calA0+"_t_stress"+t_stress+"k_l"+...
                                                k_l+"_k_a"+k_a+"_k_b"+k_b+"_w_ps"+strainRate_ps+ ...
                                                "_dsq"+deltaSq+"_k_ps"+k_ps+"_k_lp"+k_lp+...
                                                "_d_flag"+d_flag+"_taur"+tau_s+"_sm"+sm;
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
                                                innerCellShape_sd = load(innerStr);
                                                woundProperties_sd = load(woundPropertiesStr);
                                                cellID = load(innerAndBulkCellIDStr);
                                                voidArea_sd(voidArea_sd == 1e10) = NaN;
                                                innerShapes_sd = innerCellShape_sd;
                                                assert(~isempty(innerShapes_sd));
                                            catch
                                                disp('Did not load file: '+pipeline_dir+fileheader);
                                                numGoodSeeds = numGoodSeeds - 1;
                                                continue;
                                            end
                                            meanInnerShapes_sd = mean(innerShapes_sd(:,2:end),2, 'omitnan');
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
                                            %voidArea(:,4) = voidArea(:,4) + voidArea_sd(:,4);
                                            healingTime = healingTime + woundProperties_sd(1);
                                            rosetteNumber = rosetteNumber + woundProperties_sd(2);

                                            innerShapes_sd_no_nan = innerShapes_sd(end, 2:end);
                                            innerShapes_sd_no_nan = innerShapes_sd_no_nan(~isnan(innerShapes_sd_no_nan));

                                            if (length(innerShapes_sd_no_nan) > 0)
                                                % get the distribution quantities of the rosette cells in their last frame
                                                [stdRosetteShapesEnd_sd, meanRosetteShapesEnd_sd] = std(innerShapes_sd_no_nan);
                                                [stdRosetteShapesInitial_sd, meanRosetteShapesInitial_sd] = std(innerShapes_sd(1,2:end),0,2, 'omitnan');
                                            else
                                                % wound did not close in this case, so just take the largest 6
                                                % shapes for the end frame, and the average for the  initial frame
                                                [stdRosetteShapesEnd_sd, meanRosetteShapesEnd_sd] = std(maxk(bulkCellShape_sd(end,2:end),6));
                                                [stdRosetteShapesInitial_sd, meanRosetteShapesInitial_sd] = std(bulkCellShape_sd(1, 2:end));
                                            end

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
                                            stdMeanShapeArr(end+1,:) = [stdRosetteShapesEnd_sd meanRosetteShapesEnd_sd];
                                            stdMeanShapeInitialArr(end+1,:) = [stdRosetteShapesInitial_sd meanRosetteShapesInitial_sd];
                                        end

                                        voidArea(:,2) = voidArea(:,2) / numGoodSeeds;
                                        healingTime = healingTime / numGoodSeeds;
                                        rosetteNumber = rosetteNumber / numGoodSeeds;
                                        % get the minimum/maximum shape of each
                                        % cell, and take the mean over the
                                        % cells
                                        minInnerShapes = mean(min(innerShapeArr, [], 2),'omitnan');
                                        maxInnerShapes = mean(max(innerShapeArr, [], 2),'omitnan');

                                        meanRosetteInitialShape = mean(stdMeanShapeInitialArr(:,2));
                                        % to get average std, avg the
                                        % variances and take the sqrt
                                        stdRosetteInitialShape = sqrt(mean(stdMeanShapeInitialArr(:,1).^2));
                                        meanRosetteShape = mean(stdMeanShapeArr(:,2));
                                        stdRosetteShape = sqrt(mean(stdMeanShapeArr(:,1).^2));
                                        [maxVoidArea, argmaxVoidArea ]= max(voidArea(:,2));
                                        [~,timeFivePercent] = min(abs(voidArea(:,2) - 0.05*maxVoidArea));
                                        arealVelocity = (1-0.05) * maxVoidArea / (voidArea(timeFivePercent,1) - voidArea(argmaxVoidArea,1));
                                        % fprintf('%f %f %f %f %f\n', maxVoidArea, 0.05*maxVoidArea, timeFivePercent, voidArea(timeFivePercent,1), arealVelocity)

                                        if (showPlots)
                                            % plot area vs time 
                                            figure(pm1_ind) 
                                            colorList = ['r', 'k'];
                                            % if plotting reduced number of
                                            % parameters, then it's for a
                                            % production run and we're
                                            % going to also plot fits and
                                            % moving averages, and color
                                            % coordinate
                                            if (isProductionRun)
                                                figure(pm2_ind)
                                                % want (tau_1, B1) and (tau_2, B1). 
                                                % plot void area vs time
                                                %plot(voidArea(:,1)*timeConvert(pm1_ind), voidArea(:,2)*cellArea(pm1_ind),...
                                                %    'linewidth',5, 'Color', colorList(pm1_ind),'DisplayName', displayStr)
                                                n = 100;
                                                skipInt = length(voidArea(:,1))/n; % keep n points to plot 
                                                %scatter(voidArea(1:skipInt:end,1)*timeConvert(pm1_ind), voidArea(1:skipInt:end,2)*cellArea(pm1_ind),...
                                                %    20,
                                                %    colorList(pm1_ind),
                                                %    "^") % wound area in
                                                %    micron^2
                                                scatter(voidArea(1:skipInt:end,1)*timeConvert(pm1_ind), voidArea(1:skipInt:end,2)/(voidArea(1,2)),...
                                                    30, colorList(pm1_ind), "^") % wound area in area fraction
                                                legend off
                                                xlabel('Time (min)','Interpreter', 'tex','fontsize', 24);
                                                ylabel('Wound area fraction','Interpreter', 'tex','fontsize', 24);
                                                box on
                                                ax = gca;
                                                ax.TickLength = [0.025 0.025];
                                                ax.LineWidth = 1;
                                                fontsize(gcf, 20, "points")
                                                %xlim([0 inf])
                                                %ylim([0 inf])
                                                xlim([0 255]) % wing disc time range (largest time range)
                                                ylim([0 1]) % area fraction range
                                                pbaspect([1 1.3 1])

                                                % plot void area vs time,
                                                % experimental fit
                                                F = @(x,xdata)x(1)*exp(-x(2)*xdata) + x(3)*exp(-x(4)*xdata);
                                                if (pm1_ind == 1) % embryo
                                                    fit_params = [-2.7831 0.2966 3.6255 0.1172];
                                                    tlist = linspace(min(voidArea(:,1)*timeConvert(pm1_ind)), max(voidArea(:,1)*timeConvert(pm1_ind)), 20);
                                                    plot(tlist, F(fit_params, tlist), '-', 'Color', colorList(pm1_ind), 'linewidth', 2)
                                                elseif (pm1_ind == 2) % wing disc
                                                    fit_params = [0.9420 0.0690 0.3425 0.0107];
                                                    tlist = linspace(min(voidArea(:,1)*timeConvert(pm1_ind)), max(voidArea(:,1)*timeConvert(pm1_ind)), 20);
                                                    plot(tlist, F(fit_params, tlist), '-', 'Color', colorList(pm1_ind), 'linewidth', 2)
                                                end
                                                % plot shape vs time
                                               
                                                figure(pm2_ind+2);
                            
                                                skipInt = length(meanInnerShapes)/n;
                                                meanInnerShapes = mean(innerShapeArr,'omitnan');
                                                %plot(timeInnerShapes(1:skipInt:end)*timeConvert(pm1_ind),meanInnerShapes(1:skipInt:end),...
                                                %    '-','linewidth',3, 'Color', colorList(pm1_ind),'DisplayName', displayStr)
                                                scatter(timeInnerShapes(1:skipInt:end)*timeConvert(pm1_ind),meanInnerShapes(1:skipInt:end),...
                                                     30, colorList(pm1_ind), "^") 
                                                xlabel('Time (min)','Interpreter','tex','fontsize', 24);
                                                ylabel('$\mathcal{A}$','Interpreter','latex','fontsize', 24);
                                                %legend('location','southeast','fontsize', 6)
                                                legend off
                                                box on
                                                ax = gca;
                                                ax.TickLength = [0.025 0.025];
                                                ax.LineWidth = 1;
                                                fontsize(gcf, 20, "points")
                                                %xlim([0 inf])
                                                %ylim([-inf inf])
                                                pbaspect([1 1.2 1])

                                                % look for files called
                                                % shapeAll.mat, which
                                                % contain the fits of the
                                                % experimental data
                                                if (pm1_ind == 1)
                                                    load("shapeAllEmbryo.mat");
                                                elseif (pm1_ind == 2)
                                                    load("shapeAllWing.mat");
                                                end
                                                y=shapeAll(:,2);
                                                t=shapeAll(:,1);
                                                %sma = movmean(y,20); % simple moving average
                                                sma = distanceSMA(t, y, 10.0);
                                                skipInt = length(sma)/n;
                                                %plot(t,movmean(y,20),'-','linewidth',2,'Color', colorList(pm1_ind))
                                                %scatter(t(1:skipInt:end), sma(1:skipInt:end), 5, colorList(pm1_ind), 'filled')
                                                %scatter(t, sma, 10, colorList(pm1_ind), 'filled')
                                                plot(t, sma, '-','Color', colorList(pm1_ind), 'linewidth', 3)
                                                xlim([0 255]) % wing disc time range
                                                
                                                % save plot data in order
                                                % to manually overlay results varying in tau_s for
                                                % production plots 
                                                save(array_output_dir+'woundHealingPlotData/shape_pm1_'+pm1(pm1_ind)+'_pm2_'+pm2(pm2_ind)+'taus_'+tau_s+'.mat', 'meanInnerShapes')
                                                save(array_output_dir+'woundHealingPlotData/shape_time_pm1_'+pm1(pm1_ind)+'_pm2_'+pm2(pm2_ind)+'taus_'+tau_s+'.mat', 'timeInnerShapes')
                                                save(array_output_dir+'woundHealingPlotData/area_pm1_'+pm1(pm1_ind)+'_pm2_'+pm2(pm2_ind)+'taus_'+tau_s+'.mat', 'voidArea')
                                            
                                                if (false) % execute this code when producing plot for wound healing paper
                                                    % taus is 0, default
                                                    % value. but I want to show some results %where taus is nonzero as well on top
                                                    paramsToOverlay = "_pm1_"+"4915.2"+"_pm2_"+"16.0"+"taus_"+"?"; 
                                                    figure(4)
                                                    load(array_output_dir+'woundHealingPlotData/shape'+paramsToOverlay+'.mat');
                                                    load(array_output_dir+'woundHealingPlotData/shape_time'+paramsToOverlay+'.mat');
                                                    skipInt = length(meanInnerShapes)/n;
                                                    scatter(timeInnerShapes(1:skipInt:end)*timeConvert(pm1_ind),meanInnerShapes(1:skipInt:end),...
                                                     30, colorList(pm1_ind), "square")
                                                    fontsize(gcf, 14, "points")
                                                    yticks([1.2 1.4 1.6 1.8])

                                                    figure(2)
                                                    load(array_output_dir+'woundHealingPlotData/area'+paramsToOverlay+'.mat');
                                                    n = 100;
                                                    skipInt = length(voidArea(:,1))/n; % keep n points to plot 
                                                    scatter(voidArea(1:skipInt:end,1)*timeConvert(pm1_ind), voidArea(1:skipInt:end,2)/(voidArea(1,2)),...
                                                    30, colorList(pm1_ind), "square") % wound area in area fraction
                                                    fontsize(gcf, 14, "points")
                                                    yticks([0 0.5 1])
                                                end
                                            else
                                                plot(voidArea(:,1), voidArea(:,2), 'linewidth', 4, 'DisplayName', displayStr)
                                                xlabel('$t/\tau$','Interpreter','latex','fontsize', 24);
                                                ylabel('Area','Interpreter','latex','fontsize', 24);
                                            end

                                            %ylim([0 inf])
                                            legend('location','northeast','fontsize', 8)
                                            saveas(gcf, array_output_dir+"voidArea/"+"calA0"+calA0+"_"+pm1_str+"_"+pm2_str+"_"+pm1_ind+".png")
                                        
                                            % plot shape vs time
                                            figure(pm1_ind+max(5,length(pm1)))
                                            meanInnerShapes = mean(innerShapeArr,'omitnan');
                                            plot(timeInnerShapes,meanInnerShapes,...
                                                '-','linewidth',3,'DisplayName', displayStr)
                                            xlabel('Time','Interpreter','latex','fontsize', 24);
                                            ylabel('$\mathcal{A}$','Interpreter','latex','fontsize', 24);
                                            legend('location','southeast','fontsize', 6)
                                            ylim([1.0,inf])
                                            saveas(gcf, array_output_dir+"innerShapes/"+"calA0"+calA0+"_"+pm1_str+"_"+pm2_str+"_"+pm1_ind+".png")
                                        end
                                        %heatmap1(shapeii,j) = max(meanInnerShapes)/min(meanInnerShapes);
                                        %heatmap2(shapeii,j) = max(innerShapes_sd(:,1));
                                        %heatmap3(shapeii,j) = woundProperties_sd(2);

                                        heatmap1(pm1_ind,pm2_ind) = maxInnerShapes/minInnerShapes;
                                        %heatmap1(pm1_ind,pm2_ind) = ;
                                        %heatmap2(i,j) = max(innerShapes_sd(:,1));
                                        heatmap2(pm1_ind,pm2_ind) = healingTime;
                                        heatmap3(pm1_ind,pm2_ind) = rosetteNumber;
                                        heatmap4(pm1_ind, pm2_ind) = meanRosetteShape;
                                        heatmap4_std(pm1_ind, pm2_ind) = stdRosetteShape;

                                        heatmap5(pm1_ind, pm2_ind) = meanRosetteInitialShape;
                                        heatmap5_std(pm1_ind, pm2_ind) = stdRosetteInitialShape;

                                        heatmap6(pm1_ind, pm2_ind) = arealVelocity;

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
             "PS"+~isCrawling+"-"+"calA0"+calA0+"_"+pm1_str+"-"+pm2_str+".eps", 'epsc')

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
             "PS"+~isCrawling+"-"+"calA0"+calA0+"_"+pm1_str+"-"+pm2_str+".eps", 'epsc')

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
             "PS"+~isCrawling+"-"+"calA0"+calA0+"_"+pm1_str+"-"+pm2_str+".eps", 'epsc')


figure(heatmap_fig_num + numHeatmaps + 3)
% use the usual heatmap format for rosetteCellsFinalShape in order to save a heatmap object
h = heatmap(pm2, pm1, heatmap4);
h.YLabel = pm1_str;
h.XLabel = pm2_str;
h.Title = '$\langle\mathcal{A}_{rosette}\rangle$';
h.NodeChildren(3).XAxis.Label.Interpreter = 'latex';
h.NodeChildren(3).YAxis.Label.Interpreter = 'latex';
h.NodeChildren(3).Title.Interpreter = 'latex';
fontsize(gca, scale=1.5)
colormap parula
save(heatmap_filepath+"/rosetteCellsFinalShape"+ ...
                 "PS"+~isCrawling+"-"+"calA0"+calA0+"_"+pm1_str+"-"+pm2_str+".mat", 'h')
clf;
for i=1:2
    % hack - not sure why I need to reset figure to get right fontsize, but
    % not going to try to fix it now. loop iterator i runs the same code
    % twice, which gets around this fontsize problem.

    figure(heatmap_fig_num + numHeatmaps + 3)
    %convert matrices heatmap4, stdevs
    clabel = arrayfun(@(x,y){sprintf('%0.2f +/- %0.2f',x,y)}, heatmap4, heatmap4_std);
    
    heatmap_custom(heatmap4, cellstr(pm2), cellstr(pm1), clabel,'Colorbar',true,'FontSize', 6, 'TickTexInterpreter', 1);
    ylabel(pm1_str, 'interpreter', 'latex');
    xlabel(pm2_str, 'interpreter', 'latex');
    title('$\langle\mathcal{A}_{rosette}\rangle$', 'Interpreter','latex');
    fontsize(gca, scale=1.5)
    colormap parula
    saveas(gcf, heatmap_filepath+"/rosetteCellsFinalShape"+ ...
                 "PS"+~isCrawling+"-"+"calA0"+calA0+"_"+pm1_str+"-"+pm2_str+".eps", 'epsc')
end

for i=1:2
    figure(heatmap_fig_num + numHeatmaps + 4)
    %convert matrices heatmap5, stdevs
    clabel = arrayfun(@(x,y){sprintf('%0.2f +/- %0.2f',x,y)}, heatmap5, heatmap5_std);
    
    heatmap_custom(heatmap5, cellstr(pm2), cellstr(pm1), clabel,'Colorbar',true,'FontSize', 6, 'TickTexInterpreter', 1);
    ylabel(pm1_str, 'interpreter', 'latex');
    xlabel(pm2_str, 'interpreter', 'latex');
    title('$\langle\mathcal{A}_{rosette}(t=0)\rangle$', 'Interpreter','latex');
    fontsize(gca, scale=1.5)
    colormap parula
    saveas(gcf, heatmap_filepath+"/rosetteCellsInitialShape"+ ...
                 "PS"+~isCrawling+"-"+"calA0"+calA0+"_"+pm1_str+"-"+pm2_str+".eps", 'epsc')
end

figure(heatmap_fig_num + numHeatmaps + 5)
h = heatmap(pm2, pm1, heatmap6);
h.YLabel = pm1_str;
h.XLabel = pm2_str;

%h.YLabel = '$\mathcal{A}_0$';
h.Title = 'Healing speed (areal velocity)';
h.NodeChildren(3).XAxis.Label.Interpreter = 'latex';
h.NodeChildren(3).YAxis.Label.Interpreter = 'latex';
h.NodeChildren(3).Title.Interpreter = 'latex';
fontsize(gca, scale=1.5)
colormap parula
saveas(gcf, heatmap_filepath+"/arealVelocity"+ ...
             "PS"+~isCrawling+"-"+"calA0"+calA0+"_"+pm1_str+"-"+pm2_str+".eps", 'epsc')
save(heatmap_filepath+"/arealVelocity"+ ...
             "PS"+~isCrawling+"-"+"calA0"+calA0+"_"+pm1_str+"-"+pm2_str+".mat", 'h')

function movingMean = distanceSMA(t, y, deltaT) 
    % calculate the moving average of points y corresponding to the
    % distance deltaT in units of t
    movingMean = zeros(size(y));
    for ii=1:length(y)
        dataInRange = y(abs(t - t(ii)) < deltaT);
        %y(t - t(ii) < deltaT)
        %t(ii)
        %t - t(ii) < deltaT
        movingMean(ii) = mean(dataInRange);
    end
end