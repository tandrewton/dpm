addpath('/Users/AndrewTon/Documents/YalePhD/projects/dpm/bash')
addpath('C:\Users\atata\projects\dpm\bash')
addpath('/Users/AndrewTon/Documents/YalePhD/projects/dpm/matlab_funcs')
addpath('C:\Users\atata\projects\dpm\matlab_funcs')

% first identify all the strings I'll need to run over

% simulation parameters go here

% I need to identify strings that might look like the following:
% ablate_A01.10_k_a1.0_w_ps0.001_dsq0.0_k_ps4.0_k_lp4.0_t_lp1.0_d_flag3.0_bd1_sm1...
% /ablate_A01.10_k_a1.0_w_ps0.001_dsq0.0_k_ps4.0_k_lp4.0_t_lp1.0_d_flag3.0_bd1_sm1...
% _N60_Dur400_att0.2_sd1_sd1_sd1.pos

runType = "ablate";
N="60";
%ndelete="6";
calA0="1.10";
strainRate_ps="0.001";
%deltaSq = "2.0";
k_a = "1.0";
k_ps = "4.0"; %purse-string spring constant
k_lp = "4.0"; %lamellipodia spring constant
%smooth = "1";
tau_lp = "1.0"; %lamellipodia lifetime
%d_flag = "0.0"; %lamellipodia max length
%boundaryType = "0"; 
%att="0.2";
B="1.0";
Dr0="0.5";
boolCIL="0";
Duration="400";

numSeeds = 4;
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

% set up plotting windows
numPlots = 2; % area, shape
numPlotTypes = 4; % 
for i=1:numPlots
    for j=1:numPlotTypes
        figure(j+numPlotTypes*(i-1)); clf; hold on;    
        %bd 0 sm 0 = 1
        %bd 1 sm 0 = 2
        %bd 0 sm 1 = 3
        %bd 1 sm 1 = 4
    end
end

%figure(998); clf; hold on; % for debugging
%figure(999); clf; hold on; % for debugging

sm_arr = ["0" "1"];
att_arr = [0.1 0.2];
boundary_array = ["0" "1"];
color_array = ["red" "blue" "black"];
style_array = [":" "-"];
% for now, color is activity and style is attraction.

for i=1:length(boundary_array)
    boundaryType = boundary_array(i);
    for j=1:length(att_arr)
        att = att_arr(j);
        for k=1:length(sm_arr)
            smooth = sm_arr(k);
            for l=1:3
                if (l == 1) % C
                    deltaSq = "0.0";
                    d_flag = "3.0";
                elseif (l == 2) % P
                    deltaSq = "2.0";
                    d_flag = "0.0";
                elseif (l == 3) % CP
                    deltaSq = "2.0";
                    d_flag = "3.0";
                end
                % might also be looping over m=numSeeds to accumulate some results
                voidArea = zeros(0,2);
                meanInnerShapes = zeros(0,1);
                timeInnerShapes = zeros(0,1);
                timestep = 0; % determine on the fly to pad with zeros
                for m=1:numSeeds
                    % construct filenames to find the right simulation
                    bd = boundary_array(i);
                    att = att_arr(j);
                    sm = sm_arr(k);
                    seed = m;
                    run_name =runType+"_A0"+calA0+"_k_a"+k_a+"_w_ps"+strainRate_ps+ ...
                        "_dsq"+deltaSq+"_k_ps"+k_ps+"_k_lp"+k_lp+...
                        "_t_lp"+tau_lp+"_d_flag"+d_flag+"_bd"+boundaryType+"_sm"+smooth;
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
                    if (length(voidArea) < length(voidArea_sd))
                        voidArea(length(voidArea_sd),:) = 0;
                        meanInnerShapes(length(meanInnerShapes_sd),:) = 0;
                        voidArea(:,1) = voidArea_sd(:,1); %extend time column
                        timeInnerShapes = timeInnerShapes_sd; %extend time column
                    elseif (length(voidArea_sd) < length(voidArea))
                        voidArea_sd(length(voidArea),:) = 0;
                        meanInnerShapes_sd(length(meanInnerShapes),:) = 0;
                        voidArea_sd(:,1) = voidArea(:,1); %extend time column
                    end
                    % otherwise they have exactly the same length, so don't
                    % adjust lengths

                    voidArea(:,2) = voidArea(:,2) + voidArea_sd(:,2);
                    meanInnerShapes = meanInnerShapes + meanInnerShapes_sd;
                end

                voidArea = voidArea / numSeeds;
                meanInnerShapes = meanInnerShapes / numSeeds;

                % plot area vs time
                if (boundaryType == "0" && smooth == "0")
                    figure(1)
                elseif (boundaryType == "1" && smooth == "0")
                    figure(2)
                elseif (boundaryType == "0" && smooth == "1")
                    figure(3)
                elseif (boundaryType == "1" && smooth == "1")
                    figure(4)
                end

                
                plot(voidArea(:,1), voidArea(:,2), 'linewidth', 4,...
                    'Color', color_array(l),...
                    'LineStyle', style_array(j))
                xlabel('$t/\tau$','Interpreter','latex','fontsize', 24);
                ylabel('Area','Interpreter','latex','fontsize', 24);
                %set(gca,'Yscale','log')

                % plot shape vs time
                if (boundaryType == "0" && smooth == "0")
                    figure(1+numPlotTypes)
                elseif (boundaryType == "1" && smooth == "0")
                    figure(2+numPlotTypes)
                elseif (boundaryType == "0" && smooth == "1")
                    figure(3+numPlotTypes)
                elseif (boundaryType == "1" && smooth == "1")
                    figure(4+numPlotTypes)
                end
               
                % cellID row = [ci inInitialWoundNeighbors inFinalWoundNeighbors]
                % we want to access inFinalWoundNeighbors of bulkCellShape
                % bulkCellShape row = [time shape(0) shape(1) ... shape(NCELLS)]
                %innerShapes = bulkCellShape(:,[1; cellID(:,3)]==1);
                %outerShapes = bulkCellShape(:,[1; ~cellID(:,3)]==1);
                %plot(innerShapes(:,1), nanmean(innerShapes(:,2:end),2),  ...
                %    'linewidth', 4, 'DisplayName', "inner shapes",...
                %     'Color', color_array(l),'LineStyle', style_array(j))
                %plot(timeAndOuterShapes(:,1),nanmean(timeAndOuterShapes(:,2:end),2), ...
                %    'linewidth', 4, 'DisplayName', "bulk shapes")
                plot(timeInnerShapes, meanInnerShapes, 'linewidth', 4, ...
                    'DisplayName', "inner shapes", 'Color', ...
                    color_array(l),'LineStyle', style_array(j))
                xlabel('$t/\tau$','Interpreter','latex','fontsize', 24);
                ylabel('Shape','Interpreter','latex','fontsize', 24);
                %legend('location','northwest','fontsize', 14)
            end
        end
    end
end

%cd(output_dir)
% save files here
