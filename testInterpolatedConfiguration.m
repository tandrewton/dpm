isTestData = true; %uncomment if using test data
addpath('/Users/AndrewTon/Documents/YalePhD/projects/dpm/bash')
addpath('/Users/AndrewTon/Documents/YalePhD/projects/dpm/matlab_funcs')

%CHANGE THESE PARAMETERS                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   NEEDED

runType = "ablate";
%N="40";
ndelete="10";
calA0="1.10";
%strainRate_ps="0.01";
%deltaSq = "2.0";
k_ps = "1.0"; %purse-string spring constant
k_lp = "2.0"; %lamellipodia spring constant
%tau_lp = "1.0"; %lamellipodia lifetime
%d_flag = "0.0"; %lamellipodia max length
prate = "0.00"; %perimeter relaxation rate
%att="0.2";
B="1.0";
Dr0="0.5";
boolCIL="0";
Duration="800";
FSKIP = 1;

etaStr = " ";
startSeed = 1;
max_seed = 1;
makeAMovie = 1; %if makeAMovie is 0, then plot every frame separately
set(0,'DefaultFigureWindowStyle','docked')
showPeriodicImages = 0;

showverts = 0;
showBoundaries = 0;
showArea = 1;
showQuiver = 0;
walls = 0;
%disable showVoid if using printConfig on its own, outside of
%dampedNVE/dampedNP0 routines
showGlobalIndex = 0;
showVoid = 0;
showVoidBlack = 0; % print void in larger black circles to see easier
showCornersOrEdges = 0;
showPurseString = 1;
showShapeHistogram = 0;
 
%PC directory
pc_dir = "/Users/AndrewTon/Documents/YalePhD/projects/dpm/";
%pipeline is the location of data generated during simulations
subdir_pipeline = pc_dir + "pipeline/cells/"+runType+"/";

%output is location of results of this postprocessing
subdir_output = pc_dir + "output/cells/"+runType+"/";
mkdir(subdir_pipeline);
mkdir(subdir_output);


%txt = 'N = '+N+', NV = '+NV+', calA_o='+calA+', att='+att+', B='+B;
txt='test';

fnum = 1;
figure(13), clf, hold on, box on;
for seed = startSeed:max_seed
    if (isTestData)
        run_name = runType+txt;     
        pipeline_dir =  subdir_pipeline + run_name + "/";
        output_dir = subdir_output + run_name + "/";
        mkdir(pipeline_dir)
        mkdir(output_dir)
        fileheader=run_name + "_seed" + seed;
        nvestr = pc_dir+'test.pos';
        energystr = pc_dir+'test.energy';
        stressstr = pc_dir+'test.stress';
        boundaryStr = pc_dir+'test.void';
        edgeStr = pc_dir+'test.edge';
        purseStr = pc_dir+'test.purseString';
        voidAreaStr = pc_dir+'test.voidArea';
    else
        run_name =runType+"_calA0"+calA0+"_strainRate_ps"+strainRate_ps+ ...
            "_deltaSq"+deltaSq+"_k_ps"+k_ps+"_k_lp"+k_lp+...
            "_tau_lp"+tau_lp+"_d_flag"+d_flag+"_prate"+prate;
        pipeline_dir =  subdir_pipeline + run_name + "/";
        output_dir = subdir_output + run_name + "/";
        mkdir(pipeline_dir)
        mkdir(output_dir)
        fileheader=run_name +"_NCELLS"+N+"_Duration"+Duration+"_att"+att+"_startseed"+ ...
            startSeed+"_endseed"+max_seed+"_seed"+seed;
        nvestr = pipeline_dir+fileheader+'.pos';
        energystr = pipeline_dir+fileheader+'.energy';
        stressstr = pipeline_dir+fileheader+'.stress';
        boundaryStr = pipeline_dir+fileheader+ ".void";
        edgeStr = pipeline_dir+fileheader+ '.edge';
        purseStr = pipeline_dir+fileheader+ '.purseString';
        voidAreaStr = pipeline_dir+fileheader+ '.voidArea';
    end