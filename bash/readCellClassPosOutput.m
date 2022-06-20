function [cellTrajectoryData, cell_count] = readDPMClassPosOutput(fstr)
%% FUNCTION to read in cell trajectory data from cell sims given file string
% -- NOTE: ACCOUNTS FOR ARBITRARY (AND CHANGING) NUMBER OF CELLS IN EACH FRAME

% print info to console
finfo = dir(fstr);
fprintf('Reading in %s\n',finfo.name);
fprintf('File size = %f MB\n',finfo.bytes/1e6)

% open file stream
fid = fopen(fstr);

% read in sim details from first frame
NCELLS      = textscan(fid,'NUMCL %f',1,'HeaderLines',1);   NCELLS = NCELLS{1};
phi0        = textscan(fid,'PACKF %f',1);                   phi0 = phi0{1};
fline       = fgetl(fid);
Ltmp        = textscan(fid,'BOXSZ %f %f',1);
fline       = fgetl(fid);
stresstmp   = textscan(fid,'STRSS %f %f %f %f',1);
fline       = fgetl(fid);
timetmp        = textscan(fid,'TIME %f',1); 
fline       = fgetl(fid);
S = fileread(fstr);
t = regexp(S, 'NUMCL(\s+)(\d+)', 'tokens');
%cell_count(i) counts the number of cells in frame i
cell_count = str2num(char(cellfun(@(x) x(1,2), t)));


% cells to save 
NFRAMES = 1e6;
nv      = zeros(NFRAMES,NCELLS);
vrad    = cell(NFRAMES,NCELLS);
xpos    = cell(NFRAMES,NCELLS);
ypos    = cell(NFRAMES,NCELLS);
zc      = zeros(NFRAMES,NCELLS);
zv      = zeros(NFRAMES,NCELLS);
calA0   = zeros(NFRAMES,NCELLS);
a0      = zeros(NFRAMES,NCELLS);
l0      = zeros(NFRAMES,NCELLS);
area   = zeros(NFRAMES,NCELLS);
perimeter = zeros(NFRAMES, NCELLS);
time    = zeros(NFRAMES,1);
phi     = zeros(NFRAMES, 1);
cellStressXX = zeros(NFRAMES, 1);
cellStressYY = zeros(NFRAMES, 1);
cellStressXY = zeros(NFRAMES, 1);
cellStressTrace = zeros(NFRAMES, 1);
cellShapeStressXX = zeros(NFRAMES, 1);
cellShapeStressYY = zeros(NFRAMES, 1);
cellShapeStressXY = zeros(NFRAMES, 1);
cellU = zeros(NFRAMES, 1);
L       = zeros(NFRAMES,2);
stress  = zeros(NFRAMES, 3);

% number of frames found
nf = 1;

% loop over frames, read in data
while ~feof(fid)
    % get box length
    L(nf,1) = Ltmp{1};
    L(nf,2) = Ltmp{2};
    NCELLS = cell_count(nf);
    % get info about deformable particle
    for nn = 1:NCELLS
        % get cell pos and asphericity
        cInfoTmp        = textscan(fid,'CINFO %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f',1); 
        fline           = fgetl(fid);     % goes to next line in file
        NVTMP           = cInfoTmp{1};
        nv(nf,nn)       = NVTMP;
        zc(nf,nn)       = cInfoTmp{2};
        zv(nf,nn)       = cInfoTmp{3};
        a0(nf,nn)       = cInfoTmp{4};
        area(nf,nn)    = cInfoTmp{5};
        perimeter(nf,nn) = cInfoTmp{6};
        psi(nf,nn)      = cInfoTmp{7};
        cellStressXX(nf,nn)    = cInfoTmp{8};
        cellStressYY(nf,nn)    = cInfoTmp{9};
        cellStressXY(nf,nn)    = cInfoTmp{10};
        cellStressTrace(nf,nn) = cInfoTmp{11};
        cellShapeStressXX(nf,nn)    = cInfoTmp{12};
        cellShapeStressYY(nf,nn)    = cInfoTmp{13};
        cellShapeStressXY(nf,nn)    = cInfoTmp{14};
        cellU(nf,nn)                = cInfoTmp{15};
        % compute l0 from calA0 and a0       
        calA0(nf,nn)    = perimeter(nf,nn)^2/(4*pi*area(nf,nn));
        %l0(nf,nn)       = sqrt(4.0*pi*calA0(nf,nn)*a0(nf,nn))/nv(nf,nn);
        
        % get vertex positions
        vPosTmp = textscan(fid,'VINFO %*f %*f %f %f %f %f %f %f %f %f %f %f %f %f',NVTMP); 
        % note : %*f means ignore that entry.
        fline = fgetl(fid);     % goes to next line in file

        % parse data      
        xposTmp = vPosTmp{1};
        yposTmp = vPosTmp{2};
        giTmp = vPosTmp{3};
        vradTmp = vPosTmp{4};
        l0radTmp = vPosTmp{5};
        %also 6 which is going to be ignored by me for now.
        vStressXXTmp = vPosTmp{7};
        vStressYYTmp = vPosTmp{8};
        vStressXYTmp = vPosTmp{9};
        vShapeStressXXTmp = vPosTmp{10};
        vShapeStressYYTmp = vPosTmp{11};
        vShapeStressXYTmp = vPosTmp{12};
        
        % save in cell
        xpos{nf,nn} = xposTmp;
        ypos{nf,nn} = yposTmp;
        gi{nf,nn} = giTmp;
        vrad{nf,nn} = vradTmp;
        l0rad{nf,nn} = l0radTmp;
        vStressXX{nf,nn} = vStressXXTmp;
        vStressYY{nf,nn} = vStressYYTmp;
        vStressXY{nf,nn} = vStressXYTmp;
        vShapeStressXX{nf,nn} = vShapeStressXXTmp;
        vShapeStressYY{nf,nn} = vShapeStressYYTmp;
        vShapeStressXY{nf,nn} = vShapeStressXYTmp;
    end
    % increment frame count
    nf = nf + 1;
    
    % read in/discard trash information
    % NOTE: IF HOPPER GEOMETRY/NUMBER OF PARTICLES/ANYTHING ELSE CHANGES
    % DURING SIMULATION, THIS IS WHERE YOU WOULD CHECK AND MAKE ADJUSTMENTS
    
    % get ENDFR string, check that read is correct
    fline = fgetl(fid);
    endFrameStr = sscanf(fline,'%s');
    
    if ~strcmp(endFrameStr,'ENDFR')
        error('Miscounted number of particles in .pos file, ending data read in');
    end
    
    % start new frame information
    if ~feof(fid)
        % get NEWFR line
        fline       = fgetl(fid);
        newFrameStr = sscanf(fline,'%s %*s');
        if ~strcmp(newFrameStr,'NEWFR')
            disp(newFrameStr)
            disp(nf)
            disp(fline)
            error('NEWFR not encountered when expected when heading to new frame...check line counting. ending.');
        end
        
        % read in sim details from first frame
        NCELLStmp       = textscan(fid,'NUMCL %f');
        
        % read in packing fraction
        phi0            = textscan(fid,'PACKF %f',1);                   
        %phi(nf)         = phi0{1};
        phi0            = phi0{1};
        fline = fgetl(fid);
        
        % update box size
        Ltmp            = textscan(fid,'BOXSZ %f %f',1);
        L(nf,1)         = Ltmp{1};
        L(nf,2)         = Ltmp{2};
        fline = fgetl(fid);
        
        % update stress
        stresstmp       = textscan(fid,'STRSS %f %f %f',1);
        stress(nf,1)         = stresstmp{1};
        stress(nf,2)         = stresstmp{2};
        stress(nf,3)         = stresstmp{3};
        fline = fgetl(fid);
        
        % update time
        timetmp         = textscan(fid,'TIME %f',1);
        time(nf,1)      = timetmp{1};
        fline = fgetl(fid);
    end
end

% delete extra information
if (nf < NFRAMES)
    NFRAMES = nf-1;
    vrad(nf:end,:) = [];
    xpos(nf:end,:) = [];
    ypos(nf:end,:) = [];
    gi(nf:end,:) = [];
    l0(nf:end,:) = [];
    vStressXX(nf:end,:) = [];
    vStressYY(nf:end,:) = [];
    vStressXY(nf:end,:) = [];
    vShapeStressXX(nf:end,:) = [];
    vShapeStressYY(nf:end,:) = [];
    vShapeStressXY(nf:end,:) = [];
    nv(nf:end,:) = [];
    zc(nf:end,:) = [];
    zv(nf:end,:) = [];
    a0(nf:end,:) = [];
    calA0(nf:end,:) = [];
    L(nf:end,:) = [];
    psi(nf:end,:) = [];
    time(nf:end,:)   = [];
    %phi(nf:end) = [];
    area(nf:end,:) = [];
    perimeter(nf:end,:) = [];
    stress(nf:end,:) = [];
    cellStressXX(nf:end,:)    = [];
    cellStressYY(nf:end,:)    = [];
    cellStressXY(nf:end,:)    = [];
    cellStressTrace(nf:end,:) = [];
    cellShapeStressXX(nf:end,:)    = [];
    cellShapeStressYY(nf:end,:)    = [];
    cellShapeStressXY(nf:end,:)    = [];
    cellU(nf:end,:)                = [];
end

% close position file
fclose(fid);

% store cell pos data into struct
cellTrajectoryData              = struct('NFRAMES',NFRAMES,'NCELLS',NCELLS);
cellTrajectoryData.phi0         = phi0;
cellTrajectoryData.vrad         = vrad;
cellTrajectoryData.xpos         = xpos;
cellTrajectoryData.ypos         = ypos;
cellTrajectoryData.vStressXX    = vStressXX;
cellTrajectoryData.vStressYY    = vStressYY;
cellTrajectoryData.vStressXY    = vStressXY;
cellTrajectoryData.vShapeStressXX    = vShapeStressXX;
cellTrajectoryData.vShapeStressYY    = vShapeStressYY;
cellTrajectoryData.vShapeStressXY    = vShapeStressXY;
cellTrajectoryData.gi           = gi;
cellTrajectoryData.nv           = nv;
cellTrajectoryData.zc           = zc;
cellTrajectoryData.zv           = zv;
cellTrajectoryData.calA0        = calA0;
cellTrajectoryData.a0           = a0;
cellTrajectoryData.l0           = l0rad;
cellTrajectoryData.L            = L;
cellTrajectoryData.stress       = stress;
cellTrajectoryData.area         = area;
cellTrajectoryData.perimeter    = perimeter;
cellTrajectoryData.psi          = psi;
cellTrajectoryData.time         = time;
cellTrajectoryData.cellStressXX    = cellStressXX;
cellTrajectoryData.cellStressYY    = cellStressYY;
cellTrajectoryData.cellStressXY    = cellStressXY;
cellTrajectoryData.cellStressTrace = cellStressTrace;
cellTrajectoryData.cellShapeStressXX    = cellShapeStressXX;
cellTrajectoryData.cellShapeStressYY    = cellShapeStressYY;
cellTrajectoryData.cellShapeStressXY    = cellShapeStressXY;
cellTrajectoryData.cellU                = cellU;

end