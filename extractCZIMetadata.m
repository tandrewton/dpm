% download bio-formats for matlab, add to matlab path
data = bfopen('Embryo_data/250118/Part4.czi');
%% extract table values
omeMeta = data{1,4}; % read data as standardized OME format
planeCount = omeMeta.getPlaneCount(0);
sizeT = omeMeta.getPixelsSizeT(0).getValue(); % extract # time values
listOfTimes = zeros(sizeT,1); 
for ii=1:sizeT
    deltaT = omeMeta.getPlaneDeltaT(0,ii-1); % extract time value and store
    listOfTimes(ii) = deltaT.value();
end
listOfTimes = unique(listOfTimes);