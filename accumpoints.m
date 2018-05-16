function [finalArray, sX, sY] = accumpoints(data, xRange, yRange, kernel, step, fun, varargin)
% ACCUMPOINTS Construct an array by accumulation the value of discrete
% points on a grid. ACCUMPOINTS allows to apply a custom function and a
% sliding window
%
% [A, xCoor, yCoor] = ACCUMPOINTS(data, xRange, yRange, kernel, step, fun, varargin)
%
% Required inputs:
% - data:         m x 3 matrix where colums 1-2 are x-y coordinates and column 3 is the variable
% - xRange:       1 x 2 vector containing the min-max extent along x
% - yRange:       1 x 2 vector containing the min-max extent along y
% - kernel:       Size of the analysis window
% - step:         Sliding step of the window
% - fun:          Anonymous function to apply
%
% Output:
% - A:            Accumulated array
% - xCoor:        Coordinates of A along x
% - yCoor:        Coordinates of A along y
%
% Optional name-value pairs of inputs:
% To output geotif (needs mapping toolbox)
% - 'outName':    Path/name to output geotiff
% - 'CoordRef'    EPSG coordinate reference system code for the coordinates of the data (positive integer). See CoordRefSysCode
%
% Note: The function uses parfor. To disable it, just change the parfor
%       occurrence to a regular for.
%
% See also accumarray, blockproc

% Test if kernel is odd
if mod(kernel,2) == 0 || kernel == 1
    error('The kernel size should be odd and greater than 1')
end

% Test the number of subarrays
if floor(kernel/step) ~= kernel/step
   error('The ratio kernel/step should be an integer');
end

% Test extent
if min(data(:,1))<xRange(1) || max(data(:,1))>xRange(2) || min(data(:,2))<yRange(1) || max(data(:,2))>yRange(2)
    warning('Some points are outside of the specified range');
end

% Test optional input arguments
for i=1:length(varargin)
    if      strcmp(varargin{i}, 'outName'), outName = varargin{i+1};
    elseif  strcmp(varargin{i}, 'CoordRef'), CoordRef = varargin{i+1};   
    end
end

nbArray = kernel/step;   % Number of subarrays

% Define final matrix
border  = (kernel-1)/2;  % Size of the border that will be lost
sX      = (xRange(1)+border):step:(xRange(2)-border);   % x coordinates
sY      = (yRange(1)+border):step:(yRange(2)-border);   % y coordinates
finalArray = zeros(numel(sY), numel(sX));               % Final matrix

% Temporary storage
tempStorArray = cell(nbArray,nbArray);
tempStorIdx   = cell(nbArray,nbArray);

%% Main processing loop
parfor iY = 1:nbArray	% Number of arrays along Y
    for iX = 1:nbArray	% Number of arrays along X     
        xVec = xRange(1)+step*(iX-1):kernel:xRange(2);  % x bins of sub array
        yVec = yRange(1)+step*(iY-1):kernel:yRange(2);  % y bins of sub array

        % Number of bins
        nxBins = length(xVec);
        nyBins = length(yVec);

        % Bin the data (check speed for large datasets)
        [~,~,ii] = histcounts(data(:,1), xVec);
        [~,~,jj] = histcounts(data(:,2), yVec);
        % Remove empty bins and points that fall outside of the bin of the sub array
        % Note: bins a defined as min >= x > max
        idx = (ii>0 & jj>0)...
            & (data(:,1)<xVec(end) & data(:,2)<yVec(end))...
            & (data(:,1)>=xVec(1) & data(:,2)>=yVec(1));
        
        % Accumulates sub array
        tempArray = accumarray([jj(idx) ii(idx)], data(idx,3), [nyBins-1 nxBins-1], fun, NaN);
        
        % Indices of sub array on the final matrix
        xIdx = repmat(xVec(1:end-1)-xRange(1)+step, size(tempArray,1), 1) ;
        yIdx = repmat((yVec(1:end-1)-yRange(1)+step)', 1, size(tempArray,2));
        
        % Fill the temporary storage
        tempStorArray{iY,iX} = tempArray;    % Temporary array
        tempStorIdx{iY,iX}   = sub2ind([numel(sY) numel(sX)], yIdx, xIdx);   % Indices of temp array on final matrix
    end
end

%% Post processing
% Fill final array outside of the main loop because of the use of parfor
for iY = 1:nbArray	% Number of arrays along Y
    for iX = 1:nbArray	% Number of arrays along X
        % Fill the final matrix
        finalArray(tempStorIdx{iY,iX}) = tempStorArray{iY,iX};
    end
end

% Write to geotiff
if exist('outName', 'var') && exist('CoordRef', 'var')  
    if license('test', 'map_toolbox')
        R = maprefpostings([sX(1) sX(end)], [sY(1) sY(end)], step, step);
        geotiffwrite(outName, finalArray, R, 'CoordRefSysCode', CoordRef);
    else
        warning('The mapping toolbox is not available, no geotiff was written')
    end
end
