function [wPt,iPt] = caltag(I,datafile,debug)
%CALTAG Detect corner points in image of self-identifying chequerboard
%	[G,P] = CALTAG(I,datafile) takes a greyscale image I and filename "data" of a .mat
%	file describing
%	the parameters of a self-identifying chequerboard pattern. The function
%	"caltag_generate" can produce an appropriate data file
%
% Note that saddle_finder will crash is you attempt to find a saddle in a region
% of constant intensity. So you can't try to mask out saddles by painting blobs
% of solid colour over the corners.
%
% Note: it will run much slower if the input image dimensions are prime, because
% the adaptivethresh function does an FFT across each row and column. For best
% speed, dimensions with small prime factors are best.
%
% Conventions: pixel centres lie at integer coordinates. Top left pixel is (1,1)
% unless you choose the the C-style output in which case it's (0,0)
%
% return G = 2xN matrix of grid points, P = 2xN matrix of image points
% The export_caltag.m function will write out the points to file that can then
% be calibrated with OpenCV using the wrapper at https://github.com/brada/PyCVcalib
%
%	Class Support
%	-------------
%	I must be 2D numeric, real, nonsparse.
%	G and P are 2xN real matrices of class double.
%
%	For algorithm details, see
%	Atcheson, B., Heide, F., Heidrich, W. "CALTag: High Precision Fiducial
%	Markers for Camera Calibration", VMV 2010
%   http://www.cs.ubc.ca/labs/imager/tr/2010/Atcheson_VMV2010_CALTag/
%
%   Licence:
%   CALTag is free to use/modify/distribute for noncommercial use.
%   I'd appreciate an email to let me know where it's being used, but
%   that's not required. Contact atcheson at cs dot ubc dot ca for
%   commercial licencing. There is some 3rd party code from Peter Kovesi
%   and Jean-Yves Bouguet that is subject to its own licence.



%% sanity checks
if nargin ~= 3
	error( 'Three inputs required' );
end
if ischar( I )
    if ~exist( I, 'file' )
        error( 'File "%s" does not exist', I );
    else
        I = imread( I );
    end
end    
if ~ismember( ndims(I), [2,3] )
	error( 'I must be two- or three-dimensional' );
end
if ~isnumeric( I )
	error( 'I must be numeric' );
end
if ~isreal( I )
	error( 'I must be real' );
end
if issparse( I )
	error( 'I cannot be sparse' );
end
if ~exist( datafile, 'file' )
	error( 'Cannot find .mat file: %s', datafile );
end
% rather put the deploy logic into main.m
%if isdeployed
%    if ~ismember( debug, {'true','false'} )
%        error( 'debug must be "true" or "false"' );
%    end
%    debug = str2num( debug );
%else
    if ~islogical( debug )
        error ('debug must be boolean' );
    end
%end
if ~exist( 'cornerfinder_saddle_point', 'file' )
    error( 'Please first install the calibration toolbox: http://www.vision.caltech.edu/bouguetj/calib_doc/' );
end
if ~exist( 'hnormalise.m', 'file' ) ...
   || ~exist( 'homography2d.m', 'file' ) ...
   || ~exist( 'homoTrans.m', 'file' ) ...
   || ~exist( 'iscolinear.m', 'file' ) ...
   || ~exist( 'normalise2dpts.m', 'file' ) ...
   || ~exist( 'ransac.m', 'file' ) ...
   || ~exist( 'ransacfithomography.m', 'file' )
    disp( 'homography2d, homoTrans, iscolinear, normalise2dpts, ransac, ransacfithomography' );
    error( 'Please download the above files from http://www.csse.uwa.edu.au/~pk/Research/MatlabFns/' );
end
   



%% default output values
wPt = [];
iPt = [];


%% load and check datafile
load( datafile );
expectedVars = { 'resPattern', 'resMarker', 'resCode', 'idBits', ...
                 'layout', 'ID', 'CODE', 'scale', ...
                 'minHamming' };
if ~all( ismember(expectedVars,who) )
	error( 'Datafile does not contain required variables' );
end
isInt = @(x) ~isempty(x) && isequal( x, uint32(x) );
isValid = @(x) numel(x)==2 && isInt(x) && min(x)>0 && max(x)<32;
if ~all( cellfun(isValid,{resPattern,resMarker,resCode}) )
	error( 'Invalid resolution specified in datafile' );
end
if any( resCode > resMarker )
	error( 'Code resolution must be less than marker resolution' );
end
if ~isInt( idBits ) || ~isscalar( idBits )
	error( 'idBits must be a numeric scalar' );
end
if idBits > prod( resCode )
	error( 'idBits must be less than code' );
end
idBits = uint32( idBits );
if ~ismember( layout, [1,2] )
    error( 'Invalid layout' );
end
if ~all( cellfun(isInt,{ID,CODE}) )
    error( 'ID and CODE must be non-empty integers' );
end


% ugly hack to flip res to the order expected in code below. would have to go
% through entire program carefully to make this all consistent
resPattern = resPattern([2,1]);


%% derived parameters
nMarkers = prod( resPattern );
%[nPatternRows,nPatternCols] = splitvec( resPattern );
% crcBits needs to be double or else rightshift breaks
crcBits = double( prod( resCode ) - idBits );
% at minimum we expect 2x2 pixels per code dot
minMarkerArea  = prod( resMarker ) * 4;
minDotEdgePixels = 2 * 4;
minCodeArea    = prod( resCode ) * 4;
maxRegionArea  = size(I,1)*size(I,2) / 8;
minRegionArea  = minMarkerArea - minCodeArea;
minEulerNumber = -8;
maxEulerNumber = 0;


% mode: 1=old algorithm, 2=newer, better algorithm (not doing
% adaptive thresholding). This lets you easily try out new image processing
% operation sequences. 3 = just testing, probably doesn't work
mode = 2; % [brad 10 Dec 2010] change from default 2 to 3, for alex's dataset



%% initial image processing
if debug
    disp( 'Converting image...' );
end
if ndims( I ) == 3
	I = rgb2gray( I );
end
I = im2single( I );
[imgHeight,imgWidth] = size( I );
if max( [factor(imgHeight),factor(imgWidth)] ) > 21
	disp( 'Warning: image dimensions have large prime factors' );
	disp( 'adaptivethresh is fastest on image dimensions with small factors' );
end
normMatrix = [ 2/imgHeight, 0, -1; 0, 2/imgWidth, -1; 0 0 1 ];
if debug
    disp( 'Adaptive thresholding...' );
end
T = adaptivethresh( I );
% [brad] bugfix 10 Dec 2010, Speckle noise can cause isolated bright pixels, which
% cause lots of tiny regions, artificially decreasing the euler number. So we
% use a majority filter to eliminate this noise.
T = bwmorph( T, 'majority' );
% [philippe] bugfix 3 Dec 2012
T = single( T );
% /bugfix
if debug
    disp( 'Finding regions...' );
end

if layout == 1
    switch mode
        case 1
            kernel = fspecial( 'sobel' ) / 8;
            gx = imfilter( I, kernel', 'replicate' );
            gy = imfilter( I, kernel,  'replicate' );
            mag = sqrt( gx.^2 + gy.^2 ); 
            fsize = sqrt( minMarkerArea ) * min( resPattern ) / 2;
            E = adaptivethresh( mag, fsize, 5 );
        case 2
            % Note: sometimes the corners of the markers are broken by the
            % automatic default threshold selection Matlab makes. Manually
            % specifying a threshold here can fix things. Currently I don't know
            % how to choose this value automatically. Ideally we would use a LOG
            % edge filter here (or Canny) but they are slow and produce many
            % spurious pixels that are hard to filter out.
            %E = edge( I, 'sobel', 0.03, 'nothinning' );
            E = edge( I, 'sobel', 'nothinning' );
            E = bwmorph( E, 'bridge' );
            E = bwmorph( E, 'majority' );
        case 3
            %E = edge( I, 'log' );  
            % 10 Dec 2010, to get Alex's low contrast, blurry dataset to
            % work:
            E = bwmorph( T, 'thicken' ) - T;
            E = bwmorph( E, 'bridge' );
        otherwise
            error( 'Unsupported image processing mode' );
    end
    E  = bwmorph( E, 'thin', inf );
    E  = bwmorph( E, 'clean' );
    % remove very small isolated lines (noise that could be falsely treated
    % as holes inside the markers)
    CC = bwconncomp( E, 8 );
    nPixels = cellfun( @(x) numel(x), CC.PixelIdxList );
    bad = find( nPixels < minDotEdgePixels );
    CC = filtercc( CC, bad );
    for i = 1:CC.NumObjects
        E(CC.PixelIdxList{i}) = 0;
    end
    E = bwmorph( E, 'close' );
elseif layout == 2
    % here: T == adaptivethresh(I)
    E = T;
    E = imclose( E, strel('disk',2) );   
end

if debug
    disp( 'Filtering regions...' );
end
% [brad] bugfix 18 Oct 2010
% changed "CC = bwconncomp( ~E, 4 );" to 
% "CC = bwconncomp( bwmorph(~E,'erode'), 4 );" because when the edges are fuzzy
% it turns out that the thinned edge region created above has many spurs. These
% spurs are only 8-connected, so when we compute regions here with
% 4-connectivity, many pixels in those spurs count as independent "holes" with
% regard to the Euler number. We don't want that - rather a spur connecting
% to a full region should just be single hole. Note that we can't just switch to
% 8-connectivity for regions here, because then neighbouring regions with a
% single-pixel line separating them would be considered connected.
%CC = bwconncomp( ~E, 4 );
CC = bwconncomp( bwmorph(~E,'erode'), 4 );
%/bugfix
R = regionprops( CC, 'Area' );
good = find( ([R.Area]>minRegionArea) & ...
             ([R.Area]<maxRegionArea) );
if isempty( good )
    wPt = [];
    iPt = [];
    disp( 'No regions with valid areas were detected' );
    return
end
CC = filtercc( CC, good );
R = regionprops( CC, 'EulerNumber' );
good = find( ([R.EulerNumber]>=minEulerNumber) & ...
             ([R.EulerNumber]<=maxEulerNumber) );
if isempty( good )
    wPt = [];
    iPt = [];
    disp( 'No regions with valid Euler numbers were detected' );
    return
end
CC = filtercc( CC, good );


% todo: need more aggressive filtering on regions. gravel/carpet scenes
% contain tons of little blobs that take ages to run quadfit on, so we need
% to add some sort of convexity or solidity or eccentricity threshold too?
% although it's not mentioned in the paper...


% debug
if debug
    hDebug1 = figure();
    imshow( label2rgb(labelmatrix(CC)) );
    impixelinfo;
end

% early exit if not enough potential regions
if CC.NumObjects < prod(resPattern)/8
    wPt = [];
    iPt = [];
    disp( 'Too few regions detected' );
    return
end


% now loop over all regions, get their outline pixels, fit quad...
R = regionprops( CC, 'BoundingBox', 'FilledImage' );
warning off stats:kmeans:EmptyCluster
warning off stats:kmeans:EmptyClusterAllReps
for i = 1:length( R )
    [isq,cnr,cnr0] = fitquad( R(i).BoundingBox, R(i).FilledImage, layout );
    R(i).isQuad = isq;
    if isq       
        R(i).corners = cnr;
        if layout == 2
            R(i).corners_0 = cnr0;
        end
    end
end
warning on stats:kmeans:EmptyCluster
warning on stats:kmeans:EmptyClusterAllReps

if nnz( [R.isQuad] ) == 0
    disp( 'Could not fit quads to any of the regions' );
    wPt = [];
    iPt = [];
    return;
end

% cluster the nearby corners together to get initial guesses for saddles
Rq = R([R.isQuad]);
Rqc = [Rq(:).corners];
cutThresh = sqrt( minMarkerArea );
% fudge factor for layout 2, since the rotated corners are too far away from the
% actual point
if layout == 2
    cutThresh = cutThresh * 4;
end
saddlesidx = clusterdata( Rqc', 'criterion','distance', ...
                          'cutoff',cutThresh );
nSaddles = max( saddlesidx );                      
saddles_0 = zeros( 2, nSaddles );
for i = 1:nSaddles
    saddles_0(:,i) = mean( Rqc(:,saddlesidx==i), 2 );
end
% find the saddles more precisely, first with large mask then with smaller
warning off MATLAB:nearlySingularMatrix
% note: using 5x5 then 3x3 works for very sharp images, and best when you expect
% small squares. but for out-of-focus images, with large squares, a larger
% window is needed for the saddlefinder. need to figure out if large window
% adversely affects accuracy.
[saddles] = cornerfinder_saddle_point( flipud(saddles_0), I, 11,11 ); %for sharp images: 5x5
[saddles,good] = cornerfinder_saddle_point( saddles, I, 9,9 ); %for sharp images: 3x3
warning on MATLAB:nearlySingularMatrix
saddles = flipud( saddles(:,good) );
nSaddles = nnz( good );
if nSaddles == 0
    wPt = [];
    iPt = [];
    disp( 'No saddles detected' );
    return
end
 

if debug
    hDebug2 = figure;
    isQuad = [R(:).isQuad];
    imshow( ismember( labelmatrix(CC), find(isQuad) ) );
    hold on
    corners = [R(:).corners];
    plot( corners(2,:), corners(1,:), 'r+' );
    plot( saddles_0(2,:), saddles_0(1,:), 'gO' );
    plot( saddles(2,:), saddles(1,:), 'g+' );
end

% find the closest saddles for each corner of the quad
% would have been better to keep track of the indices into Rqc, saddles_0
% and saddlesabove to get these directly without computing distances
dist2Saddles = @(p) sum( sqrt((saddles-repmat(p,[1,nSaddles])).^2) );
% define a unit square and sample points inside it (homogeneous coords)
unitSquare = [ 0 1 1 0; 0 0 1 1; 1 1 1 1 ];
resBorder = (resMarker-resCode)/2;
tl = 0.5./resMarker + resBorder./resMarker;
br = (resBorder+resCode)./resMarker - 0.5./resMarker;
Sy = linspace( tl(1), br(1), resCode(1) );
Sx = linspace( tl(2), br(2), resCode(2) );
[Gy,Gx] = ndgrid( Sy, Sx );
S = [ Gy(:)'; Gx(:)'; ones(1,numel(Gx)) ];
for i = 1:length(R)
    if R(i).isQuad
        [d1,idx1] = min( dist2Saddles(R(i).corners(:,1)) );
        [d2,idx2] = min( dist2Saddles(R(i).corners(:,2)) );
        [d3,idx3] = min( dist2Saddles(R(i).corners(:,3)) );
        [d4,idx4] = min( dist2Saddles(R(i).corners(:,4)) );
        if all( [d1,d2,d3,d4] < 10 ) ...
           && length( unique([idx1,idx2,idx3,idx4]) ) == 4
            R(i).saddles = saddles(:,[idx1,idx2,idx3,idx4]);
            % get homography mapping the unit square to this one
            quadSquare = [ R(i).saddles; 1 1 1 1 ];
            R(i).H = homography2d( unitSquare, quadSquare );
            R(i).HS = homoTrans( R(i).H, S );
            % rotate sampling grid if using the rotated layout (2)
            if layout == 2
                centroid = mean( R(i).HS, 2 );
                nPts = size( R(i).HS, 2 );
                centroid2pts = R(i).HS - repmat( centroid, [1,nPts] );
                tformMatrix = [0.5 -0.5 0; 0.5 0.5 0; 0 0 1];
                centroid2pts = tformMatrix * centroid2pts;
                R(i).HS = centroid2pts + repmat( centroid, [1,nPts] );
            end
            % sample the code from the image
            R(i).code = uint32( reshape( interp2(T,R(i).HS(2,:), ...
                                                   R(i).HS(1,:), ...
                                                   '*linear'), resCode ) );
        else
            R(i).isQuad = false;
        end
    end
end
R = R([R.isQuad]);
if isempty( R )
    wPt = [];
    iPt = [];
    disp( 'No quads detected' );
    return
end
        
if debug
    hDebug3 = figure;
    imshow( I );
    hold on;
    saddles = [R(:).saddles];
    corners = [R(:).corners];
    plot( saddles(2,:), saddles(1,:), 'r+','MarkerSize',10 );
    plot( saddles_0(2,:), saddles_0(1,:), 'rO', 'MarkerSize',5 );
    plot( corners(2,:), corners(1,:), 'r*','MarkerSize',5 );
    samples = [R(:).HS];
    plot( samples(2,:), samples(1,:), 'g+' );
    if layout == 2
        corners_0 = [R(:).corners_0];
        plot( corners_0(2,:), corners_0(1,:), 'rd', 'MarkerSize',5 );       % it appears that in layout2, the rotated points can end up 
                                                                            % quite far from the actual saddle, leading to them being
                                                                            % missed unless you get 4 nearby points that give a good
                                                                            % average, does the viewing angle affect this?
                                                                      
    end
end
    

% for converting from binary to decimal
pows = uint32( 2.^(idBits+crcBits-1:-1:0) );
b2d = @(x) sum( reshape(x',1,[]).*pows, 'native' );
isvalid=@(x)crc_validate(b2d(x),idBits,crcBits)&&ismember(b2d(x),CODE);
% check each possible orientation of the code to see if exactly one of
% them is valid, and if so record the angle to the origin corner of the
% marker from the centre
% the [y,x] coordinates in R(i).corners are in [row,column] order
idmask = leftshift( 2^idBits-1, crcBits );
for i = 1:length(R)
    c = R(i).code;
    validity = [ isvalid(c),           isvalid(rot90(c,-1)), ...
                 isvalid(rot90(c,-2)), isvalid(rot90(c,-3)) ];
    R(i).isValid = (nnz(validity) == 1);
    if R(i).isValid
        % the '1' in validity indicates the top left corner of the
        % marker, which is at (x=0,y=1) in marker coordinates. Topleft
        % is also (arbitrarily) the 4th point in our anticlockwise
        % ordering
        % shift corner/saddle array so origin point is first
        shift = find( validity ) - 1;
        R(i).code = rot90( R(i).code, -shift );
        data = b2d( R(i).code );
        R(i).id = rightshift( bitand(data,idmask), crcBits );
        R(i).crc = bitand( data, 2^crcBits-1 );
        R(i).saddles = circshift( R(i).saddles,[0,-shift-1] );
        R(i).corners = circshift( R(i).corners,[0,-shift-1] );   
        origin = R(i).saddles(:,1);
        second = R(i).saddles(:,2);
        vec = second - origin;
        R(i).orientation = cart2pol( vec(2), vec(1) );
        if debug
            pts = flipud( [ origin, origin+vec/2 ] );                
            line( pts(1,:), pts(2,:), 'LineWidth',3, 'Color','r' );
        end
    end
end
R = R([R.isValid]);
   

% early exit if no valid quads found
if isempty(R)
    wPt = [];
    iPt = [];
    disp( 'No valid quads detected' );
    return;
end


% Bugfix 6/June/2011 [Nassir] - solve problem whereby perfectly axis-aligned image
% could not be detected because some markers had orientation epsilon and
% others had orientation 2*pi-epsilon, so they were all very different
% from the median and ended up all being rejected as outliers. Solution
% is to clamp values near to the -pi/pi crossover point to one side
newOrientation = [R.orientation] + pi;
wrapAround = abs( newOrientation - 2*pi ) < 0.001;
for i = 1:length(R)
    if wrapAround(i)
        R(i).orientation = -pi;
    end
end


% filter out the regions with orientations more than 30 degrees different
% from the median orientation
% can't use mean because epsilon and 2pi-epsilon should be considered very
% close together but their mean is pi
[mox,moy] = pol2cart( median([R.orientation]), 1 );
medianOrientation = [mox;moy];
[ox,oy] = pol2cart( [R.orientation], 1 );
angleDiffs = real( acos([ox;oy]'*medianOrientation) ) * 180/pi; % dot products
R = R(angleDiffs<30);


% find the mappings from grid (world) to image coordinates for each corner
% of each marker
for i = 1:length(R)
    [found,loc] = ismember( R(i).id, ID );
    if found
        [row,col] = ind2sub( size(ID), loc );
        % flip to postscript origin@bottomleft coordinate system
        row = resPattern(2) - row;
        col = col - 1;
        if debug
            % change sample markers to indicate this is a valid marker
            samples = [R(i).HS];
            plot( samples(2,:), samples(1,:), 'g*' );
        end
        % record the point correspondences for this marker
        % [brad] [20 Dec 2010] bugfix: support scales other than 1.0 in
        % generated pattern (generate pattern script changed too)
        %R(i).wPt = [row col; row col+1; row+1 col+1; row+1 col]';
        R(i).wPt = [row col; row col+1; row+1 col+1; row+1 col]' * scale;
        R(i).iPt = R(i).saddles;            
    else
        R(i).isValid = false;
    end
end
R = R([R.isValid]);    
wPt = [R.wPt]';
iPt = [R.iPt]';
[wPt,idx] = unique( wPt, 'rows' );
iPt = iPt(idx,:);

if isempty( iPt ) || isempty( wPt )
    disp( 'No valid codes detected' );
    return;
end



% estimate the amount of radial distortion
obj_fun = @(r) distortion_metric(wPt,iPt,resPattern,normMatrix,r);
options = optimset( 'Display','off' );
radialDist = lsqnonlin( obj_fun, 0, [],[], options );
if debug
    disp( 'Estimated radial distortion coefficient:' );
    disp( radialDist );
end

% 1 May 2012:
% Chen Xing discovered bug here - radial distortion estimator (or undistorter)
% function is bad. For lens with moderate distortion, where caltag grid appears
% only on extreme edge of image, it gives poor estimate, which throws off rest
% of calibration and causes other weird errors to appear. Solution is to fix the
% radial distortion estimation/undistortion process, or for now, just disable it
if debug
    disp('FIXME: radialDist');
end
radialDist = 0.0;

% now look for any points we may have missed
% undistort detected points, then try to fit a homography
iPtUndist = radial_distort( iPt, normMatrix, 1-1/(1+radialDist) );    %-radialDist );
%if debug
%    plot( iPtUndist(:,2), iPtUndist(:,1), 'c*' );
%end
H = ransacfithomography( wPt', iPtUndist', 0.1 );
% NxM squares imply (N+1)x(M+1) grid points
resGrid = fliplr( resPattern + 1 );

% [brad] [20 Dec 2010] bugfix: support scales other than 1.0. Note that
% there could be a potential problem here if the scale is not a "nice"
% number, ie if both it and its reciprocal can be easily stored within
% numerical precision. If not, we may have to start rounding before calling
% ind2sub.
%idxFound = sub2ind( resGrid, wPt(:,1)+1, wPt(:,2)+1 );
%idxMissing = setdiff( 1:prod(resPattern+1), idxFound' );
%[mrow,mcol] = ind2sub( resGrid, idxMissing );
%mrow = mrow - 1;
%mcol = mcol - 1;
idxFound = sub2ind( resGrid, round(wPt(:,1)/scale)+1, round(wPt(:,2)/scale)+1 );
idxMissing = setdiff( 1:prod(resPattern+1), idxFound' );
[mrow,mcol] = ind2sub( resGrid, idxMissing );
mrow = (mrow-1) * scale;
mcol = (mcol-1) * scale;
% /end of bugfix


trialPoints = homoTrans( H, [mrow;mcol;ones(1,length(mrow))] ); %??? (rol,row)?
trialPoints = trialPoints(1:2,:);
% apply the radial distortion estimate
trialPoints = radial_distort( trialPoints', normMatrix, -radialDist )'; % radialDist+radialDist^2+radialDist^3+radialDist^4 )';
% exclude points that lie outside the image boundaries
% saddlefinder default windowsize is 5, we just pad a bit extra to be safe
borderPadding = 5 + 3;
bad = any( trialPoints < borderPadding ) | ...
      trialPoints(2,:) > size(I,2)-borderPadding | ...
      trialPoints(1,:) > size(I,1)-borderPadding;
trialPoints = trialPoints(:,~bad);
%[brad] from Anika's dataset: bugfix September 2010
mrow = mrow(~bad);
mcol = mcol(~bad);
%[/brad]
if debug
   plot( trialPoints(2,:), trialPoints(1,:), 'm*', 'MarkerSize',5 );
end
% find the saddles more precisely
% choose the radius to be about half the distance from the marker corner to the
% nearest corner of code dots
%rad = 0.75 * norm( H(1:2,1:2) * (resBorder./resMarker)' );
rad = floor(0.5 * norm( H(1:2,1:2) * (resBorder./resMarker)' ));
warning off MATLAB:nearlySingularMatrix
[trialPoints] = cornerfinder_saddle_point( flipud(trialPoints), I, rad,rad );
[trialPoints,good] = cornerfinder_saddle_point( trialPoints, I, 9,9 ); %before trying Anika's dataset: 3x3
warning on MATLAB:nearlySingularMatrix
% reject the trial points where the saddlefinder didn't converge
iPtMissed = flipud( trialPoints(:,good) )';
wPtMissed = [mrow(good)', mcol(good)'];


% reject any points we found so far that don't really looks like saddles (could
% happen if occluder lies near a corner point and saddle finder converges to
% something on the occluder's edge).
nPoints = size( iPt, 1 );
valid = true( nPoints, 1 );
for i = 1:nPoints
    valid(i) = validate_point( I, iPt(i,1), iPt(i,2), rad );
end
%uncomment this to see which saddles get rejected due to this test
if debug
    for i = 1:length( valid )
        if ~valid(i)
            plot( iPt(i,2), iPt(i,1), 'yO', 'MarkerSize',15 );
        end
    end
end
wPt = wPt(valid,:);
iPt = iPt(valid,:);



% validate saddles in regions around trial points

% examine the saddle points we did detect (so are confident about) and use
% the neighbourhoods of those points as training data to find the beta
% distribution that best describes those regions (depending on the angle of
% the plane in the image, there could be more 'white' around a saddle than
% 'black'. Note that roughly half the saddles will be b/w/b/w while the other
% half will be w/b/w/b. So if each quadrant doesn't take up exactly one quarter
% of the square sample patch (due to perspective warping etc) then you could get
% more black pixels in half the points and more white in the other half, which
% would confuse things when comparing distributions. So we invert the sample
% values from half of them (based on their world coordinates to determine
% whether even or odd).
% NOTE: should really sample in a disc shaped region around the point, but it's
% must simpler to code a square sampling grid
[sampleCols,sampleRows] = meshgrid( -rad:rad, -rad:rad );
nPoints = size( iPt, 1 );
beta = zeros( nPoints, 2 );
samples = cell( nPoints, 1 );
for i = 1:nPoints
    xi = iPt(i,2) + sampleCols - 0.5;
    yi = iPt(i,1) + sampleRows + 0.5;
% uncomment to see locations sampled around each point
%     if debug
%         plot( xi, yi, 'c+' );
%     end
    zi = interp2( I, xi, yi, '*linear' );
    zi = zi(~isnan(zi)); % in case of extrapolating outside image border
    % stretch patch to full 0..1 dynamic range
    zi = imadjust( zi, [min(zi(:)),max(zi(:))], [0.01,0.99] );
    % invert half the patches
    if mod( sum(wPt(i,:)), 2 )
        zi = 1 - zi;
    end
    %imwrite( zi, sprintf('betanhood_%2d.png',i) );
    beta(i,:) = betafit( double(zi(:)) );
    %samples{i} = zi(:);
end
% reject outliers (points we thought were saddles but actually aren't - this can
% happen when an occluder just slightly blocks a corner so we still detect the
% slightly warped quad, and then the location of the true corner will be
% detected somewhere at random on the edge of the occluder)
beta = abs( beta(:,1) - beta(:,2) );
beta_median = median( beta );
beta_std = std( beta );
outliers = (beta > 0.05) & (abs( beta - beta_median ) > (3 * beta_std));

% bhattacharyya method doesn't really work. going back to beta dist. perhaps i
% should check stddev on all alpha/beta/abs(alpha-beta) params?
%sampleHist = hist( vertcat(samples{:}), round(rad) );
%bhatDist = zeros( nPoints, 1 );
%for i = 1:nPoints
%    h = hist( samples{i}, round(rad) );
%    bhatDist(i) = bhattacharyya( sampleHist, h );
%end
%outliers = abs(bhatDist-mean(bhatDist)) > 3*std(bhatDist);

% uncomment this to see which saddles get rejected due to this test
if debug
    for i = 1:length( outliers )
        if outliers(i)
            plot( iPt(i,2), iPt(i,1), 'yO', 'MarkerSize',15, 'LineWidth',3 );
        end
    end
end
wPt = wPt(~outliers,:);
iPt = iPt(~outliers,:);
beta = beta(~outliers,:);
beta_median = median( beta );
beta_std = std( beta );


% now check the distributions around the potential points, and see if they
% are similar to the one found above
[sampleCols,sampleRows] = meshgrid( -rad:rad, -rad:rad );
nPoints = size( iPtMissed, 1 );
valid = true( nPoints, 1 );
for i = 1:nPoints
    xi = iPtMissed(i,2) + sampleCols - 0.5;
    yi = iPtMissed(i,1) + sampleRows + 0.5;
    zi = interp2( I, xi, yi, '*linear' );
    zi = zi(~isnan(zi)); % in case of extrapolating outside image border
    zi = imadjust( zi, [min(zi(:)),max(zi(:))], [0.01,0.99] );
    % invert half the patches
    if mod( sum(wPtMissed(i,:)), 2 )
        zi = 1 - zi;
    end
    beta = betafit( double(zi(:)) ); 
    underOne = all( beta < 1.0 );
    beta = abs( beta(1) - beta(2) );
    valid(i) = underOne & (beta < 0.2) & (abs( beta - beta_median ) < (6 * beta_std));
    
    %debug: always return true to disable the test
    %valid(i) = true;
    
    % now that the distribution is correct, do a second test to check if
    % opposite points of a circle sampled around the trial point have
    % matching intensities. alas we can't simply check to see if all opposite
    % points are below some small threshold, because sample points that happen
    % to lie near the black/white boundary could easily be very different from
    % their opposite peers. so instead we check that _most_ of the points are
    % below a small threshold, and _all_ of them are below a larger threshold
    if valid(i)
        valid(i) = validate_point( I, iPtMissed(i,1), iPtMissed(i,2), rad );
    end
    
end


if debug
    plot( iPt(:,2), iPt(:,1), 'rO', 'MarkerSize',10,'LineWidth',3 );
end


if debug
    iPtMissedBad = iPtMissed(~valid,:);
    iPtMissedGood = iPtMissed(valid,:);
    plot( iPtMissedGood(:,2), iPtMissedGood(:,1), 'm+', 'MarkerSize',10 );
    plot( iPtMissedGood(:,2), iPtMissedGood(:,1), 'mO', 'MarkerSize',10, 'LineWidth',3 );
end


% take only the valid missed points
wPt = [wPt; wPtMissed(valid,:)];
iPt = [iPt; iPtMissed(valid,:)];
[wPt,idx] = sortrows( wPt, [1,2] );
iPt = iPt(idx,:);



% by default the topleft of the image is the origin, but this can cause
% strange y-axis flips in the calibrated 3D world, so we can choose to make
% the bottom left point the origin, making the image the standard
% upper-right positive quadrant
bottomleftorigin = false;
if bottomleftorigin
    iPt(:,1) = size(I,1) - iPt(:,1);
end



% matlab coordinates are [row,col] with (1,1) being the centre of the top
% left pixel, so images run from (0.5,0.5) to (M+0.5,N+0.5). C coordinates
% are [col,row] with (0,0) being the centre of the top left pixel
output_C_coordinates = false;
if output_C_coordinates
    iPt = fliplr( iPt ) - 1;
end




% output to XML string
% switching from matlab's row/col to imagespace x,y coordinate order
%xmlCorr = sprintf( '<c u="%f" v="%f" x="%f" y="%f" z="0" />\n', ...
%                   [iPt(:,2), iPt(:,1), wPt(:,2), wPt(:,1)]' );
%disp( xmlCorr );               





% colour coding:
% red * = initial seed points
% red + = initial saddle points
% red circles = saddles found easily
% magenta * = seed points for missed points, found via ransac homography&undist
% magenta + = saddle points for missed points that converged in saddlefinder
% magenta circles = missed points that converged
% yellow circles (if uncommented) = points rejected due to filters



% simple function allowing one to write "[x,y,...] = splitvec([1,2,...])"
% the built-in function "deal" appears not to support this, and simply
% calling "num2cell" inline doesn't work either
function varargout = splitvec( v )
	if ~isnumeric( v ) || ~isvector( v )
		error( 'Can only split simple numeric vectors' );
	end
	varargout = num2cell( v );

    
% remove some elements from a connected-components data structure
% given a boolean mask vector this will return a modified connected
% components structure where the unflagged (0) entries are removed
% mask may alternatively be an index vector (eg [2,12,13,207] ) of values
% to keep
function CC = filtercc( CC, mask )
    if ~isvector( mask )
        error( 'Mask must be a vector' );
    end
    CC.PixelIdxList = CC.PixelIdxList(mask);
    if islogical( mask )
        CC.NumObjects = nnz( mask );
    else
        CC.NumObjects = numel( mask );
    end
    
    
% old code: %gradients = ( circshift(perim,-1) - circshift(perim,1) ) / 2;
% mask: binary image representing the filled region (tightly cropped)
% isQuad: boolean, true if the region is approximately a parallelogram
% corners: 2x4 matrix of corner coordinates
% corners_0: 2x4 matrix of original marker corners, before rotation in
%            layout2
function [isQuad,corners,corners_0] = fitquad( bbox, mask, layout )
    % default return values
    isQuad = false;
    corners = [];
    corners_0 = [];
    % remove small spurs in edge image (holes in mask) that mess things up
    % faster to do it only on potential regions than on the whole image
    mask = padarray( mask, [3,3] );
    mask = imclose( mask, strel('disk',2) );
    mask = mask(4:end-3,4:end-3);
    % trace the perimeter
    seedPoint = find( mask, 1 );
    [sy,sx] = ind2sub( size(mask), seedPoint );  
    perim = bwtraceboundary( mask, [sy,sx], 'NE', 8, Inf, 'clockwise' );
    % using gaussian-smoothed 1st derivative (central differences)
    % precomputed kernel is conv([-1,0,1]/2,fspecial('gaussian',[3,1],1))
    %kernel = [-0.1370;-0.2259;0;0.2259;0.1370];
    % precomputed kernel is conv([-1,0,1]/2,fspecial('gaussian',[5,1],5/3))
    %kernel = [-0.0668;-0.1146;-0.0704;0;0.0704;0.1146;0.0668];
    %bugfix: too small kernel allows small holes and spurs in perimeter to be
    %treated as separate gradient clusters. therefore we need to smooth by a
    %kernel roughly on the order of the size of the smallest quad edge
    kernelSize = floor( min(size(mask))/3 );
    kernel = conv([-1,0,1]/2,fspecial('gaussian',[kernelSize,1],kernelSize/6));
    gradients = imfilter( perim, kernel, 'circular' );
    % this could throw an empytycluster warning if it can't find 4 clusters
    % we disable the warning message in the caller, and check the returned
    % K afterwards to see if it worked
    %[clusteridx,clustermeans] = kmeans( gradients, 4, 'replicates',1 );     % maybe reduce repl for performance reasons?
        
    %debug: the above kmeans clustering works mostly, but is quite slow and can
    % sometimes miss the true clusters because of its random initialisation. 
    % Below we are trying to manually set seed locations to avoid having to do
    % multiple random trial of kmeans. because we know the perimeter
    % istraversed in order, we can choose equispaced gradients as the seeds
    seedIdx1 = floor( linspace( 1, size(gradients,1), 5 ) );
    seedIdx2 = seedIdx1 + floor( diff(seedIdx1(1:2))/2 );
    seeds = cat( 3, gradients(seedIdx1(1:4),:), gradients(seedIdx2(1:4),:) );
    [clusteridx,clustermeans] = kmeans( gradients, 4, 'start',seeds, 'emptyaction','drop' );
    
    nClusters = sum( isfinite( clustermeans(:,1) ) );
    if nClusters ~= 4
        return;
    end    
    % initial assignments of points to lines
    points = { perim(clusteridx==1,:), perim(clusteridx==2,:), ...
               perim(clusteridx==3,:), perim(clusteridx==4,:) };        
    % check to see if any one cluster has too few points
    nPointsPerCluster = cellfun( @length, points );  
    if min(nPointsPerCluster) < max(nPointsPerCluster)/16
        return;
    end          
    % fit a line through the points in each cluster       
    lines  = cellfun( @fitline, points, 'UniformOutput',false );   
    % sometimes little spurs along the edges have the same gradient as
    % other sides of the quad, so they belong to same cluster, but we
    % really don't want to include those points in the linefit. So use
    % Lloyd's Algorithm to fit the lines
    %point2line = @(P,L) abs( dot( repmat(L.n,[size(P,1),1]), ...             % does binding this function cost much time?
    %                              repmat(L.p,[size(P,1),1])-P, 2 ) );        % changed to script function
    maxIter = 6;
    iter = 0;
    isConverged = false; 
    isDegenerate = false;
    %minPtsPrClstr = length( gradients ) / 4 * (0.75);
    minPtsPrClstr = length( gradients ) / 16 * (0.75);
    while ~isDegenerate() && iter<maxIter && ~isConverged 
        % uncomment to check the quad fitting function...
        %debugPlotQuadFit( mask, points, lines );
        distances = cellfun( @point2line, {perim,perim,perim,perim}, ...
                             lines, 'UniformOutput',false );
        [mindist,minidx] = min( cell2mat(distances), [], 2 );
        newPoints = { perim(minidx==1,:), perim(minidx==2,:), ...
                      perim(minidx==3,:), perim(minidx==4,:) };      
        iter = iter + 1;
        isConverged = isequal( points, newPoints ); 
        %isDegenerate = any( cellfun(@numel,newPoints) < 2*minPtsPrClstr );
        isDegenerate = any( cellfun(@numel,newPoints) < minPtsPrClstr );
        if ~isConverged   
            points = newPoints;
            lines  = cellfun( @fitline, points, 'UniformOutput',false );      
        end
    end
    if ~isConverged
        return;
    end   
    % figure out which line is parallel to the first
    dotprods = [0,0,0,0];
    for i = 2:4
        dotprods(i) = abs( dot( lines{1}.d, lines{i}.d ) );
    end
    [dotprods,idx] = sort( dotprods );
    %if dotprods(4) < 0.75 || any( dotprods(2:3)>0.5 )
    %    return;
    %end
    if dotprods(4) < 0.75 || abs(diff(dotprods(2:3)))>0.2
        return;
    end
    lines = lines(idx);        
    % get the line intersection points
    corners = [ intersectLines( lines{1}, lines{2} ),...
                intersectLines( lines{2}, lines{4} ),...
                intersectLines( lines{4}, lines{3} ),...
                intersectLines( lines{3}, lines{1} ) ];   
    % check if any corners points lie well outside the boundingbox    
    midToCorner = corners - repmat( size(mask)'/2, [1,4] );
    [theta,rho] = cart2pol( midToCorner(1,:), midToCorner(2,:) );
    if any( rho > 1.1*sqrt(sum(size(mask).^2))/2 )
        return;
    end    
    % sort the corners into counterclockwise order
    [theta,idx] = sort( theta );
    %rho = rho(idx); % don't need rho again
    corners = corners(:,idx);
    % if using the rotated layout (2) then we need to transform these
    % corners so that they lie at the bowtie centres
    if layout == 2
        corners_0 = corners + repmat( bbox([2,1])', [1,4] );
        centroid = mean( corners, 2 );
        centroid2corners = corners - repmat( centroid, [1,4] );
        % 45 degree rotation + sqrt(2) uniform scale matrix
        tformMatrix = [1 1; -1 1];
        centroid2corners = tformMatrix * centroid2corners;
        corners = centroid2corners + repmat( centroid, [1,4] );
    end
    % convert corners into global image coordinates
    corners = corners + repmat( bbox([2,1])', [1,4] );
    isQuad = true;
    
    


% least square line fit goes through the mean of the 2D points with
% direction being the eigenvector corresponding to the largest eigenvalue
% of the covarience matrix
% points is Nx2
% dir is a unit vector
function line = fitline( points )
    points = points(~isnan(points(:,1)),:);
    if size( points, 2 ) < 2
        line = [];
        return;
    end
    point = mean( points );
    covmatrix = cov( points );
    [evecs,evals] = eig( covmatrix );
    [maxeval,maxidx] = max( diag(evals) );
    dir = evecs(:,maxidx);
    normal = [0,1;-1,0] * dir;
    line = struct( 'p',point, 'd',dir', 'n',normal' );
    
    
% return the 2D intersection point (column vec) of two lines
% assuming they are not colinear
% column vec: [vertical coord; horizontal coord]
function p = intersectLines( line1, line2 )
    alpha = [line1.d',line2.d'] \ (line2.p'-line1.p');
    p = line1.p' + alpha(1) * line1.d';

   
% orthogonal distances from Nx2 array of points to a line
function y = point2line( P, L ) 
    y = abs( dot( repmat(L.n,[size(P,1),1]), ...
                  repmat(L.p,[size(P,1),1])-P, 2 ) );
  
              
% plot debug into to see if lloyds algorithm is converging
% region is the binary image corresponding to 'FilledImage'
% lines and points are 4-element cell arrays
function debugPlotQuadFit( region, points, lines )
    nLines = length( lines );
    assert( nLines == length(points) );
    colours = 'rgbycmk';
    region = bwperim( region );
    imshow( region * 0.2 );
    hold on;
    dbgPlotPoints = @(p,c) plot( p(:,2), p(:,1), strcat(c,'O') );
    k = max( size(region) );
    dbgPlotLine = @(l,c) line( [ l.p(2)-k*l.d(2), l.p(2)+k*l.d(2) ], ...
                               [ l.p(1)-k*l.d(1), l.p(1)+k*l.d(1) ], ...
                               'Color',c );
    for i = 1:nLines
        dbgPlotPoints( points{i}, colours(i) );
        dbgPlotLine( lines{i}, colours(i) );
    end
    hold off;
    
    
 function y = crc_encode( data, dBits, cBits )
     if data > 2^dBits
         disp( 'Warning: data exceeds permissible bit depth' );
         data = bitand( data, 2^dBits-1 );
     end
     % from wikipedia CRC page (normal representation, decimal)
     % zeros represent invalid generator
     generators = [1,0,0,3,21,3,9,7,0,563,901,2063,0,0,17817,4129];
     poly = bitor( generators(cBits), 2^cBits );
     crc = leftshift( data, cBits );
     div = leftshift( poly, dBits );
     while div >= poly
         if bitxor( crc, div ) < crc
             crc = bitxor( crc, div );
         end
         div = div / 2;
     end
     y = bitor( leftshift(data,cBits), crc );
     
     
function y = crc_validate( payload, idBits, crcBits )
    % bitshifting in matlab is very counterintuitive, so we work with
    % doubles here instead of integers
    idmask = leftshift( 2^idBits-1, crcBits );
    data = rightshift( bitand(payload,idmask), crcBits );
    %data = floor( rightshift(double(payload),double(crcBits)) );
    y = (crc_encode(data,idBits,crcBits) == payload);
     
  
% these function are duplicated (and modified) here from the fixed-point
% toolbox because the mcc compiler does not support fixed point toolboxes
% in deployed applications. These functions don't really work nicely when
% bits slide off the edge, so the calling code takes care to zero out those
% bits first
function y = leftshift(x, n)
    y = x * 2.^(n);

function y = rightshift(x, n)
    y = x * 2.^(-n);
    
% the matlab stats toolbox has a betafit function that works nicely, but it
% uses fminsearch to find the alpha/beta parameters. This gives results
% that do not correspond to the parameter estimate formulae I found on the
% web. Those estimates are computed by this function, which is much faster.
% I'm not sure yet whether these are accurate enough to discriminate saddle
% points, but you can easily switch back to the matlab function by
% commenting out this function. Update: this one seems to work
function phat = betafit( x )
    m = mean( x(:) );
    v = var( x(:) );
    alpha = m*( m*(1-m)/v - 1 );
    beta = (1-m)*( m*(1-m)/v - 1 );
    phat = [alpha, beta];
        

    
% apply radial distortion to image points. Only considering the first order
% term, since it probably dominates the distortion, and we won't really have
% enough data to get a stable estimate of higher order distortion coefficients
function iPt = radial_distort( iPt, normMatrix, r )
    iPt(:,3) = 1;
    iPt = (normMatrix * iPt');
    iPt = iPt(1:2,:)';
    mag = sqrt( sum(iPt.^2,2) );
    scale = r * mag;
    iPt(:,1) = iPt(:,1) + iPt(:,1) .* scale;
    iPt(:,2) = iPt(:,2) + iPt(:,2) .* scale;
    iPt(:,3) = 1;
    iPt = normMatrix \ iPt';
    iPt = iPt(1:2,:)';
    
    
% given the point correspondences, return a measure of the radial
% distortion by fitting straight lines though the rows and columns of the
% grid and summing up the distances of points in those rows/cols to the
% appropriate lines
function y = distortion_metric( wPt, iPt, resPattern, normMatrix, r )
    
    if length( iPt ) < 8
        y = 0;
        return;
    end
    
    % apply radial shift by r    
    iPt = radial_distort( iPt, normMatrix, r );
    
    % now compute colinearity error metric for each point by fitting a straight
    % line through all the points in each grid row and column, then finding the
    % distance from each of those points to that line.
    dists = cell( sum(resPattern+1), 1 );
    for i = 0:resPattern(2)
        pts = iPt(wPt(:,1)==i,:);
        if size( pts, 1 ) > 2
            L = fitline( pts );
            dists{i+1} = point2line( pts, L );
        end
    end
    for j = 0:resPattern(1)
        pts = iPt(wPt(:,2)==j,:);
        if size( pts, 1 ) > 2
            L = fitline( pts );
            dists{i+j+2} = point2line( pts, L );
        end
    end
    y = vertcat( dists{:} );
    y = sqrt( y );
 
    
%function check_3rd_party_code( filename, url )
%    if ~exist( filename, 'file' )
%        sprintf( 'Warning: %s not found (3rd party code)\n', filename );
%        sprintf( 'Would you like to download it from %s now? (1/0)', url );
%        if input( ' > ' )
%            wget( url );
%        end
%    end


% sample opposing points on a small circle around a point, and return true if
% most of the opposing points have similar values to each other
function v = validate_point( I, row, col, radius )
     % get sample points on a circle
    nCircSamples = ceil( 2 * pi * radius );
    [sx,sy] = pol2cart( linspace(0,2*pi,nCircSamples+1), radius );
    % no point in sampling at 2*pi as well as 0
    % [bugfix 24 Jan 2012: remove incorrect 0.5 shift (Ofri)]
    sx = sx(1:end-1) + col;%- 0.5;
    sy = sy(1:end-1) + row;%+ 0.5;
    zi = interp2( I, sx, sy, '*linear' );
    %zi = imadjust( zi, [min(zi(:)),max(zi(:))], [0,1] );
%    zi = [ zi(1:nCircSamples/2); zi(nCircSamples/2+1:end) ];
    %v = all( abs(diff(zi)) < 0.1 ); % this threshold is too aggressive
    % rather we say that all of the points should be within loose threshold,
    % while just most of them should satisfy a tighter threshold. Because the
    % points sampled on the boundary between black/white can easily not match
    % their partner
    %  delta = abs( diff(zi) );
    %  v = all(delta<0.6) && (nnz(delta<0.3) > 0.5*(nCircSamples/2));
    % now switching to using binary image instead, so we just want most of the
    % opposing points to be equal (some lying very near the black/white boundary
    % might not be)
%    similar = zi(1,:) == zi(2,:);
%    v = nnz( similar ) > floor(0.75 * nCircSamples/2);
    % differtent test
    v = false;
    zi = zi(~isnan(zi));
    if ~isempty( zi )
        kernel = fspecial( 'gaussian', [round(nCircSamples/4),1], nCircSamples/4/6 );
        zi_smooth = imfilter( zi, kernel, 'circular' );
        zi_smooth = zi_smooth > mean( zi_smooth );
        nFlips = sum( abs( diff( [zi_smooth,zi_smooth(1)] ) ) );
        v = nFlips == 4;
    end
    
    %debug: always return true to disable this test
    %v = true;


% distance metric between two probability distributions (1d vectors)
function y = bhattacharyya( p, q )
    p = p(:);
    q = q(:);
    p = p / sum( p );
    q = q / sum( q );
    y = -log( sum( sqrt(p.*q) ) );
    
