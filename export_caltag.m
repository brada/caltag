function export_caltag( imagefilepattern, datafile, destfile )
%EXPORT_CALTAG Run Caltag on set of images and write results to HDF5 file
%	EXPORT_CALTAG(imagefilepattern,datafile,destfile)
%   args:
%   - imagefilepattern (eg. 'frame_*.png')
%   - datafile (eg. 'pattern_14x9.mat')
%   - destfile (eg. 'caltagpoints.h5')
%
% Contributed by Felix Heide, 3 Oct 2011
% Modified by Brad Atcheson, 25 Oct 2011

debug = false;
[srcpath, ~] = fileparts( imagefilepattern );
imagefiles = dir( imagefilepattern );

%Results
WPT = [];
IPT = [];
POINTCOUNT = [];
                                                                
%Now run Caltag on the images
for i = 1:length(imagefiles)     

    name = fullfile( srcpath, imagefiles(i).name );
    disp( ['Running CALTag on ', name] );
    [wPt, iPt] = caltag( name, datafile, debug );
    disp( [' found ', num2str(size(iPt,1)), ' points'] );

    %Concatenate WPT and IPT
    WPT = [WPT;wPt];
    IPT = [IPT;iPt];
    POINTCOUNT = [POINTCOUNT;size(iPt,1)];
    
    if i == 1
        I = imread( name );
        IMGSIZE = [size(I, 2), size(I, 1)];
    end
        
end

%Add third dimension to WPT
WPT(:,3) = 0;

%Convert to C-style coordinates
%That is [col,row] with (0,0) being the centre of the top left pixel
IPT = fliplr( IPT ) - 1;

%Save result
hdf5write( destfile, '/WPT', WPT, '/IPT', IPT, ...
            '/POINTCOUNT', POINTCOUNT, '/IMGSIZE', IMGSIZE, ...
            '/IMGFILES', {imagefiles.name} );
disp( 'Done' );
