function main(imageFiles, dataFile, outputFile)
% eg: main('test/*.png', 'test/output.mat', 'test/points.h5')


% 'dir' doesn't work because you don't get full paths 
[~,files] = system( ['ls ', imageFiles] );
files = textscan( files, '%s' );
files = files{1};


nFiles = length( files );

if nFiles == 0
    error( 'No input files found' );
end
  
 
for i = 1:nFiles
    file = files{i};
    disp( file );
    [wPt,iPt] = caltag( file, dataFile, false );
    
    dset = ['/', file, '/world'];
    h5create( outputFile, dset, size(wPt), 'ChunkSize',[16,2], 'Deflate',9, 'Shuffle',true );
    h5write( outputFile, dset, wPt );
     
    dset = ['/', file, '/image'];
    h5create( outputFile, dset, size(iPt), 'ChunkSize',[16,2], 'Deflate',9, 'Shuffle',true );
    h5write( outputFile, dset, iPt );
end
