% adaptive thresholding based on code by Peter Kovesi at
% http://www.csse.uwa.edu.au/~pk/Research/MatlabFns/
% Modified by author to run faster by exploting seperability of Gaussian kernel
% Note: I believe this function can be improved by making t somehow dependent
% upon the dynamic range of the input image. For very low-contrast images, the t
% value needs to be reduced. As a workaround, we first stretch the image using
% imadjust to hopefully provide a consistent contrast.
function Y = adaptivethresh( I, k, t )

% default arguments
if nargin < 2
    k = round( length(I)/20 );
end
if nargin < 3
    t = 15;
end

% note: should be removed for speed if we can figure out how to make t dependent
% upon contrast
I = imadjust( I );

% large area blur
g = single( fspecial( 'gaussian', [size(I,1),1], k ) );
G = abs( fft(g) );
F = fft( I );
for col = 1:size(F,2)
    F(:,col) = F(:,col).*G;
end
F = fft( real( ifft(F) )' );
g = single( fspecial( 'gaussian', [size(I,2),1], k ) );
G = abs( fft(g) );
for col = 1:size(F,2)
    F(:,col) = F(:,col).*G;
end
F = real( ifft(F) )';

% threshold
Y = I > F*(1-t/100);
