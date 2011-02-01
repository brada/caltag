
function mySaddleFinder

% see www.dip.ee.uct.ac.za/publications/theses/MScArne.pdf  (Appendix1)

% example code...
if ndims(I) == 3
    I = rgb2gray( I );
end;
x = [779.3; 1680.3]; % initial saddle guess [y (vert); x (horiz)]
winSize = 17;
assert( mod(winSize,2)==1 );
halfw = (winSize-1)/2;
window = im2double( I(x(1)-halfw:x(1)+halfw, x(2)-halfw:x(2)+halfw) );