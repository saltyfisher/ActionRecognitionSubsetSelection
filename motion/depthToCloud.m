function [pcloud] = depthToCloud(depth,xstep,ystep)
topleft = [1 1];
depth= double(depth);
[xsize ysize]=size(depth);
k=1;
pcloud = zeros(xsize,ysize,3);
for i=1:xsize
    for j=1:ysize
        [ pcloud(i,j,1) pcloud(i,j,2) pcloud(i,j,3) ] = myDepth2RealWorldCoordinate( i, j, depth(i,j),xstep,ystep );
    end
end





