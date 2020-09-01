function [ pWorldX, pWorldY, pWorldZ ] = myDepth2RealWorldCoordinate( depthX, depthY, depthZ,resolutionX,resolutionY )
horizontalFov=58.5;%degree
verticalFov=45.6;%degree
xzFactor = tan(horizontalFov /2 *(pi/180)) * 2;
yzFactor = tan(verticalFov /2 *(pi/180)) * 2;
normalizedX = depthX / resolutionX - 0.5;
normalizedY = 0.5 - depthY / resolutionY;
pWorldX = normalizedX * depthZ * xzFactor;
pWorldY = normalizedY * depthZ * yzFactor;
pWorldZ = depthZ;
end

