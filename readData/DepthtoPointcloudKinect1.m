function P=DepthtoPointcloudKinect1(depth)
    r = size(depth,1);
    c = size(depth,2);
    
    camera_factor = 1;
    camera_cx = 325.5;
    camera_cy = 253.5;
    camera_fx = 518.0;
    camera_fy = 519.0;
    
    P = zeros(r*c,3);
    for i = 1:r
        for j = 1:c
            z = double(depth(i,j)) / camera_factor;
            x = (j - camera_cx) * z / camera_fx;
            y = (i - camera_cy) * z / camera_fy;

            P((i-1)*r+j,:) = [x,y,z];
        end
    end
end