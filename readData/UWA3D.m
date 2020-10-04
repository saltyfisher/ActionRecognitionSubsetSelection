clear;clc;

actor_id = {'01','02','03','04','05','06','07','08','09','10'};
activity_label = {'01','02','03','04','05','06','07','08','09','10','11','12','13','14','15',...
                    '16','17','18','19','20','21','22','23','24','25','26','27','28','29','30'};

i = 1;
j = 1;
depth_list = dir(['E:\KangLi\Datasets\UWA3D_Single_View_depth\UWA3D_Single_View_depth\s' actor_id{i} '_a' activity_label{j} '_e01\*.png']);

allpt = [];
for i = 1:length(depth_list)
    depth_img = imread(fullfile(depth_list(i).folder,depth_list(i).name));
    pt = DepthtoPointcloudKinect1(depth_img);
    tmp = pt;
    tmp(sum(tmp,2)==0,:) = [];
    tmp(:,[1,3]) = tmp(:,[3,1]);
    tmp(:,[2,3]) = tmp(:,[3,2]);
    tmp(:,3) = tmp(:,3)*-1;
    allpt = [allpt;tmp];
%     pcshow(tmp);
%     xlabel('X');ylabel('Y');zlabel('Z');
%     pause(1);
end
pcshow(allpt);
xlabel('X');ylabel('Y');zlabel('Z');