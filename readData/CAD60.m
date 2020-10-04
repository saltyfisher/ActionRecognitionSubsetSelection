clear;clc;

actor_id = ['1','2','3','4'];
activity_label = {{'0512172825','still'};{'0512171649','talking on the phone'};
                  {'0512175502','writing on whiteboard'};{'0512173312','drinking water'};
                  {'0512164800','rinsing mouth with water'};{'0512164529','brushing teeth'};
                  {'0512165243','wearing contact lenses'};{'0512165327','wearing contact lenses'};
                  {'0512174513','talking on couch'};{'0512174643','relaxing on couch'};
                  {'0512171207','cooking (chopping)'};{'0512171444','cooking (stirring)'};
                  {'0512173520','opening pill container'};{'0512173548','opening pill container'};
                  {'0512173623','opening pill container'};{'0512170134','working on computer'};
                  {'0512174930','random'}};

i = 1;
j = 4;
j = activity_label{j};
depth_list = dir(['E:\KangLi\Datasets\CAD60\data\data' actor_id(i) '\' j{1} '\Depth*.png']);
sample_rate = 1:10:length(depth_list);

depth_list = depth_list(sample_rate);
allpt = [];
for i = 1:length(depth_list)
    depth_img = imread(fullfile(depth_list(i).folder,depth_list(i).name));
    pt = DepthtoPointcloudKinect1(depth_img);
    tmp = pt;
    tmp(sum(tmp,2)==0,:) = [];
    tmp(:,[1,3]) = tmp(:,[3,1]);
    tmp(:,[2,3]) = tmp(:,[3,2]);
    tmp(:,3) = tmp(:,3)*-1;
    allpt = [allpt;pt];
%     pcshow(tmp);
%     xlabel('X');ylabel('Y');zlabel('Z');
%     pause(1);
end
pcshow(allpt);
xlabel('X');ylabel('Y');zlabel('Z');