clear;close all;
action_id = '2';
pt_list = dir(['./PointCloud//a' action_id '*.mat']);
vd_list = dir(['./Depth//a' action_id '*.mat']);
tf_list = dir(['./tHOPC3Dfeatures1/a' action_id '*.mat']);
sf_list = dir(['./sHOPC3Dfeatures/a' action_id '*.mat']);
tf = load(fullfile(tf_list(1).folder,tf_list(1).name));
sf = load(fullfile(sf_list(1).folder,sf_list(1).name));

sample_id = 1;
pointcloud = load(fullfile(pt_list(1).folder, pt_list(sample_id).name));
videoCell = load(fullfile(vd_list(1).folder, vd_list(sample_id).name));
videoCell = videoCell.depth;
allpcloud = pointcloud.allpcloud;
video = zeros(240, 320, size(videoCell, 2));
for j = 1:size(videoCell, 2)
    video(:, :, j) = videoCell{j};
end
allpcloud = permute(allpcloud, [1,2,4,3]);
pointcloud = reshape(allpcloud,[],3);

sf = sf.Q;
tf = tf.Q;

t1 = sum((sf - tf).^2, ndims(sf));
t2 = sum(sf + tf, ndims(sf)) + 1e-5;
fea = t1./t2;

% f_pt = pointCloud(pointcloud);
% f_pt = pcdenoise(f_pt);
% pcshow(f_pt);

tps=3;%%constant temporal scale for MSRAction3D Dataset, same value for our UWA3D dataset
Redious=140;%%%constant spatial scale for MSRAction3D Dataset. This dataset captured by time-of-flight camera
%Redious=140;%%%constant spatial scale for UWA3D Dataset. This dataset captured by Kinect camera. So the Redious value is different from MSRAction3D Dataset.
fVector=[];
xstep=40; ystep=32; zstep=3;%%size of each video is 240*160*f, xstem denotes the number of pixels along X dimension. So, 240/40=6, 160/32=5,
%%i.e. We divide each depth video into 6*5*3 spatio-temporal cells along
%%X,Y, and T dimensions.
patchsize=10;%%to reduce the time complexity, we limit the number of points around each point.

f=size(video,3);
c=size(video,2);
r=size(video,1);

%tps = floor(f/2);
% PC1x=zeros(r+patchsize,c+patchsize,f+tps);
% PC1y=zeros(r+patchsize,c+patchsize,f+tps);
% PC1z=zeros(r+patchsize,c+patchsize,f+tps);
% PCV1=zeros(r+patchsize,c+patchsize,f+tps);
% PC2x=zeros(r+patchsize,c+patchsize,f+tps);
% PC2y=zeros(r+patchsize,c+patchsize,f+tps);
% PC2z=zeros(r+patchsize,c+patchsize,f+tps);
% PCV2=zeros(r+patchsize,c+patchsize,f+tps);
% PC3x=zeros(r+patchsize,c+patchsize,f+tps);
% PC3y=zeros(r+patchsize,c+patchsize,f+tps);
% PC3z=zeros(r+patchsize,c+patchsize,f+tps);
% PCV3=zeros(r+patchsize,c+patchsize,f+tps);

video = padarray(video,[patchsize patchsize]);
allpcloud = padarray(allpcloud,[patchsize patchsize]);
outVideo=zeros(size(video,1),size(video,2),size(video,3)+2*tps,3);
outVideo(:,:,tps+1:f+tps,:)=allpcloud;
while true
    frameno = randi([tps+1,f+tps],1);
    colno = randi([patchsize+1,c+patchsize],1);
    rowno = randi([patchsize+1,r+patchsize],1);

    pt = outVideo(rowno,colno,frameno,:);
    pt = pt(:);
    if pt(1)<-200 || pt(1)>-100 || pt(2)<0 || pt(2)>50 || pt(3)<350 || pt(3)>450
        continue
    end
    if pt(1)~=0 || pt(2)~=0 || pt(3)~=0 >0 && fea(rowno,colno,frameno)>0
        break;
    end
end
Q = {};
%初始化起始点,决策点集,起点点集
ind0 = [rowno,colno,frameno];
% ind0 = [123 220 18];
p0 = outVideo(ind0(1),ind0(2),ind0(3),:);
p0 = p0(:)';
decision_set = [];
start_set = [];
start_set = [start_set; ind0];
local_start_set = [];
t = 0;
orientation0 = [];
local_motion.pointcloud = fullfile(pt_list(1).folder, pt_list(sample_id).name);
tmp = reshape(outVideo,[],3);
tmp(find(tmp(:,3)==0),:) = [];
f_pt = pointCloud(tmp);
pcshow(f_pt);
hold on;
pcshow(p0,[1,0,0],'MarkerSize',400);
hold off;
flag = zeros(1,(r+2*patchsize)*(c+2*patchsize)*(f+2*tps));
flag(sum(tmp,2)==0) = 1;
figure;

ang1 = 45;
ang2 = 45;
while true
    local_start_set = [];
    for i = 1:size(start_set, 1)
        rowno = start_set(i,1);
        colno = start_set(i,2);
        frameno = start_set(i,3);
        
        pt=outVideo(rowno,colno,frameno,:);
        pt=pt(:);
%         start_point = pointCloud(pt');
%         start_point.Color = uint8([0,0,0]);
%         hold on;
%         pcshow(start_point,'MarkerSize',50);
%         hold off;
        local_motion.start_ind = start_set(i,:);
        local_motion.start_point = pt';
        cubex=outVideo(rowno-patchsize+1:rowno+patchsize,colno-patchsize+1:colno+patchsize,frameno-tps+1:frameno+tps,1);
        cubey=outVideo(rowno-patchsize+1:rowno+patchsize,colno-patchsize+1:colno+patchsize,frameno-tps+1:frameno+tps,2);
        cubez=outVideo(rowno-patchsize+1:rowno+patchsize,colno-patchsize+1:colno+patchsize,frameno-tps+1:frameno+tps,3);
        points(:,1)=cubex(:);
        points(:,2)=cubey(:);
        points(:,3)=cubez(:);
        outside_point = sum(points,2);

        dataTrans = [points(:,1)-pt(1) points(:,2)-pt(2) points(:,3)-pt(3)];
        dist = sqrt(sum(dataTrans.^2,2));
        ind = find((dist(outside_point~=0)<Redious));
        dataLocal = dataTrans(ind,:);
        if isempty(ind)
            continue
        end

        avg=sum(dataLocal)./size(dataLocal,1);
        dataLocal=dataLocal-repmat(avg,size(dataLocal,1),1);
        CM=dataLocal'*dataLocal;
        CM=CM./size(dataLocal,1);
        [V,D]=eig(CM);
        [B,ix] = sort(diag(D),'descend');
        orientation = V(:,ix(1));
        local_motion.orientation = orientation';
        if t == 0
            orientation0 = orientation;
        end

        %找到与中心点相距r/2的点
        ind = dist(outside_point~=0)<Redious;
        [ia,ix] = sort(dist(ind),'descend');
        ind = ind(ix(1:ceil(size(ind,1)*0.1),:));
        if isempty(ind)
            continue;
        end
        %找到与当前中心点同向的点
        tmp = repmat(orientation',size(dataTrans,1),1);
        tmp = sum(tmp.*dataTrans,2)./sqrt(sum(tmp.^2,2))./sqrt(sum(dataTrans.^2,2));
        tmp(~ind) = 0;
        tmp = rad2deg(acos(tmp));
        ind1 = tmp>0 & tmp<ang1;
        %找到与起始点同向的点
        tmp = repmat(orientation0',size(dataTrans,1),1);
        tmp = sum(tmp.*dataTrans,2)./sqrt(sum(tmp.^2,2))./sqrt(sum(dataTrans.^2,2));
        tmp(~ind) = 0;
        tmp = rad2deg(acos(tmp));
        ind2 = tmp>0 & tmp<ang2;
        %找到同时满足要求的点
        ind3 = (ind1&ind2)==1;
        local_motion.valid_point = [];
        if sum(ind3) == 0
            continue;
        end
        %映射回全局坐标
        tmp = find(ind3);
        [row,col,frm] = ind2sub([patchsize*2,patchsize*2,tps*2],tmp);
        %舍弃已经选择过的点
        row = rowno-patchsize+row;
        col = colno-patchsize+col;
        frm = frameno-tps+frm;
        tmp_ind = sub2ind([r+2*patchsize,c+2*patchsize,f+2*tps],row,col,frm);
        a = find(flag(tmp_ind)==0);
        if isempty(a)
            continue;
        end
        flag(tmp_ind(a)) = 1;
        [row,col,frm] = ind2sub([r+patchsize*2,c+patchsize*2,f+tps*2],tmp_ind(a));
        b = [row,col,frm];
        b(b(:,1)<patchsize | b(:,1)>(r),:) = [];
        b(b(:,2)<patchsize | b(:,2)>(c),:) = [];
        b(b(:,3)<tps | b(:,3)>(f),:) = [];
        local_start_set = [local_start_set;b];
        b = sub2ind([r+patchsize*2,c+patchsize*2,f+tps*2],b(:,1),b(:,2),b(:,3));
        valid_point = reshape(outVideo,[],3);
        local_motion.valid_point = valid_point(b,:);
%         hold on;
%         pcshow(pointCloud(local_motion.valid_point));
%         pcshow(pt',[1,0,0],'MarkerSize',400);
%         hold off;
        Q = [Q, {local_motion}];
    end
%     valid_point = [];
%     start_point = [];
%     orientation = [];
%     for j = 1:length(Q)
%         valid_point = [valid_point;Q{j}.valid_point];
%         start_point = [start_point;Q{j}.start_point];
%         orientation = [orientation;Q{j}.orientation];
%     end
%     hold on;
%     cmatrix = ones(size(start_point,1),3).*[1,0,0];
%     pcshow(start_point,cmatrix,'MarkerSize',400);
%     pcshow(valid_point);
%     pcshow(p0,[0,0,0],'MarkerSize',400);
%     quiver3(start_point(:,1),start_point(:,2),start_point(:,3),orientation(:,1),orientation(:,2),orientation(:,3),'-b','LineWidth',0.2,'MarkerSize',20);
%     hold off;
    size(local_start_set)
    if isempty(local_start_set)
        break
    end
    start_set = local_start_set;
    t = t+1;
    size(Q),t
end



% hold off;
save('Q.mat','Q');
outVideo=[];
clear('temp', 'points','dataLocal','cubex','cubey','cubez');