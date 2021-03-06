clear;
close;
warning off;

action_id = {'01','02','03','04','05','06','07','08','09','10','11','12','13','14','15','16','17','18','19','20'};
subject_num = 10;
repeat_num = 3;
action_cube = cell(1, length(action_id));
for a = 1:length(action_id)
    depth_path = ['E:\KangLi\HOPC Matlab Code\HOPC Matlab Code\Depth\a' action_id{a} '*.mat'];
    depth_list = dir(depth_path);
    pointcloud_path = ['E:\KangLi\HOPC Matlab Code\HOPC Matlab Code\PointCloud\a' action_id{a} '*.mat'];
    pointcloud_list = dir(pointcloud_path);
    sfeature_path = ['E:\KangLi\HOPC Matlab Code\HOPC Matlab Code\sHOPC3Dfeatures\a' action_id{a} '*.mat'];
    sfeature_list = dir(sfeature_path);
    tfeature_path = ['E:\KangLi\HOPC Matlab Code\HOPC Matlab Code\tHOPC3Dfeatures1\a' action_id{a} '*.mat'];
    tfeature_list = dir(tfeature_path);
    x_max=-inf;x_min=inf;y_max=-inf;y_min=inf;
    alldepth = cell(1, length(depth_list));
    %%生成决策空间
    for i = 1:length(depth_list)
        load(fullfile(depth_list(i).folder, depth_list(i).name));
        depth = cell2mat(depth);
        depth = reshape(depth, 240, 320, []);
        z_max = max(max(max(depth)));
        z_min = min(min(depth(depth~=0)));
        [r,c,~] = find(depth~=0);
        c = rem((c-1),320)+1;
        y_min=min(r);x_min=min(c);y_max=max(r);x_max=max(c);
        
        sfeature = load(fullfile(sfeature_list(i).folder, sfeature_list(i).name));
        tfeature = load(fullfile(tfeature_list(i).folder, tfeature_list(i).name));
        sf = sfeature.Q;tf = tfeature.Q;

        t1 = sum((sf - tf).^2, ndims(sf));
        t2 = sum(sf + tf, ndims(sf)) + 1e-5;
        t3 = t1./t2;
        t3 = floor(t3);
        sf = floor(sf);
        t3 = t3(y_min:y_max, x_min:x_max,:,:);
        t3 = (t3-min(t3(:)))./(max(t3(:))-min(t3(:))-1);
        gamma_ = 0.02;
        t3 = t3 - gamma_;
        t3(t3<0) = 0;
        sf = sf(y_min:y_max, x_min:x_max,:,:);
        depth = depth(y_min:y_max, x_min:x_max, :);
        %//TODO:可以计算向量各个维度的统计量，加权求和
        sf = sum(sf, ndims(sf));
        %%转换到决策空间中
        %按时间顺序记录决策空间中的索引
        [h, w, f] = size(depth);
        depth = depth-z_min;
        %前f维是原函数变量，后f维是约束函数变量
        input_tmp = zeros(h, w, z_max-z_min, 2*f);
        [r,c,v] = find(depth>=0);
        for l = 1:length(v)
            idy = rem((c(l)-1),w)+1;
            idt = fix((c(l)-1)/w)+1;
            input_tmp(r(l),idy,v(l),idt+f) = t3(r(l),idy,idt);
            input_tmp(r(l),idy,v(l),idt) = sf(r(l),idy,idt);
        end
        input_var = zeros(h,w,z_max-z_min,f+1);
        input_var(:,:,:,1) = sum(input_tmp(:,:,:,f+1:end), ndims(input_tmp));
        input_var(:,:,:,2:end) = input_tmp(:,:,:,1:f);
        input_var = reshape(input_var, h*w*(z_max-z_min), []);
        clear input_tmp;
        %%搜索关键点
        EAMC(input_var(:,2:end),input_var(:,1),2000);
    end
%     cubexyz = [x_min y_min z_min;x_max y_max z_max];
%     action_cube(a) = {cubexyz};
end

function output = EAMC(fsinput, csinput, B)
%myFun - Description
%
% Syntax: output = EAMC(input)
%
% Long description
%%初始化参数
    validpos = find(csinput~=0);

    n = length(validpos(:));
    X=zeros(n+1, n);Y=zeros(n+1, n);
    %Z=zeros(n+1, n);W=zeros(n+1, n);    
    population = zeros(1, n);
    %f(x), c(x), |x|, g(x)
    Xfitness=zeros(n+1, 4);Yfitness=zeros(n+1, 4);
    %Zfitness=zeros(n+1, 4);Wfitness=zeros(n+1, 4);
    %Wfitness(:, 2) = inf;
    %f(x), c(x), |x|, g(x)
    offSpringFit = zeros(1, 4);
    xysame=zeros(1, n+1);
    %zwsame=zeros(1, n+1);
    xysame(1) = 1;
    %zwsame(1) = 1;
    popSize = 1;t = 0;iter1 = 1;
    T = ceil(n*n*10);kn = n;
    figure;
    subplot(4,1,1);
    h1 = animatedline;
    title('原函数值');
    subplot(4,1,2);
    h2 = animatedline;
    title('约束值');
    subplot(4,1,3);
    h3 = animatedline;
    title('符合要求的个体总数');
    subplot(4,1,4);
    h4 = animatedline;
    title('代理函数值');
    tt = 0;
    ri = [];
    resultIndex = -1;
    stable = 0;
%%迭代更新种群
    while t<T
        if iter1 == kn
            iter1 = 0;
            maxValue = -inf;
            for p=1:n+1
                if Yfitness(p, 2)<=B && Yfitness(p, 1)>maxValue
                    maxValue = Yfitness(p, 1);
                    if resultIndex == p
                        stable = stable+1;
                    end
                    resultIndex = p;
                end
            end
            if stable = kn // 4
                break
            end
            Yfitness(resultIndex, :), popSize, resultIndex
            ri = [ri, resultIndex];
            addpoints(h1,tt,Yfitness(resultIndex,1));
            addpoints(h2,tt,Yfitness(resultIndex,2));
            addpoints(h3,tt,popSize);
            addpoints(h4,tt,Yfitness(resultIndex,4));
            drawnow limitrate;
            tt = tt+1;
        end
        iter1 = iter1 + 1;
        s = population(randi(popSize), :);
        offSpring = mutation(s);    
        offSpringFit(1, 1) = FS(offSpring, csinput, validpos);
        offSpringFit(1, 2) = CS(offSpring, csinput, validpos);
        offSpringFit(1, 3) = sum(offSpring(1, :));
        offSpringFit(1, 4) = GS(B, 1, offSpringFit);
        %indice 记录上一次生成的种群
        indice = offSpringFit(1, 3);
        if offSpringFit(1, 3) < 1
            t = t+1;
            continue;
        end
        isadd1 = 0;
        isadd2 = 0;
        %Xfitness代理函数主导的解，X记录此类解，Yfitness原函数主导的解，Y记录此类解
        if offSpringFit(1, 2) <= B
            if offSpringFit(1, 4) >= Xfitness(indice, 4)
                X(indice, :) = offSpring;
                Xfitness(indice, :) = offSpringFit;
                isadd1 = 1;
            end
            if offSpringFit(1, 1) >= Yfitness(indice, 1)
                Y(indice, :) = offSpring;
                Yfitness(indice, :) = offSpringFit;
                isadd2 = 1;
            end
            if isadd1+isadd2 == 2
                xysame(indice) = 1;
            else
                if isadd1+isadd2 == 1
                    xysame(indice) = 0;
                end
            end
        end
        tempSize = 1;
        %计算种群大小
        for i = 2:n+1
            if Xfitness(i, 3) > 0
                %XY种群被同一个后代更新
                if Yfitness(i, 3)>0 && xysame(i)==1
                    tempSize = tempSize+1;
                end
                %XY种群被不同后代更新
                if Yfitness(i, 3)>0 && xysame(i)==0
                    tempSize = tempSize+2;
                end
                if Yfitness(i, 3) == 0
                    tempSize = tempSize+1;
                end
            else
                if Yfitness(i, 3) > 0
                    tempSize = tempSize+1;
                end
            end
        end
        if popSize ~= tempSize
            population = zeros(tempSize, n);
        end
        popSize = tempSize;
        j = 2;
        %更新种群
        for i = 2:n+1
            if Xfitness(i, 3) > 0
                if Yfitness(i, 3) > 0 && xysame(i) == 1
                    population(j, :) = X(i, :);
                    j = j+1;
                end
                if Yfitness(i, 3) > 0 && xysame(i) == 0
                    population(j, :) = X(i, :);
                    j = j+1;
                    population(j, :) = Y(i, :);
                    j = j+1;
                end
                if Yfitness(i, 3) > 0
                    population(j, :) = X(i, :);
                end
            else
                if Yfitness(i, 3) > 0
                    population(j, :) = Y(i, :);
                    j = j+1;
                end
            end
        end
        t = t+1;
%         resultIndex = -1;
%         maxValue = -inf;
%         for p = 1:n+1
%             if Yfitness(p, 2) <= B && Yfitness(p, 1) > maxValue
%                 maxValue = Yfitness(p, 1);
%                 resultIndex = p;
%             end
%         end
%         Yfitness(resultIndex, :), popSize, t
    end
end
%%变异操作
function Y = mutation(X)
    n = length(X);
    rand_rate = 1.0 / n * 2000;
    change = binornd(1, rand_rate, 1, n);
    Y = abs(X-change);
end

function Y = FS(offspring, fsinput, validpos)
%FS - Description
%value of objective function
% Syntax: Y = FS(X)
% Long description
    %SelectedPointFeature
    %CM Covariance Matrix
    [~, frms] = size(fsinput);
    selectedPoint = logical(offspring);
    Y = sum(fsinput(selectedPoint));
%     nums = sum(selectedPoint);
%     sameIndex = [];
%     count = [];
%     if nums > 0
%         for i = 1:frms
%             a = fsinput(validpos(selectedPoint), i);
%             hasSame = false;
%             j = 1;
%             for k = 1:length(sameIndex)
%                 if sum(a == fsinput(validpos(selectedPoint), sameIndex(k))) == nums
%                     count(j) = count(j)+1;
%                     hasSame = true;
%                 end
%                 j = j+1;
%             end
%             if hasSame == false
%                 sameIndex = [sameIndex, i];
%                 count = [count, 1];
%             end
%         end
%     end
%     Y = 0.0;
%     for i = 1:length(count)
%         prob = 1.0*count(i)/frms;
%         Y = Y - prob*log2(prob);
%     end
end

function Y = CS(offspring, csinput, validpos)
%CS - Description
%value of cost function
% SynYutpCSmXCS)
%
% Long descrivalue of cost functionption
    selectedPoint = logical(offspring);
    selectedPoint = csinput(validpos(selectedPoint));
    if sum(selectedPoint<0) > 0
        Y = 0;
    else
        Y = sum(offspring);
    end
end
%%代理函数值
function Y = GS(B, alpha, offSpringFit)
%GS - Description
%value of surrogate function
% SynYutpCSmXCS)
%
% Long descrivalue of cost functionption
    if offSpringFit(1, 3) >= 1
        Y = 1.0*offSpringFit(1, 1)/(1.0-(1.0/exp(alpha*offSpringFit(1, 2)/B)));
    else
        Y = 0;
    end
end