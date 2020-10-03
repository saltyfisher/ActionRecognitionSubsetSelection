clear;
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
    tfeature_path = ['E:\KangLi\HOPC Matlab Code\HOPC Matlab Code\tHOPC3Dfeatures\a' action_id{a} '*.mat'];
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
        y_min=min(r);x_min=min(c);y_max=max(r);x_max=mod(max(c), 320);
        
        sfeature = load(fullfile(sfeature_list(i).folder, sfeature_list(i).name));
        tfeature = load(fullfile(tfeature_list(i).folder, tfeature_list(i).name));
        sfeature = sfeature.Q;tfeature = tfeature.Q;
        sf = reshape(sfeature, [], 30);
        tf = reshape(tfeature, [], 30);

        t1 = sum((sf - tf).^2, 2);
        t2 = sum(sf + tf, 2) + 1e-5;
        t3 = t1./t2 - 1;
        t3 = reshape(t3, size(t1));
        t3 = floor(t3);
        sf = floor(sf);
        
        [h, w, f] = size(depth);
        temp = cat(length(size(sf)), sf, t3);
        temp = reshape(temp, h*w, f, []);
        [r, c, l] = find(depth~=0);
        r=round(r);c=round(c);l=round(l);
%         c = mod(c, w);
        id(j) = {[[r,c] depth(r,c,l)-z_min]};
        input_var = zeros(h, w, z_max-z_min, f, size(temp, ndims(temp)));
        for k = 1:length(depth)
            input_var(:,:,depth(r,c,k),k,:) = temp(r,c,k,:);
        end
        input_var = input_var(y_min:y_max, x:min_x_max, :, :, :);
        alldepth = {depth};
    end
    cubexyz = [x_min y_min z_min;x_max y_max z_max];
    action_cube(a) = {cubexyz};
    %%转换到决策空间中
    %按时间顺序记录决策空间中的索引
    id = cell(1, length(depth));
    for j = 1:length(alldepth)
        depth = alldepth{j};
        depth = cell2mat(depth);
        depth = reshape(depth, 240, 320, []);
        
        sfeature = load(fullfile(sfeature_list(j).folder, sfeature_list(j).name));
        tfeature = load(fullfile(tfeature_list(j).folder, tfeature_list(j).name));
        sfeature = sfeature.Q;tfeature = tfeature.Q;
        sf = reshape(sfeature, [], 30);
        tf = reshape(tfeature, [], 30);

        t1 = sum((sf - tf).^2, 2);
        t2 = sum(sf + tf, 2) + 1e-5;
        t3 = t1./t2 - 1;
        t3 = reshape(t3, size(t1));
        t3 = floor(t3);
        sf = floor(sf);
        
        [h, w, f] = size(depth);
        temp = cat(length(size(sf)), sf, t3);
        temp = reshape(temp, h*w, f, []);
        [r, c, l] = find(depth~=0);
        id(j) = {[[r,c] depth([r,c,l])-z_min]};
        input_var = zeros(h, w, z_max-z_min, f, size(temp, ndims(temp)));
        for k = 1:length(depth)
            input_var(:,:,depth(r,c,k),k,:) = temp(r,c,k,:);
        end
        
    end
end

function output = EAMC(fsinput, csinput, B)
%myFun - Description
%
% Syntax: output = EAMC(input)
%
% Long description
%%初始化参数
    g = 10;
    validpos = int(find(csinput>g));

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
    T = ceil(n*n*10);kn = n*n;
%%迭代更新种群
    while t<T
        if iter1 == kn
            iter1 = 0;
            resultIndex = -1;
            maxValue = -inf;
            for p=1:n+1
                if Yfitness(p, 2)<=B && Yfitness(p, 1)>maxValue
                    maxValue = Yfitness(p, 1);
                    resultIndex = p;
                end
            end
            Yfitness(resultIndex, :), popSize
        end
        iter1 = iter1 + 1;
        s = population(randi(popSize), :);
        offSpring = mutation(s);    
        offSpringFit(1, 1) = FS(offSpring, fsinput, validpos);
        offSpringFit(1, 2) = CS(offSpring, csinput, validpos);
        offSpringFit(1, 3) = sum(offSpring(1, :));
        offSpringFit(1, 4) = GS(B, 1.0, offSpringFit, fsinput, validpos);
        %indice 记录上一次生成的种群
        indice = int(offSpringFit(1, 3));
        if offSpringFit(1, 3) < 1
            t = t+1;
            continue;
        end
        isadd1 = 0;
        isadd2 = 0;
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
        %count the population size, 0也是种群之一
        tempSize = 1;
        for i = 2:n+1
            if Xfitness(i, 3) > 0
                if Yfitness(i, 3)>0 && xysame(i)==1
                    tempSize = tempSize+1;
                end
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
        %融合X，Y，Z，W
        for i = 2:n+1
            if Xfitness(i, 3) > 0
                if Yfitness(i, 3) > 0 && xysame(i) == 1
                    population(j, :) = X(i, :);
                    j = j+1;
                end
                if Yfitness(o, 3) > 0 && xysame(i) == 0
                    population(j, :) = X(i, :);
                    j = j+1;
                    population(j, :) = Y(i, :)
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
    resultIndex = -1;
    maxValue = -inf;
    for p = 1:n+1
        if Yfitness(p, 2) <= B && Yfitness(p, 1) > maxValue
            maxValue = Yfitness(p, 1);
            resultIndex = p;
        end
    end
    Yfitness(resultIndex, :), popSize
    end
end
%%变异操作
function Y = mutation(X)
    n = length(X);
    rand_rate = 1.0 / n;
    change = int8(binornd(1, rand_rate, 1, n));
    Y = abs(X-change);
end

function Y = FS(offspring, fsinput, validpos)
%FS - Description
%value of objective function
% Syntax: Y = FS(X)
% Long description
    %SelectedPointFeature
    %CM Covariance Matrix
    [~, frms] = size(input);
    selectedPoint = logical(offspring);
    nums = sum(selectedPoint);
    sameIndex = [];
    count = [];
    if nums > 0
        for i = 1:frms
            a = fsinput(validpos(selectedPoint));
            hasSame = false;
            j = 0;
            for k = 1:length(sameIndex)
                if sum(a == fsinput(validpos(selectedPoint), sameIndex(k))) == nums
                    count(j) = count(j)+1;
                    hasSame = true;
                end
                j = j+1;
            end
            if hasSame == false
                sameIndex = [sameIndex, i];
                count = [count, 1];
            end
        end
    end
    for i = 1:length(count)
        prob = 1.0*count(i)/frms;
        Y = Y - prob*log2(prob);
    end
end

function Y = CS(offspring, csinput, validpos)
%CS - Description
%value of cost function
% SynYutpCSmXCS)
%
% Long descrivalue of cost functionption
    selectedPoint = csinput(validpos);
    Y = sum(selectedPoint);

end
%%代理函数值
function Y = GS(B, alpha, offSpringFit, input)
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