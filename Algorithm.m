%covariance matrix
warning off;
slist = dir('E:\KangLi\HOPC Matlab Code\HOPC Matlab Code\sHOPC3Dfeatures\*.mat');
tlist = dir('E:\KangLi\HOPC Matlab Code\HOPC Matlab Code\tHOPC3Dfeatures\*.mat');

for i = 1:length(tlist)
    filename = [slist(i).folder, '\', slist(i).name];
    sfeature = load(filename);
    sf = reshape(sfeature.Q, [], 30);
    filename = [tlist(i).folder, '\', tlist(i).name];
    tfeature = load(filename);
    tf = reshape(tfeature.Q, [], 30);
    
    t1 = sum((sf - tf).^2, 2);
    t2 = sum(sf + tf, 2) + 1e-5;
    t3 = t1./t2;
    
    STK = EAMC(t3, sf, 200);
end
function output = EAMC(fsinput, csinput, B)
%myFun - Description
%
% Syntax: output = EAMC(input)
%
% Long description
%%初始化参数
    g = 10;
    validpos = find(csinput>g);

    n = length(validpos);
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
            if Xfitness(i, 3) > 0:
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
                if Yfitness(i, 3) > 0:
                    population(j, :) = Y(i, :);
                    j = j+1;
                end
            end
        end
        t = t+1;
    resultIndex = -1;
    maxValue = -inf;
    for p = 1:n+1
        if Yfitness(p, 2) <= B and Yfitness(p, 1) > maxValue
            maxValue = Yfitness(p, 1);
            resultIndex = p;
        end
    end
    Yfitness(resultIndex, :), popSize
    end
end

function Y = mutation(X)
    n = length(X);
    rand_rate = 1.0 / n;
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
    [~, frms] = size(input);
    selectedPoint = logical(offspring);
    nums = sum(selectedPoint);
    sameIndex = [];
    count = [];
    tempSum = 0.0;
    if nums > 0
        for i = 1:frms
            a = fsinput(validpos(selectedPoint));
            hasSame = false;
            j = 0;
            for k = 1:length(sameIndex)
                if sum(a == fsinput(validpos(selectedPoint), sameIndex(k))) == length
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
        tempSum = tempSum - prob*log2(prob);
    end
end

function Y = CS(offspring, input)
%CS - Description
%value of cost function
% SynYutpCSmXCS)
%
% Long descrivalue of cost functionption
    SPFeature = input(logical(offspring));
    SPFeature(SPFeature<0) = 0;

    Y = sum(SPFeature);
end

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