%covariance matrix
global CM;

FeatureSet = [];
function output = EAMC(input, B)
%myFun - Description
%
% Syntax: output = EAMC(input)
%
% Long description
%%初始化参数
    n = length(input);
    X=zeros(n+1, n);Y=zeros(n+1, n);Z=zeros(n+1, n);W=zeros(n+1, n);    
    population = zeros(1, n);
    %f(x), c(x), |x|, g(x)
    Xfitness=zeros(n+1, 4);Yfitness=zeros(n+1, 4);Zfitness=zeros(n+1, 4);Wfitness=zeros(n+1, 4);
    Wfitness(:, 2) = inf;
    %f(x), c(x), |x|, g(x)
    offSpringFit = zeros(1, 4);
    xysame=zeros(1, n+1);zwsame=zeros(1, n+1);
    xysame(1) = 1;zwsame(1) = 1;
    popSize = 1;t = 0;iter1 = 1;
    T = ceil(n*n*10);kn = n*n;
%%迭代更新种群
    while t<T
        if iter1 == kn
            iter1 = 0;
            resultIndex = -1;
            maxValue = -inf;
            for p=1:n+1
                if Yfitness(p, 2)<=B and Yfitness(p, 1)>maxValue
                    maxValue = Yfitness(p, 1);
                    resultIndex = p;
                end
            end
            Yfitness(resultIndex, :), popSize
        end
        iter1 = iter1 + 1;
        s = population(randi(popSize), :);
        offSpring = mutation(s);    
        offSpringFit(1, 1) = FS(offSpring, input);
        offSpringFit(1, 2) = CS(offSpring, input);
        offSpringFit(1, 3) = sum(offSpring(1, :));
        offSpringFit(1, 4) = GS(B, 1.0, offSpringFit, input);
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
                if Yfitness(i, 3)>0 and xysame(i)==1
                    tempSize = tempSize+1;
                end
                if Yfitness(i, 3)>0 and xysame(i)==0
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
                if Yfitness(i, 3) > 0 and xysame(i) == 1
                    population(j, :) = X(i, :);
                    j = j+1;
                end
                if Yfitness(o, 3) > 0 and xysame(i) == 0
                    population(j, :) = X(i, :);
                    j = j+1;
                    population(j, :) = Y(i, :):
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

function Y = FS(offspring, input)
%FS - Description
%value of objective function
% Syntax: Y = FS(X)
% Long description
    %SelectedPointFeature
    %CM Covariance Matrix
    SPFeature = logical(offspring);
    SPCM = CM(SPFeature, :);
    SPCM = SPCM(:, SPFeature);
    %UnselectedPointFeature
    USPFeature = ~logical(offspring);
    USPCM = CM(USPFeature, :);
    USPCM = USPCM(:, USPFeature);

    %Entropy
    [n, attrNum] = size(input);
    SPH = mvnpdf(input, zeros(1, sum(offspring)), SPCM);
    USPH = mvnpdf(input, zeros(1, n-sum(offspring)), USPCM);
    % VH = mvnpdf(input, zeros(1, n), CM);
    % // TODO: 会爆内存
    SPH = sum(SPH.*log(SPH));
    USPH = sum(USPH.*log(USPH));
    % VH = sum(VH.*log(VH));

    Y = SPH+USPH;
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