classdef SubsetSelection < handle
    properties
        depth
        pointcloud
        
        fsinput
        csinput
        frms
        %ÄŁĐÍłŹ˛ÎĘý
        B;k
        alpha
        %˝řťŻËăˇ¨ąäÁż
        X;Y;population;Xfitness;Yfitness
    end
    methods
        function obj = SubsetSelection(B,alpha,k)
           %alphaĘÇ´úŔíşŻĘý˛ÎĘý
            obj.B = B;
            obj.k = k;
            obj.alpha = alpha;
        end
        function GenerateInput(obj,depth,pointcloud,sfeature,tfeature)
            z_max = max(max(max(depth)));
            z_min = min(min(depth(depth~=0)));
            [r,c,~] = find(depth~=0);
            c = rem((c-1),320)+1;
            y_min=min(r);x_min=min(c);y_max=max(r);x_max=max(c);

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
            %//TODO:żÉŇÔźĆËăĎňÁż¸÷¸öÎŹśČľÄÍłźĆÁżŁŹźÓČ¨ÇóşÍ
            sf = sum(sf, ndims(sf));
            %%×Şťťľ˝žö˛ßżŐźäÖĐ
            %°´ĘąźäËłĐňźÇÂźžö˛ßżŐźäÖĐľÄË÷Ňý
            [h, w, f] = size(depth);
            depth = depth-z_min;
            %Ç°fÎŹĘÇÔ­şŻĘýąäÁżŁŹşófÎŹĘÇÔźĘřşŻĘýąäÁż
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
            obj.fsinput = input_var(:,2:end);
            obj.csinput = input_var(:,1);
        end
        
        function Y = mutation(obj,X)
            n = length(X);
            rand_rate = 1.0 / n * obj.k;
            change = binornd(1, rand_rate, 1, n);
            Y = abs(X-change);
        end
        %%ČŤžÖËŃË÷
        %kąíĘž×ÓźŻľÄšćÄŁÉĎ˝ç?
        function Y = GlobalSearch(obj,pos,k,T,validpos)
            pos = mutation(pos,k);
            t=0;
            %vźÇÂźłőĘź˝âŁŹuźÇÂźąťŃĄÖĐľÄłőĘź˝â
            %W¸ůžÝśÔÄżąęşŻĘýľÄšąĎ×´Ó´óľ˝ĐĄźÇÂź˝â
            v = pos;
            u = zeros(size(v));
            W = zeros(2, sum(pos(:)));
            pos_ind = find(W==1);
            origin_fs = FS(pos,validpos);
            %łőĘźťŻW
            for i = 1:sum(pos(:))
                W(1,i) =  pos_ind(i);
                W(2,i) = origin_fs-FS(pos((pos_ind(i))=0),validpos);
            end
            while t<T
                %ŃĄłöŇŞ¸üĐÂľÄ˝â
                [~,p] = max(W(2,:));
                p = W(1,p);
                u(p) = 1;

            end
        end
        %%ĺą?é¨ćç´˘çŽćł?
        function Y = LocalSearch(obj)
            1234;
        end
        function Y = FS(obj,offspring, validpos)
            %FS - Description
            %value of objective function
            % Syntax: Y = FS(X)
            % Long description
                %SelectedPointFeature
                %CM Covariance Matrix
            [~, obj.frms] = size(obj.fsinput);
            selectedPoint = logical(offspring);
            Y = sum(obj.fsinput(selectedPoint));
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
        function Y = CS(obj,offspring,validpos)
        %CS - Description
        %value of cost function
        % SynYutpCSmXCS)
        %
        % Long descrivalue of cost functionption
            selectedPoint = logical(offspring);
            selectedPoint = obj.csinput(validpos(selectedPoint));
            if sum(selectedPoint<0) > 0
                Y = 0;
            else
                Y = sum(offspring);
            end
        end
        %%´úŔíşŻĘý?
        function Y = GS(obj,offSpringFit)
        %GS - Description
        %value of surrogate function
        % SynYutpCSmXCS)
        %
        % Long descrivalue of cost functionption
            if offSpringFit(1, 3) >= 1
                Y = 1.0*offSpringFit(1, 1)/(1.0-(1.0/exp(obj.alpha*offSpringFit(1, 2)/obj.B)));
            else
                Y = 0;
            end
        end
        function output = EAMC(obj)
            %myFun - Description
            %
            % Syntax: output = EAMC(input)
            %
            % Long description
            %%ĺĺ§ĺĺć?
            validpos = find(obj.csinput~=0);
        
            n = length(validpos(:));
            obj.X=zeros(n+1, n);obj.Y=zeros(n+1, n);
            %Z=zeros(n+1, n);W=zeros(n+1, n);    
            obj.population = zeros(1, n);
            %f(x), c(x), |x|, g(x)
            obj.Xfitness=zeros(n+1, 4);obj.Yfitness=zeros(n+1, 4);
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
            title('Ô­şŻĘýÖľ');
            subplot(4,1,2);
            h2 = animatedline;
            title('ÔźĘřÖľ');
            subplot(4,1,3);
            h3 = animatedline;
            title('ˇűşĎŇŞÇóľÄ¸öĚĺ×ÜĘý');
            subplot(4,1,4);
            h4 = animatedline;
            title('´úŔíşŻĘýÖľ');
            tt = 0;
            ri = [];
            resultIndex = -1;
            maxValue = -inf;
            count = 0;
            %%čż­äťŁć´ć°ç§çž¤
            while t<T
                if iter1 == kn
                    iter1 = 0;
                    stable = true;
                    for p=1:n+1
                        if obj.Yfitness(p, 2)<=obj.B && obj.Yfitness(p, 1)>maxValue
                            maxValue = obj.Yfitness(p, 1);
                            resultIndex = p;
                            stable = false;
                            count = 0;
                        end
                    end
                    if stable == true
                        count = count+1;
                    end
                    if (count-kn/4)>10
                        break
                    end
                    obj.Yfitness(resultIndex, :), popSize, resultIndex
                    ri = [ri, resultIndex];
                    addpoints(h1,tt,obj.Yfitness(resultIndex,1));
                    addpoints(h2,tt,obj.Yfitness(resultIndex,2));
                    addpoints(h3,tt,popSize);
                    addpoints(h4,tt,obj.Yfitness(resultIndex,4));
                    drawnow limitrate;
                    tt = tt+1;
                end
                iter1 = iter1 + 1;
                s = obj.population(randi(popSize), :);
                offSpring = obj.mutation(s);    
                offSpringFit(1, 1) = obj.FS(offSpring, obj.csinput, validpos);
                offSpringFit(1, 2) = obj.CS(offSpring, obj.csinput, validpos);
                offSpringFit(1, 3) = sum(offSpring(1, :));
                offSpringFit(1, 4) = obj.GS(obj.B, 1, offSpringFit);
                %indice čŽ°ĺ˝ä¸ä¸ćŹĄçćçç§çž¤
                indice = offSpringFit(1, 3);
                if offSpringFit(1, 3) < 1
                    t = t+1;
                    continue;
                end
                isadd1 = 0;
                isadd2 = 0;
                %XfitnessäťŁçĺ˝ć°ä¸ťĺŻźçč§ŁďźXčŽ°ĺ˝ć­¤çąťč§ŁďźYfitnessĺĺ˝ć°ä¸ťĺŻźçč§ŁďźYčŽ°ĺ˝ć­¤çąťč§?
                if offSpringFit(1, 2) <= obj.B
                    if offSpringFit(1, 4) >= obj.Xfitness(indice, 4)
                        obj.X(indice, :) = offSpring;
                        obj.Xfitness(indice, :) = offSpringFit;
                        isadd1 = 1;
                    end
                    if offSpringFit(1, 1) >= obj.Yfitness(indice, 1)
                        obj.Y(indice, :) = offSpring;
                        obj.Yfitness(indice, :) = offSpringFit;
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
                %čŽĄçŽç§çž¤ĺ¤§ĺ°
                for i = 2:n+1
                    if obj.Xfitness(i, 3) > 0
                        %XYç§çž¤č˘Ťĺä¸?ä¸ŞĺäťŁć´ć?
                        if obj.Yfitness(i, 3)>0 && xysame(i)==1
                            tempSize = tempSize+1;
                        end
                        %XYç§çž¤č˘Ťä¸ĺĺäťŁć´ć?
                        if obj.Yfitness(i, 3)>0 && xysame(i)==0
                            tempSize = tempSize+2;
                        end
                        if obj.Yfitness(i, 3) == 0
                            tempSize = tempSize+1;
                        end
                    else
                        if obj.Yfitness(i, 3) > 0
                            tempSize = tempSize+1;
                        end
                    end
                end
                if popSize ~= tempSize
                    obj.population = zeros(tempSize, n);
                end
                popSize = tempSize;
                j = 2;
                %ć´ć°ç§çž¤
                for i = 2:n+1
                    if obj.Xfitness(i, 3) > 0
                        if obj.Yfitness(i, 3) > 0 && xysame(i) == 1
                            obj.population(j, :) = obj.X(i, :);
                            j = j+1;
                        end
                        if obj.Yfitness(i, 3) > 0 && xysame(i) == 0
                            obj.population(j, :) = obj.X(i, :);
                            j = j+1;
                            obj.population(j, :) = obj.Y(i, :);
                            j = j+1;
                        end
                        if obj.Yfitness(i, 3) > 0
                            obj.population(j, :) = obj.X(i, :);
                        end
                    else
                        if obj.Yfitness(i, 3) > 0
                            obj.population(j, :) = obj.Y(i, :);
                            j = j+1;
                        end
                    end
                end
                t = t+1;
            end
        end
    end
end