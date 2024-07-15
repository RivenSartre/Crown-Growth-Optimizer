function [Best_pos,Best_fitness,Convergence_curve]=CGO(N,Max_iter,lb,ub,dim,fobj)

    if(max(size(ub)) == 1)
       ub = ub.*ones(1,dim);
       lb = lb.*ones(1,dim);  
    end
    
    %% initialization X
    X = initialization_CGO(N,dim,ub,lb);      % initialization
    
    Best_pos = zeros(1,dim);                  % Best position
    Best_fitness = inf;                       % Best fitness
    Objective_values = zeros(1,size(X,1));    % initialize fitness
    
    Convergence_curve = [];                   % curve
    N_half = floor(N*0.5);                    
                            
    for i = 1:size(X,1)
        Objective_values(1,i) = fobj(X(i,:));
    end
    
    [~,idx1] = sort(Objective_values);
    first_best = X(idx1(1),:);
    second_best = X(idx1(2),:);
    third_best = X(idx1(3),:);
    sum1 = 0;
    for i = 1:N_half
        sum1 = sum1+X(idx1(i),:);
    end
    half_best_mean = sum1/N_half;
    
    Pool = [];
    Pool(1,:) = first_best;
    Pool(2,:) = second_best;
    Pool(3,:) = third_best;
    Pool(4,:) = half_best_mean;
    Pool(5,:) = mean(X);
    
    Convergence_curve(1) = fobj( Pool(1,:));
    
    for i = 1:N
        index(i) = i;
    end
    
    Ns = floor(N*0.5);          % Sprouting broach
    Ng = N-Ns;                    % Growing broach
    b  = 0.8;
    Dis = (ub(1)-lb(1)) * 0.001;
    alpha = 0.2;
    M(1) = 1;
    distances = zeros(N,1);

    X_temp = zeros(N,dim);

    %% main loop
    l = 2; 
    while l <= Max_iter
        GR = randn(N,dim);          
        M(l) = 1-0.85/( 1 + exp(-10*b*(2*l/(Max_iter)- 1)) );
    
        %% X_centroid
        for j = 1:dim
            sum1 = 0;
            for i = 1:N
                sum1 = sum1+X(i,j);
            end
            X_centroid(j) = sum1/N;
        end
        
        %% growing and Sprouting
        
        index_s = randperm(N,Ns);
        index_g = setdiff(index,index_s);
        
        % Sprouting
        for i = 1:Ns
            r1 = rand;
            k1 = randperm(5,1);
            for j = 1:size(X,2) 
                X(index_s(i),j) = Pool(k1,j)+GR(index_s(i),j)*(r1*(Best_pos(j)-X(index_s(i),j))+(1-r1)*(X_centroid(j)-X(index_s(i),j)));
            end
        end
        
        % Growing
        if Ng >= 1
            for i = 1:Ng
                r2 = 2*rand-1;
                for j = 1:size(X,2) 
                    X_temp(index_g(i),j) = X(index_g(i),j)+M(l)*GR(index_g(i),j)*(r2*(Best_pos(j)-X(index_g(i),j))+(1-r2)*(X_centroid(j)-X(index_g(i),j)));
                end
            end

            for k=1:Ng
                vec=zeros(N,1);
                distances(k) = sqrt(sum((X_temp(index_g(i),j) - X(index_g(k), :)).^2));
                if distances(k) < Dis
                    vec(k) = X(index_g(k), j) - X_temp(index_g(i),j);
                end
                Freq = sum(vec);
            end
            X(index_g(i),j) =  X_temp(index_g(i),j) + Freq*alpha;

        end
        
    
        %% Check to see if the solution is beyond the scope of the search and bring it back
        for i = 1:size(X,1)
            for j = 1:dim
                if X(i,j) > ub(j)
                    X(i,j) = ub(j);
                end
                if X(i,j) < lb(j)
                    X(i,j) = lb(j);
                end
            end
            
            % Calculate the Objective values
            Objective_values(1,i) = fobj(X(i,:));
            % If there is a better solution, update the Best_pos
            if Objective_values(1,i) < Best_fitness
                Best_pos = X(i,:);
                Best_fitness = Objective_values(1,i);
            end
        end
        
        if Ns < N && mod(l,Max_iter/100)==0
            Ns = Ns+1;
            Ng = Ng-1;
            [~, sortedIndices] = sort(Objective_values);
            indices = sortedIndices(1:floor(N*0.382));
            for k = 1:length(indices)
                X(indices(k),:) = initialization_CGO(1,dim,ub,lb);      % initialization
            end
        end  
       
        %% update the pool
        [~,idx1] = sort(Objective_values);
        second_best = X(idx1(2),:);
        third_best = X(idx1(3),:);
        half = X(idx1(N_half),:);
        % sum1 = 0;
        % for i = 1:N_half
        %     sum1 = sum1+X(idx1(i),:);
        % end
        Pool(1,:) = Best_pos;
        Pool(2,:) = second_best;
        Pool(3,:) = third_best;
        Pool(4,:) = half;
        Pool(5,:) = mean(X);

        Convergence_curve(l) = Best_fitness;
        l = l+1;
    
    end

end
