function [Best_pos,Best_fitness,Convergence_curve]=CGO(N,Max_iter,lb,ub,dim,fobj)

    if(max(size(ub)) == 1)
       ub = ub.*ones(1,dim);
       lb = lb.*ones(1,dim);  
    end
    
    %% initialization X
    X = initialization_CGO(N,dim,ub,lb);      % initialization, Eq.1
    
    Best_fitness = inf;                       % Best fitness
    Objective_values = zeros(1,size(X,1));    % initialize fitness
    
    Convergence_curve = [];                   % Convergence_curve
    N_half = floor(N*0.5);                    
                            
    for i = 1:size(X,1)
        Objective_values(1,i) = fobj(X(i,:));
    end
    
    
    % Elite pool, Eq.9
    [~,idx1] = sort(Objective_values);
    Best_pos = X(idx1(1),:);
    second_best = X(idx1(2),:);
    third_best = X(idx1(3),:);
    
    sum1 = 0;
    for i = 1:N_half
        sum1 = sum1+X(idx1(i),:);
    end
    half_best_mean = sum1/N_half;           % Half centroid position, Eq.10

    Pool = [];                              % Elite pool, Eq.9
    Pool(1,:) = Best_pos;                   % Best position
    Pool(2,:) = second_best;                % 2nd best position
    Pool(3,:) = third_best;                 % 3rd best position
    Pool(4,:) = half_best_mean;             % Half centroid position
    Pool(5,:) = mean(X);                    % Centroid position
    
    Convergence_curve(1) = fobj( Pool(1,:));
    
    for i = 1:N
        index(i) = i;
    end
    
    % Parameters
    Ns = floor(N*0.5);                      % Sprouting broach
    Ng = N-Ns;                              % Growing broach
    b  = 0.8;
    Dis = (ub(1)-lb(1)) * 0.01;
    alpha = 2;
    V(1) = 1;
    tcut = Max_iter/N/2;

    X_temp = zeros(N,dim);

    %% main loop
    l = 2; 
    while l <= Max_iter
        GR = randn(N,dim);                  % Gaussian random number, Eq.4
        V(l) = 1-0.85/( 1 + exp(-10*b*(2*l/(Max_iter)- 1)) );  % Branch growth velocity, Eq.2
    
        %% X_centroid
        X_centroid = mean(X);               % Eq.6

        %% growing and Sprouting
        index_s = randperm(N,Ns);
        index_g = setdiff(index,index_s);
        
        % Sprouting
        for i = 1:Ns
            r1 = rand;
            k1 = randperm(5,1);
            for j = 1:size(X,2) 
                % Sprouting vectors, Eq.11
                VecS1 = Best_pos(j) - X(index_s(i),j);
                VecS2 = X(index_s(i),j) - X_centroid(j);
                % Sprouting stage, Eq.12
                X(index_s(i),j) = Pool(k1,j)+GR(index_s(i),j)*( r1*VecS1 + (1-r1)*VecS2 );
            end
        end
        
        % Growing
        if Ng >= 1
            for i = 1:Ng
                r2 = rand;
                for j = 1:size(X,2) 
                    % Growing vectors, Eq.3
                    VecG1 = Best_pos(j) - X(index_g(i),j);
                    VecG2 = X(index_g(i),j) - X_centroid(j);
                    % Growing stage, Eq.5
                    X_temp(index_g(i),j) = X(index_g(i),j)+V(l)*GR(index_g(i),j)*( r2*VecG1 + (1-r2)*VecG2 );
                end
            end
            
            % Calculate distance and repulsion, Eq.7
            vec=zeros(N,dim);
            distances=zeros(N,1);
            for k=1:Ng
                for i=1:N
                    distances(i) = sqrt(    sum( (X_temp(index_g(k), :) - X(i, :) ).^2  )     );
                    if distances(k) < Dis
                        vec(i,:) = X_temp(index_g(k), :) - X(i, :);
                    end
                end
                Freq = sum(vec);
                X(index_g(k),:) = X_temp(index_g(k),:) + Freq*alpha;   % Eq.8
            end
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
        
        if Ns < N && mod(l,tcut)==0
            Ns = Ns+1;
            Ng = Ng-1;

            [~, sortedIndices] = sort(Objective_values,'descend');
            indices = sortedIndices(1:floor(N*0.382));
            for k = 1:length(indices)
                X(indices(k),:) = initialization_CGO(1,dim,ub,lb);      % initialization
            end
        end  
       
        %% update the pool
        [~,idx1] = sort(Objective_values);
        second_best = X(idx1(2),:);
        third_best = X(idx1(3),:);
    
        sum1 = 0;
        for i = 1:N_half
            sum1 = sum1+X(idx1(i),:);
        end
        half_best_mean = sum1/N_half;
        
        Pool(1,:) = Best_pos;
        Pool(2,:) = second_best;
        Pool(3,:) = third_best;
        Pool(4,:) = half_best_mean;
        Pool(5,:) = mean(X);

        Best_fitness = Objective_values(1,idx1(1));
        Convergence_curve(l) = Best_fitness;
        l = l+1;
    
    end

end
