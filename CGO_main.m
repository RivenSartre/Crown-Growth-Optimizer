clear all 
close all
clc

format longG
format compact

%% initialization
N=10;               % agents
Max_iter=1000;       % Max iteration
dim = 10;            % dimension

lb=-100;             % lower bound
ub=100;              % upper bound
Max_test=20;

Result=zeros(30,4);


%% testing CGO on Function 1 to 30 in CEC2017
for j=1:30
    Function_name=j;    % number of Function
    fobj = @(x) cec17_func(x',Function_name);
    
    WaitbarInter = Max_test / 100;  
    tic
    h = waitbar( 0, ['completed:0%   running...time:', num2str(toc),'   Function:', num2str(j)]);
    
    for i=1:Max_test
        [Best_pos(i,:),Best_fitness(i),CGO_curve(i,:)]=CGO(N,Max_iter,lb,ub,dim,fobj);  % run
    
        if mod(i, WaitbarInter) == 0
            waitbar(i / Max_test, h, ['completed:' num2str(i / Max_test * 100) ...
            '%   running...time:', num2str(toc),'/',num2str(toc/(i / Max_test)), '   Function:', num2str(j)])
        end
    end
    close(h)
    Result(j,1)=min(Best_fitness);
    Result(j,2)=mean(Best_fitness);
    Result(j,3)=max(Best_fitness);
    Result(j,4)=std(Best_fitness);
end


%% testing CGO on Function 1 to 48 in CEC2020-RW
Max_test=48;

WaitbarInter = Max_test / 100;  % 一个和进度条有关的参数
tic
a=1;
h = waitbar( 0, ['completed:0%   running...time:', num2str(toc),'   Function:', num2str(a)]);

for k = 1:Max_test
    [lb,ub,dim,fobj] = Problem_Define(k);
    [~,CGO_fitness(k,:),CGO_curve(k,:)] = CGO(N,Max_iter,lb,ub,dim,fobj);  
        
    if mod(k, WaitbarInter) == 0
            waitbar(k / Max_test, h, ['completed:' num2str(k / Max_test * 100) ...
            '%   running...time:', num2str(toc),'/',num2str(toc/(k / Max_test)), '    Function:', num2str(k)])
    end
end

close(h);


