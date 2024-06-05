function Positions=initialization_CGO(SearchAgents_no,dim,ub,lb)
    Boundary_no= size(ub,2); % Number of borders
    % If the boundaries of all variables are equal and the user enters a single number for ub and lb
    if Boundary_no==1
        Positions=rand(SearchAgents_no,dim).*(ub-lb)+lb;
    end
    % If each variable has a different lb and ub
    if Boundary_no>1 % If the number of columns on the upper boundary is greater than 1
        for i=1:dim
            ub_i=ub(i);
            lb_i=lb(i);
            Positions(:,i)=rand(SearchAgents_no,1).*(ub_i-lb_i)+lb_i;
        end
    end
end