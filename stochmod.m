function [yr, tr] = stochmod(p, numsims, x)
totals=cell(1,numsims);  % Stores total infected tomatoes at each time point
p=abs(p); % Model parameters
c=cell(1,numsims); % Stores all of the time points

Output_val = zeros(numsims, length(x));

% Loop for 10,000 simulation
for q=1:numsims
 
    
X = zeros(4,1); % state vector 

X(1) = 800-(p(3));      % initial number of susceptible beets
X(2) =150;              % initial number of susceptible toms
X(3) = (p(3));          % initial number of infective beet
X(4)= 0 ;               % initial number of Infective toms

 
V = [-1 0; 0 -1; 1 0; 0 1];   % Transition matrix 
time=8; % Initial time
c{q}=0;
totals{q}=0; 

tfinal = 70;    % total time (days)

 % Extracting number of infected tomatoes as desired time points 
    while time < tfinal
        a(1) = p(1)*X(1)*X(3); % Rate of infection in beets
        a(2) =p(2)*X(2)*X(3); % Rate of infection in toms
        asum = sum(a);

        j = find(rand<cumsum(a/asum), 1 );
        tau = abs(log(1/rand)/asum);  % time step 

        X = X+V(:,j);
        time = time+tau;  
        
    c{q}=[c{q}, time];  % updating time vector
    
    totals{q}=[totals{q}, X(4)];  % Updating state vector  
        
    end
    
    if tau==inf
            X=[0,0,0,0];
            error('Tau is infinite')
    end

for kk = 1:length(x)
    in=find(abs(c{q}-x(kk))==min(abs(c{q}-x(kk))));
    timestamp = c{q};
    output = totals{q};
    if  (timestamp(in) <= x(kk))
        Output_val(q, kk) = output(in);
    else
        Output_val (q, kk) = output(in-1); 
    end
end


end

yr=Output_val; % returning number of infected tomatos
tr=x;  % returning timestamps
