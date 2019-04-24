function [ E,I,Err ] = State_to_Class( a,b,ProbStop, N,kB,kL,TB,TL,ticker, maxAge, points)

%State_to_Class Calculates by MC integration expected number of children in
%the age class with endpoints (a,b) in a household of size N at ticker
%point ticker. I gives the individual terms which when summed give E,
%whilst Err is the error (variance) involved in summing
%   Interval_Prob works out probability of getting age intervals from y_i
%   to y_N after we have integrated the other stuff, which goes to 1.
%   Summing each probability of being in interval gives us expected number.

Nc = N-2; % Nc gives the number of children in the household
if Nc==0
    E = 0; I = 0; Err = 0; % No kids, everything is zero
else
    Nc_max = find(ProbStop==1,1)-2;
    if Nc>(ticker<=kB+kL+kB)*(Nc_max-(ticker>kB+kL))
        E = 0; I = 0; Err = 0; % If we have more children than allowed
        disp(['Demographic state impossible under parameters, #children=' num2str(Nc) ]);
    elseif ticker>kB+kL+kB
        E = 0; I = 0; Err = 0; % If all children should have left
        disp('Children gone');
    else
        I = zeros(1,Nc); % Probability that each child is in class
        Err = zeros(1,Nc); % Error
        
        for i=1:Nc-1 % Cycle over all but youngest child
            V = maxAge^(Nc-i+1); % Size of region of integration
            samples = zeros(1,points);
            for p=1:points
                Y_sample = rand(1,Nc-i+1)*maxAge; % Set of age intervals
                if (Y_sample(1)<b-sum(Y_sample(2:end))) && (Y_sample(1)>a-sum(Y_sample(2:end)))
                    P = IntervalProb(Y_sample,ProbStop,kB,kL,kB,TB,TL,ticker);
                    samples(p) = P;
                end
            end
            I(i) = V*mean(samples);
            Err(i) = V*var(samples);
        end
        
        % Now do youngest child
        V = maxAge;
        samples = zeros(1,points);
        for p=1:points
            Y_sample = rand(1,1)*maxAge;
            if Y_sample<b && Y_sample>a
                P = IntervalProb(Y_sample,ProbStop,kB,kL,kB,TB,TL,ticker);
                samples(p) = P;
            end
        end
        I(Nc) = V*mean(samples);
        Err(Nc) = V*var(samples);
        
        E = sum(I);
    end
end
end

function [P] = IntervalProb(Y,P_Stop,kB,kL,TB,TL,ticker)
% This function calculates the probability of getting a given set of age
% intervals at a specific ticker point

Nc = length(Y); % Number of children
N_max = find(P_Stop==1,1)-2; % Max number of CHILDREN

p_Erl = 1; % Store joint probability of A1,...,AN-1 (Erlangs)
lambda_B = kB/TB; % Rate of Erlang events

for i=1:Nc-1 % This loop covers all but youngest individual
    p_Erl = p_Erl*lambda_B*exp(-lambda_B*Y(i));
    for j=1:kB-1
        p_Erl = Y(i)*lambda_B*p_Erl;
        p_Erl = p_Erl/(kB-j);
    end
end

% The procedure for the last individual depends on the ticker value
if ticker<kB+1 % Waiting for another birth
    p_last = gampdf(Y(Nc),ticker,1/lambda_B); % Erlang up to current ticker
elseif ticker<kB+kL+1 % Waiting for first one to leave
    p_last = gampdf(Y(Nc),ticker-kB,TL/kL); % Erlang again
elseif ticker<kB+kL+kB+1 % Waiting for others to leave
    p_last = 0;
    if N_max>1
        Size_Is = zeros(1,N_max-Nc);
        for j=1:N_max-Nc % Nc+j gives total number of children to be born
            Size_Is(j)=P_Stop(Nc+j+2); % Probability that total number of kids was N_max-(i-j)
            for l=0:Nc+j-1
                Size_Is(j)=Size_Is(j)*(1-P_Stop(l+2));
            end
            lambda_L1 = (kL/TL)*ones(1,kL);
            lambda_L2 = (kB/TB)*ones(1,ticker-kB-kL+kB*(j-1));
            lambda = [lambda_L1 lambda_L2];
            p_last = p_last + Size_Is(j)*HypoExpPdf(Y(Nc),lambda);
        end
        p_last = p_last/sum(Size_Is);
    end
else
    p_last = 0;
end

P = p_Erl*p_last;

end

function [P_Hyp] = HypoExpPdf(x,lambda)
% Hypoexponential distribution - see wikipedia etc
range = length(x);
n = length(lambda);

alpha = zeros(1,n);
alpha(1) = 1;

index = 1:n;
Theta = sparse(index,index,-lambda,n,n);
Theta = Theta + sparse(index(1:end-1),index(1:end-1)+1,lambda(1:end-1),n,n);

for l=range:-1:1
    xTheta = x(l)*Theta;
    xTheta = expm(xTheta);
    P_Hyp(l) = -alpha*xTheta*Theta*ones(n,1);
end

P_Hyp = full(P_Hyp);

end