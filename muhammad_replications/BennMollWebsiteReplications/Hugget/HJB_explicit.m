clear all; clc;

tic;

s = 2; %CRRA utility with parameter s
r = 0.03; %interest rate
rho = 0.05; %discount rate
z1 = .1;
z2 = .2;
z = [z1,z2];
la1 = 0.02; %lambda_1
la2 = 0.03; %lambda_2
la = [la1,la2];

% state space for assets 
I = 100; %number of grid points for assets
amin = -.02; %minimum asset level
amax = 2; %maximum asset level
a = linspace(amin, amax, I)'; %asset grid
da = a(2) - a(1); %asset grid spacing
aa = [a, a]; %asset grid for both statesz
zz = ones(I,1)*z;   %state grid for poisson shocks, a column of I for each state
maxit= 200;
crit = 10^(-6);
% set up matrices for consumption, value function and savings policy
c = zeros(I,2); %consumption matrix
v = zeros(I,2); %value function matrix
s = zeros(I,2); %savings policy matrix
% set up forward and backward difference matrices 
dVf = zeros(I,2); %forward difference matrix
dVb = zeros(I,2); %backward difference matrix

% initialize value function to utility at the current income and wealth level
v(:,1) = ((z1*aa(:,1)).^(1-s)/(1-s))/rho; %value function for state 1
v(:,2) = ((z2*aa(:,2)).^(1-s)/(1-s))/rho; %value function for state 2

% iterate on the value function until convergence
for it = 1:maxit
    % calculate forward and backward differences
    dVf(1:I-1,1) = (v(2:I,1) - v(1:I-1,1))/da; %forward difference for state 1
    dVf(1:I-1,2) = (v(2:I,2) - v(1:I-1,2))/da; %forward difference for state 2
    dVb(2:I,1) = (v(2:I,1) - v(1:I-1,1))/da; %backward difference for state 1
    dVb(2:I,2) = (v(2:I,2) - v(1:I-1,2))/da; %backward difference for state 2
    
    % calculate consumption and savings policy with both forward and backward differences
    cf = (dVf)^(-1/s); %forward consumption
    cb = (dVb)^(-1/s); %backward consumption
    sf = zz+r*aa - cf; %forward savings policy
    sb = zz+r*aa - cb; %backward savings policy

    % at ss when savings policy is 0, we have SS consunption and derivative of value function by 
    c0 = zz+r*aa; %consumption at steady state
    dV0 = (c0).^(-s); %derivative of value function at steady state
    % calculate derivative using upwind scheme 
    If = sf > 0 ; %indicator for forward savings policy or positive drift 
    Ib = sb < 0 ; %indicator for backward savings policy or negative drift
    I0 = 1-If - Ib; %indicator for no drift

    
    v = v_new; %update value function
end