% Supply Chain Modeling & Analysis
% Primary Author - Ajinkya Bharambe

clear all; clc
% Given
% For all give, 1 represents Diesel, 2 represents Gas engine oil,
% 3 represents natural gas engine oil
% Diesel Engine Oil - DE495 
k1 = 7500;          % Setup cost
p1 = 2;             % Unit stockout cost 
d1 = 1000000;       % Annual demand  
sd1 = 450000;       % standard deviation
var_d1 = sd1^2;     % Variance     
Leadt = 1;          % lead time
leadt = 1/12;
h1 = 1;             % Holding cost

% Gas Engine oil - GE275
k2 = 7500;          % Setup cost
p2 = 10;            % Unit stockout cost   
d2 = 400000;        % Annual demand  
sd2 = 140000;       % standard deviation
var_d2 = sd2^2;     % Variance     
Leadt = 1;          % lead time
leadt = 1/12;
h2 = 5;             % Holding cost

% Natural Gas Engine oil - NG645
k3 = 7500;          % Setup cost
p3 = 15;            % Unit stockout cost  
d3 = 500000;        % Annual demand  
sd3 = 65000;        % standard deviation
var_d3 = sd3^2;     % Variance     
Leadt = 1;          % lead time
leadt = 1/12;
h3 = 5;             % Holding cost
%% Question 1
[EOQ1,RI1,RL1] = eoq(k1,d1,h1,leadt);
[EOQ2,RI2,RL2] = eoq(k2,d2,h2,leadt);
[EOQ3,RI3,RL3] = eoq(k3,d3,h3,leadt);

fprintf('|---------------------|\n')
fprintf('      Question 1\n')
fprintf('|---------------------|\n')
Product = {'Diesel Engine Oil- DE495';'Gas Engine Oil- GE275';'Natural Gas Engine Oil- NG645'};
EOQ = [EOQ1;EOQ2;EOQ3];
Reorder_Int = [RI1;RI2;RI3];
Reorder_Level = [RL1;RL2;RL3];
T = table(Product,EOQ,Reorder_Int,Reorder_Level);
disp(T)

%% Question 2
fprintf('|---------------------|\n')
fprintf('      Question 2\n')
fprintf('|---------------------|\n')
[mew1,var1,std1] = mvs(d1,var_d1,Leadt);
[mew2,var2,std2] = mvs(d2,var_d2,Leadt);
[mew3,var3,std3] = mvs(d3,var_d3,Leadt);
fprintf('Following are the parameters during lead time---->\n')
Mean = [mew1;mew2;mew3];
Variance = [var_d1;var_d2;var_d3];
St_Deviation = [std1;std2;std3];
T1 = table(Product,Mean,Variance,St_Deviation);disp(T1)

fprintf('\nQ,R Iterations for Diesel Engine Oil---->\n')
[Q1,R1,n_R1,z1,in1,F_R1] = qr1(EOQ1,k1,d1,h1,std1,mew1,p1);
fprintf('\nQ,R Iterations for Gas Engine Oil---->\n')
[Q2,R2,n_R2,z2,in2,F_R2] = qr1(EOQ2,k2,d2,h2,std2,mew2,p2);
fprintf('\nQ,R Iterations for Natural Gas Engine Oil---->\n')
[Q3,R3,n_R3,z3,in3,F_R3] = qr1(EOQ3,k3,d3,h3,std3,mew3,p3);

fprintf('\nSummarize the above results----->\n');
Iterations = [in1;in2;in3];
Q = [Q1(in1);Q2(in2);Q3(in3)];
R = [R1(in1);R2(in2);R3(in3)];
T3 = table(Product,Iterations,Q,R); disp(T3)

%% Problem 6
% Special Engine Oil - IEE534 
kn = 7500;           % Setup cost
pn = 10;             % Unit stockout cost 
dn = d1+d2+d3;       % Annual demand  
var_dn = sd1^2+sd2^2+sd3^2;     % Variance     
sdn = sqrt(var_dn);
Leadt = 1;          % lead time
leadt = 1/12;
hn = 5;             % Holding cost

[EOQn,RIn,RLn] = eoq(kn,dn,hn,leadt);

fprintf('|---------------------|\n')
fprintf('      Question 6\n')
fprintf('|---------------------|\n')
[mewn,varn,stdn] = mvs(dn,var_dn,Leadt);
fprintf('\nQ,R Iterations for Special Engine Oil---->\n')
[Qn,Rn,n_Rn,zn,inn,F_Rn] = qr1(EOQn,kn,dn,hn,stdn,mewn,pn);

%% Functions required
function [EOQ,RI,RL] = eoq(k,d,h,leadt)
EOQ = sqrt(2*k*d/h);
RI = EOQ/d;
RL = EOQ*(leadt/RI);
end

function [mew,var,std] = mvs(d,var_d,Leadt)
mew = d*Leadt/12;         % Mean of demand with lead time
var = var_d*Leadt/12;     % Variance of demand with lead time
std = sqrt(var);          % Standard deviation of demand with lead time

%fprintf('Mean during lead time = %f\n', mew);
%fprintf('Variance during lead time is %f\n', var);
%fprintf('The standard deviation during lead time = %f\n', std);
end
function [q,r,n_R,z,in,F_R] = qr1(EOQ,k,d,h,std,mew,p)

F_R = 1-((EOQ*h)/(p*d));
z = norminv(F_R); 
fprintf('<--------Iteration 1-------->\n')
R_1 = std*z + mew;                          % Reorder point
L = (normpdf(z) - (z*cdf('Normal',-z,0,1)));
n_R = L*std;
Q_1 = sqrt(2*d*(k+(p*n_R))/h);
fprintf('Z = %f\n',z)
fprintf('Q = %f, R = %f.\n',Q_1,R_1);

Q1(1) = Q_1;
R1(1) = R_1;
n = 100;  % Just some random number so that we can run iterations for (Q,R)
for i = 2:n
    fprintf('<--------Iteration %d-------->\n',i)
    F_R = 1-((Q1(i-1)*h)/(p*d));
    z = norminv(F_R);
    R = std*z + mew;
    L = (normpdf(z) - (z*cdf('Normal',-z,0,1)));
    n_R = L*std;
    Q = sqrt(2*d*(k+(p*n_R))/h);
    fprintf('Z = %f\n',z)
    fprintf('Q = %f, R = %f.\n',Q,R);
    Q1(i) = Q; R1(i) = R;
    q(i) = Q; r(i) = R;

    if round(Q) == round(Q1(i-1)) && round(R) == round(R1(i-1))
        in=i;
        break;
    end
end

end

