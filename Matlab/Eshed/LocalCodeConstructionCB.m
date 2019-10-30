%% code parameters
kappa = 13;gamma_l=3; gamma = 6; gamma_c = gamma - gamma_l;
%from Tanner's paper 
p = 157; 
% erased circulants
nu = 12; a = floor(a/gamma_l); b = nu - a; 

%% create codes
% Unbalnced partition
B=ones(gamma,kappa); %CC are irelevant
B(gamma_l+1,end-nu+1:end)=0;
H = CBlift(B,p);
localRow = gamma_c*p+1; HL_1 = H(localRow:end,:);

% Balanced partition
B=ones(gamma,kappa);endj= 0;
for ii = 1:gamma_l
    startj = endj+1; endj = startj+a-1; 
    B(gamma_l+ii,startj:endj)= 0;
end
B = fliplr(B);
H = CBlift(B,p); HL_2 = H(localRow:end,:);