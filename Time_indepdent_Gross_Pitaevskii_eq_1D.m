function [] = Time_indepdent_Gross_Pitaevskii_eq_1D 
%%%
%%% --- 1D time-independent Gross-Pitaevskii equation for Bose-Einstein condesation at T = 0 K
%
% [-0.5*nabla^2_xx + V_ext + V_mf]*psi(x) = mu*psi(i);
% mu is chemical potential per particle
% V_mf = lambda*|psi|^2 is mean-field potential
% V_ext = 0.5*x^2 is external harmonic potenial 
% E = mu - 0.5*< |lambda*|psi|^2| > is energy per particle
% 
% Uses: function legDC2 - differentation matrix based on Legendre-Gauss-Lobatto nodes
%%% ---
%
% References:
%          1.Ts. Tsogbayar, Ts. Banzragch, & Kh. Tsookhuu, Proceeding of
%           International Conference of Applied Science and Engineering -
%           2021, Ulaanbaatar, Mongolia 
%          2. Ts.Banzragch & Ch.Khurelbaatar, Booklet of Young Researchers' Meeting - Khureltogoot-2017, Mongolia
%          3. Yu-Ren S., Guang-Hui W., and Cong-Bo L., et al., (2012) Chin. Phys. Lett. 29, 110302
%           
%
% Written by Tsogbayar Tsednee (PhD), Institute of Physics and Technology,
%                                     Mongolian Academy of Sciences, Ulaanbaatar, Mongolia 
% Contact: tsog215@gmail.com
% Date: October 22, 2021
%
format short
clear;
%
clc;
N = 256; a = -10.0; b = -a; % grid parameters & you may change them.
itermax = 300; tol = 10^(-8); % computational parameters & you may change them.
%
eig_num = 1.; % eigenvalue index & you may change it
%
alf = 0.50; % convergence parameter & you may change it 
%
lambda = -1.0; % lambda > 0 is for repulsive potential; 
               % lambda < 0 is for attractive potential
%
%%%
[x,wx,D]=legDC2(N,a,b);
wx = wx(2:N); x = x(2:N); 
D2 = (2/(b-a))^2*D^2;
D2 = D2(2:N,2:N);         % differentation matrix taking into account the boundary condition
%
wf_1s = exp(-0.5*x.^2);%./(pi)^(1/4); % initial trial function 
sm = sum(wx.*wf_1s.*wf_1s);
n_c = 1./sqrt(abs(sm));      % The normalization constant
wf_1s = n_c*wf_1s;
%%%
V_old = lambda*wf_1s.^2; 
%
V_mf = sum(wx.*wf_1s.*V_old.*wf_1s); % expectation value of the mean-field potential
%
V_sq = wf_1s.^2;
%
for iter = 1:itermax
%    
%    iter 
    V_old_sq = V_sq;
    V_new = lambda*V_old_sq;
    V_mf = sum(wx.*V_old_sq.*V_new);
%    
    V_n = alf*V_new + (1-alf)*V_old;
    V_old = V_n;
%%% The Hamiltonian
    H_ham = -0.5*D2 + 0.5*diag(x.^2) + diag(V_n);
%
    [Vec,En] = eig(H_ham);                                     % Eigenvalue problem
    En = diag(En);
    [foo, ij] = sort(En);
    En = En(ij);
%    [En(1),En(2),En(3),En(4),En(5)];
%
    Vec = Vec(:,ij);                       % The unnormalized eigenfunctions
    V1 = Vec(:,eig_num);                         % The unnormalized eigenfunction for the ground state, 
                                       % here 1 is for the ground state, 2,3... are corresponded 
                                       % for the excited states
%V1 = [0.,;V1,;0.];
    sm = sum(wx.*V1.*V1);    
    %
    n_c = 1./sqrt(abs(sm));      % The normalization constant
    V1 = n_c*V1;
    %
    V_ext = sum(wx.*V1.*(0.5*x.^2).*V1);
    %
    En_total = En(eig_num) - 0.5*V_mf;
    
    V_nsq = V1.*V1;
    
    cor = rms(V_nsq - V_old_sq);     
%
    V_sq = V_nsq;
%
    if (cor < tol)
        break
    end  
    
[iter, V_mf, V_ext, En(eig_num), En_total  ]  ;

end
%
[lambda, iter, V_mf, V_ext, En(eig_num), En_total  ]  
%
plot(x,V1) % plot of psi(x)
%
%
end
%%%
% [lambda, iter, V_mf, V_ext, En(eig_num), En_total  ]
%
% for repulsive potential (lambda > 0) & alf = 0.2
%  1.0000   41.0000    0.3609    0.2996    0.8699    0.6895
%  5.0000   34.0000    1.3908    0.4842    2.0115    1.3161
% 10.0000   34.0000    2.3202    0.6835    3.1072    1.9471
%%%%
% for attactive potential (lambda < 0) & alf = 0.5
% -1.0000   25.0000   -0.4510    0.2005    0.0627    0.2882
% -5.0000   57.0000   -4.3856    0.0582   -3.1728   -0.9801


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xi,w,D]=legDC2(N,a,b)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% legDc.m
%
% Computes the Legendre differentiation matrix with collocation at the 
% Legendre-Gauss-Lobatto nodes.
%
% Reference: 
%   C. Canuto, M. Y. Hussaini, A. Quarteroni, T. A. Tang, "Spectral Methods
%   in Fluid Dynamics," Section 2.3. Springer-Verlag 1987
%
% Written by Greg von Winckel - 05/26/2004
% Contact: gregvw@chtm.unm.edu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Truncation + 1
N1=N+1;

% CGL nodes
xc=cos(pi*(0:N)/N)';

% Uniform nodes
xu=linspace(-1,1,N1)';

% Make a close first guess to reduce iterations
if N<3
    x=xc;
else
    x=xc+sin(pi*xu)./(4*N);
end

% The Legendre Vandermonde Matrix
P=zeros(N1,N1);

% Compute P_(N) using the recursion relation
% Compute its first and second derivatives and 
% update x using the Newton-Raphson method.

xold=2;
while max(abs(x-xold))>eps

    xold=x;
        
    P(:,1)=1;    P(:,2)=x;
    
    for k=2:N
        P(:,k+1)=( (2*k-1)*x.*P(:,k)-(k-1)*P(:,k-1) )/k;
    end
     
    x=xold-( x.*P(:,N1)-P(:,N) )./( N1*P(:,N1) );
end

X=repmat(x,1,N1);
Xdiff=X-X'+eye(N1);

L=repmat(P(:,N1),1,N1);
L(1:(N1+1):N1*N1)=1;
D=(L./(Xdiff.*L'));
D(1:(N1+1):N1*N1)=0;
D(1)=(N1*N)/4;
D(N1*N1)=-(N1*N)/4;

% Linear map from[-1,1] to [a,b]
xi=(a*(1-x)+b*(1+x))/2;      % added by Tsogbayar Tsednee

% Compute the weights
w=(b-a)./(N*N1*P(:,N1).^2); % added by Tsogbayar Tsednee
return
end
