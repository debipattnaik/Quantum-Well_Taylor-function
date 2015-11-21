%program to find the lowest energy states of E of 1 D time independent
%Schrodinger equation for potential specified by user.
%Debi prasad Pattnaik, Uni of Nottingham,UK

%parameters for solving
clc; close all
L=5; % length of the interval
N=1000;
x=linspace(-L,L,N)';
dx=x(2)-x(1);
% declaration of the potential
%U=1/2*10*x.^(2); % quadriatic harmonic oscillator potential
U=10*(-1+1/2*x.^(2)-1/24*x.^(4)+1/720*x.^(6))
% replace this potential with any potential of your choice. :)

% solving here
e=ones(N,1);
lap=spdiags([e -2*e e],[-1 0 1],N,N)/dx^2;
hbar=1;
m=1;
H=-1/2*(hbar^2/m)*lap +spdiags(U,0,N,N);
nmodes=3;options.disp=0;
[V,E]=eigs(H,nmodes,'sa',options); % find eigen values
[E,ind]=sort(diag(E));
V=V(:,ind);
Usc=U*max(abs(V(:)))/max(abs(U));
plot(x,V,x,Usc,'--k');
lgnd_str=[repmat('Energy:',nmodes,1),num2str(E)];
legend(lgnd_str);
