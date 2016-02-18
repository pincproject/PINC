% Testing the probability of particles leaving.

close all
clear all

% Number of grid points is Ng. Using normalized length L=Ng, and assigns
% uniform position in this interval, and velocity between +/- 1.

for p=8:25;

    No(p) = 2^p;
    N = No(p);
    Ng = 64;

    pos = Ng*rand(N,1);
    vel = 2*rand(N,1)-1;

    pos = pos+vel;

    left = numel(find(pos<0))+numel(find(pos>Ng));

    fraction(p) = 100*left/N;

end

analytical = 100/(2*Ng)*ones(size(No))


semilogx(No,fraction);
hold all
semilogx(No,analytical)
grid on