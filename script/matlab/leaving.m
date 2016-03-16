% Testing the probability of particles leaving.

close all
clear all

% Number of grid points is Ng. Using normalized length L=Ng, and assigns
% uniform position in this interval, and velocity between +/- 1.

Ngx = 32;
Ngy = 16;

for p=10:25;

    No(p) = 2^p;
    N = No(p);


    x = Ngx*rand(N,1);
    y = Ngy*rand(N,1);
    vx = 2*rand(N,1)-1;
    vy = 2*rand(N,1)-1;

    x = x+vx;
    y = y+vy;

    leftx = (x<0)+(x>Ngx);
    lefty = (y<0)+(y>Ngy);
    left_corner = sum(leftx.*lefty);
    left  = logical(leftx+lefty);
%     for i=1:N
%         if(left(i) || 0)
%             disp([num2str(i) ': (' num2str(x(i)) ',' num2str(y(i)) ') leftx: ' num2str(leftx(i)) ', lefty: ' num2str(lefty(i))]);
%         end
%     end
    left = sum(left);
    left_edge = left-left_corner;

    frac_corner(p) = 100*left_corner/N;
    frac(p) = 100*left/N;
    frac_edge(p) = 100*left_edge/N;

end

ana_corner = 1/((2*Ngx)*(2*Ngy))
ana = 0.5/Ngx+0.5/Ngy-0.25/(Ngx*Ngy);
ana_edge = 0.5/Ngx+0.5/Ngy-0.5/(Ngx*Ngy);

ana_corner = ana_corner*100*ones(size(No))
ana = ana*100*ones(size(No))
ana_edge = ana_edge*100*ones(size(No))



semilogx(No,frac_corner);
hold all
grid on
semilogx(No,ana_corner);

figure
semilogx(No,frac);
hold all
grid on
semilogx(No,ana);

figure
semilogx(No,frac_edge);
hold all
grid on
semilogx(No,ana_edge);