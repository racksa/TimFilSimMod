close all
clear

N = 6;

% Use spiral seeding as initial condition:
X = zeros(3,N);

for n=1:N
    
    X(3,n) = 2*(n-1)/(N-1) - 1;
    
    r = sqrt(1 - X(3,n)*X(3,n));
    
    if (n==1) || (n==N)
        
        phi = 0;
        
    else
        
        phi = phi + 3.6/(sqrt(N)*r);
        
    end
    
    X(1,n) = r*cos(phi);
    X(2,n) = r*sin(phi);
    
end

h1 = figure;

plot3(X(1,:),X(2,:),X(3,:),'k.');
axis equal;
axis([-1 1 -1 1 -1 1]);
title('Initial seeding');

% Begin the repulsive-potential simulation:
f0 = 5;
dt = 1.0;
Xold = 10^4 + X;

h2 = figure;

iter = 0;

while max(abs(X-Xold), [], 'all') > 10^-4
    
    F = zeros(3,N);
    Xold = X;
    
    for n=1:N-1
        
        Xn = X(:,n);
        
        for m=n+1:N
            
            d = Xn - X(:,m);
            dnorm = norm(d);
            d = d/dnorm;
            
            F(:,n) = F(:,n) + (2-dnorm)*f0*d;
            F(:,m) = F(:,m) - (2-dnorm)*f0*d;
            
        end
    end
    
    X = X + dt*F;
    
    for n=1:N
        
        X(:,n) = X(:,n)/norm(X(:,n));
        
    end
    
    if mod(iter,10)==0
        
        figure(h2);
        plot3(X(1,:),X(2,:),X(3,:),'k.');
        axis equal;
        axis([-1 1 -1 1 -1 1]);
        pause(0.00001);
        
    end
    
    iter = iter + 1;
    
end

figure(h2);
plot3(X(1,:),X(2,:),X(3,:),'k.');
axis equal;
axis([-1 1 -1 1 -1 1]);
title('Equilibrium seeding');

dlmwrite(sprintf('sphere%i.seed', N), X', 'delimiter', ' ');