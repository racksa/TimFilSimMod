clear

cd data/

D = load('plane_bent_in_one_direction_sync_data.dat');

hold on;

for i=1:size(D,1)
    
    L_over_R = D(i,2);
    N1_over_N2 = 10/D(i,1);
    direction = D(i,3);
    
    if direction==1 % Straight direction
        
        text(L_over_R,N1_over_N2,num2str(1),'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
        
    elseif direction==2 % Bent direction
        
        text(L_over_R,N1_over_N2,num2str(2),'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
        
    elseif direction==3
        
        plot(L_over_R,N1_over_N2,'ko');
        
    end
    
end

% Bounding curves

N2 = 10./(0:0.05:2);
val = 2*pi./(1+N2);

plot(-val,0:0.05:2,'r-');
plot(min(val,1),0:0.05:2,'r-');

set(gca,'FontSize',24,'FontName','Times');
xlabel('$L/R_2$','Interpreter','latex');
ylabel('$N_1 / N_2$','Interpreter','latex');
axis([-0.5 0.5 0.5 2]);
set(gca,'Yscale','log');

hold off;

cd ..