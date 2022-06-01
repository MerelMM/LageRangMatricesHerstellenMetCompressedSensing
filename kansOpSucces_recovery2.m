drempel = 10e-5;
aantalMetingen = 20;
aantalKeer = 20;
x = zeros(200,1);
y = zeros(200,1);
z = zeros(200,1)
k=0;
m=0;
for i=1:aantalMetingen
    for j=1:aantalKeer
        k=k+1;
        x(k)=i;
        [y(k),z(k)] = recover2(5,4,2,i);
    end
end

%plot kans op succes y
g=0;
k=0;
succesRecovery = zeros(20,1);
succesRecoveryZ = zeros(20,1);
for j=1:20
    g=g+1;
    succes =0;
    succesZ = 0;
    for i = 1:20
        k=k+1;
        if (y(k)<=drempel)
            succes = succes + 1;
        end
        if(z(k)<=drempel)
            succesZ = succesZ+1;
        end
    end
    succesRecovery(g) = succes/aantalKeer;
    succesRecoveryZ(g) = succesZ/aantalKeer;
end
x=1:1:20;
fig1 = figure(1)
plot(x,succesRecovery,'--ks')
hold on
plot(x,succesRecoveryZ, '--rs')
title('kans op succes (afwijking <= 1e-5) matrixrecovery 5x4-matrix van rang 2 met fminunc en gradient door Matlab afhankelijk van factorisatie matrix z')
legend('z=AB', 'z=(A/10^6)*(10^6*B)')
xlabel('aantal metingen')
ylabel('succesvol hersteld [%]')
grid on
