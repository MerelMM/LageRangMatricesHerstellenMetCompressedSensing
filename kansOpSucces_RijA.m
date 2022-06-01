
%recoverMetGradient
% met kolom 1 *10^-3 -6 -9
% kijk ook hoeveel iteraties er nodig zijn
conditieV = 1e-9;
drempel = 10e-5;
aantalMetingen = 20;
aantalKeer = 10;
x = zeros(200,1);
y = zeros(200,1);
k=0;
m=0;
z= zeros(200,1);
rijen = 6;
kolommen = 5;
rang = 3;
for i=1:aantalMetingen
    for j=1:aantalKeer
        k=k+1;
        x(k)=i;
        [y(k), z(k)] = conditie(rijen, kolommen, rang, i);
    end
end

%%gemiddelde
g=0;
k=0;
succesRecovery = zeros(20,1);
succesRecoverz = zeros(20,1);
for(j=1:aantalMetingen)
    g=g+1;
    succes =0;
    succesz = 0;
    for i = 1:aantalKeer
        k=k+1;
        if (y(k)<=drempel)
            succes = succes + 1;
        end
         if (z(k)<=drempel)
            succesz = succesz + 1;
        end
    end
    gemiddeldI
    succesRecovery(g) = succes/aantalKeer;
    succesRecoverz(g) = succesz / aantalKeer
end
x= 1:1:20;
fig2= figure(2)
plot(x,succesRecoverz, '--ks')
title('goed conditie')
xlabel('aaltal metingen')
ylabel ('aantal iteraties')
grid on

fig1 = figure(1)
plot(x,succesRecovery,'--ks')
title('slechte conditie')
xlabel('aantal metingen')
ylabel('succesvol hersteld [%]')
grid on
