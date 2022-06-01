%willekeurig beginpunt en willekeurige matrix van afmetingen rijen x
%kolommen: kans dat resultaat gevonden?
%rijen kolommen rang aantalMetingen

%i: aantal metingen
%j: verschillende startscenario's
% 1 keer met A en B random, 2de keer met kleine 1 en grote B maar product
% blijft gelijk
x = zeros(200,1);
y = zeros(200,1);
z = zeros(200,1);
k=0;
m=0;
for(i=1:20)
    for(j=1:20)
        k=k+1;
        x(k)=i;
        [y(k),z(k)] = recover2(5,4,2,i);
    end
end
%eerste 5 metingen zijn allemaal vor-or 1 meting, 2de 5 voor 2 meting ..;
fig1 = figure(1)
plot(x,y,'bs')
hold on
plot(x,z,'rs')
legend('z=AB', 'z=(A/10^6)*(10^6*B)')
title('afwijking benadering 5x4-matrix van rang 2 met fminunc en gradient door Matlab, startpunt z met A en B random gegenereerd')
xlabel('aantal metingen')
ylabel('norm(matrix - result)/norm(matrix)')
grid on

g=0;
k=0;
yGemiddeld = zeros(20,1);
zGemiddeld = zeros(20,1);
for(j=1:20)
    g=g+1;
    gemiddeldey =0;
    gemiddeldez=0;
    for i = 1:20
        k=k+1;
        gemiddeldey = gemiddeldey + y(k)
        gemiddeldez = gemiddeldez + z(k)
    end
    yGemiddeld(g) = gemiddeldey/20;
    zGemiddeld(g) = gemiddeldez/20;
end
x=1:1:20;
fig2 = figure(2)
plot(x,yGemiddeld,'bs')
hold on
plot(x,zGemiddeld,'rs')
legend('z=AB', 'z=(A/10^6)*(10^6*B)')
title('gemiddelde afwijking benadering 5x4-matrix van rang 2 met fminunc en gradient door Matlab, startpunt z met A en B random gegenereerd')
xlabel('aantal metingen')
ylabel('norm(matrix - result)/norm(matrix)')
grid on