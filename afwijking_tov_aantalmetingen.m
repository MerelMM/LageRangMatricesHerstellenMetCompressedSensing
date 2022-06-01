%willekeurig beginpunt en willekeurige matrix van afmetingen rijen x
%kolommen: kans dat resultaat gevonden?
%rijen kolommen rang aantalMetingen

%i: aantal metingen
%j: verschillende startscenario's
x = zeros(200,1);
y = zeros(200,1);
k=0;
m=0;
for(i=1:20)
    for(j=1:20)
        k=k+1;
        x(k)=i;
        y(k) = recover(5,4,2,i);
    end
end
%eerste 5 metingen zijn allemaal vor-or 1 meting, 2de 5 voor 2 meting ..;
fig1 = figure(1)
plot(x,y,'ks')
title('afwijking matrixrecovery 5x4-matrix van rang 2 met fminunc en gradient door Matlab ')
xlabel('aantal metingen')
ylabel('norm(matrix - result)/norm(matrix)')
grid on

%%gemiddelde
g=0;
k=0;
yGemiddeld = zeros(20,1);
for(j=1:20)
    g=g+1;
    gemiddeldey =0;
    for i = 1:20
        k=k+1;
        gemiddeldey = gemiddeldey + y(k)
    end
    yGemiddeld(g) = gemiddeldey/20;
end
x=1:1:20;
fig2 = figure(2)
plot(x,yGemiddeld,'bs')
title('gemiddelde afwijking matrixrecovery 5x4-matrix van rang 2 met fminunc en gradient door Matlab')
xlabel('aantal metingen')
ylabel('norm(matrix - result)/norm(matrix)')
grid on