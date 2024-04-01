function main()
    clear all
    close all
    clc

    Xmax = 499;
    Xmin = 0;
    H = input('Combien d''harmoniques voulez-vous?');
    a = input('Choisissez votre amplitude');
    f = input('Entrez la frequence du signal');
    fe= input('Indiquez la frequence d''echantillonnage');
    fin=input('Donnez la fin de votre vecteur temps');
    t = 0:1/fe:fin;

    [TabConcat, TabArche, TabReconstr, reconstruction_arche] = calculateSignals(H, a, f, t);

    plotResults(TabConcat, TabArche, TabReconstr, reconstruction_arche, Xmin, Xmax);
end

function [TabConcat, TabArche, TabReconstr, reconstruction_arche] = calculateSignals(H, a, f, t)
    TabConcat = zeros(2*H, 1);
    TabArche = zeros(2*H, length(t));

    for i = 1:2*H
        CoeffArche = calculateXArche(i);
        TabConcat(i,1) = CoeffArche;

        HarmoArche = cos(2*pi*i*f*t);
        TabHarmos(i,:) = HarmoArche;
        TabArche(i,:) = CoeffArche*HarmoArche;
    end

    TabReconstr = cumsum(TabArche);
    reconstruction_arche = TabReconstr(2*H,:);
end

function CoeffArche = calculateXArche(i)
    CoeffArche = -(1)/(i^2*pi)*(-1)^i;

end

function plotResults(TabConcat, TabArche, TabReconstr, reconstruction_arche, Xmin, Xmax)
    TabAbsolu = abs(TabConcat(find(abs(TabConcat)>=3.898171832519376e-17)));
    
    figure(1)
    hold on
    stem(abs(TabAbsolu))
    plot(abs(TabAbsolu),'r--')
    title('Les harmoniques du signal arche')
    xlabel('Coefficent du signal arche')
    ylabel('Amplitude')

    figure(2)
    plot(TabReconstr')
    title('Représentation 2D des harmoniques du signal arche')
    xlabel('Temps(ms)')
    ylabel('Amplitude')
    axis([Xmin Xmax min(reconstruction_arche)-0.1 max(reconstruction_arche)+0.1])

    figure(3)
    plot(reconstruction_arche)
    title('Signal arche reconstruit grâce à la série de Fourier réelle ')
    xlabel('Temps (ms)')
    ylabel('Amplitude')
    axis([Xmin Xmax min(reconstruction_arche)-0.1 max(reconstruction_arche)+0.1])

    figure(4)
    mesh(TabReconstr)
    title('Visualisation 3D de la décomposition en série de Fourier réelle du signal arche')
    xlabel('Temps (ms)')
    ylabel('Le rang des coefficients')
    zlabel('Amplitude')

    figure(5)
    stem(abs(fft(reconstruction_arche)))
    axis([Xmin 60 0 160])
    title('Estimation de la TF du signal arche grâce à l''algorithme FFT')
    xlabel('Coefficient du signal arche')
    ylabel('Amplitude')

    figure(6)
    plot(TabArche')
    title('Représentation des harmoniques pondérés')
    xlabel('Fréquence')
    ylabel('Amplitude')
    axis([Xmin Xmax+1 -0.7 0.7])

    figure(7)
    mesh(TabArche)
    title('Rerésentation 3D de la décroissance en 1/n du signal arche')
    xlabel('Le nombre d''harmoniques')
end
