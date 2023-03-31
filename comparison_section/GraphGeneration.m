Times = readmatrix('Times.csv'); %We saved our Sets variable as 'Times.csv' - this pulls that into this code

%We did do a longer run - up to 450x450, instead of 150x150, however with
%fixed alpha = 2

DeepTimes = readmatrix('DeepTimes.csv');

%Now we want to split our datasets according to: - Method (so we can graph
%them separately), and by Alpha

%we have 8 values of alpha, and 3 methods
% if we include 1 column as an empty space separator between alphas
% and one column for our dimension, that makes... 33 columns

DataForPlot = zeros(141,33);


for i1 = 1:1128
    Dimension = Times(i1,1);
    VarAlpha = Times(i1,2);
    VarRank = (VarAlpha + 1)/2;
    Col1 = (4*VarRank)-1;
    Col3 = (4*VarRank)+1;

    DataRow = Dimension - 9;
    DataForPlot(DataRow,1) = Dimension;
    DataForPlot(DataRow,Col1:Col3) = Times(i1,3:5);
end

%Now DataForPlot is in the form:
%Dim 0 Method1 Method2 Method3  0 Method1 Method2 Method3
%etc.
1+1;

tiledlayout(2,2)
for i2 = 0:2
    %3 methods, 3 setups

    nexttile
    hold on
    for i3 = 1:2:15
        VarRank_1 = (i3 + 1)/2;
        Column = (4*VarRank_1)-1 + i2;
        plot(DataForPlot(:,1),DataForPlot(:,Column))


    end
    hold off
    if i2 == 0
        X = ' GMRES';%sets correct method label
    elseif i2 == 1
        X = ' Conjugate Gradient';
    elseif i2 == 2
        X = ' Biconjugate Gradient';
    end
    %plots
    X = strcat('Time for the ',X,' method to complete with varied Alpha and Dimensions');
    title(X)
    xlabel('Dimensions of the X by X matrix')
    ylabel('Time for method to complete (seconds)')
    legend({'1','3','5','7','9','11','13','15'},'Location','northeastoutside')
end

%Now we're going to plot our last
nexttile
hold on
plot(DeepTimes(:,1),DeepTimes(:,3))
plot(DeepTimes(:,1),DeepTimes(:,4))
plot(DeepTimes(:,1),DeepTimes(:,5))
xlabel('Dimensions of the X by X matrix')
ylabel('Time for method to complete (seconds)')
legend({'GMRES','Conjugate Gradient','Biconjugate Gradient'},'Location','northeastoutside')
title('Time for various methods to complete, with constant alpha = 2')
hold off

