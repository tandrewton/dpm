figure(); hold on;
A = load('energyAndForceTest.txt');
numIts = size(A,1);
plot(0:numIts-1, A(:,1), 'DisplayName', "force x")
plot(0:numIts-1, A(:,2), 'DisplayName', "force y")
plot(0:numIts-1, A(:,3), 'DisplayName', "energy")

xlabel("frame number")
ylabel("force or energy")

numDirections = 2;
for i=1:numDirections
    xline(i*numIts/numDirections - 1, 'HandleVisibility','off')
end

legend('Location','east')