energy = load("energy.test");
timestep = [1:1e2];

figure(1), clf, hold on, box on;
plot(timestep, energy(:,3)+energy(:,4));