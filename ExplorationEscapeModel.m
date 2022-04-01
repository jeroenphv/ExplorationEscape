clear all
clc


%% INPUTS TO THE MODEL

% internal state of the animal
fearBias = exp(-1);
rewardSensitivity = exp(5.5);
anxiety = 80;

% adversity level of the environment
numberOfPredators = 5; %number of predators

% other
boolFig = false; %set true if you want to visualize model -- will severely slow down model



%% Set other parameters of model and pre-allocate vectors

%parameters of model
rangePredators = [100 100 1600 1600]; %coordinates where predators can roam freely
numberOfFakePredators = 20-numberOfPredators; %number of predator-like (safe) stimuli
rangeFakePredators = [100 100 1600 1600]; %coordinates where predator-like (safe) stimuli can roam freely
speedOfPredators = 7.5; %speed of the predators in meter / min
speedOfAgent = 5; %speed of agent
totalminutes = 60*24*365*3; %maximum survival in minutes; set to 3 years


%generate path of predators
cellCoordinatesPredators = cell(numberOfPredators,1);
for i = 1: numberOfPredators
    cellCoordinatesPredators{i} = generatePredatorCoordinates(totalminutes, speedOfPredators, rangePredators);
end

%generate path of predator-like (safe) stimuli
cellCoordinatesFakePredators = cell(numberOfFakePredators,1);
for i = 1: numberOfFakePredators
    cellCoordinatesFakePredators{i} = generatePredatorCoordinates(totalminutes, speedOfPredators, rangeFakePredators);
end


%initialize model parameters
nestLocation = [25 25];
agentLocation = nestLocation; %inital location of agent
energy = 100;
fitness = 100;
boolGo = 0;
boolAttacked = 0;
boolEsc = 0;
fitnessDecay = 100*cos((1:totalminutes)/(totalminutes/(0.5*pi))); %decay of fitness with age
fear_delayed = 0;
foragingLocation = [100+round(200*rand) 100+round(200*rand)]; %initial foraging location
k = 0; %k is number of minutes passed
muPred = 20;
sigmaPred = 25;
muFakePred = 38;
sigmaFakePred = 20;
eta = 0.95; %forgetting of danger is 5% per minute
maxDangerLevel = rand(numberOfFakePredators,1)*60; %max danger level of predator-like (safe) stimuli is 60% of that of actual predators

%pre-allocate vectors to save data
vecEnergy = nan(1,totalminutes);
vecFitness = nan(1,totalminutes);
matLocation = nan(totalminutes,2);
vecValueForage = nan(1,totalminutes);
vecValueFood = nan(1,totalminutes);
vecFear = nan(1,totalminutes);
totalattacks = 0;
totalfeedingbouts = 0;
totalminutesExpl = 0;
totalexits = 0;
totalescapes = 0;

if boolFig
    figure %open figure
end


while fitness > 0 && k < totalminutes %loop through minutes
    
    k = k+1; %add a minute
    
    %retrieve position of predators in this minute
    matScat = [];
    for l = 1:numberOfPredators
        matScat = [matScat; cellCoordinatesPredators{l}(k,:)];
    end
    
    %retrieve position of predator-like (safe) stimuli
    matScatFake = [];
    for l = 1:numberOfFakePredators
        matScatFake = [matScatFake; cellCoordinatesFakePredators{l}(k,:)];
    end
    
    %save data in vectors
    matLocation(k,:) = agentLocation; %save agent location
    vecEnergy(k) = energy; %save energy level
    vecFitness(k) = fitness; %save fitness level
    
    %calculate objetive danger level arising from predators
    if numberOfPredators > 0
        vecDistanceToPredators = sqrt((agentLocation(1)-matScat(:,1)).^2 + (agentLocation(2)-matScat(:,2)).^2); %pythagoras
        vecTotalDangerPred = 100*(exp(-0.5 * ((vecDistanceToPredators - muPred)./sigmaPred).^2) ./ (sqrt(2*pi) .* sigmaPred))/(exp(-0.5 * ((muPred - muPred)./sigmaPred).^2) ./ (sqrt(2*pi) .* sigmaPred)); %feed through gaussians
    else vecTotalDangerPred = [];
    end
    
    %calculate objective danger level arising from predator-like (safe)
    %stimuli
    vecDistanceToFakePredators = sqrt((agentLocation(1)-matScatFake(:,1)).^2 + (agentLocation(2)-matScatFake(:,2)).^2); %pythagoras
    vecTotalDangerFakePred = maxDangerLevel.*(exp(-0.5 * ((vecDistanceToFakePredators - muFakePred)./sigmaFakePred).^2) ./ (sqrt(2*pi) .* sigmaFakePred))/(exp(-0.5 * ((muFakePred - muFakePred)./sigmaFakePred).^2) ./ (sqrt(2*pi) .* sigmaFakePred)); %feed through gaussians
    
    %calculate subjective danger level arising from all stimuli together
    fear = 100*sum([vecTotalDangerPred/100; vecTotalDangerFakePred/100].^(fearBias));
    
    % determine value of food
    valueFood = -log((energy/100))*rewardSensitivity;
    
    %slowly forget danger (with 5%  per minute) when danger has passed
    if fear > fear_delayed
        fear_delayed = fear;
    else fear_delayed = fear_delayed*eta;
    end
    if fear_delayed < 0.05
        fear_delayed = 0; %beause of the 'forgetting' function it never actually becomes 0, set to 0 when reaching < 0.05
    end
    
    %compute net value of foraging and determine behavior of agent
    valueForage = valueFood - fear_delayed - anxiety;
    if valueForage > 0 && agentLocation(1) == nestLocation(1) %if value outweighs dangers and agent is in home location
        boolGo = 1; %go forage
        totalexits = totalexits + 1;
    elseif valueForage < 0
        boolGo = 0; %reatreat back to nest otherwise
        
        %determine if this was an escape response, and count it as such
        if k > 1 && vecEnergy(k) < vecEnergy(k-1) && fear_delayed > 0 && boolEsc == 0 && agentLocation(1) > nestLocation(1) && agentLocation(1) ~= foragingLocation(1)
            boolEsc = 1;
            totalescapes = totalescapes + 1;
            fprintf('[minute %.0f] Escape from predator or predator-like stimulus\n', k)
        end
        
    end
    
    %save data in vector
    vecValueForage(k) = valueForage;
    vecValueFood(k) = valueFood;
    vecFear(k) = fear_delayed;
    
    %set new foraging location if agent is not on its way to food zone
    if boolGo == 0
        foragingLocation = [100+round(200*rand) 100+round(200*rand)];
    end
    
    %change position of animal
    if boolGo == 1 %go foraging
        distance = sqrt((agentLocation(1)-foragingLocation(1))^2 + (agentLocation(2)-foragingLocation(2))^2); %distance to food
        
        if distance < speedOfAgent*(fitness/100) %note that speed of agent declines linearly with fitness
            agentLocation = foragingLocation; %it reaches the foraging location
        else %move towards foraging location
            steps = distance/(speedOfAgent*(fitness/100));
            xChange = (foragingLocation(1)-agentLocation(1))/steps;
            yChange = (foragingLocation(2)-agentLocation(2))/steps;
            agentLocation(1) = agentLocation(1) + xChange;
            agentLocation(2) = agentLocation(2) + yChange;
        end
        
    elseif boolGo == 0 && agentLocation(1) > nestLocation(1) %go back to nest
        distance = sqrt((agentLocation(1)-nestLocation(1))^2 + (agentLocation(2)-nestLocation(1))^2); %distance to nest
        
        if distance < speedOfAgent*(fitness/100) %note that speed of agent declines linearly with fitness
            agentLocation = [nestLocation(1) nestLocation(1)]; %it reaches nest
        else %move towards the nest
            steps = distance/(speedOfAgent*(fitness/100)); %note that speed of agent declines linearly with fitness
            xChange = (nestLocation(1)-agentLocation(1))/steps;
            yChange = (nestLocation(2)-agentLocation(2))/steps;
            agentLocation(1) = agentLocation(1) + xChange;
            agentLocation(2) = agentLocation(2) + yChange;
        end
    elseif boolGo == 0 && agentLocation(1) < nestLocation(1) % remain in nest
        agentLocation = [nestLocation(1) nestLocation(1)];
    end
    
    
    %adjust energy levels
    if agentLocation(1) == foragingLocation(1) &&  agentLocation(2) == foragingLocation(2) %agent is feeding
        energy = energy+30;
        totalfeedingbouts = totalfeedingbouts + 1;
        fprintf('[minute %.0f] Feeding bout\n', k)
    elseif agentLocation(1) == nestLocation(1) %agent is stationary in nest
        energy = energy - 100/2880; %baseline energy loss - animals can live 2-3 days without food before they die (2880 minutes = 2 days)
    else %if animal is not in nest
        energy = energy - 0.05; %extra energy loss when agent is moving
        energy = energy - 100/2880; %baseline energy loss
    end
    
    %adjust energy levels so they stay within [0, 100]
    if energy > 100
        energy = 100;
    elseif energy < 0
        energy = 0;
    end
    
    %calculate fitness
    if energy*fitnessDecay(k)/100 > fitness
        fitness = fitness + (energy*fitnessDecay(k)/100-fitness)*0.00001; %restoration of fitness when energy > fitness
    elseif energy*fitnessDecay(k)/100 < fitness && energy > 0 %deplection of fitness when energy < fitness
        fitness = fitness + (energy*fitnessDecay(k)/100-fitness)*0.001;                                                                                                                                                                                                                                                                                                                  ;
    elseif energy*fitnessDecay(k)/100 < fitness && energy == 0 %rapid decline of fitness when energy = 0
        fitness = fitness - 0.05;
        
        if fitness < 0
            fitness = 0; %fitness can not be negative
        end
        
    end
    
    %adjust fitness to max fitness in age
    if fitness > fitnessDecay(k)
        fitness = fitnessDecay(k);
    end
    
    
    %predator attacks
    if boolAttacked == 0 %if the animal is not currently attacked
        for s = 1 : numberOfPredators %see if any of the predators is around
            if abs(matScat(s,1)-agentLocation(1)) < muPred && abs(matScat(s,2)-agentLocation(2)) < muPred && agentLocation(1) ~= nestLocation(1);
                boolAttacked = 2; %agent is deemed attack if predator enters square area of 20 units (equal to muPred, see fig 1A of paper) around agent
            end
        end
    end
    
    if boolAttacked == 2 % agent is attacked
        fitness = fitness*0.3333; %decrease fitness
        boolAttacked = 1; %makes sure animal can't be attacked until back at nest
        boolGo = 0; %go back to nest
        fprintf('\n[minute %.0f] Agent attacked', k); pause(0.25); fprintf('.'); pause(0.25); fprintf('.'); pause(0.25); fprintf('.'); pause(0.25); fprintf('\n\n')
        totalattacks = totalattacks+1; %count attack
    elseif boolAttacked == 1 && agentLocation(1) == nestLocation(1) && agentLocation(2) == nestLocation(2)
        boolAttacked = 0; %if agent gets back to the nest, reset attack boolean
    end
    
    %count exploration minutes
    if agentLocation(1) > nestLocation(1)
        totalminutesExpl = totalminutesExpl + 1;
    elseif agentLocation(1) == nestLocation(1)
        boolEsc = 0; %reset boolean for escape
    end
    
    
    %% Visualize model if boolFig is true
    if boolFig
        if k == 1
            a = subplot(4,3,1);
            if length(matScat) > 0
                b1 = scatter(matScat(:,1), matScat(:,2), 'r');
            end
            hold on
            b2 = scatter(matScatFake(:,1), matScatFake(:,2), [], [1 0.5 0]);
            c = scatter(agentLocation(1), agentLocation(2), 'k', 'MarkerEdgeColor' , 'k');
            hold off
            xlim([0 1500])
            ylim([0 1500])
            d = text(1200,200,sprintf('Minute %.0f', k));
            rectangle('Position', [0 0 50 50], 'EdgeColor', 'g')
            rectangle('Position', [100 100 200 200], 'EdgeColor', 'k')
            
            subplot(4,3,4:6)
            e = plot(vecEnergy);
            hold on
            f = plot(vecFitness, 'g');
            ylim([0 100])
            xlim([0 200])
            set(gca,'XTick', [])
            ylabel(sprintf('\\color{blue}Energy\n\\color{green}Fitness'))
            xlabel('Last 200 minutes')
            
            subplot(4,3,7:9)
            g = plot(valueForage, 'r');
            hold on
            line([0 200], [0 0], 'Color', 'g');
            i = plot(fear_delayed, 'Color', [1 0.5 0]);
            xlim([0 200])
            set(gca,'XTick', [])
            ylabel(sprintf('\\color{red}V_f_o_r_a_g_i_n_g\n\\color{green}Threshold'))
            xlabel('Last 200 minutes')
            
            subplot(4,3,10:12)
            i = plot(fear_delayed, 'Color', [1 0.5 0]);
            xlim([0 200])
            set(gca,'XTick', [])
            ylabel(sprintf('\\color[rgb]{1, 0.5, 0}V_p_r_e_d_a_t_o_r'))
            xlabel('Last 200 minutes')
            
            drawnow
        else
            if length(matScat) > 0
                b1.XData = matScat(:,1);
                b1.YData = matScat(:,2);
            end
            b2.XData = matScatFake(:,1);
            b2.YData = matScatFake(:,2);
            c.XData = agentLocation(1);
            c.YData = agentLocation(2);
            d.String = sprintf('Minute %.0f', k);
            if k < 200
                e.XData = 1:length(vecEnergy);
                e.YData = vecEnergy;
                f.XData = 1:length(vecFitness);
                f.YData = vecFitness;
                g.XData = 1:length(vecValueForage);
                g.YData = vecValueForage;
                i.XData = 1:length(vecFear);
                i.YData = vecFear;
            else
                e.XData = 1:200;
                e.YData = vecEnergy(k-199:k);
                f.XData = 1:200;
                f.YData = vecFitness(k-199:k);
                g.XData = 1:200;
                g.YData = vecValueForage(k-199:k);
                i.XData = 1:200;
                i.YData = vecFear(k-199:k);
                
            end
            drawnow
        end
        
        if fitness == 0
            a.Color = [1 0 0];
        end
        
    end
    
    
end

%remove nan from data vectors
vecEnergy = vecEnergy(~isnan(vecEnergy));
vecFitness = vecEnergy(~isnan(vecEnergy));
matLocation = matLocation(find(~isnan(matLocation(:,1))),:);
vecValueFood = vecValueFood(~isnan(vecValueFood));
vecFear = vecFear(~isnan(vecFear));
vecValueForage = vecValueForage(~isnan(vecValueForage));

%print summary
fprintf('\n********************************************************************************\n')
fprintf('\nAgent died after %.2f days\n', k/(60*24))
fprintf('- Total predator attacks: %.0f\n', totalattacks)
fprintf('- Total feeding bouts: %.0f\n\n', totalfeedingbouts)
fprintf('Behavioral parameters:\n')
fprintf('- Fraction time exploring: %.2f\n', totalminutesExpl/k)
fprintf('- Avg nest exits per day: %.2f\n', totalexits/(k/(60*24)))
fprintf('- Fraction nest exits that were escaped: %.2f\n\n', totalescapes/totalexits)
fprintf('Vectors with data:\n')
fprintf('- vecEnergy: contains energy level of agent per minute\n')
fprintf('- vecFitness: contains fitness level of agent per minute\n')
fprintf('- matLocation: contains coordinates (x,y location) of agent per minute\n')
fprintf('- vecValueFood: value of food per minute\n')
fprintf('- vecFear: subjective fear level per minute\n')
fprintf('- vecValueForage: net value of foraging per minute [value food - anxiety - fear]\n')
fprintf('\n********************************************************************************\n')

%% functions

function coordinates = generatePredatorCoordinates(totalminutes, speed, range);

coordinates = nan(totalminutes,2);

%initiate position animal coordinates
xpredold = range(1)+round((range(3)-range(1))*rand);
ypredold = range(2)+round((range(4)-range(2))*rand);
count = 0;

while count < totalminutes
    
    
    %new coordinates
    xprednew = range(1)+round((range(3)-range(1))*rand);
    yprednew = range(2)+round((range(4)-range(2))*rand);
    
    %distance
    dist = sqrt((xpredold-xprednew)^2+(ypredold-yprednew)^2);
    steps = dist/speed;
    stepsize = dist/steps;
    
    for j = 1 : steps
        count = count + 1;
        newcoors = [xpredold+j*(xprednew-xpredold)/steps, ypredold+j*(yprednew-ypredold)/steps];
        coordinates(count,:) = newcoors;
    end
    
    xpredold = xprednew;
    ypredold = yprednew;
    
    
end

coordinates = coordinates(1:totalminutes,:);

end
