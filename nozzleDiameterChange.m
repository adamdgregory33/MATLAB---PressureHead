clc
clear
clf('reset')
cla reset
close all;
pump = [0	,0,	22;
2,	0.002,	17.43;
3.4,	0.0034,	11.73;
4.17,	0.00417,	7.77;
4.61,	0.00461,	5.31;
5.19,	0.00519,	1.69;
5.41,	0.00541,	0.29];


%Pump Data, in form: Flow rate L/s, Flow rate m^3/s, Head

pumpEquation = polyfit(pump(:,2),pump(:,3),3);
%generate equation for pump head vs flow rate

xValsPump = [0:0.0001:0.0055];
%get x values between 0 and 5.5 L/s or 0.0055 m^3/s
yValsPump = polyval(pumpEquation,xValsPump);

figure(1)
plot(pump(:,2),pump(:,3),xValsPump,yValsPump);
axis([0 0.006 0 30])
hold on


torquePower = csvread('rpmvstorquevspower.csv');
%imports torque data in format rpm, torque, power

rpmTorqueEquation = polyfit(torquePower(:,1),torquePower(:,2),1);
%Get an equation for torque vs rpm all in SI units

rpmPowerEquation = polyfit(torquePower(:,1),torquePower(:,3),1);
xValsRpm = linspace(0,4000);
%get an equation for power vs Rpm all in SI units

yValsRpmTorque = polyval(rpmTorqueEquation,xValsRpm);
yValsRpmPower = polyval(rpmPowerEquation,xValsRpm);


%generate the y values using this equation

reynoldsNum = 0;


initialDiameter = 0.018;
%inital diameter in metres
decrement = 0.0005;
finalDiameter = initialDiameter - decrement;
finalDiameter = 0.014;
%final diameter in metres

cd = 0.95;
%pipe losses

gravity = 9.81;
%gtravity

bucketAngle = (141/180)*pi;
%angle water is turned through

density = 1000;
%density of water

radius = 0.1;
%Radius of water wheel


size = size(pump);
rows = size(1);
%Gets the number of rows of data in the 
finished = false;

operatingPoint = zeros(2);

yValsFinal = zeros(rows,1);
%create a variable with enough slots 


counter = 0;
maxPower = 0;
efficiency = 0;

while finished == false
    
    counter = counter+1;
    finalArea = (pi * (finalDiameter^2))/4;
    disp("Nozzle Diameter "+finalDiameter);
    
    figure(5)
    title('Final Area')
    plot(counter,finalArea,'ro');
    hold on
    
    
    beta = finalDiameter/initialDiameter;
    disp("Beta :"+beta);
    
    deltaH = zeros(rows,1);
    for i =rows:-1:1
        deltaH(i) = (1/(2*gravity))*(1-(beta^4))*(pump(i,2)/(cd*finalArea))^2;
        %calculate the pressure head for each corresponding pump point
    end
    
    pressureEquation = polyfit(pump(:,2),deltaH(:),3);
    tempYVals = polyval(pressureEquation,xValsPump);
    %Generate the trendline and the points from this trendline
    
    figure(1)
    plot(xValsPump,tempYVals);
    hold on
    operatingPoint = InterX([xValsPump;yValsPump],[xValsPump;tempYVals]);
    disp("Flow Rate m3/s: "+operatingPoint(1));
    title('Pressure vs flow rate')
    %Find the intercept of the pressureEquation and the pump equation
    
    tempVelocity = operatingPoint(1) / finalArea;
    %velocity of the output jet is calculated
    
    disp("Velocity+ "+tempVelocity);
    
    figure(6)
    title('velocity')
    plot(counter,tempVelocity,'ro');
    hold on
    
    massFlow = tempVelocity * density * finalArea;
    %Mass flow rate at the end of the outlet (rho * u * A)
    disp("MassFlow+ "+ massFlow);
    
    torque =(1-cos(bucketAngle))* massFlow * radius * tempVelocity;
    figure(7)
    title('torque')
    plot(counter,torque,'ro');
    hold on
    
    %torque equation given in tutorial: (1-cos theta) * m\dot * Radius *
    %Velocity bucket
    
    rpm = (60/(2*pi))*(tempVelocity/radius);
    
    
    %RPM = (v/r) * (60/(2*pi) for omega in RPM
    disp("RPM :"+rpm);
    
    linearTorque = polyfit([0,rpm],[torque,0],1);   
    tempYValsTorque = polyval(linearTorque,xValsRpm);
    %create a line connecting stall torque and no load for the water wheel
    
    
    rpmTorqueCrosspoint = InterX([xValsRpm;yValsRpmTorque],[xValsRpm;tempYValsTorque]);
    %find the point where the linear torque of wheel crosses the motor
    %curve to find operating point
    
    if(isempty(rpmTorqueCrosspoint) == 0)
        figure(3)
        plot(xValsRpm,yValsRpmTorque,xValsRpm,tempYValsTorque);
        hold on
        title('Torque vs rpm')
        axis([0 2000 0 15])

        tempLinePower = zeros(100,2);
        tempLinePower(:,1) =  rpmTorqueCrosspoint(1);
        tempLinePower(:,2) = linspace(0,max(yValsRpmPower));
        %create a vertical line of the RPM that 

        rpmPowerCrosspoint = InterX([tempLinePower(:,1).';tempLinePower(:,2).'],[xValsRpm;yValsRpmPower]);
        %find where the motor curve and this line meet


        if(isempty(rpmPowerCrosspoint) ==0 && (finalDiameter - decrement) > 0 )
            
            tempPower = rpmPowerCrosspoint(2);
            figure(9)
            plot(counter,tempPower,'ro')
            title('Power')
            hold on
            
            
            efficiency = tempPower/(0.5*massFlow*tempVelocity^2);
            reynoldsNum = (tempVelocity * finalDiameter) / (1.0035*10^-6);
            figure(8)
            plot(counter,reynoldsNum,'ro')
            title('reynolds num')
            hold on
            
            %{
            disp(reynoldsNum);
            disp(efficiency);
            disp(tempPower);
            %}
            finalDiameter = finalDiameter - decrement;
            figure(2)
            title('power vs rpm')
            axis([0 2000 0 max(yValsRpmPower)])
            plot(tempLinePower(:,1),tempLinePower(:,2),xValsRpm,yValsRpmPower);
            hold on
            finished = true;

        else
            
            yValsFinal = tempYVals;
            disp("No power cross point");
            finished = true;
        end
        
    else
        finished = true;
        yValsFinal = tempYVals;
        disp("No torque cross point");

    end
    disp(" ");
end

hold off

%disp(finalDiameter);

