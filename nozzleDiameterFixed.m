clc
clear
clf('reset')
cla reset
close all;

%This one basically finds the gear ratio that provides the highest RPM

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



torquePower = csvread('rpmvstorquevspower.csv');
%imports torque data in format rpm, torque, power

rpmTorqueEquation = polyfit(torquePower(:,1),torquePower(:,2),1);
%Get an equation for torque vs rpm all in SI units

rpmPowerEquation = polyfit(torquePower(:,1),torquePower(:,3),2);
xValsRpm = linspace(0,4000);
%get an equation for power vs Rpm all in SI units

yValsRpmTorque = polyval(rpmTorqueEquation,xValsRpm);
yValsRpmPower = polyval(rpmPowerEquation,xValsRpm);
%generate the y values using this equation

reynoldsNum = 0;




initialDiameter = 0.018;
%inital diameter in metres
finalDiameter = 0.014;
%final diameter in metres

cd = 1;
gravity = 9.81;
bucketAngle = (141/180)*pi;
density = 1000;
radius = 0.092;
counter = 0;
maxPower = 0;
efficiency = 0;


gearRatio = 44/30;

targetMotorRPM = 4000;

size = size(pump);
rows = size(1);
%Gets the number of rows of data in the 

operatingPoint = zeros(2);

yValsFinal = zeros(rows,1);
%create a variable with enough slots 


counter = counter+1;
finalArea = (pi * (finalDiameter^2))/4;
disp("Nozzle Diameter "+finalDiameter);


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




operatingPoint = InterX([xValsPump;yValsPump],[xValsPump;tempYVals]);
disp("Flow Rate m3/s: "+operatingPoint(1));

figure(1)
plot(xValsPump,yValsPump,xValsPump,tempYVals);
title("Figure 1: Hydrostatic Head vs Flow rate")
xlabel("Flow Rate (m^3/s)");
ylabel("Hydrostatic Head (m)");
legend("Pump Operating Curve","Nozzle Operating Curve")
axis([0 0.006 0 30])
grid on;

%Find the intercept of the pressureEquation and the pump equation

tempVelocity = operatingPoint(1) / finalArea;
%velocity of the output jet is calculated

disp("Velocity+ "+tempVelocity);


massFlow = tempVelocity * density * finalArea;
%Mass flow rate at the end of the outlet (rho * u * A)
disp("MassFlow+ "+ massFlow);

torque =((1-cos(bucketAngle))* massFlow * radius * tempVelocity)/gearRatio;
disp("Stall Torque : "+torque);

%torque equation given in tutorial: (1-cos theta) * m\dot * Radius *
%Velocity bucket

rpm = gearRatio*(60/(2*pi))*(tempVelocity/radius);


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
    legend("Motor Torque RPM curve", "Turbine Operting Curve");
    hold on;
    xlabel("RPM");
    ylabel("Torque (Nm)");
    title('Figure 3: Torque vs RPM of turbine');
    axis([0 2000 0 15])
    grid on;

    tempLinePower = zeros(100,2);
    tempLinePower(:,1) =  rpmTorqueCrosspoint(1);
    tempLinePower(:,2) = linspace(0,max(yValsRpmPower));
    %create a vertical line of the RPM that 

    rpmPowerCrosspoint = InterX([tempLinePower(:,1).';tempLinePower(:,2).'],[xValsRpm;yValsRpmPower]);
    %find where the motor curve and this line meet


    if(isempty(rpmPowerCrosspoint) ==0)

        tempPower = rpmPowerCrosspoint(2);
        disp("Power = "+tempPower);


        efficiency = tempPower/(0.5*massFlow*tempVelocity^2);
        reynoldsNum = (tempVelocity * finalDiameter) / (1.0035*10^-6);

        figure(2)

        plot(tempLinePower(:,1),tempLinePower(:,2),xValsRpm,yValsRpmPower);
        title('Figure 4: Power Output vs RPM');
        xlabel('RPM')
        ylabel('Power (W)')
        legend("Turbine RPM Value","Motor Operating Curve");
        axis([0 2000 0 1000]);
        grid on;
        
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


