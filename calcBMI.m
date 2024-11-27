function output = calcBMI(prob,wl)

load('dataPandey2024.mat')

%% calculate new BMI distribution for giving prob obese adults GLP inducing wl weight loss
p1=p0;

[ d, ix30 ] = min( abs( x-30 ) );
p1(ix30:end)=(1-prob).*p0(ix30:end);

for i=1:length(x)
    if x(i)>=(1-wl)*30
        [ d, jTemp ] = min( abs( x-x(i)/(1-wl) ) );
        p1(i)=p1(i)+prob*p0(jTemp)/(1-wl);
    end
end

if abs(trapz(x,p1)-1)>.01
    display('ERROR: something went wrong with p1');
end

%% calculate BMI percentages in 7 standard categories

S1=S0;
for i=1:length(BMIthresholds)
    [ d, ixTemp ] = min( abs( x-BMIthresholds(i)) );
    S1(i)=trapz(x(ixTemp:end),p1(ixTemp:end));
end

class1=class0;
class1(1)=1-S1(1);
class1(end)=S1(end);
for i=2:length(class0)-1
    class1(i)=S1(i-1)-S1(i);
end

if abs(sum(class1)-1)>.01
    display('ERROR: something went wrong with classes');
end

%% estimating annual lives saved

deaths1=P*dot(mu0,class1);

livesSaved=deaths0-deaths1;

%% output

output=[S1(4);S1(5);S1(6);livesSaved];

end