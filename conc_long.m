    %Anil Cengiz - Oct 19,2024
    %This function generates LONG time courses of drug concentration based on the PK
    %parameters of a drug 
    
    function [time,conc]=conc_long(par_struct,tau,max_time)
        
        %For Semaglutide 
        if par_struct.name=="semaglutide" 
        ka=par_struct.ka;
        V=par_struct.V;
        ke=0.0488/V; %clearance/elimination rate (1/hour)
        F=par_struct.F;
        weekly_dose=par_struct.weekly_dose;

        n_dose=floor(max_time/tau);
        dosing_times=0:tau:n_dose*tau; %Dosing times
        time=0:1:(n_dose*tau+tau);
        c=zeros(size(time)) ; %concentration holder
        g=zeros(size(time)); %drug amount holder 
        D=weekly_dose*ones(1,length(dosing_times));
    
    D(1:4)=0.25;
    D(5:8)=0.5;
    D(9:12)=1;
    D(13:16)=1.7;
    D(17:end)=2.4;
        
        %Loop over dosing times 
        for j=1:length(dosing_times)
            i=dosing_times(j)+1;
            %Update drug amount at absorption site
            g(i)=g(i)+D(j)*F;
            y0=[g(i),c(i)];    
    
            %Flow ODE to find g and c during 1 dosing cycle 
            tspan=0:1:tau;
            [t,y] = ode45(@(t,y) odefcn_semag(t,y,V,ka,ke), tspan, y0);
            g(i:i+tau)=y(:,1);
            c(i:i+tau)=y(:,2);
        end
    
        %Create a vector that has the timecourse for the whole drug amount (at absorption site) and drug concentration 
        %conc=[c(1:(end-1))];
        %time=time(1:end-1);
        conc=c;
        
        
        %For Tirzepatide 
        elseif par_struct.name=="tirzepatide"
        ka=par_struct.ka;
        Vc=par_struct.Vc;
        Vp=par_struct.Vp;
        cl=par_struct.cl;
        q=par_struct.q;
        F=par_struct.F;
        weekly_dose=par_struct.weekly_dose;

        n_dose=floor(max_time/tau);
        dosing_times=0:tau:n_dose*tau; %Dosing times
        time=0:1:(n_dose*tau+tau)    ;    
        c=zeros(size(time)) ; %holder for drug concentration in the central compt.
        cp=zeros(size(time)) ; %holder for drug concentration in the peripheral compt.
        g=zeros(size(time)); %holder for drug concentration at absorption site
        D=weekly_dose;
        %Tirzepatide uses a dose escalation scheme in which, dose is increased
        %by 2.5 mg every 4 weeks (equivalently, every 4 doses)
    
        
    
        %Loop over dosing times 
        for j=1:length(dosing_times)
            i=dosing_times(j)+1;
        
            %Update drug at absorption site
            g(i)=g(i)+D*F;
            y0=[g(i),c(i),cp(i)];
            %Flow ODE to find g and c and cp during 1 dosing cycle 
            tspan=0:1:tau;
            [~,y] = ode45(@(t,y) odefcn_tirz(t,y,Vc,Vp,cl,q,ka), tspan, y0);
        
            g(i:i+tau)=y(:,1); %drug amount
            c(i:i+tau)=y(:,2); %central conc
            cp(i:i+tau)=y(:,3); %peri. conc
        end
        
            %Create a vector that has the timecourse for the whole drug amount
            %(at absorption site) and drug concentration at central and
            %peripheral compartments
       % conc=[c(1:(end-1))];
       conc=c;
      % time=time(1:end-1);

        end
        
        
        function dydt = odefcn_tirz(~,y,Vc,Vp,cl,q,ka)
        %function dydt = odefcn(t,y,Vc,Vp,cl_time,cl_vec,q,ka)
        %ind=find(t>=cl_time*7*24,1,'last');
        %cl=cl_vec(ind);
        
          dydt = zeros(3,1);
          dydt(1)=-ka*y(1);
          dydt(2) = ka*y(1)/Vc - cl*y(2)/Vc - q/Vc*y(2) + q/Vp*y(3);
          dydt(3)=  q/Vc*y(2) - q/Vp*y(3);
        end
        
        function dydt = odefcn_semag(~,y,V,ka,ke)
          dydt = zeros(2,1);
          dydt(1)=-ka*y(1);
          dydt(2) = ka*y(1)/V - ke*y(2);
        end
        end