function obj = Flow_HA1(SWCr, SWCs, a,n, Ks, l, SWCfc, ...
    D1, D2, S, R, l_s, l_d, ln_s, ln_d, ...
     Hleaf, RH, T, Va, day)


% Alternate the unit of plant resistance, 10 min cm-2
Ra_r1=R.Ra_r1*10^4/600;
Ra_r1_2=R.Ra_r1_2*10^4/600;
Ra_r2=R.Ra_r2*10^4/600;
Ra_s1=R.Ra_s1*10^4/600;
Ra_s2=R.Ra_s2*10^4/600;
Ra_b1=R.Ra_b1*10^4/600;
Ra_b2=R.Ra_b2*10^4/600;
Rc_r=R.Rc_r*10^4/600;
Rc_r_2=R.Rc_r2*10^4/600;
Rc_s=R.Rc_s*10^4/600;
Rr_s=R.Rr_s*10^4/600;
Rr_s_d=R.Rr_s*10^4/600*2;
Rr_d=R.Rr_d*10^4/600;

% Set the initial SWC as the saturated SWC
obj.SWC1(1, 1)=SWCs;
obj.SWC2(1, 1)=SWCs;

% Iterative simulations over time
      for i=1:1:144*day % each time iteration is 10 min and there are 144 iterations in a day
          
          % Soil calculation
          obj.Ks1(i, 1) = swc2ks(obj.SWC1(i, 1), SWCr, SWCs, n, Ks, l); % soil hydraulic conductivity of shallow layer, cm day-1
          obj.Ks2(i, 1) = swc2ks(obj.SWC2(i, 1), SWCr, SWCs, n, Ks, l); % soil hydraulic conductivity of deep layer, cm day-1
          obj.Rs1(i, 1)=1/obj.Ks1(i, 1); % soil hydraulic resistivity of shallow layer, day cm-1
          obj.Rs2(i, 1)=1/obj.Ks2(i, 1); % soil hydraulic resistivity of deep layer, day cm-1
          obj.H1(i, 1)=-swc2h(obj.SWC1(i, 1), SWCr, SWCs, a, n); % soil head pressure of shallow layer, cm
          obj.H2(i, 1)=-swc2h(obj.SWC2(i, 1), SWCr, SWCs, a, n); % soil head pressure of deep layer, cm
          obj.Rrs1(i, 1)=obj.Rs1(i, 1)*ln_s/2/pi/l_s*24*6; % soil resistance of shallow root, day cm-2
          obj.Rrs1_2(i, 1)=obj.Rs2(i, 1)*ln_d/2/pi/l_s*2*24*6; % soil resistance of deep root from shallow root, day cm-2
          obj.Rrs2(i, 1)= obj.Rs2(i, 1)*ln_d/2/pi/l_d*24*6; % soil resistance of deep root, day cm-2
          
          % Water flow between soil layers
          qs1=(obj.H2(i, 1)-obj.H1(i, 1)-D1/2-D2/2)/(obj.Rs1(i, 1)*D1/2+obj.Rs2(i, 1)*D2/2); % water flow from shallow layer to deep layer, cm day-1
          obj.qs1(i, 1)=qs1/24/6; % alternate the unit, cm 10min-1
          qs2=(0-obj.H2(i, 1)-D2/2)/(obj.Rs2(i, 1)*D2/2); % water flow from deep layer to ground water table, cm day-1
          obj.qs2(i, 1)=qs2/24/6; % alternate the unit, cm 10min-1
          
          % Calculate soil evaporation
          obj.Esoil(i, 1)=Esoil(obj.SWC1(i, 1), SWCr, SWCfc, Va(i, 1), RH(i, 1), T(i, 1));
          
          %Solve the water flow of network
          syms q1 q2 q3 p1 p2 p3; % cm3 10min-1
          eqns = [-q1*(2*obj.Rrs1(i, 1)+2*Rr_s+Ra_r1)+p1*Rc_r+q2*(2*obj.Rrs1(i, 1)+2*Rr_s+Ra_r1)==0,...
               -(q1+p1)*(Ra_r1+Ra_s1/2)+p3*Rc_s+(q3-p2)*(Ra_s2/2+Ra_r2/2)-p2*Rc_r_2+(q2-p1)*(Ra_r1+Ra_r1_2/2)-p1*Rc_r==0,...
               -q3*(obj.Rrs2(i, 1)+Rr_d+Ra_r2/2)-p2*Rc_r_2-(q2-p1+p2)*(obj.Rrs1_2(i, 1)+Rr_s_d+Ra_r1_2/2)==0,...
               (q3-p2-p3)*(Ra_s2/2+Ra_b2)-p3*Rc_s-(q1+p1+p3)*(Ra_s1/2+Ra_b1)==0,...
              (obj.H2(i, 1)-D1-D2/2)-q3*(obj.Rrs2(i, 1)+Rr_d+Ra_r2/2)-(q3-p2)*(Ra_r2/2+Ra_s2/2)-(q3-p2-p3)*(Ra_s2/2+Ra_b2)-Hleaf(i, 1)==0,...
              (obj.H2(i, 1)-D1-D2/2)+(q2-p1+p2)*(obj.Rrs1_2(i, 1)+Rr_s_d+Ra_r1_2/2)+(q2-p1)*(Ra_r1+Ra_r1_2/2)+q2*(2*obj.Rrs1(i, 1)+2*Rr_s+Ra_r1)-(obj.H1(i, 1)-D1/2)==0,...
              (obj.H1(i, 1)-D1/2)-q1*(2*obj.Rrs1(i, 1)+2*Rr_s+Ra_r1)-(q1+p1)*(Ra_r1+Ra_s1/2)-(q1+p1+p3)*(Ra_s1/2+Ra_b1)-Hleaf(i, 1)==0];
          Q=solve(eqns);
          obj.q1(i, 1)=double(Q.q1); obj.q2(i, 1)=double(Q.q2); obj.q3(i, 1)=double(Q.q3);
          obj.p1(i, 1)=double(Q.p1); obj.p2(i, 1)=double(Q.p2); obj.p3(i, 1)=double(Q.p3);

              
          % root surface head pressure (water potential), cm
          obj.H_int_1(i, 1)=obj.H1(i, 1)-obj.q1(i, 1)*2*obj.Rrs1(i, 1);
          obj.H_int_2(i, 1)=obj.H1(i, 1)-obj.q2(i, 1)*2*obj.Rrs1(i, 1);
          
          % root xylem head pressure (water potential), cm
          obj.H_root_1(i, 1)=obj.H1(i, 1)-obj.q1(i, 1)*(2*obj.Rrs1(i, 1)+2*Rr_s+Ra_r1);
          obj.H_root_2(i, 1)=obj.H1(i, 1)-obj.q2(i, 1)*(2*obj.Rrs1(i, 1)+2*Rr_s+2*Ra_r1);
          
          % solve the soil water content in next iteration
          obj.SWC1(i+1, 1)=obj.SWC1(i, 1)+obj.qs1(i, 1)/D1-obj.Esoil(i, 1)/D1-(obj.q1(i, 1)+obj.q2(i, 1))/D1/S;
          obj.SWC2(i+1, 1)=obj.SWC2(i, 1)-obj.qs1(i, 1)/D2+obj.qs2(i, 1)/D2-obj.q3(i, 1)/D2/S;
          
          
          i=i+1
       end
end


