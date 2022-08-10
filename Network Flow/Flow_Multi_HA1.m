function obj = Flow_Multi_HA1(SWCr1, SWCs1, a1, n1, Ks1,...
    SWCr2, SWCs2, a2, n2, Ks2, l, SWCfc, ...
    D1, D2, S, R, l_s, l_d, ln_s, ln_d, ...
     Hleaf, RH, T, Va, day)

 % The meaning of all codes is similar to Flow_HA1.m
 
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

obj.SWC1(1, 1)=SWCs1;
obj.SWC2(1, 1)=SWCs2;

      for i=1:1:144*day
          
          obj.Ks1(i, 1) = swc2ks(obj.SWC1(i, 1), SWCr1, SWCs1, n1, Ks1, l);
          obj.Ks2(i, 1) = swc2ks(obj.SWC2(i, 1), SWCr2, SWCs2, n2, Ks2, l);
          obj.Rs1(i, 1)=1/obj.Ks1(i, 1);
          obj.Rs2(i, 1)=1/obj.Ks2(i, 1);
          obj.H1(i, 1)=-swc2h(obj.SWC1(i, 1), SWCr1, SWCs1, a1, n1);
          obj.H2(i, 1)=-swc2h(obj.SWC2(i, 1), SWCr2, SWCs2, a2, n2);
          obj.Rrs1(i, 1)=obj.Rs1(i, 1)*ln_s/2/pi/l_s*24*6;
          obj.Rrs1_2(i, 1)=obj.Rs2(i, 1)*ln_d/2/pi/l_d*2*24*6;
          obj.Rrs2(i, 1)= obj.Rs2(i, 1)*ln_d/2/pi/l_d*24*6;
          
          qs1=(obj.H2(i, 1)-obj.H1(i, 1)-D1/2-D2/2)/(obj.Rs1(i, 1)*D1/2+obj.Rs2(i, 1)*D2/2);
          obj.qs1(i, 1)=qs1/24/6;
          qs2=(0-obj.H2(i, 1)-D2/2)/(obj.Rs2(i, 1)*D2/2);
          obj.qs2(i, 1)=qs2/24/6;

          obj.Esoil(i, 1)=Esoil(obj.SWC1(i, 1), SWCr1, SWCfc, Va(i, 1), RH(i, 1), T(i, 1));
          
          syms q1 q2 q3 p1 p2 p3;                                                          
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

          obj.H_int_1(i, 1)=obj.H1(i, 1)-obj.q1(i, 1)*2*obj.Rrs1(i, 1);
          obj.H_int_2(i, 1)=obj.H1(i, 1)-obj.q2(i, 1)*2*obj.Rrs1(i, 1);
          
          obj.H_root_1(i, 1)=obj.H1(i, 1)-obj.q1(i, 1)*(2*obj.Rrs1(i, 1)+2*Rr_s+Ra_r1)-(obj.q1(i, 1)+obj.p1(i, 1))*Ra_r1;
          obj.H_root_2(i, 1)=obj.H1(i, 1)-obj.q2(i, 1)*(2*obj.Rrs1(i, 1)+2*Rr_s+2*Ra_r1)-(obj.q2(i, 1)-obj.p1(i, 1))*Ra_r1;
          
          obj.SWC1(i+1, 1)=obj.SWC1(i, 1)+obj.qs1(i, 1)/D1-obj.Esoil(i, 1)/D1-(obj.q1(i, 1)+obj.q2(i, 1))/D1/S;
          obj.SWC2(i+1, 1)=obj.SWC2(i, 1)-obj.qs1(i, 1)/D2+obj.qs2(i, 1)/D2-obj.q3(i, 1)/D2/S;
          
          
          i=i+1
       end
end


