classdef bh_class
    %BH_CLASS BH routines
    %   Detailed explanation goes here
    
    properties
        Extrap_Spline_B
        Extrap_Spline_pB
    end
    
    methods
        function obj = bh_class(material_ID)
            %BH_CLASS preparing the BH data
            %   Detailed explanation goes here
            if      material_ID==1
%                 [H,B] = textscan('BH_Cogent_M800_50A_50Hz.txt','%f');
%                 Bsat=(B(end-1)+B(end))./2; % Saturation field, no tangent further
%                 
%                 Bmin=0.1379;
%                     Interpolation
%                 mu0 = 4*pi*1e-7;
%                 SB  = spline(B,H);
%                 x0  = Bmod;
%                 x0(x0>Bsat) = Bsat;
%                 x0(x0<Bmin) = Bmin;
            elseif  material_ID==2
                fileID = fopen('BH_Cogent_M270_35A_50Hz.txt');
                out = textscan(fileID,'%f %f');  H = out{1}; B = out{2};
                fclose(fileID);
                mu0 = 4*pi*1e-7;
                mu_init = (B(1)-B(2))/((H(1)-H(2))*mu0);
                mu_end  = (B(end-1)-B(end))/((H(end-1)-H(end))*mu0);
%                 ik = B>0 & B<0.9;
%                 B(ik) = []; H(ik) = [];
                htest = H; btest = B;
                bad_points = true(numel(btest));
                
                while sum(bad_points)>0
                    SBtest     = spline(btest,[1/(mu_init*mu0); htest; 1/(mu_end*mu0)]);
                    SBextest   = fnxtr(SBtest,2);
                    fpptest    = fnder(SBextest,2);
                    bad_points = ppval(fpptest,btest)<0; 
                    bad_points(1) = false;
                    btest(bad_points) = []; htest(bad_points) = [];
                mu_init = (btest(1)-btest(2))/((htest(1)-htest(2))*mu0);
                mu_end  = (btest(end-1)-btest(end))/((htest(end-1)-htest(end))*mu0);
                end
                
                
                
                SB   = spline(btest,[1/(mu_init*mu0); htest; 1/(mu_end*mu0)]);
                SBex = fnxtr(SB,2);
                fp  = fnder(SBex,1);
                
            else
                warning ('Undefined material option - Please import new material data as explained in description');
            end
            
            % Assigning the classes data
            obj.Extrap_Spline_B  = SBex;
            obj.Extrap_Spline_pB = fp;
            
        end
        
        function [Br, mur] = BHtool(obj,Bmod)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            f   = @(x) ppval(obj.Extrap_Spline_B,x);
            
            x0  = Bmod;
            df  = ppval(obj.Extrap_Spline_pB,x0);
            % Output
            mu0 = 4*pi*1e-7;
            mur = 1./df./mu0;
            Br  = x0 - f(x0)./df;
        end
    end
end

