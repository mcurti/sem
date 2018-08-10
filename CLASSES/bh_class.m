classdef bh_class
    %BH_CLASS BH routines
    %   Detailed explanation goes here
    
    properties
        Extrap_Spline_B
        Extrap_Spline_pB
        Extrap_Spline_ppB
        Extrap_Spline_n2
        Extrap_Spline_pn2
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
                mu0 = 4*pi*1e-7; %Js = 2; mur = 1000;
%                 bh = @(H) mu0*H+2*Js/pi.*atan((pi*(mur-1)*mu0*H)./(2*Js));
%                 B = bh(H);
                mu_init = (B(1)-B(2))/((H(1)-H(2))*mu0);
                mu_end  = (B(end-1)-B(end))/((H(end-1)-H(end))*mu0);
%                 ik = B>0 & B<0.9;
%                 B(ik) = []; H(ik) = [];

                % Compute geometrical parameters
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
                
                
                % Preparing the splines
                % H(B)
                SB   = spline(btest,[1/(mu_init*mu0); htest; 1/(mu_end*mu0)]);
                SBex = fnxtr(SB,2);
                
                % H'(B)
                fp   = fnder(SBex,1);
                
                % H''(B)
                fpp  = fnder(SBex,2);
                
                % nu(B)
%                 nu = htest./btest; nu(1) = 1/(mu_init*mu0);
                
%                 Snu = spline
                
                
%                 htest = H; btest = B;
%                 bad_points = true(numel(btest));
                nutest = htest./btest; nutest(1) = 1/(mu_init*mu0);
                nu_init = (btest(2)^2-btest(3)^2)/((nutest(2)-nutest(3)));
                nu_end  = (btest(end-1)^2-btest(end)^2)/((nutest(end-1)-nutest(end)));
                
%                 while sum(bad_points)>2
%                     SNB2test     = spline(btest(2:end).^2,[1/nu_init; nutest(2:end); 1/nu_end]);
%                     SNB2extest   = fnxtr(SNB2test,2);
%                     fpptest    = fnder(SNB2extest,2);
%                     bad_points = ppval(fpptest,btest.^2)<0;
% %                     bad_points(1) = false;
%                     btest(bad_points) = []; nutest(bad_points) = [];
%                     nu_init = (btest(2)^2-btest(3)^2)/((nutest(2)-nutest(3)));
%                     nu_end  = (btest(end-1)^2-btest(end)^2)/((nutest(end-1)-nutest(end)));
% %                     figure(1); clf; plot(B.^2,ppval(fpptest,B.^2))
%                 end
                 nu2b = ppval(SBex,0:0.01:50)./(0:0.01:50); nu2b(1) = nu2b(2);
                SNB2 = spline((0:0.01:50).^2,nu2b);
%                 SNB2ex = fnxtr(SNB2,2);
                fpn  = fnder(SNB2,1);
                
            else
                warning ('Undefined material option - Please import new material data as explained in description');
            end
            
            % Assigning the classes data
            obj.Extrap_Spline_B   = SBex;
            obj.Extrap_Spline_pB  = fp;
            obj.Extrap_Spline_ppB = fpp;
            obj.Extrap_Spline_n2  = SNB2;
            obj.Extrap_Spline_pn2 = fpn;
            
        end
        
        function [Br, dmur, mur] = BHtool(obj,Bmod)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            f   = @(x) ppval(obj.Extrap_Spline_B,x);
            
            x0  = Bmod;
            df  = @(x) ppval(obj.Extrap_Spline_pB,x);
            % Output
            mu0 = 4*pi*1e-7;
            dmur = 1./df(x0)./mu0;
            nu = ppval(obj.Extrap_Spline_n2,x0.^2);
            
            mur  = 1./nu./mu0;
            Br  = x0 - f(x0)./df(x0);
        end
        function [dnu ]  = dnu(obj,Bmod)
            %METHOD2 Summary of this method goes here
            %   Detailed explanation goes here
%             f   = @(x) ppval(obj.Extrap_Spline_B,x);
            
            x0  = Bmod.^2;
            df  = @(x) ppval(obj.Extrap_Spline_pn2,x);
            % Output
            mu0 = 4*pi*1e-7;
            dnu = df(x0).*mu0;
        end
    end
end

