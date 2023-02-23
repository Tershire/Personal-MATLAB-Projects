function [] = Multi_DOF_Vib ()
close all; clear all; clc

gifON = false;
gifON = true;

%% PLOT SETTINGS
clr1 = [0      0.4470 0.7410];
clr2 = [0.8500 0.3250 0.0980];
clr3 = [0.9290 0.6940 0.1250];
clr4 = [0.4940 0.1840 0.5560];

szMarker = 32;

%% SYSTEM PARAMETERS
    % Initial Conditions
        yo    = [3; 0; 0];
        yDOTo = [0; 0; 0];
        
    % M, B, K Values
        mS = [1 1 1];          % [m1, m2, m3, ...] 
        bS = [1 1 1 1];        % [b1, b2, b3, b4, ...]
        kS = [3 1 1 3];        % [k1, k2, k3, k4, ...]

        N = length(mS);
        M = zeros(N);
        B = zeros(N);
        K = zeros(N);
        for n = 1:N     
            M(n, n)   = mS(n);

            B(n, n)   = bS(n) + bS(n+1);
            B(n, n+1) = -bS(n+1);
            B(n+1, n) = -bS(n+1);        
            B = B(1:N, 1:N);

            K(n, n)   = kS(n) + kS(n+1);
            K(n, n+1) = -kS(n+1);
            K(n+1, n) = -kS(n+1);        
            K = K(1:N, 1:N);
        end
        M;
        B;
        K;
    
    % Steady-State Distance between Masses
        dist = 9;
    
%% MAIN (FREE VIBRATION: B = 0)

    syms w y(t) q(t)

    % ---------------------------------------------------------------------------------------------------------
    % Solve (K - w.^2.*M)v = 0 to Get w and v
        A = M^-1 * K;
        
        [V, Osq] = eig(A);     % V: Eigen Vectors, Osq: Eigen Values (Square Matrix of w.^2 s)
        
        O = sqrt(Osq);
        wS = diag(O)
     
    % ---------------------------------------------------------------------------------------------------------    
    % Coordinate Conversion (y --> q)
        % m_ = v'*M*v
        % v_ = v./sqrt(m_)
        % V_ = [v_1, v_2, ...]
        
        m_S = [];
        V_ = [];
        for n = 1:N
            v = V(:, n); 
            
            m_ = v'*M*v;
            m_S = [m_S, m_];
            
            v_ = v./sqrt(m_);
            V_ = [V_, v_];
        end
        m_S;
        V_;
        
            % CHECK!
                if V_'*M*V_ == eye(N)
                    disp("V_'*M*V_ == eye(N), Check!")
                end
                if V_'*K*V_ == Osq
                    disp("V_'*M*V_ == O_sq,   Check!")
                end
        
        % Apply Initial Conditions
            qo    = V_^-1 * yo;
            qDOTo = V_^-1 * yDOTo;

        % Solve DE of q(t) Coordinate
        syms q1(t) q2(t) q3(t)
        
            q = [q1; q2; q3];
            Dq = diff(q);
            
            ode = diff(q, t, 2) + Osq * q == 0;

            % Solve the matrix equation using dsolve. Simplify the solution by using the simplify function.
            [q1Sol(t), q2Sol(t), q3Sol(t)] = dsolve(ode);
            
            % The constants C1 and C2 appear because no conditions are specified.
            % Solve the system with the initial conditions u(0) = 2 and v(0) = –1. When specifying equations in matrix form, you must specify initial conditions in matrix form too. dsolve finds values for the constants that satisfy these conditions.

            condo    =  q(0) == qo;
            condDOTo = Dq(0) == qDOTo;
            
            Conds = [condo condDOTo];
            
%             [q1Sol(t), q2Sol(t), q3Sol(t)] = dsolve(ode)
            [q1Sol(t), q2Sol(t), q3Sol(t)] = dsolve(ode, Conds);

            q1Sol(t) = simplify(q1Sol(t));
            q2Sol(t) = simplify(q2Sol(t));
            q3Sol(t) = simplify(q3Sol(t));
            
            q = [q1Sol(t); q2Sol(t); q3Sol(t)];
    
    % ---------------------------------------------------------------------------------------------------------    
    % Coordinate Conversion (q --> y)
        % y = V_*q; 
        y = V_ * q;
                
%% ANIMATE!!    
    nameFileGIF = ['Multi_DOF_3' ' yo= ' num2str(yo') ' yDOTo= ' num2str(yDOTo') ' mS= ' num2str(mS) ' kS= ' num2str(kS) ' .gif'];

    % Apply Steady-State Distance Between Masses
    for n = 1:N
        yn = y(n);
        y(n) = yn + n.*dist;
    end

    % Time Range
        if gifON
            tRange = linspace(0, 60, 300);    % for GIF 
        else
            tRange = linspace(0, 30, 1200);
        end

    % -------------------------------------------------------------------------------------
    % MAIN
        hold on
        p = plot(yo, 0, 'o', 'MarkerSize', szMarker, 'MarkerFaceColor', 'k');
        hold off
        axis([0, +(N+1).*dist, -0.1, +0.1])
        % axis tight manual

        k = 1;
        for tNow = tRange
            yATtNow = double(subs(y, t, tNow));

            set(p, {'XData'}, num2cell(yATtNow));

            drawnow

            if gifON
            % @@@ ANIMATED GIF @@@
                % Capture the plot as an image 
                  frame = getframe(gcf); 
                  im = frame2im(frame); 
                  [imind, cm] = rgb2ind(im, 256); 

                % Write to the GIF File 
                  if k == 1 
                      imwrite(imind, cm, nameFileGIF, 'gif', 'Loopcount', inf); 
                  else 
                      imwrite(imind, cm, nameFileGIF, 'gif', 'DelayTime', 0.0, 'WriteMode', 'append'); 
                  end 
                  k = k+1;
            end
        end
         
end