%-----------------------------------------------------------------------
%    software    :   DotX Control Software
%    author      :   DotX Control Solutions BV
%                    Netherlands
%                    www.dotxcontrol.com
%                    info@dotxcontrol.com
%-----------------------------------------------------------------------
%   PID Controller 
%-----------------------------------------------------------------------
%  PID     = dotxPidController(Kp, Ti, Td, Ts)      : Creates a dotx PID Controller Object
%  [u, du] = PID.pidNext(error)                     : Get the PID output for the next step
% [tf_Pid] = PID.getTransferFunction()              : Returns the transfer function of the PID
%-----------------------------------------------------------------------

classdef dotxPidController < handle    
    properties
        aciStatus = 0                 % Initialisation Status
        
        du_min    = -100000      % Constraints
        du_max    = +100000      % Constraints
        u_min     = -100000      % Constraints
        u_max     = +100000      % Constraints

        Kp        = 1
        Ti        = 99999
        Td        = 0
        Ts        = 1
        u_prev    = 0                 % Previous Output of the PIDs
        i         = 0;
        c
        a

        error_f1    = 0;    % Internal State of the PID
        error_f2_m1 = 0;    % Internal State of the PID
    end
    
    % --------------------------------------------------------------------%
    % Constructor
    % --------------------------------------------------------------------%    
    methods
        function o = dotxPidController(Kp, Ti, Td, Ts)
            if nargin > 0
                o.Kp = Kp;
                o.Td = Td;
                o.Ti = Ti;        
                o.Ts = Ts;
                
                % Set Filter coefficients
                o.a = 0.1;
                o.c = exp(-o.Ts/(o.a*o.Td));
            end
        end
    end
    
    % --------------------------------------------------------------------%
    % Methods
    % --------------------------------------------------------------------%    
    methods
        
        % ----------------------------------------------------------------%
        % PID Next function
        % ----------------------------------------------------------------%    
        function [u, du] = pidNext(o, error)
            % computes discrete-time PID control action, with
            % PID in velocity form. Let k be the discrete-time, i.e. k = 0, 1, 2, ...
            % INPUTS:
            % er = control error
            % aciStatus = 0 or 1. If 0, then initialise, else, continue
            %
            % OUTPUTS:
            % du = u(k+1) - u(k)
            % xnext(1) = er(k
            error_f1 = o.error_f1;
            
            % Below is an implementation of a digital lead-lag filter:
            % tfLeadLag = tf([o.Td 1],[o.a*o.Td + 1]);
            if(o.Td > 0)
                error_f1 =   o.c*error_f1  + (1.0-o.c)    *error;
                error_f2   = 1.0/o.a*error + (1.0-1.0/o.a)*error_f1;
            else
                error_f2 = error;
            end
            
            dError           = error_f2 - o.error_f2_m1;

            % PID Equation
            du               = o.Kp* ( (dError) + (o.Ts/o.Ti)*o.error_f2_m1 );

            o.error_f2_m1 = error_f2;
            o.error_f1    = error_f1;
            
            %--------------------------------------------
            % Apply speed constraints
            %--------------------------------------------
            du = min(o.du_max*o.Ts,du);
            du = max(o.du_min*o.Ts,du); 

            u  = o.u_prev + du;

            %--------------------------------------------
            % Apply absolute constraints
            %--------------------------------------------
            u  = min(u,o.u_max); 
            u  = max(u,o.u_min);

            du = u - o.u_prev; %recompute Delta u 

            o.u_prev = u;
            %--------------------------------------------
            % Update PID state for next call
            %--------------------------------------------
            %o.x(1) = Error;
            %o.x(2) = errorFilt;
        end
        
        
        % ----------------------------------------------------------------%
        % Get PID Transfer Function
        % ----------------------------------------------------------------%    
        function [tf_Pid] = getTransferFunction(o)                        
            tf_Lead     = tf([o.Td 1],[o.a*o.Td 1]);           
            tf_Integral = tf([o.Ti 1],[o.Ti 0]);
            tf_Pid      = o.Kp * (tf_Lead * tf_Integral);                  
        end
        
        
        % ----------------------------------------------------------------%
        % Autotune #1: Relay Feedback
        % ----------------------------------------------------------------%    
        % TBD
        
        % ----------------------------------------------------------------%
        % Autotune #2: SIMC
        % ----------------------------------------------------------------%    
        function [Kp, Ti, Td] = autoTuneSIMC( o, Model, ControllerType )
            if nargin == 1
                ControllerType = 'PI'               
            end
                        
            Kp = 0;
            Ti = 0;
            Td = 0;
            
            Model.Tdt = Model.Tdt + o.Ts;
            % In case of a PI Controller
            if( strcmp(ControllerType, 'PI' ) )
                Kp = 0.35 *  Model.Tp / (  Model.Tdt * Model.K ) ;
                Ti = min(8*Model.Tdt, Model.Tp);
                Td = 0;
            end
            
            % In case of a PID Controller
            if( strcmp(ControllerType, 'PID' ))
                Kp = 0.5 * Model.Tp / (Model.Tdt*Model.K);
                Ti = min( 8*Model.Tdt, Model.Tp);
                Td = 0.5*Model.Tdt;
            end
            
            o.Kp = Kp;
            o.Ti = Ti;
            o.Td = Td;
            
        end
            
        
        
        
    end
    
end

