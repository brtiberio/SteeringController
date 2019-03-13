% Copyright (c) 2019 J. Sequeira
%                    Bruno Tiberio
%                    F. Ferreira da Silva
%
% Permission is hereby granted, free of charge, to any person obtaining a
% copy of this software and associated documentation files (the
% "Software"), to deal in the Software without restriction, including
% without limitation the rights to use, copy, modify, merge, publish,
% distribute, sublicense, and/or sell copies of the Software, and to permit
% persons to whom the Software is furnished to do so, subject to the
% following conditions:
%
% The above copyright notice and this permission notice shall be included
% in all copies or substantial portions of the Software.
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
% OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
% MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN
% NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
% DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
% OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
% USE OR OTHER DEALINGS IN THE SOFTWARE.
%
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% CAR TRANSVERSAL CONTROL - STEERING ANGLE CONTROL

classdef SteeringController < handle
    properties
        %------------------------------------------------------------------
        % Constants section
        %------------------------------------------------------------------
        
        % look ahead distance
        path_look_ahead = 3;
        % look ahead number of points;
        samples_look_ahead;
        % portion of path to track
        window_start = [];
        window_end = [];
        % minimum distance from car to path
        min_distance;
        
        %------------------------------------------------------------------
        % state vector containing pose of vehicle
        %------------------------------------------------------------------
        
        % current position of car [x y] in world frame 
        car_position = [];
        % direction of car in world frame (psi)
        car_pose=[];
        % current path index
        path_index = [];
        % sampling time
        sampling_time = 0.01;
        % last wheel angle
        last_wheel_angle = 0;
        % wheel angle velociy hard limiter
        ws_limiter = deg2rad(45); %0.3; 10.9;
        % max steering angle
        max_angle = deg2rad(28);
        % selection of method.
        method = 'pure pursuit';
        % reference path points
        map_points = [];
        % length of reference path
        num_points = [];
        %------------------------------------------------------------------
        % gains related to Sequeira 
        %------------------------------------------------------------------
        Kp_dist = 0.01;
        % wheel velocity proportional gain;
        Kp_ws = 0.01;
        Kd_ws = 0.1;
        %------------------------------------------------------------------
        % Gains related to pure pursuit
        % If equal to zero, lateral distance contribution is ignored.
        %------------------------------------------------------------------
        % pure pursuit gains and cumulative sum
        Kp = 0;
        Ki = 0;
        cumsum_dist = 0;
        
        
        wheel_base = 2.2;
        last_angle_diff = 0;
        ws_smooth = 0;
        
        % steering velociy hard limiter
        rate_limiter = deg2rad(12);
        
        
    end
    properties(SetAccess = private)
        watchdog_timer;
        isUpdated=0;
    end
    
    methods
        %% Public methods
        function obj = SteeringController(varargin)
            % Create a controller for steering wheel angle using a
            % reference trajectory and current position of the car.
            
            valid_methods ={'sequeira', 'pure pursuit','stanley'};
            
            args = inputParser();
            args.addParameter('SampleTime',0.01,...
                @(x)validateattributes(x,{'numeric'},...
                {'positive','real','scalar'}));
            args.addParameter('maxAngle',deg2rad(28),...
                @(x)validateattributes(x,{'numeric'},...
                {'positive','real','scalar'}));
            args.addParameter('method','sequeira',...
                @(x)~isempty(validatestring(x,valid_methods)));
            args.addParameter('wheelBase',2.2,...
                @(x)validateattributes(x,{'numeric'},...
                {'positive','real','scalar'}));
            args.addParameter('rateLimiter',12,...
                @(x)validateattributes(x,{'numeric'},...
                {'positive','real','scalar'}));
            
            args.parse(varargin{:});
            obj.sampling_time = args.Results.SampleTime;
            obj.max_angle = args.Results.maxAngle;
            obj.method = args.Results.method;
            obj.wheel_base = args.Results.wheelBase;
            obj.rate_limiter = args.Results.rateLimiter;
            
            % TODO implement code for watchdog timer
        end
        
        function [varargout] = clip_angle(obj, angle, ws_smooth)
            %
            % Verify if the hard limits for max angle are atained
            % If so, limit it to plus minus max angle
            % Args:
            %     angle:
            %     ws_smooth:
            % Return:
            %   :angle: wheel angle in radians
            %   :ws_smooth: wheel speed in radians/s
            switch obj.method
                case 'sequeira'
                    if (abs(angle) > obj.max_angle)
                        angle = sign(angle)*obj.max_angle;
                        % do not try to go over hardlimiter
                        ws_smooth = 0;
                    end
                    varargout = {angle, ws_smooth};
                case 'pure pursuit'
                    if (abs(angle) > obj.max_angle)
                        angle = sign(angle)*obj.max_angle;
                    end
                    varargout = {angle};
                otherwise
                    % TODO
            end
            
        end
        function wheel_speed = clip_wheel_speed(obj, speed)
            %
            % clip speed between max values.
            % Args:
            %     speed: velocity of wheel speed in radians.
            % Return: 
            %     speed in radians.
            %
            wheel_speed= sign(speed)*min(abs(speed), obj.ws_limiter);
        end
        
        
        
        function load_map(obj, map)
            %
            % Load map points as a matrix
            % Args:
            %     map: A matrix containing X and Y positions as column
            %          vector
            %
            obj.map_points= map;
            obj.num_points= length(obj.map_points);
        end
        function load_car(obj, x_start, y_start, pose_start)
            %
            % Load car location and orientation in world frame. 
            %
            % Args:
            %     x_start: Initial position in X coordinate.
            %     y_start: Initial position in Y coordinate.
            %     pose_start: Initial orientation in radians relative to
            %                 world frame.
            %
            %
            obj.car_pose=pose_start;
            obj.car_position = [x_start y_start];
        end
        
        function find_min(obj)
            %--------------------------------------------------------------------------
            % Find minimum lateral distance between initial position and reference
            % trajectory
            %--------------------------------------------------------------------------
            if isempty(obj.path_index)
                % if it is the first time here, find the closest point in
                % reference trajectory.
                [min_distance, index] = min(vecnorm(obj.map_points-repmat(obj.car_position,obj.num_points,1),2,2));
                obj.path_index = index;
                cumsum = 0;
                index = 1;
                % find the number of points required to fullfill the look
                % ahead distance. 
                while cumsum< obj.path_look_ahead && obj.path_index+index < obj.num_points
                    cumsum = cumsum+norm(obj.map_points(obj.path_index+index,:) - obj.map_points(obj.path_index+index-1,:));
                    index = index +1;
                end
                obj.samples_look_ahead = index-1;
            else
                % after first point being discovered, search only for a
                % limited window. If index is greater than one, slide
                % window to front.
                [min_distance, index] = min(vecnorm(obj.map_points(obj.window_start:obj.window_end,:)-repmat(obj.car_position,obj.window_end-obj.window_start+1,1),2,2));
                obj.path_index = obj.path_index+index-1;
            end
            % get window for look ahead distance
            obj.window_start = obj.path_index;
            obj.window_end = min(obj.path_index+obj.samples_look_ahead, obj.num_points);
            obj.min_distance = min_distance;
            
        end
        
        function [varargout] = update(obj)
            switch obj.method
                case 'sequeira'
                     % for easy understanding
                     psi = obj.car_pose;
                     obj.find_min();
                     %------------------------------------------------------
                     % get vector from car to look ahead
                     % get angle between look ahead and car
                     %------------------------------------------------------
                     direction_look_ahead = obj.map_points(obj.window_end,:)- obj.car_position;
                     angle_look_ahead = atan2(direction_look_ahead(2),direction_look_ahead(1));
                     % get difference from two angles (vector look ahead and
                     % car direction).
                     angle_diff = (angle_look_ahead-psi);
                     angle_diff = obj.normalize_angle(angle_diff);
                     
                     
                     % compute derivative term
                     dot_angle_diff = (angle_diff - obj.last_angle_diff) / obj.sampling_time;
                     
                     dist_to_look_ahead = norm(direction_look_ahead);
                     % steering wheel speed
                     ws = obj.Kp_dist*dist_to_look_ahead + obj.Kd_ws*dot_angle_diff + obj.Kp_ws*angle_diff;
                     ws_smooth = 0.9*obj.ws_smooth + 0.1*ws;
                     % limit max velocity
                     ws_smooth = obj.clip_wheel_speed(ws_smooth);
                     new_wheel_angle = obj.last_wheel_angle + obj.sampling_time*ws_smooth;
                     % add to be changed!
                     % can not try to force past the limits
                     [new_wheel_angle, ws_smooth] = obj.clip_angle(new_wheel_angle, ws_smooth);
                     obj.ws_smooth = ws_smooth;
                     obj.last_wheel_angle = new_wheel_angle;
                     varargout = {new_wheel_angle, ws_smooth};
                case 'stanley'
                    % todo
                case 'pure pursuit'
                    psi = obj.car_pose; % angle of the car in world frame
                    obj.find_min();     % find closest point of trajectory
                    obj.cumsum_dist = obj.cumsum_dist + obj.min_distance*obj.sampling_time;
                    %------------------------------------------------------
                    % get vector from car to look ahead
                    % get angle between look ahead and car
                    %------------------------------------------------------
                    direction_look_ahead = obj.map_points(obj.window_end,:)- obj.car_position;
                    angle_look_ahead = atan2(direction_look_ahead(2),direction_look_ahead(1));
                    % get difference from two angles (vector look ahead and
                    % car direction).
                    angle_diff = (angle_look_ahead-psi); 
                    angle_diff = obj.normalize_angle(angle_diff);
                    % distance of look ahead vector
                    norm_look_ahead = norm(direction_look_ahead, 2);
                    
                    % curvature of circle to reach look ahead point
                    curvature = 2*sin(angle_diff)/norm_look_ahead;
                    % angle required to adjust car curvature to path
                    % see ref: http://dx.doi.org/10.4218/etrij.15.0114.0123
                    curvature_angle = atan(curvature*obj.wheel_base);
                    
                    
                    % get angle based in lateral distance of rear wheel
                    dist_angle = sign(angle_diff)*(obj.Kp*obj.min_distance +obj.Ki*obj.cumsum_dist);
                    
                    new_wheel_angle = obj.clip_angle(curvature_angle+dist_angle);
                    % calculate rate of change to avoid destroying the
                    % mechanical parts. Assume same rate for "falling" and
                    % "rising".
                    rate = (new_wheel_angle-obj.last_wheel_angle)/obj.sampling_time;
                    if abs(rate) > obj.rate_limiter
                        new_wheel_angle = sign(rate)*obj.sampling_time*obj.rate_limiter+obj.last_wheel_angle; 
                    end
                    obj.last_wheel_angle = new_wheel_angle;
                    varargout = {new_wheel_angle};    
                otherwise
                    error(strcat('Unknown method: ',obj.method));
            end
        end
    end
    methods(Access = private)
        function enableWatchdog(obj)
            % TODO
            obj.watchdog_timer.start();
        end
        function disableWatchdog(obj)
            % TODO
            obj.watchdog_timer.stop();
        end
        
    end
    methods(Static)
        function angle = normalize_angle(angle)
            %
            % Normalize an angle to [-pi, pi].
            % Args:
            %     angle: angle in radians to be limited
            % Return:
            %     Angle in radian in [-pi, pi]
            %
            while angle > pi
                angle = angle - 2 * pi;
            end
            while angle < -pi
                angle = angle + 2 * pi;
            end
        end
    end
end

