#include <fstream>
#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "json.hpp"
#include "spline.h"

using namespace std;
using json = nlohmann::json;

/**
  Define some reusable data structures for points in different coordinate systems.
*/

/**
  A 2-dimensional point in the Cartesian coordinate system.
*/
typedef struct xyPoint {
  double x;
  double y;
} xyPoint;

/**
  A 2-dimensional point in the Frenet coordinate system.
*/
typedef struct frPoint {
  double s;
  double d;
} frPoint;

/**
  When solving the path planning problem, using a state machine is a common approach.
  This contains the actions the car will use to perform path planning.
*/
enum State {
  /** The car should follow the target lane with the speed limit, while keeping a safe
      distance with the car in the front. */
  KEEP_LANE,
  /** The car has an intent to move to the left lane, but it should make sure the move
      is safe. */
  PREP_LCH_LEFT,
  /** The car has an intent to move to the right lane, but it should make sure the move
      is safe. */
  PREP_LCH_RIGHT,
  /** The car should move to the left lane. */
  LCH_LEFT,
  /** The car should move to the right lane. */
  LCH_RIGHT,
  /** The car should stop. This isn't implemented nor used anywhere. */
  STOP
};

/** Allow cout-ing an enum State */
ostream& operator<<(ostream& of, State st) {
  switch (st) {
    case KEEP_LANE: return of << "KEEP_LANE"; break;
    case PREP_LCH_LEFT: return of << "PREP_LCH_LEFT"; break;
    case PREP_LCH_RIGHT: return of << "PREP_LCH_RIGHT"; break;
    case LCH_LEFT: return of << "LCH_LEFT"; break;
    case LCH_RIGHT: return of << "LCH_RIGHT"; break;
    case STOP: return of << "STOP"; break;
    default: return of; break;
  }
}

// The car's maximum safe speed is slightly below the speed limit.
constexpr double MAX_SPEED = 49.5;
// The car's maximum safe acceleration is 5 m/(s^2).
constexpr double MAX_ACCEL = 5.0;

// Some functions implemented in the project base code.

constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s) {
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.find_first_of("}");
  if (found_null != string::npos) {
    return "";
  } else if (b1 != string::npos && b2 != string::npos) {
    return s.substr(b1, b2 - b1 + 2);
  }
  return "";
}

double distance(double x1, double y1, double x2, double y2)
{
	return sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
}
int ClosestWaypoint(double x, double y, const vector<double> &maps_x, const vector<double> &maps_y)
{

	double closestLen = 100000; //large number
	int closestWaypoint = 0;

	for(int i = 0; i < maps_x.size(); i++)
	{
		double map_x = maps_x[i];
		double map_y = maps_y[i];
		double dist = distance(x,y,map_x,map_y);
		if(dist < closestLen)
		{
			closestLen = dist;
			closestWaypoint = i;
		}

	}

	return closestWaypoint;

}

int NextWaypoint(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y)
{

	int closestWaypoint = ClosestWaypoint(x,y,maps_x,maps_y);

	double map_x = maps_x[closestWaypoint];
	double map_y = maps_y[closestWaypoint];

	double heading = atan2( (map_y-y),(map_x-x) );

	double angle = abs(theta-heading);

	if(angle > pi()/4)
	{
		closestWaypoint++;
	}

	return closestWaypoint;

}

// Transform from Cartesian x,y coordinates to Frenet s,d coordinates
frPoint getFrenet(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y)
{
	int next_wp = NextWaypoint(x,y, theta, maps_x,maps_y);

	int prev_wp;
	prev_wp = next_wp-1;
	if(next_wp == 0)
	{
		prev_wp  = maps_x.size()-1;
	}

	double n_x = maps_x[next_wp]-maps_x[prev_wp];
	double n_y = maps_y[next_wp]-maps_y[prev_wp];
	double x_x = x - maps_x[prev_wp];
	double x_y = y - maps_y[prev_wp];

	// find the projection of x onto n
	double proj_norm = (x_x*n_x+x_y*n_y)/(n_x*n_x+n_y*n_y);
	double proj_x = proj_norm*n_x;
	double proj_y = proj_norm*n_y;

	double frenet_d = distance(x_x,x_y,proj_x,proj_y);

	//see if d value is positive or negative by comparing it to a center point

	double center_x = 1000-maps_x[prev_wp];
	double center_y = 2000-maps_y[prev_wp];
	double centerToPos = distance(center_x,center_y,x_x,x_y);
	double centerToRef = distance(center_x,center_y,proj_x,proj_y);

	if(centerToPos <= centerToRef)
	{
		frenet_d *= -1;
	}

	// calculate s value
	double frenet_s = 0;
	for(int i = 0; i < prev_wp; i++)
	{
		frenet_s += distance(maps_x[i],maps_y[i],maps_x[i+1],maps_y[i+1]);
	}

	frenet_s += distance(0,0,proj_x,proj_y);

	return {.s = frenet_s, .d = frenet_d};

}

// Get the D centerline coordinate for a given lane.
// The left edge of the lane will be getDForLane(d) - 2
// The right edge of the lane will be getDForLane(d) + 2
// Waypoints follow the yellow lines
// Each lane is 4 m wide, and to stay in the center, you want to add 2 m.
// The starting lane is second from left (lane = 1)
double getDForLane(int lane) {
  if (lane < 0) {
    cout << "Lane is in the opposite direction" << endl;
  }
  
  return 2 + 4 * lane;
}

// Get the lane that matches a car's Frenet coordinate D value.
int getLaneForD(double d) {
  return round(d / 4.0) - 2;
}

// Convert m/s speed to mph.
double msToMph(double ms) {
  return ms * 2.23694;
}

// Convert mph speed to m/s.
double mphToMs(double mph) {
  return mph / 2.23694;
}

// Transform world coordinates to car coordinates.
xyPoint transformToLocal(double x, double y, double car_x, double car_y, double car_theta) {
  xyPoint p;
  double dx = x - car_x;
  double dy = y - car_y;
  p.x = dx * cos(-car_theta) - dy * sin(-car_theta);
  p.y = dx * sin(-car_theta) + dy * cos(-car_theta);
  return p;
}

// Transform from Frenet s,d coordinates to Cartesian x,y
xyPoint getXY(double s, double d, const vector<double> &maps_s, const vector<double> &maps_x, const vector<double> &maps_y)
{
	int prev_wp = -1;

	while(s > maps_s[prev_wp+1] && (prev_wp < (int)(maps_s.size()-1) ))
	{
		prev_wp++;
	}

	int wp2 = (prev_wp+1)%maps_x.size();

	double heading = atan2((maps_y[wp2]-maps_y[prev_wp]),(maps_x[wp2]-maps_x[prev_wp]));
	// the x,y,s along the segment
	double seg_s = (s-maps_s[prev_wp]);

	double seg_x = maps_x[prev_wp]+seg_s*cos(heading);
	double seg_y = maps_y[prev_wp]+seg_s*sin(heading);

	double perp_heading = heading-pi()/2;

	double x = seg_x + d * cos(perp_heading);
	double y = seg_y + d * sin(perp_heading);

	return { .x = x, .y = y };

}


int main() {
  uWS::Hub h;
  
  //
  // The state of the car.
  //
  
  // The current action that the car is performing. Initially the car will keep lane.
  enum State state = KEEP_LANE;
  
  // The last PREP_LCH_X state that the car performed. This is useful later.
  enum State last_lane_change = PREP_LCH_LEFT;
  
  // The car's target lane.
  int target_lane = 1;

  // The car's speed. This will be used to generate spaced points.
  double set_speed = 0.;

  // The car's target speed. set_speed will approach this smoothly.
  double target_speed = MAX_SPEED;

  // Load up map values for waypoint's x,y,s and d normalized normal vectors
  vector<double> map_waypoints_x;
  vector<double> map_waypoints_y;
  vector<double> map_waypoints_s;
  vector<double> map_waypoints_dx;
  vector<double> map_waypoints_dy;

  // Waypoint map to read from
  string map_file_ = "../data/highway_map.csv";
  // The max s value before wrapping around the track back to 0
  double max_s = 6945.554;
  
  ifstream in_map_(map_file_.c_str(), ifstream::in);

  string line;
  while (getline(in_map_, line)) {
  	istringstream iss(line);
  	double x;
  	double y;
  	float s;
  	float d_x;
  	float d_y;
  	iss >> x;
  	iss >> y;
  	iss >> s;
  	iss >> d_x;
  	iss >> d_y;
  	map_waypoints_x.push_back(x);
  	map_waypoints_y.push_back(y);
  	map_waypoints_s.push_back(s);
  	map_waypoints_dx.push_back(d_x);
  	map_waypoints_dy.push_back(d_y);
  }
  
  h.onMessage([&map_waypoints_x,
               &map_waypoints_y,
               &map_waypoints_s,
               &map_waypoints_dx,
               &map_waypoints_dy,
               &target_lane,
               &target_speed,
               &set_speed,
               &state,
               &last_lane_change](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length, uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    //auto sdata = string(data).substr(0, length);
    //cout << sdata << endl;
    if (length && length > 2 && data[0] == '4' && data[1] == '2') {

      auto s = hasData(data);

      if (s != "") {
        auto j = json::parse(s);
        
        string event = j[0].get<string>();
        
        if (event == "telemetry") {
            // j[1] is the data JSON object
        	  // Main car's localization Data
          	double car_x = j[1]["x"];
          	double car_y = j[1]["y"];
          	double car_s = j[1]["s"];
          	double car_d = j[1]["d"];
            // Car D:
            // 0 1 2 3 4 5 6 7 8 9 10 11 12
            // =       |       |         |
          	double car_yaw = j[1]["yaw"];
          	double car_speed = j[1]["speed"];
            
            // True if we are too close to a vehicle in front.
            bool too_close = false;
            
            // Show the received telemetry data about the car.
             cout << "T x=" << car_x << " y=" << car_y << " s=" << car_s << " d=" << car_d << " yaw=" << car_yaw
                 << " speed=" << car_speed << endl;
            

          	// Previous path data given to the Planner
          	auto previous_path_x = j[1]["previous_path_x"];
          	auto previous_path_y = j[1]["previous_path_y"];
            
          	// Previous path's end s and d values
          	double end_path_s = j[1]["end_path_s"];
          	double end_path_d = j[1]["end_path_d"];
            
          	// Sensor Fusion Data, a list of all other cars on the same side of the road.
          	auto sensor_fusion = j[1]["sensor_fusion"];

          	json msgJson;

            // X coordinates of path planner points
            vector<double> ptsx;
            
            // Y coordinates of path planner points
            vector<double> ptsy;

            // X coordinate list sent to car
          	vector<double> next_x_vals;
            
            // Y coordinate list sent to car
          	vector<double> next_y_vals;
            
            xyPoint xy;
                  
            /** We append the future generated path to the already unused existing path.
              The simulator sends data on how many points it did not use from the last
              iteration of the path planner.
              This number shows how long the unused existing path is.
            */
            const int path_size = previous_path_x.size();
            
            // The reference path is the x/y point from where we start to calculate the future
            // path. This is either the car position or the last point of the previous_path.
            double ref_x = car_x;
            double ref_y = car_y;
            double ref_yaw = deg2rad(car_yaw);
            
            // We are actually working on the end of the trajectory, not the car's current location.
            // This simplifies the problem and lets us anticipate actions ahead.
            if (path_size > 0) {
              car_s = end_path_s;
              car_d = end_path_d;
            }
            
            // Keep the car at a distance from cars in the front. This is the action the car uses
            // when keeping lane, and while performing other maneuvers.
            for (int i = 0; i < sensor_fusion.size(); i++) {
              vector<double> sf_car = sensor_fusion[i];
              
              // Car data vector: [., ., ., vx, vy, s, d]
              //                   0  1  2   3   4  5  6
              double d = sf_car[6];
              
              // When one of the sensor fusion cars is in our lane:
              if (d < (getDForLane(target_lane) + 2) &&
                  d > (getDForLane(target_lane) - 2)) {
                double s = sf_car[5];
                double vx = sf_car[3];
                double vy = sf_car[4];
                double v  = sqrt(vx * vx + vy * vy);

                // extrapolate the car's S value into the future since we're working on the far end
                // of the path.
                s += path_size * .02 * v;
                
                // If the car is approximately 1.2 seconds ahead of us,
                // then take action.
                if ((s > car_s) && ((s - car_s) < 1.2 * mphToMs(car_speed))) {
                  cout << "Car in front! car_s=" << car_s << " s=" << s << " car_d=" << car_d << "d=" << d << endl;
                  // Mark that we're too close to the front of the car and can consider changing
                  // the lane.
                  too_close = true;
                  
                  // Reduce our speed to match with the car in front, but keep a +1 mph speed difference so we can
                  // edge closer to the car to find a lane changing opportunity.
                  // If we're too close to the car in front, then stop edging closer to the front car.
                  if ((s - car_s < 18)) {
                    target_speed = msToMph(v) - 5;
                  } else {
                    target_speed = msToMph(v) + 1.;
                  }
                }
              }
            }
            
            // The PREP_LCH block decided that we can change the lane, so change the lane to the left.
            // The smoother will take care of making the change look good.
            if (state == LCH_LEFT) {
              if (target_lane != 0) {
                target_lane -= 1;
                last_lane_change = PREP_LCH_LEFT;
              }
              // After the lane has been changed, move back to the KEEP_LANE state.
              state = KEEP_LANE;
            }
            
            // The PREP_LCH block decided that we can change the lane, so change the lane to the left.
            // The smoother will take care of making the change look good.
            else if (state == LCH_RIGHT) {
              if (target_lane != 2) {
                target_lane += 1;
                last_lane_change = PREP_LCH_RIGHT;
              }
              state = KEEP_LANE;
            }
            
            // Prepare to change the lane left or right.
            // Assigns the state to LCH_LEFT or LCH_RIGHT if it's safe.
            else if (state == PREP_LCH_LEFT || state == PREP_LCH_RIGHT) {
              
              // Keep the name of the state for logging purposes.
              const string sn = state == PREP_LCH_LEFT ? "Left" : "Right";
              
              // Calculate how the lane changes based on the state, either -1 or +1
              const int delta = state == PREP_LCH_LEFT ? -1 : 1;

              // If the lane change is impossible (change lane right on the rightmost lane), log and continue with
              // keeping lane.
              if ((state == PREP_LCH_LEFT && target_lane == 0) ||
                  (state == PREP_LCH_RIGHT && target_lane == 2)) {
                cout << "PREP_LCH: target lane == " << target_lane << ", no change!" << endl;
                state = KEEP_LANE;
              }
              
              // True if the car can safely swap the lane.
              bool ready_to_change = true;
              
              // Find all cars in the target lane, and if some of them are in the way, mark that it's unsafe to swap
              // the lane.
              for (int i = 0; i < sensor_fusion.size(); i++) {
                vector<double> sf_car = sensor_fusion[i];
                double s = sf_car[5];
                double vx = sf_car[3];
                double vy = sf_car[4];
                double v  = sqrt(vx * vx + vy * vy);
                
                // Extrapolate that the car will be moving linearly in frenet for as long as our path lasts.
                s += path_size * .02 * v;
                double d = sf_car[6];

                // If the car is in the target lane,
                if (d < (getDForLane(target_lane + delta) + 2) &&
                    d > (getDForLane(target_lane + delta) - 2)) {
                  
                  cout << "Car in " << sn << " lane: " << s << "; " << d << " ";
                  
                  // and the car is less than 12m ahead of or behind us, the car is in the way.
                  if (fabs(s - car_s) < 12.) {
                    cout << "IN THE WAY " << fabs(s - car_s);
                    ready_to_change = false;
                  } else {
                    cout << "OK " << fabs(s - car_s);
                  }
                  
                  cout << endl;
                }
              }
              
              
              // If it's safe to swap the lane, then change state to LCH_*.
              if (ready_to_change) {
                cout << "Go go lane " << sn << "!" << endl;
                state = state == PREP_LCH_LEFT ? LCH_LEFT : LCH_RIGHT;
              } else {
                state = KEEP_LANE;
              }
            }

            // When the lane keeper marked that we should be looking for a lane change,
            if (too_close) {
              // and we're safely in our target lane (KEEP_LANE can be in effect in the middle of a lane while still
              // settling from it)
              if (fabs(car_d - getDForLane(target_lane)) < 0.2) {
                
                // Swap to left lane if it's the only possible action
                if (target_lane == 2) {
                  state = PREP_LCH_LEFT;
                }
                // Swap to right lane if it's the only possible action
                else if (target_lane == 0) {
                  state = PREP_LCH_RIGHT;
                }
                // If we're in the middle lane, repeat the last lane change to not get stuck in a back-and-forth behind
                // two cars
                else {
                  state = last_lane_change;
                }
              }
              else {
                target_speed *= .98;
              }
            }
            // When in KEEP_LANE and no car ahead, we can safely drive max speed.
            else {
              target_speed = MAX_SPEED;
            }

            //
            // Trajectory Smoother!
            // This takes the desired target lane and speed and generates a spline that smooths the car's path and
            // generates points that follow the desired speed.
            //

            // If available, take the two last points of the previous_path to help generate a smooth spline.
            // If not available, take the car's position and extrapolate a previous point based on the car yaw.
            if (path_size < 2) {
              ptsx.push_back(car_x - cos(car_yaw));
              ptsy.push_back(car_y - sin(car_yaw));
              
              ptsx.push_back(car_x);
              ptsy.push_back(car_y);
            } else {
              const double ref_x_prev = previous_path_x[path_size - 2];
              const double ref_y_prev = previous_path_y[path_size - 2];
              
              ref_x = previous_path_x[path_size - 1];
              ref_y = previous_path_y[path_size - 1];
              
              ref_yaw = atan2(
                ref_y - ref_y_prev,
                ref_x - ref_x_prev
              );
              
              ptsx.push_back(ref_x_prev);
              ptsy.push_back(ref_y_prev);
              
              ptsx.push_back(ref_x);
              ptsy.push_back(ref_y);
            }

            // Add 3 points spaced 30 meters apart, in the desired state. This allows the spline to generate a 30m
            // curve to the desired state.
            for (int i = 1; i < 4; i++) {
              xy = getXY(car_s + i * 30, getDForLane(target_lane), map_waypoints_s, map_waypoints_x, map_waypoints_y);
              ptsx.push_back(xy.x);
              ptsy.push_back(xy.y);
            }
            
            // Transform the ptsx list to the car's coordinates.
            for (int i = 0; i < ptsx.size(); i++) {
              double dx = ptsx[i] - ref_x;
              double dy = ptsy[i] - ref_y;
              
              ptsx[i] = (dx * cos(- ref_yaw) - dy * sin(- ref_yaw));
              ptsy[i] = (dx * sin(- ref_yaw) + dy * cos(- ref_yaw));
            }
            
            // Make the spline happen
            tk::spline s;
            s.set_points(ptsx, ptsy);
            
            // Push the previous_path onto the next_path so the car still follows it.
            for (int i = 0; i < path_size; i++) {
              next_x_vals.push_back(previous_path_x[i]);
              next_y_vals.push_back(previous_path_y[i]);
            }
            
            /*
             In car-local coordinates, figure out the x & y to follow in order to follow the spline.
             
             Space each xy point out by a specified distance to make sure that when a point is attained every 0.2s,
             the speed traveled from point to point will be the target speed in m/s.
             
             
             
             y              + + + +
             ^          + +
             |      + +
             |    +                 30
             | +                  /
             |-|-|-|-|-|-|-|-|-|-|-> x
            */
            
            // Set the target x/y in car coordinates to 30 m away.
            const double target_x = 30.0;
            const double target_y = s(target_x);
            const double target_dist = sqrt(target_x * target_x + target_y * target_y);
            
            // Start creating the X/Y points to get to the target XY in car coordinates:
            double x_addon = 0.;
            
            // Create enough points to keep the next_x/y_vals exactly 50 long
            for (int i = 1; i <= (50 - path_size); i++) {
              
              // While creating the points, also increment/decrement speed to keep near the target_speed.
              // Use a speed increment that keeps us within max acceleration bounds.
              if ((target_speed - set_speed) > 0.2) {
                set_speed += MAX_ACCEL * .02;
              } else if ((target_speed - set_speed) < -0.2)  {
                set_speed -= MAX_ACCEL * .02;
              }
              
              // Find the number of points to generate in order to get 30m forward:
              // 30m / (200ms per iteration) * set speed in m/s
              double N = (target_dist / (.02 * mphToMs(set_speed)));
              
              // The next point will be 30 meters / number of points to the right
              double x_pt = x_addon + target_x / N;
              // The next y will be spline(x)
              double y_pt = s(x_pt);
              x_addon = x_pt;
              
              // Transform X and Y back to global coordinates
              double x_rel = x_pt;
              double y_rel = y_pt;
              x_pt = (x_rel * cos(ref_yaw) - y_rel * sin(ref_yaw));
              y_pt = (x_rel * sin(ref_yaw) + y_rel * cos(ref_yaw));
              
              x_pt += ref_x;
              y_pt += ref_y;
              
              next_x_vals.push_back(x_pt);
              next_y_vals.push_back(y_pt);
            }
            
            cout << "target_speed = " << target_speed << ", target_lane = " << target_lane << ", state = " << state << endl;
            
            // Push the data to the car simulator.
          	msgJson["next_x"] = next_x_vals;
          	msgJson["next_y"] = next_y_vals;

          	auto msg = "42[\"control\","+ msgJson.dump()+"]";
            
          	ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
          
        }
      } else {
        // Manual driving
        std::string msg = "42[\"manual\",{}]";
        ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
      }
    }
  });

  // We don't need this since we're not using HTTP but if it's removed the
  // program
  // doesn't compile :-(
  h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data,
                     size_t, size_t) {
    const std::string s = "<h1>Hello world!</h1>";
    if (req.getUrl().valueLength == 1) {
      res->end(s.data(), s.length());
    } else {
      // i guess this should be done more gracefully?
      res->end(nullptr, 0);
    }
  });

  h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
    std::cout << "Connected!!!" << std::endl;
  });

  h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code,
                         char *message, size_t length) {
    ws.close();
    std::cout << "Disconnected" << std::endl;
  });

  int port = 4567;
  if (h.listen(port)) {
    std::cout << "Listening to port " << port << std::endl;
  } else {
    std::cerr << "Failed to listen to port" << std::endl;
    return -1;
  }
  h.run();
}
