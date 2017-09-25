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

// for convenience
using json = nlohmann::json;

/** X/Y point */
typedef struct xyPoint {
  double x;
  double y;
} xyPoint;

/** Frenet point */
typedef struct frPoint {
  double s;
  double d;
} frPoint;

/** Car state */
enum State {
  KEEP_LANE,
  PREP_LCH_LEFT,
  PREP_LCH_RIGHT,
  LCH_LEFT,
  LCH_RIGHT,
  STOP,
  START
};

// The car's maximum allowed speed in mph
constexpr double MAX_SPEED = 49.5;
// The car's maximum allowed acceleration in m/s2
constexpr double MAX_ACCEL = 5.0;


// For converting back and forth between radians and degrees.
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


// Get the D coordinate for a given lane.
// The leftmost lane is 0
// Waypoints follow the yellow lines
// Each lane is 4 m wide, and to stay in the center, you want to add 2 m.
// The starting lane is second from left (lane = 1)
double getDForLane(int lane) {
  if (lane < 0) {
    cout << "Lane is in the opposite direction" << endl;
  }
  
  return lane * 4 + 2;
}

int getLaneForD(double d) {
  return round(d / 4.0) - 2;
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
  
  // Car state:
  
  // The car's target lane.
  int target_lane = 1;
  // The car's target speed. 50mph is the limit.
  double target_speed = 0.;

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
  
  h.onMessage([&map_waypoints_x,&map_waypoints_y,&map_waypoints_s,&map_waypoints_dx,&map_waypoints_dy,&target_lane,&target_speed](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
                     uWS::OpCode opCode) {
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
            
            // True if we are too close to the particular vehicle.
            bool too_close = false;
            
            // Telemetry!
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
            
            // int current_lane = getLaneForD(car_d);

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
                        
            const int path_size = previous_path_x.size();
            
            double ref_x = car_x;
            double ref_y = car_y;
            double ref_yaw = deg2rad(car_yaw);
            
            // Simple case of KEEP LANE:
            // Is there a car in front of us?
            // Rather... Is there a car in front of our path?
            
            if (path_size > 0) {
              car_s = end_path_s;
            }
            
            for (int i = 0; i < sensor_fusion.size(); i++) {
              vector<double> sf_car = sensor_fusion[i];
              
              // Car vector: [., ., ., vx, vy, s, d]
              //              0  1  2   3   4  5  6
              double d = sf_car[6];
              
              if (d < (2 + 4*target_lane + 2) && d > (2 + 4 * target_lane - 2)) {
                    double s = sf_car[5];
                    double vx = sf_car[3];
                    double vy = sf_car[4];
                    double v  = sqrt(vx * vx + vy * vy);
                    
                    // Extrapolate that the car will be moving linearly in frenet for as long as our path lasts.
                    s += path_size * .02 * v;
                    
                    // If the car is in front of the path and if it's too close,
                    // reduce the car's target velocity.
                    if ((s > car_s) && ((s - car_s) < 30)) {
                      cout << "Car in front! cars=" << car_s << " s=" << s << " card=" << car_d << "d=" << d << endl;
                      too_close = true;
                      
                      if (target_lane != 0) {
                        target_lane--;
                      } else if (target_lane != 2) {
                        target_lane++;
                      }
                    }
              }
            }
            
            if (too_close) {
              target_speed -= .224;
            } else if (target_speed < 49.5) {
              target_speed += .224;
            }
            
            
            // Grab some historical data for the car.
            // This will help in creating a better spline.
            // xref_prev---xref---CAR--->x--->x--->
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

            // Add 3 future points spaced 30m apart, for spline generation
            for (int i = 1; i < 4; i++) {
              xy = getXY(car_s + i * 30, getDForLane(target_lane), map_waypoints_s, map_waypoints_x, map_waypoints_y);
              ptsx.push_back(xy.x);
              ptsy.push_back(xy.y);
            }
            
            // Transform the entire current ptsx list to car-local coordinates
            for (int i = 0; i < ptsx.size(); i++) {
              /*xyPoint t = transformToLocal(ptsx[i], ptsy[i], ref_x, ref_y, ref_yaw);
              ptsx[i] = t.x;
              ptsy[i] = t.y; */
              
              double dx = ptsx[i] - ref_x;
              double dy = ptsy[i] - ref_y;
              
              ptsx[i] = (dx * cos(- ref_yaw) - dy * sin(- ref_yaw));
              ptsy[i] = (dx * sin(- ref_yaw) + dy * cos(- ref_yaw));
            }
            
            //
            // Smooth the resulting path via spline.h
            //
            
            tk::spline s;
            s.set_points(ptsx, ptsy);
            
            
            //
            // Generate the future path:
            // - Add unused historical points onto the future path
            // - Add enough future points to keep future path at 50 long
            //
            
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
            
            const double target_x = 30.0;
            const double target_y = s(target_x);
            const double target_dist = sqrt(target_x * target_x + target_y * target_y);
            
            // Keeps track of the X coordinate over the loop.
            double x_addon = 0.;
            
            for (int i = 1; i <= (50 - path_size); i++) {
              // Find the number of points to generate in order to get 30m forward:
              // 30m / (200ms per iteration) * ()(target velocity in mph / 2.24) = (target vel in m/s)).
              double N = (target_dist / (.02 * target_speed / 2.24));
              
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
            
          	msgJson["next_x"] = next_x_vals;
          	msgJson["next_y"] = next_y_vals;

          	auto msg = "42[\"control\","+ msgJson.dump()+"]";
            
          	//this_thread::sleep_for(chrono::milliseconds(1000));
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
