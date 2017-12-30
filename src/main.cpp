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
#include <sys/time.h>
#include "spline.h"
#include <stdlib.h>
#include <cmath>

using namespace std;

// for convenience
using json = nlohmann::json;

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

	double heading = atan2((map_y-y),(map_x-x));

	double angle = fabs(theta-heading);
  angle = min(2*pi() - angle, angle);

  if(angle > pi()/4)
  {
    closestWaypoint++;
  if (closestWaypoint == maps_x.size())
  {
    closestWaypoint = 0;
  }
  }

  return closestWaypoint;
}

// Transform from Cartesian x,y coordinates to Frenet s,d coordinates
vector<double> getFrenet(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y)
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

	return {frenet_s,frenet_d};

}

// Transform from Frenet s,d coordinates to Cartesian x,y
vector<double> getXY(double s, double d, const vector<double> &maps_s, const vector<double> &maps_x, const vector<double> &maps_y)
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

	double x = seg_x + d*cos(perp_heading);
	double y = seg_y + d*sin(perp_heading);

	return {x,y};

}

bool is_adjacent_lane(int lane_id, int current_lane)
{
  return lane_id + 1 == current_lane || current_lane + 1 == lane_id;
}

// The max s value before wrapping around the track back to 0
double max_s = 6945.554;
double ref_vel = 0.0;
int lane = 1;

struct AdjacentCarInfo {
  int front_car_id = -1;
  int front_car_distance = -1;
  int back_car_id = -1;
  int back_car_distance = -1;
};

struct LaneStatusInfo {
  bool can_change_lane = true;
  double current_front_car_distance = 0;
  double future_front_car_distance = 0;
};

struct CarInstruction {
  int change_to_lane = -1;
  double ref_vel = 0;
};

int get_lane_id_from_d(double d) {
  // Each lane's width is 4
  if (d < 4.0) {
    return 0;
  } else if (d < 8.0) {
    return 1;
  } else {
    return 2;
  }
}

// For each lane, fetch two cars, one in front of my car and the other following the car
// assume my car is in the lane with the position of car_s
vector<AdjacentCarInfo> get_adjacent_car_per_lane(
    const vector<vector<double>>& sensor_fusion,
    double car_s) {
  vector<AdjacentCarInfo> infos;
  infos.resize(3);
  for (int i=0; i < sensor_fusion.size(); ++i) {
    double check_car_d = sensor_fusion[i][6];
    double check_car_s = sensor_fusion[i][5];
    int lane_id = get_lane_id_from_d(check_car_d);
    if (check_car_s >= car_s) { // front
      double dist = check_car_s - car_s;
      if (infos[lane_id].front_car_id < 0 ||
          dist < infos[lane_id].front_car_distance) {
        infos[lane_id].front_car_id = i;
        infos[lane_id].front_car_distance = dist;
      }
    } else { // back
      double dist = car_s - check_car_s;
      if (infos[lane_id].back_car_id < 0 ||
          dist < infos[lane_id].back_car_distance) {
        infos[lane_id].back_car_id = i;
        infos[lane_id].back_car_distance = dist;
      }
    }
  }

  for (int i=0; i<3; ++i) {
    auto& info = infos[i];
    if (info.front_car_id < 0) {
      info.front_car_distance = 400;
    }
    if (info.back_car_id < 0) {
      info.back_car_distance = 400;
    }
  }

  cout << "adjacent_car_info: "
       << "0: front: " << infos[0].front_car_id << " " << infos[0].front_car_distance
       << ", back: " << infos[0].back_car_id << " " << infos[0].back_car_distance
       << ", 1: front: " << infos[1].front_car_id << " " << infos[1].front_car_distance
       << ", back: " << infos[1].back_car_id << " " << infos[1].back_car_distance
       << ", 2: front: " << infos[2].front_car_id << " " << infos[2].front_car_distance
       << ", back: " << infos[2].back_car_id << " " << infos[2].back_car_distance
       << endl;
  return infos;
}

vector<LaneStatusInfo> get_lane_status_infos(
    const vector<AdjacentCarInfo>& adjacent_infos,
    const vector<vector<double>>& sensor_fusion,
    double car_speed,
    double car_s) {
  vector<LaneStatusInfo> lane_status_infos;
  lane_status_infos.resize(3);
  for (int i=0; i<3; ++i) {
    const auto& adjacent_info = adjacent_infos[i];

    auto& lane_status_info = lane_status_infos[i];
    if (adjacent_info.front_car_distance < 5 || adjacent_info.back_car_distance < 5) {
      lane_status_info.can_change_lane = false;
    }

    lane_status_info.current_front_car_distance = adjacent_info.front_car_distance;

    // check front car
    if (adjacent_info.front_car_id >= 0) {
      double vx = sensor_fusion[adjacent_info.front_car_id][3];
      double vy = sensor_fusion[adjacent_info.front_car_id][4];
      double check_car_speed = sqrt(vx*vx + vy*vy);
      double check_car_s = sensor_fusion[adjacent_info.front_car_id][5];

      // future distance in 1.5 seconds
      lane_status_info.future_front_car_distance = \
        (check_car_s + check_car_speed * 1.5) - (car_s + car_speed * 1.5);
      if (lane_status_info.future_front_car_distance <= 0) {
        lane_status_info.can_change_lane = false;
      }
    } else {
      // no front car
      lane_status_info.future_front_car_distance = lane_status_info.current_front_car_distance;
    }

    // check back car
    if (adjacent_info.back_car_id >= 0) {
      double vx = sensor_fusion[adjacent_info.back_car_id][3];
      double vy = sensor_fusion[adjacent_info.back_car_id][4];
      double check_car_speed = sqrt(vx*vx + vy*vy);
      double check_car_s = sensor_fusion[adjacent_info.back_car_id][5];

      double future_back_car_distance = \
        (car_s + car_speed * 1.5) - (check_car_s + check_car_speed * 1.5);
      if (future_back_car_distance < 0) {
        lane_status_info.can_change_lane = false;
      }
    }
  }

  cout << "lane_status_infos: "
       << "0: can_change_lane: " << lane_status_infos[0].can_change_lane
       << " front_dist: " << lane_status_infos[0].current_front_car_distance
       << " future front_dist: " << lane_status_infos[0].future_front_car_distance
       << ", 1: can_change_lane: " << lane_status_infos[1].can_change_lane
       << " front_dist: " << lane_status_infos[1].current_front_car_distance
       << " future front_dist: " << lane_status_infos[1].future_front_car_distance
       << ", 2: can_change_lane: " << lane_status_infos[2].can_change_lane
       << " front_dist: " << lane_status_infos[2].current_front_car_distance
       << " future front_dist: " << lane_status_infos[2].future_front_car_distance
     << endl;
  return lane_status_infos;
}

int best_score_lane(const vector<double>& lane_scores, int current_lane) {
  // initialize best_lane_id as current_lane, so the best lane's score has to be bigger
  int best_lane_id = current_lane;
  for (int lane_id = 0; lane_id < lane_scores.size(); ++lane_id) {
    if (lane_scores[best_lane_id] < lane_scores[lane_id]) {
      best_lane_id = lane_id;
    }
  }
  return best_lane_id;
}

int next_lane_to_drive(
    const vector<LaneStatusInfo>& lane_status_infos,
    int current_lane,
    double car_speed) {

  const auto& current_lane_info = lane_status_infos[current_lane];

  if (current_lane_info.current_front_car_distance <= 5 ||
      car_speed <= 20) {
    // no enough front space to change lane
    return current_lane;
  }

  vector<double> lane_scores;
  lane_scores.resize(3);
  for (int lane_id = 0; lane_id < 3; ++lane_id) {
    lane_scores[lane_id] = lane_status_infos[lane_id].future_front_car_distance;

    if (current_lane != lane_id) {
      if (!lane_status_infos[lane_id].can_change_lane) {
        // cannot change the lane
        lane_scores[lane_id] = -1000;
      } else {
        if (is_adjacent_lane(lane_id, current_lane)) {
          lane_scores[lane_id] -= 20;
        } else {
          if (lane_status_infos[1].can_change_lane) {
            lane_scores[lane_id] -= 40;
          } else {
            // block it, as it has to pass the middle lane
            lane_scores[lane_id] = -1000;
          }

        }
      }
    }
  }

  int change_to_lane = best_score_lane(lane_scores, current_lane);
  cout << "best_score_lane: "
     << "current lane: "<< current_lane
     << ", best score lane: "<< change_to_lane
     << ", 0: score " << lane_scores[0]
     << ", 1: score: " << lane_scores[1]
     << ", 2: score: " << lane_scores[2]
     << endl;

  if (change_to_lane == current_lane ||
      lane_scores[change_to_lane] <= 0) {
    // no need to change the lane
    return current_lane;
  }

  if (!is_adjacent_lane(change_to_lane, current_lane)) {
    // has to go through the middel lane to reach it
    change_to_lane = 1;
  }

  // calculate the lane score to determine if it can change lane
  if (current_lane_info.future_front_car_distance <= 10) {
    // too close to the front car, change the lane now
    return change_to_lane;
  }

  static int previous_best_lane = -1;
  static int accumulate_count = 0;
  if (change_to_lane == previous_best_lane && accumulate_count > 30) {
    previous_best_lane = -1;
    accumulate_count = 0;
    return change_to_lane;
  }

  if (change_to_lane == previous_best_lane) {
    ++accumulate_count;
    // wait for more round of the same estimation
    return current_lane;
  }

  // wait for more round of estimation
  previous_best_lane = change_to_lane;
  accumulate_count = 1;
  return current_lane;
}

CarInstruction update_car_instruction(
    const vector<LaneStatusInfo>& lane_status_infos,
    double ref_vel,
    int current_lane) {
  CarInstruction instruction;

  // shall the car change lane
  instruction.change_to_lane = next_lane_to_drive(lane_status_infos, current_lane, ref_vel);

  if (instruction.change_to_lane == current_lane &&
      lane_status_infos[current_lane].future_front_car_distance < 0) {
    instruction.ref_vel = ref_vel - 0.224;
  } else if (instruction.change_to_lane == current_lane &&
             lane_status_infos[current_lane].future_front_car_distance < 10) {
    instruction.ref_vel = ref_vel;
  } else if (ref_vel < 49.5) {
    instruction.ref_vel = ref_vel + 0.224;
  } else {
    instruction.ref_vel = ref_vel;
  }

  cout << "update_car_instruction: "
       << "current lane: "<< current_lane
       << " next lane: " << instruction.change_to_lane
       << ", old ref_vel: "<< ref_vel
       << ", new ref_vel: " << instruction.ref_vel
       << endl;

  return instruction;
}

int main() {
  uWS::Hub h;

  // Load up map values for waypoint's x,y,s and d normalized normal vectors
  vector<double> map_waypoints_x;
  vector<double> map_waypoints_y;
  vector<double> map_waypoints_s;
  vector<double> map_waypoints_dx;
  vector<double> map_waypoints_dy;

  // Waypoint map to read from
  string map_file_ = "../data/highway_map.csv";

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


  h.onMessage([&map_waypoints_x,&map_waypoints_y,&map_waypoints_s,&map_waypoints_dx,&map_waypoints_dy](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
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
        	double car_yaw = j[1]["yaw"];
        	double car_speed = j[1]["speed"];

        	// Previous path data given to the Planner
        	auto previous_path_x = j[1]["previous_path_x"];
        	auto previous_path_y = j[1]["previous_path_y"];
        	// Previous path's end s and d values
        	double end_path_s = j[1]["end_path_s"];
        	double end_path_d = j[1]["end_path_d"];

        	// Sensor Fusion Data, a list of all other cars on the same side of the road.
        	auto sensor_fusion = j[1]["sensor_fusion"];

          int prev_size = previous_path_x.size();

          if (prev_size == 0) {
            end_path_s = car_s;
          }

          vector<AdjacentCarInfo> adjacent_infos = get_adjacent_car_per_lane(
            sensor_fusion, car_s);

          vector<LaneStatusInfo> lane_status_infos = get_lane_status_infos(
            adjacent_infos, sensor_fusion, car_speed, car_s);

          CarInstruction instruction = update_car_instruction(
            lane_status_infos, ref_vel, lane);
          if (instruction.change_to_lane >= 0) {
            lane = instruction.change_to_lane;
          }
          if (instruction.ref_vel >= 0) {
            ref_vel = instruction.ref_vel;
          }

        	json msgJson;

          car_s = end_path_s;
        	// TODO: define a path made up of (x,y) points that the car will visit sequentially every .02 seconds
          // Reference x, y, yam status
          // either refer the starting location where the car is or previous path's end points
          double ref_x = car_x;
          double ref_y = car_y;
          double ref_yaw = deg2rad(car_yaw);

          vector<double> ptsx;
          vector<double> ptsy;

          // if the previous size is almost empty, use the car's starting point
          if (prev_size < 2) {
            double prev_car_x = car_x - cos(car_yaw);
            double prev_car_y = car_y - sin(car_yaw);

            ptsx.push_back(prev_car_x);
            ptsx.push_back(car_x);
            ptsy.push_back(prev_car_y);
            ptsy.push_back(car_y);
          } else {
            // Refine the reference points as the previous car's position
            ref_x = previous_path_x[prev_size - 1];
            ref_y = previous_path_y[prev_size - 1];

            double ref_x_prev = previous_path_x[prev_size - 2];
            double ref_y_prev = previous_path_y[prev_size - 2];
            ref_yaw = atan2(ref_y - ref_y_prev, ref_x - ref_x_prev);

            ptsx.push_back(ref_x_prev);
            ptsx.push_back(ref_x);
            ptsy.push_back(ref_y_prev);
            ptsy.push_back(ref_y);
          }

          // In Frenet, add evenly 30m spaced points ahead of starting reference
          vector<double> next_wp0 = getXY(car_s + 30, (2+4*lane), map_waypoints_s, map_waypoints_x, map_waypoints_y);
          vector<double> next_wp1 = getXY(car_s + 60, (2+4*lane), map_waypoints_s, map_waypoints_x, map_waypoints_y);
          vector<double> next_wp2 = getXY(car_s + 90, (2+4*lane), map_waypoints_s, map_waypoints_x, map_waypoints_y);

          ptsx.push_back(next_wp0[0]);
          ptsx.push_back(next_wp1[0]);
          ptsx.push_back(next_wp2[0]);

          ptsy.push_back(next_wp0[1]);
          ptsy.push_back(next_wp1[1]);
          ptsy.push_back(next_wp2[1]);

          for (int i=0; i<ptsx.size(); ++i) {
            double shift_x = ptsx[i] - ref_x;
            double shift_y = ptsy[i] - ref_y;

            ptsx[i] = shift_x * cos(0 - ref_yaw) - shift_y * sin(0 - ref_yaw);
            ptsy[i] = shift_x * sin(0 - ref_yaw) + shift_y * cos(0 - ref_yaw);
          }

          // create a spline
          tk::spline s;
          s.set_points(ptsx, ptsy);

          vector<double> next_x_vals;
          vector<double> next_y_vals;

          for (int i=0; i < previous_path_x.size(); ++i) {
            next_x_vals.push_back(previous_path_x[i]);
            next_y_vals.push_back(previous_path_y[i]);
          }

          double target_x = 30.0;
          double target_y = s(target_x);
          double target_dist = sqrt(target_x * target_x + target_y * target_y);

          double x_add_on = 0;

          for (int i=0; i<50 - previous_path_x.size(); ++i) {
            double N = target_dist / (0.02 * ref_vel / 2.24);
            double x_point = x_add_on + target_x / N;
            double y_point = s(x_point);

            x_add_on = x_point;

            double x_ref = x_point;
            double y_ref = y_point;

            x_point = x_ref * cos(ref_yaw) - y_ref * sin(ref_yaw);
            y_point = x_ref * sin(ref_yaw) + y_ref * cos(ref_yaw);

            x_point += ref_x;
            y_point += ref_y;

            //cout << i << ": point " << x_point << ", " << y_point << endl;
            next_x_vals.push_back(x_point);
            next_y_vals.push_back(y_point);
          }

          // calculate how to break up spline points so we travel at our desired reference velocity

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
