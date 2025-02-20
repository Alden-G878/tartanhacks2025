#include <Arduino.h>
#include <vector>
#include <Servo.h>

#define SAIL_DIAMETER 2*0.057
#define ROPE_ANGLE M_PI/4
#define DEFAULT_DIST 0.025
#define SERVO_RAD 0.015
#define SERVO0 0
#define SERVO1 1
#define SERVO2 2
#define SERVO3 3
#define SERVO4 4
#define SERVO5 5
#define SERVO6 6
#define SERVO7 7
#define SERVO8 8

Servo s0;
Servo s1;
Servo s2;
Servo s3;
Servo s4;
Servo s5;
Servo s6;
Servo s7;
Servo s8;

double dot(std::vector<double> v1, std::vector<double> v2) {
  double res = 0.0;
  size_t len = max(v1.capacity(), v2.capacity());
  size_t v1_cap = v1.capacity();
  size_t v2_cap = v2.capacity();
  for(size_t i=0;i<len;i++) {
    double v1_a;
    double v2_a;
    if(i<v1_cap) v1_a = v1.at(i);
    else v1_a = 0;
    if(i<v2_cap) v2_a = v2.at(i);
    else v2_a = 0;
    res+=(v1_a*v2_a);
  }
  return res;
}

std::vector<double> prop(std::vector<double> v, double d) {
  size_t s = v.capacity();
  for(size_t i=0;i<s;i++) {
    v.at(i) = v.at(i) * d;
  }
  return v;
}

std::vector<double> sub(std::vector<double> v1, std::vector<double> v2) {
  size_t len = max(v1.capacity(), v2.capacity());
  std::vector<double> res(len);
  size_t v1_cap = v1.capacity();
  size_t v2_cap = v2.capacity();
  for(size_t i=0;i<len;i++) {
    double v1_a;
    double v2_a;
    if(i<v1_cap) v1_a = v1.at(i);
    else v1_a = 0;
    if(i<v2_cap) v2_a = v2.at(i);
    else v2_a = 0;
    res.at(i) = (v1_a-v2_a);
  }
  return res;
}

double size(std::vector<double> v) {
  return sqrt(dot(v,v));
}

std::vector<double> normalize(std::vector<double> v) {
  return prop(v, 1/size(v));
}

void setup() {
  Serial.begin(115200);
  Serial.setTimeout(5000);
  s0.attach(SERVO0); // front middle R
  s1.attach(SERVO1); // back middle R
  s2.attach(SERVO2);
  s3.attach(SERVO3);
  s4.attach(SERVO4);
  s5.attach(SERVO5);
  s6.attach(SERVO6);
  s7.attach(SERVO7);
  s8.attach(SERVO8);
}
void loop() {
  Serial.print("Angle of wind: ");
  String s = Serial.readString();
  double ang = atof(s.c_str());
  Serial.println();
  Serial.printf("The angle of the wind is %f\n",ang);
  Serial.print("Heading angle: ");
  s = Serial.readString();
  double h_ang = atof(s.c_str());
  Serial.println();
  Serial.printf("The angle of the heading is %f\n",h_ang);
  std::vector<double> heading(2);
  std::vector<double> wind_head(2);
  wind_head.at(0) = cos(ang*(M_PI/180));
  wind_head.at(1) = sin(ang*(M_PI/180));
  heading.at(0) = cos(h_ang*(M_PI/180));
  heading.at(1) = sin(h_ang*(M_PI/180));

  double rot_ang = -h_ang*(M_PI/180);

  double wh_t_0 = wind_head.at(0);
  double wh_t_1 = wind_head.at(1);
  double h_t_0 = heading.at(0);
  double h_t_1 = heading.at(1);
  wind_head.at(0) = cos(rot_ang)*wh_t_0 - sin(rot_ang)*wh_t_1;
  wind_head.at(1) = sin(rot_ang)*wh_t_0 + cos(rot_ang)*wh_t_1;
  heading.at(0) = cos(rot_ang)*h_t_0 - sin(rot_ang)*h_t_1;
  heading.at(1) = sin(rot_ang)*h_t_0 + cos(rot_ang)*h_t_1;

  std::vector<double> h_norm(2);
  std::vector<double> w_norm(2);
  h_norm.at(0) = heading.at(1);
  h_norm.at(1) = -heading.at(0);
  w_norm.at(0) = wind_head.at(1);
  w_norm.at(1) = -wind_head.at(0);

  std::vector<double> app = prop(w_norm, SAIL_DIAMETER/2);

  double dist = DEFAULT_DIST - app.at(1); // offset of the boom
  
  double height = SAIL_DIAMETER / 2;
  double r = sqrt(pow(height, 2) + pow(dist, 2)); 
  double def_r = sqrt(pow(height, 2) + pow(DEFAULT_DIST, 2));
  double dr = r - def_r;
  double ndr = r + dr;

  // Step 4: Compute servo rotation angle (in degrees)
  double servo_angle = dr / SERVO_RAD * (180/M_PI);
  double n_servo_angle = ndr / SERVO_RAD * (190/M_PI);
  while(1) {
    delay(500);
    Serial.printf("Write servo: %f\n",servo_angle);
    s2.write(servo_angle);
    //s1.write(n_servo_angle);
    delay(500);
    Serial.printf("Write servo: %f\n",0);
    s2.write(0);
    //s1.write(0);
  }
}