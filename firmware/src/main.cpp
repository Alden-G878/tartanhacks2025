#include <Arduino.h>
#include <vector>

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

void setup() {
  Serial.begin(115200);
}
void loop() {
  Serial.print("Angle of wind: ");
  double ang = Serial.read();
  Serial.println();
  Serial.printf("The angle of the wind is %f\n",ang);
  Serial.print("Heading angle: ");
  double h_ang = Serial.read();
  Serial.println();
  Serial.printf("The angle of the heading is %f\n",h_ang);
  std::vector<double> heading(2);
  // TODO: normalize heading
  std::vector<double> wind_head(2);
  wind_head.at(0) = cos(ang);
  wind_head.at(1) = sin(ang);
  heading.at(0) = cos(h_ang);
  heading.at(1) = sin(h_ang);
  double d = dot(heading, wind_head);
  std::vector<double> norm = prop(wind_head, d);
  norm = sub(norm, wind_head);
  Serial.printf("(x,y)=(%f,%f)\n",norm.at(0), norm.at(1));
}