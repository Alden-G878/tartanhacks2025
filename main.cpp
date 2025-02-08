#include <Arduino.h>
#include <vector>

void setup() {
    std::vector<double> *heading = new std::vector<double>(2);
    // TODO: normalize heading
    std::vector<double> *wind_head = new std::vector<double>(2);
    *wind_head.at(0) = cos(ang);
    *wind_head.at(1) = sin(ang);
    double d = std::dot(heading, wind_head);
    std::vector<double> norm = wind_head*d;
    norm -= wind_head;

}
void loop() {

}