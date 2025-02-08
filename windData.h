#ifndef WIND_DATA_H
#define WIND_DATA_H

#include <string>

std::tuple<double, double> returnWindData(double latitude, double longitude);
std::string fetchWeatherData(double latitude, double longitude);
void displayWindData(const std::string &jsonResponse);

#endif  // WIND_DATA_H
