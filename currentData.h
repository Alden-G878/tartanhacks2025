#ifndef CURRENT_DATA_H
#define CURRENT_DATA_H

#include <tuple>

std::tuple<double, double> returnData(double latitude, double longitude);
void fetchOceanCurrentData(double latitude, double longitude, double &direction, double &magnitude);

#endif  // CURRENT_DATA_H
