


#include <cmath>

#include <iostream>
#include <string>
#include <curl/curl.h>
#include <json/json.h>
#include <sstream>

size_t WriteCallback(void* contents, size_t size, size_t nmemb, std::string* output) {
    size_t total_size = size * nmemb;
    output->append((char*)contents, total_size);
    return total_size;
}

void fetchOceanCurrentData(double latitude, double longitude, double &current_speed, double &current_direction) {
    CURL* curl;
    CURLcode res;
    std::string response_string;

    std::string url = "https://marine-api.open-meteo.com/v1/marine?latitude=" 
                      + std::to_string(latitude) + "&longitude=" + std::to_string(longitude)
                      + "&hourly=ocean_current_velocity,ocean_current_direction";

    curl = curl_easy_init();
    if (curl) {
        curl_easy_setopt(curl, CURLOPT_URL, url.c_str());
        curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, WriteCallback);
        curl_easy_setopt(curl, CURLOPT_WRITEDATA, &response_string);

        res = curl_easy_perform(curl);
        curl_easy_cleanup(curl);


        if (res == CURLE_OK) {
            Json::CharReaderBuilder reader;
            Json::Value jsonData;
            std::string errors;

            std::istringstream stream(response_string);
            if (Json::parseFromStream(reader, stream, &jsonData, &errors)) {
                // Extract ocean current velocity & direction
                if (!jsonData["hourly"]["ocean_current_velocity"].empty() && !jsonData["hourly"]["ocean_current_direction"].empty()) {
                    current_speed = jsonData["hourly"]["ocean_current_velocity"][0].asDouble();
                    current_direction = jsonData["hourly"]["ocean_current_direction"][0].asDouble();
                } else {
                    std::cerr << "No ocean current data available for this location." << std::endl;
                }
            } else {
                std::cerr << "Failed to parse JSON: " << errors << std::endl;
            }
        } else {
            std::cerr << "Error fetching data: " << curl_easy_strerror(res) << std::endl;
        }
    }
}

std::tuple<double, double> returnData(double latitude, double longitude) {
    double current_speed = 0, current_direction = 0;

    fetchOceanCurrentData(latitude, longitude, current_speed, current_direction);

    std::cout << "Ocean Current Speed: " << current_speed << " m/s" << std::endl;
    std::cout << "Ocean Current Direction: " << current_direction << "°" << std::endl;

    return std::make_tuple(current_speed, current_direction);
}
