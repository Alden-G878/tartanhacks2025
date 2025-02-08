#include <iostream>
#include <curl/curl.h>
#include <json/json.h>  // JsonCpp library
#include <sstream>

using namespace std;

// API Key (Replace with your actual API key)
const string API_KEY = "193cace1120fd7161808797e8b8ccbbe";
const string BASE_URL = "http://api.openweathermap.org/data/2.5/weather";

// Callback function to write the response into a string
size_t WriteCallback(void* contents, size_t size, size_t nmemb, string* output) {
    size_t totalSize = size * nmemb;
    output->append((char*)contents, totalSize);
    return totalSize;
}

// Function to fetch weather data using latitude & longitude
string fetchWeatherData(double latitude, double longitude) {
    CURL* curl;
    CURLcode res;
    string readBuffer;

    curl = curl_easy_init();
    if (curl) {
        string url = BASE_URL + "?lat=" + to_string(latitude) + "&lon=" + to_string(longitude) +
                     "&appid=" + API_KEY + "&units=metric";  // Using metric for speed in m/s

        curl_easy_setopt(curl, CURLOPT_URL, url.c_str());
        curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, WriteCallback);
        curl_easy_setopt(curl, CURLOPT_WRITEDATA, &readBuffer);
        curl_easy_setopt(curl, CURLOPT_SSL_VERIFYPEER, 0L); // Disable SSL verification (use cautiously)

        res = curl_easy_perform(curl);
        if (res != CURLE_OK) {
            cerr << "curl_easy_perform() failed: " << curl_easy_strerror(res) << endl;
        }

        curl_easy_cleanup(curl);
    }

    return readBuffer;
}

// Function to parse and display wind speed & direction using JsonCpp
std::tuple<double, double> displayWindData(const string& jsonResponse) {
    Json::Value jsonData;
    Json::CharReaderBuilder reader;
    string errors;

    istringstream jsonStream(jsonResponse);
    if (!Json::parseFromStream(reader, jsonStream, &jsonData, &errors)) {
        cerr << "Error parsing JSON: " << errors << endl;
        return std::make_tuple(0, 0);
    }

    if (jsonData.isMember("wind")) {
        double windSpeed = jsonData["wind"]["speed"].asDouble();
        int windDirection = jsonData["wind"]["deg"].asInt();
        return std::make_tuple(windSpeed, windDirection);
    } else {
        cerr << "Wind data not available in API response!" << endl;
        return std::make_tuple(0, 0);
    }
    
}

std::tuple<double, double> returnWindData(double latitude, double longitude) {

    string response = fetchWeatherData(latitude, longitude);
    if (!response.empty()) {
        return displayWindData(response);
    } else {
        cout << "Failed to get weather data." << endl;
        return std::make_tuple(0, 0);
    }
}
