#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <json/json.h>

const double EARTH_RADIUS_KM = 6371.0;  // Earth's radius in km
const double GRID_SPACING_KM = 5.0;     // 5km grid spacing

// Structure to hold vector data
struct VectorData {
    double direction;  // in degrees
    double magnitude;  // in m/s
};

// Structure to hold a grid cell
struct GridCell {
    double latitude;
    double longitude;
    VectorData wind;
    VectorData current;
    double wind_scale;  // Each square has its own wind scale factor
    double net_velocity_x;
    double net_velocity_y;
};

// Convert degrees to radians
double degToRad(double degrees) {
    return degrees * M_PI / 180.0;
}

// Convert radians to degrees
double radToDeg(double radians) {
    return radians * 180.0 / M_PI;
}

// Compute destination coordinates given start lat/lon, bearing, and distance
void computeGridCoordinates(double startLat, double startLon, int rows, int cols, std::vector<std::vector<GridCell>> &grid) {
    for (int i = 0; i < rows; ++i) {
        double newLat = startLat + (i * GRID_SPACING_KM / EARTH_RADIUS_KM) * (180.0 / M_PI);
        for (int j = 0; j < cols; ++j) {
            double newLon = startLon + (j * GRID_SPACING_KM / EARTH_RADIUS_KM) * (180.0 / M_PI) / cos(degToRad(newLat));
            grid[i][j] = {newLat, newLon, {0, 0}, {0, 0}, 1.0, 0, 0};  // Default wind scale = 1.0
        }
    }
}

// Convert direction & magnitude to vector components
void computeVectorComponents(double direction, double magnitude, double &vx, double &vy) {
    double rad = degToRad(direction);
    vx = magnitude * cos(rad);
    vy = magnitude * sin(rad);
}

// Process wind & current data with wind scaling per grid cell
void processGrid(std::vector<std::vector<GridCell>> &grid) {
    for (auto &row : grid) {
        for (auto &cell : row) {
            double wind_vx, wind_vy, current_vx, current_vy;

            computeVectorComponents(cell.wind.direction, cell.wind.magnitude * cell.wind_scale, wind_vx, wind_vy);
            computeVectorComponents(cell.current.direction, cell.current.magnitude, current_vx, current_vy);

            cell.net_velocity_x = wind_vx + current_vx;
            cell.net_velocity_y = wind_vy + current_vy;
        }
    }
}

// Save grid data to a JSON file
void saveGridToJson(const std::vector<std::vector<GridCell>> &grid, const std::string &filename) {
    Json::Value root;
    for (const auto &row : grid) {
        for (const auto &cell : row) {
            Json::Value gridPoint;
            gridPoint["latitude"] = cell.latitude;
            gridPoint["longitude"] = cell.longitude;
            gridPoint["wind_direction"] = cell.wind.direction;
            gridPoint["wind_magnitude"] = cell.wind.magnitude;
            gridPoint["wind_scale"] = cell.wind_scale;
            gridPoint["current_direction"] = cell.current.direction;
            gridPoint["current_magnitude"] = cell.current.magnitude;
            gridPoint["net_velocity_x"] = cell.net_velocity_x;
            gridPoint["net_velocity_y"] = cell.net_velocity_y;
            root.append(gridPoint);
        }
    }
    std::ofstream file(filename);
    file << root;
    file.close();
}

// Main function
int main() {
    double startLat, startLon, endLat, endLon;
    
    std::cout << "Enter start latitude and longitude: ";
    std::cin >> startLat >> startLon;
    std::cout << "Enter end latitude and longitude: ";
    std::cin >> endLat >> endLon;

    // Calculate number of rows and columns
    int rows = std::round((endLat - startLat) * EARTH_RADIUS_KM / GRID_SPACING_KM);
    int cols = std::round((endLon - startLon) * EARTH_RADIUS_KM / GRID_SPACING_KM);

    std::vector<std::vector<GridCell>> grid(rows, std::vector<GridCell>(cols));

    // Compute grid coordinates
    computeGridCoordinates(startLat, startLon, rows, cols, grid);

    // Simulated data input (Replace with actual API data)
    for (auto &row : grid) {
        for (auto &cell : row) {
            cell.wind = {rand() % 360, (rand() % 10) + 1};     // Random wind direction and speed
            cell.current = {rand() % 360, (rand() % 5) + 1};   // Random current direction and speed
            cell.wind_scale = 0.5 + static_cast<double>(rand()) / RAND_MAX;  // Random scaling factor between 0.5 - 1.5
        }
    }

    // Process grid data
    processGrid(grid);

    // Save to JSON file
    saveGridToJson(grid, "grid_output.json");

    std::cout << "Grid data saved to 'grid_output.json'.\n";
    return 0;
}
