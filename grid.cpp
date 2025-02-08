#include <iostream>
#include <vector>
#include <cmath>
#include <math.h>
#include <fstream>
//#include <jsoncpp/json/json.h>
#include "currentData.h"
#include "windData.h"

#define MAXV 1000

using namespace std;

class EdgeNode{
    public:
        int key;
        double weight;
        EdgeNode *next;
        EdgeNode(int, double);
        EdgeNode(int key, double weight);
};

EdgeNode::EdgeNode(int key, double weight){
    this->key = key;
    this->weight = weight;
    this->next = NULL;
}

class Graph{
    bool directed;
    public:
        EdgeNode *edges[MAXV + 1];
        Graph(bool);
        ~Graph();
        void insert_edge(int, int, double, bool);
        void print();
};

Graph::Graph(bool directed){
    this->directed = directed;
    for(int i = 1; i < (MAXV + 1); i ++){
        this->edges[i] = NULL;
    }
}

Graph::~Graph(){
    //to do
}

void Graph::insert_edge(int x, int y, double weight, bool directed){
    if(x > 0 && x < (MAXV + 1) && y > 0 && y < (MAXV + 1)){
        EdgeNode *edge = new EdgeNode(y, weight);
        edge->next = this->edges[x];
        this->edges[x] = edge;
        if(!directed){
            insert_edge(y, x, weight, true);
        }
    }
}

void Graph::print(){
    for(int v = 1; v < (MAXV + 1); v ++){
        if(this->edges[v] != NULL){
            cout << "Vertex " << v << " has neighbors: " << endl;
            EdgeNode *curr = this->edges[v];
            while(curr != NULL){
                cout << curr->key << endl;
                curr = curr->next;
            }
        }
    }
}

void init_vars(bool discovered[], int distance[], int parent[]){
    for(int i = 1; i < (MAXV + 1); i ++){
        discovered[i] = false;
        distance[i] = std::numeric_limits<int>::max();
        parent[i] = -1;
    }
}

void dijkstra_shortest_path(Graph *g, int parent[], int distance[], int start){

    bool discovered[MAXV + 1];
    EdgeNode *curr;
    int v_curr;
    int v_neighbor;
    int weight;
    int smallest_dist;

    init_vars(discovered, distance, parent);

    distance[start] = 0;
    v_curr = start;

    while(discovered[v_curr] == false){

        discovered[v_curr] = true;
        curr = g->edges[v_curr];

        while(curr != NULL){

            v_neighbor = curr->key;
            weight = curr->weight;

            if((distance[v_curr] + weight) < distance[v_neighbor]){
                distance[v_neighbor] = distance[v_curr] + weight;
                parent[v_neighbor] = v_curr;
            }
            curr = curr->next;
        }

        //set the next current vertex to the vertex with the smallest distance
        smallest_dist = std::numeric_limits<int>::max();
        for(int i = 1; i < (MAXV + 1); i ++){
            if(!discovered[i] && (distance[i] < smallest_dist)){
                v_curr = i;
                smallest_dist = distance[i];
            }
        }
    }
}

void print_shortest_path(int v, int parent[]){
    if(v > 0 && v < (MAXV + 1) && parent[v] != -1){
        cout << parent[v] << " ";
        print_shortest_path(parent[v], parent);
    }
}

void print_distances(int start, int distance[]){
    for(int i = 1; i < (MAXV + 1); i ++){
        if(distance[i] != std::numeric_limits<int>::max()){
            cout << "Shortest distance from " << start << " to " << i << " is: " << distance[i] << endl;
        }
    }
}

const double EARTH_RADIUS_KM = 6371.0;  // Earth's radius in km
const double GRID_SPACING_KM = 5.0;     // 5km grid spacing

// Structure to hold vector data
struct VectorData {
    double direction;  // in degrees
    double magnitude;  // in m/s

    VectorData(double dir = 0, double mag = 0) : direction(dir), magnitude(mag) {}
};

// Structure to hold a grid cell
struct GridCell {
    double latitude;
    double longitude;
    VectorData wind;
    VectorData current;
    double up_conn;
    double down_conn;
    double left_conn;
    double right_conn;
    double wind_scale;  // Each square has its own wind scale factor
    double dot;
    double heading;
    
    GridCell(double lat = 0, double lon = 0, VectorData windData = {0, 0}, VectorData currentData = {0, 0},
             double scale = 1.0, double dotProd = 0, double head = 0)
        : latitude(lat), longitude(lon), wind(windData), current(currentData),
          wind_scale(scale), dot(dotProd), heading(head) {}
    
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
void computeGridCoordinates(double startLat, double startLon, int rows, int cols, double heading, std::vector<std::vector<GridCell>> &grid) {
    for (int i = 0; i < rows; ++i) {
        double newLat = startLat + (i * GRID_SPACING_KM / EARTH_RADIUS_KM) * (180.0 / M_PI);
        for (int j = 0; j < cols; ++j) {
            double newLon = startLon + (j * GRID_SPACING_KM / EARTH_RADIUS_KM) * (180.0 / M_PI) / cos(degToRad(newLat));
            grid[i][j] = GridCell{newLat, newLon, VectorData{0, 0}, VectorData{0, 0}, 1.0, 0, heading};  // Default wind scale = 1.0
        }
    }
}

// Convert direction & magnitude to vector components
void computeVectorComponents(double direction, double magnitude, double *vx, double *vy) {
    double rad = degToRad(direction);
    *vx = magnitude * cos(rad);
    *vy = magnitude * sin(rad);
}

// Process wind & current data with wind scaling per grid cell
void processGrid(std::vector<std::vector<GridCell>> &grid) {
    for (size_t i = 0; i<grid.size(); i++) {
        std::vector<GridCell> *row = &grid.at(i);
        for (size_t j = 0; j<row->size(); j++) {
            GridCell *cell = &row->at(j);
            double wind_vx, wind_vy, current_vx, current_vy;

            if (abs(cell->heading-cell->wind.direction)<90){
                cell->wind_scale = 0;
            }

            computeVectorComponents(cell->wind.direction, cell->wind.magnitude * cell->wind_scale, &wind_vx, &wind_vy);
            computeVectorComponents(cell->current.direction, cell->current.magnitude, &current_vx, &current_vy);

            //cell.dot = wind_vx*current_vx + wind_vy*current_vy;
            
            int col_up = j-1;
            int col_down = j+1;
            int row_l = i-1;
            int row_r = i+1;
            if(col_up>=0) {
                double wind_dot = wind_vy;
                double cur_dot = current_vy;
                double res = cur_dot;
                if(wind_dot>=0) res+=wind_dot;
                cell->up_conn = abs(res);
            } else cell->up_conn = -1.0;
            if(col_down < grid.at(0).capacity()) {
                double wind_dot = wind_vy;
                double cur_dot = current_vy;
                double res = cur_dot;
                if(wind_dot<=0) res+=wind_dot;
                cell->down_conn = abs(res);
            } else cell->down_conn = -1.0;
            if(row_l>=0) {
                double wind_dot = wind_vx;
                double cur_dot = current_vx;
                double res = cur_dot;
                if(wind_dot<=0) res+=wind_dot;
                cell->left_conn = abs(res);
            } else cell->left_conn = -1.0;
            if(row_r<grid.capacity()) {
                double wind_dot = wind_vx;
                double cur_dot = current_vx;
                double res = cur_dot;
                if(wind_dot>=0) res+=wind_dot;
                cell->right_conn = abs(res);
            } else cell->right_conn = -1.0;

        }
    }
}


// Main function
int main() {
    double startLat, startLon, endLat, endLon, heading;
    
    std::cout << "Enter start latitude and longitude: ";
    std::cin >> startLat >> startLon;
    std::cout << "Enter end latitude and longitude: ";
    std::cin >> endLat >> endLon;
    std::cout << "Enter the heading: ";
    std::cin >> heading;


    // Calculate number of rows and columns
    int rows = std::round((endLat - startLat) * EARTH_RADIUS_KM / GRID_SPACING_KM);
    int cols = std::round((endLon - startLon) * EARTH_RADIUS_KM / GRID_SPACING_KM);

    std::vector<std::vector<GridCell>> grid(rows, std::vector<GridCell>(cols));

    // Compute grid coordinates
    std::cout << "Starting to compute grid inputs" << std::endl;
    computeGridCoordinates(startLat, startLon, rows, cols, heading, grid);
    std::cout << "Finished computing grid inputs" << std::endl;


    // Simulated data input (Replace with actual API data)
    for (auto &row : grid) {
        for (auto &cell : row) {

            //cell.wind = {std::get<0>(returnWindData(cell.latitude, cell.longitude)),std::get<1>(returnWindData(cell.latitude, cell.longitude)) };     //INSERT RYANS API THING
            //cell.current = {std::get<0>(returnData(cell.latitude, cell.longitude)),std::get<1>(returnData(cell.latitude, cell.longitude)) };   // Random current direction and speed
            cell.wind = {0,0};
            cell.current = {0,0};
        }
    }
    std::cout << "Finished sim data" << std::endl;
    // Process grid data
    processGrid(grid);
    std::cout << "Finished processing grid" << std::endl;

    Graph *g = new Graph(false);

    std::cout << "Entering graph stuff" << std::endl;

    // turn grid into graph
    for(size_t i=0;i<grid.size();i++) {
        std::vector<GridCell> *row = &grid.at(i);
        std::cout << "Row idk" << std::endl;
        for(size_t j=0;j<row->size();j++) {
            printf("%i,%i\n",i,j);
            GridCell *cell = &row->at(j);
            size_t ind = cols*j+i;
            if(cell->up_conn!=-1.0) g->insert_edge(ind, cols*(j-1)+i, cell->up_conn, false);
            if(cell->down_conn!=-1.0) g->insert_edge(ind, cols*(j+1)+i, cell->down_conn, false);
            if(cell->left_conn!=-1.0) g->insert_edge(ind, cols*j+i-1, cell->left_conn, false);
            if(cell->right_conn!=-1.0) g->insert_edge(ind, cols*j+i+1, cell->right_conn, false);
        }
    }
    std::cout << "Finished processing graph" << std::endl;
    // Save to JSON file
    //saveGridToJson(grid, "grid_output.json");

    std::cout << "Grid data saved to 'grid_output.json'.\n";
    return 0;
}
