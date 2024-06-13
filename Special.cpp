#include <vector>
#include <algorithm>
#include <map>
#include <cmath>
#include "grafoPMC.h"
#include "./LibGrafos/alg_grafoPMC.h"
#include "./LibGrafos/matriz.h"
#include "./LibGrafos/particion.h"


/*A distributor of a drink distribution company has to visit all its customers every day. But, at the beginning of their workday, he does not know how much drinks each client has to serve, so he cannot plan an optimal route to visit them all.

Therefore, our delivery man decides to carry out the following strategy:
- The truck starts from the warehouse with the maximum load allowed for your client more
next.
- The distributor downloads the beverage boxes that the client asks for. If you don't have enough
Boxes in the truck, delivers all the ones. This client will finish being served in
some other time throughout the day, when the cast strategy leads to
Distributor to him.

After serving a client:
- If drinks are left in the truck, the delivery man consults its navigation system
GPS based to know the route that takes you to your closest client
pending to be served.
- If there are no drinks in the truck, return to the warehouse along the shortest road and
Again load the truck completely.

After loading the truck, the dealer consults its navigation system and goes for
The shortest way to visit the client pending being closer.
Implements a subprogram that calculates and returns the total distance traveled in a
day by our delivery man, from the following:
1. GRAFO represented by cost matrix with the distances of the direct roads
between customers and between them and the central.
2. Maximum truck capacity (number of beverage boxes).
3. We will assume that there is an int request () function that returns the number of boxes that
The client is left to serve the deliveryman.*/
#pragma region Ejercicio1
int request(size_t i,std::vector<int> &clientBoxes){ return clientBoxes[i]; }

template <typename tValue>
tValue delivery(const GrafoP<tValue> &distances, const int totalCapacity, const size_t central, std::vector<int> &clientBoxes){
    int actualCapacity = totalCapacity, totalBoxes = 0;
    tValue totalDistance,minPath = distances.INFINITO;
    matriz<typename GrafoP<tValue>::vertice> matrizAux;
    GrafoP<tValue> minCost = Floyd(distances,matrizAux);
    size_t aux,actualPosition = central;

    for(size_t i = 0; i < distances.numVert(); i++) totalBoxes += request(i);

    while(totalBoxes > 0){
        minPath = distances.INFINITO;
        for(size_t i = 0; i < minCost.numVert(); i++){
            if(i != actualPosition && i != central && minPath > minCost[actualPosition][i]){
                minPath = minCost[actualPosition][i];
                aux = i;
            }
        }

        if(actualCapacity == 0){
            actualCapacity = totalCapacity;
            actualPosition = central;
            totalDistance += minCost[aux][actualPosition];
        }else if(request(aux,clientBoxes) > 0){
            if(request(aux,clientBoxes) >= totalCapacity){
                clientBoxes[aux] -= totalCapacity;
                totalCapacity = actualCapacity;
                totalDistance += minCost[actualPosition][aux];
                actualPosition = central;
                totalDistance += minCost[aux][actualPosition];
            }else{
                totalCapacity -= request(aux,clientBoxes);
                clientBoxes[aux] = 0;
                totalDistance += minCost[actualPosition][aux];
                actualPosition = aux;
            }
        }
    }

    return totalDistance;
}
#pragma endregion


/*We will model narnia as a nxm square matrix. The
Faun movements will be modernized with the movements of a chess horse. In other words, each movement of the faun must be a chess horse movement.
narnia is reached through its entrance, box (0.0), and one leaves through a single exit, in the box (N-1, M-1).
It would be a fairly easy problem, but narnia is a country full of dangers, in particular if you are a faun.
To begin with, the locals have put within narnia a series of traps in certain boxes, so that if you go through them you die.
But not happy with that, the inhabitants of narnia have Hired knights, which are also placed in strategic boxes. 
In this case, the faun does not death if it falls into one of them, but its dies in case of falling into any of the boxes that surround the gentleman,
 between 3 and 8, depending on their position, since its sword has length 1.
The problem asks us two things, the first to know if the faun can safely make the path between the entrance and exit of narnia, ç
and, if so, what would be the minimum number of jumps necessary to achieve it.*/
#pragma region Ejercicio2
struct cells{
    size_t x,y;
};

template <typename tValue>
tValue Faun(const std::vector<cells> &traps, const std::vector<cells> &knights, const size_t N, const size_t M){
    GrafoP<tValue> narnia(N*M);

    for(size_t i = 0; i < narnia.numVert(); i++){
        for(size_t j = 0; j < narnia.numVert(); j++){
            narnia[i][j] = 0;
        }
    }

    for(size_t i = 0; i < traps.size(); i++) narnia[traps[i].x][traps[i].y] = narnia.INFINITO;
    for(size_t i = 0; i < knights.size(); i++){
        for(size_t i = -1; i < 1; i++){
            for(size_t j = -1; j < 1; j++){
                if(0 <= knights[i].x-i < narnia.numVert() && 0 <= knights[i].y-j < narnia.numVert()){
                    narnia[knights[i].x-i][knights[i].y-j] = narnia.INFINITO;
                }
            }
        }    
    }

    for(size_t i = 0; i < narnia.numVert(); i++){
        for(size_t j = 0; j < narnia.numVert(); j++){
            if(0 <= narnia[i+2][j+1] < narnia.numVert() && narnia[i][j] != narnia.INFINITO) narnia[i+2][j+1] = 1;
            if(0 <= narnia[i-2][j+1] < narnia.numVert() && narnia[i][j] != narnia.INFINITO) narnia[i-2][j+1] = 1;
            if(0 <= narnia[i+2][j-1] < narnia.numVert() && narnia[i][j] != narnia.INFINITO) narnia[i+2][j-1] = 1;
            if(0 <= narnia[i-2][j-1] < narnia.numVert() && narnia[i][j] != narnia.INFINITO) narnia[i-2][j-1] = 1;
            if(0 <= narnia[i+1][j+2] < narnia.numVert() && narnia[i][j] != narnia.INFINITO) narnia[i+1][j+2] = 1;
            if(0 <= narnia[i-1][j+2] < narnia.numVert() && narnia[i][j] != narnia.INFINITO) narnia[i-1][j+2] = 1;
            if(0 <= narnia[i+1][j-2] < narnia.numVert() && narnia[i][j] != narnia.INFINITO) narnia[i+1][j-2] = 1;
            if(0 <= narnia[i-1][j-2] < narnia.numVert() && narnia[i][j] != narnia.INFINITO) narnia[i-1][j-2] = 1;
        }
    }

    typename GrafoP<tValue>::vertice origin = 0;
    std::vector<typename GrafoP<tValue>::vertice> vect;
    std::vector<tValue> minCost = Dijkstra(narnia,origin,vect);

    return minCost[minCost.size()-1];
}
#pragma endregion


/*There is a Labyrinth of NxNxN boxes from which the entrance and exit boxes are known.
If you are in a box you can only move (in the best case) in the next six directions (up, below, right, left, inside, outside).
On the other hand, among some of the boxes there is a stone that prevents the engine towards it. Implement of the oni subprogram:
-n (Labyrinth dimension),
-The list of boxes that have a stone,
-The entrance box, and
-The exit box,
Calculate the cost of the shortest path to go from the entrance to the exit and its length.
*Note: Define data molerons, prototype of the operations of the TADS and graph algorithms.*/
#pragma region Ejercicio3
struct cells3{
    int x,y,z;
};

size_t cells3ToNodo(cells3 c, size_t N){ return c.x + c.y*N + c.z*(N*N); }

cells3 NodoToCells3(size_t v, size_t N){
    cells3 aux;

    aux.x = v % N;
    aux.y = (v/N) % N;
    aux.z = v / (N*N);

    return aux;
}

bool checkCoords(size_t i, size_t j, size_t N){
    cells3 auxI = NodoToCells3(i,N);
    cells3 auxJ = NodoToCells3(j,N);
    
    if(auxI.z == auxJ.z) return std::abs((auxI.x-auxJ.x) + (auxI.y-auxJ.y)) == 1;
    else return  (auxI.x == auxJ.x) && (auxI.y == auxJ.y) && std::abs(auxI.z-auxJ.z) == 1;
}

template <typename tValue>
tValue shortPathLabyrinthN3(const size_t N, const std::vector<cells3> &boxes, const size_t entrance, const size_t exit){
    GrafoP<tValue> labyrinth(pow(N,3));

    for(size_t i = 0; i < labyrinth.numVert(); i++){
        for(size_t j = 0; j < labyrinth.numVert(); j++){
            labyrinth[i][j] = 1;
        }
    }

    for(size_t i = 0; i < labyrinth.numVert(); i++){
        for(size_t j = 0; j < boxes.size(); j++){
        labyrinth[i][cells3ToNodo(boxes[j],N)] = labyrinth.INFINITO;
        labyrinth[cells3ToNodo(boxes[j],N)][i] = labyrinth.INFINITO;
        }
    }

    std::vector<typename GrafoP<tValue>::vertice> vect;
    std::vector<tValue> minCost = Dijkstra(labyrinth,entrance,vect);

    return minCost[exit];
}
#pragma endregion