#include <vector>
#include <algorithm>
#include <map>
#include <cmath>
#include "grafoPMC.h"
#include "./LibGrafos/alg_grafoPMC.h"
#include "./LibGrafos/matriz.h"
#include "./LibGrafos/particion.h"


#pragma region Ejercicio1
template <typename tValue>
struct coordinates{
    size_t x;
    size_t y;
};

template <typename tValue>
struct timbuktuSol{
    matriz<tValue> cost;
    Particion partition;
};

template <typename tValue>
timbuktuSol<tValue> Timbuktu(const GrafoP<tValue> &adjacency, const std::vector<coordinates<tValue>> &cities){
    timbuktuSol<tValue> timbuktu;
    matriz<tValue> costTimbuktu(adjacency.numVert());
    Particion partition(adjacency.numVert());

    for(size_t i = 0; i < adjacency.numVert(); i++){
        for(size_t j = 0; j < adjacency.numVert(); j++){
            if(adjacency[i][j]){
                unsigned auxI = partition.encontrar(i);
                unsigned auxJ = partition.encontrar(j);

                if(auxI != auxJ) partition.unir(auxI,auxJ);
            }
        }
    }

    for(size_t i = 0; i < adjacency.numVert(); i++){
        for(size_t j = 0; j < adjacency.numVert(); j++){
            if(adjacency[i][j]) costTimbuktu[i][j] = sqrt(pow(cities[i].x-cities[j].x,2)+pow(cities[i].y-cities[j].y,2));
        }
    }

    matriz<typename GrafoP<tValue>::vertice> aux;
    timbuktu.cost = Floyd(costTimbuktu,aux);
    timbuktu.partition = partition;

    return timbuktu;
}
#pragma endregion


#pragma region Ejercicio2
template <typename tValue>
struct linePlane{ 
    matriz<tValue> cost;
    std::vector<coordinates2<tValue>> cities;
};

struct coordinates2{ size_t x,y; };

template <typename tValue>
linePlane<tValue> linePlaneTimbuktu(const std::vector<coordinates2> &cities, const GrafoP<tValue> &adjacency){
    linePlane<tValue> line;
    Particion partition(adjacency.numVert());

    for(size_t i = 0; i < adjacency.numVert(); i++){
        for(size_t j = 0; j < adjacency.numVert(); j++){
            if(adjacency[i][j]){
                unsigned auxI = partition.encontrar(i);
                unsigned auxJ = partition.encontrar(j);
                
                if(auxI != auxJ) partition.unir(auxI,auxJ);
            }
        }
    }

    matriz<tValue> cost(adjacency.numVert());

    for(size_t i = 0; i < adjacency.numVert(); i++){
        for(size_t j = 0; j < adjacency.numVert(); j++){
            if(adjacency[i][j]) cost[i][j] = sqrt(pow(cities[i].x-cities[j].x,2)+pow(cities[i].y-cities[j].y,2));

            if(partition.encontrar(i) != partition.encontrar(j)){
                cost[i][j] = sqrt(pow(cities[i].x-cities[j].x,2)+pow(cities[i].y-cities[j].y,2));
                line.cities.push_back(coordinates2(i,j));
            }
        }
    }

    matriz<typename GrafoP<tValue>::vertice> vecAux;
    line.cost = Floyd(cost,vecAux);

    return line;
}
#pragma endregion


#pragma region Ejercicio3
template <typename tValue>
GrafoP<tValue> maximumSpanningTree(const GrafoP<tValue> &grafo){
    GrafoP<tValue> grafoKruskall(grafo.numVert()),aux(grafo.numVert());

    for(size_t i = 0; i < grafo.numVert(); i++){
        for(size_t j = 0; j < grafo.numVert(); j++){
            aux[i][j] = -grafo[i][j];
        }
    }

    grafoKruskall = Kruskall(aux);

    for(size_t i = 0; i < grafo.numVert(); i++){
        for(size_t j = 0; j < grafo.numVert(); j++){
            grafoKruskall[i][j] = -grafoKruskall[i][j];
        }
    }

    return grafoKruskall;
}
#pragma endregion


#pragma region Ejercicio4
template <typename tValue>
GrafoP<tValue> valleJerte(const GrafoP<tValue> &grafo){
    GrafoP<tValue> grafoKruskall(grafo.numVert()),aux(grafo.numVert());

    for(size_t i = 0; i < grafo.numVert(); i++){
        for(size_t j = 0; j < grafo.numVert(); j++){
            aux[i][j] = -grafo[i][j];
        }
    }

    grafoKruskall = Kruskall(aux);

    for(size_t i = 0; i < grafo.numVert(); i++){
        for(size_t j = 0; j < grafo.numVert(); j++){
            grafoKruskall[i][j] = -grafoKruskall[i][j];
        }
    }

    return grafoKruskall;
}
#pragma endregion

#pragma region Ejercicio5
template <typename tValue>
tValue RETEUNI3(const GrafoP<tValue> &grafo){
    tValue minimumLenght;
    GrafoP<tValue> cost = Kruskall(grafo);

    for(size_t i = 0; i < cost.numVert(); i++){
        for(size_t j = 0; j < cost.numVert(); j++){
            if(cost[i][j] != cost.INFINITO) minimumLenght += cost[i][j];
        }
    }

    return minimumLenght;
}
#pragma endregion