#include <vector>
#include <algorithm>
#include <map>
#include "grafoPMC.h"
#include "./LibGrafos/alg_grafoPMC.h"
#include "./LibGrafos/matriz.h"
#include "./LibGrafos/particion.h"


#pragma region Ejercicio1
template <typename tValue>
tValue amazingCharacter(const GrafoP<tValue> &grafo, typename GrafoP<tValue>::vertice &origin, typename GrafoP<tValue>::vertice &destiny){
    tValue value = 0;
    GrafoP<tValue> grafoAux = grafo;
    matriz<typename GrafoP<tValue>::vertice> matrizAux;

    for(size_t i = 0; i < grafoAux.numVert(); i++){
        for(size_t j = 0; j < grafoAux.numVert(); j++){
            grafoAux[i][j] = -grafoAux[i][j]; // o -1 * grafoAux[i][j] if tValue isnt an arithmetic type
        }
    }

    matriz<tValue> matrizCostes = Floyd(grafoAux,matrizAux);

    for(size_t i = 0; i < grafoAux.numVert(); i++){
        for(size_t j = 0; j < grafoAux.numVert(); j++){
            grafoAux[i][j] = -grafoAux[i][j]; // o -1 * grafoAux[i][j] if tValue isnt an arithmetic type
            if(value < grafoAux[i][j]){
                value = grafoAux[i][j];
                origin = i;
                destiny = j;
            }
        }
    }

    return value;
}
#pragma endregion


#pragma region Ejercicio2
typedef struct coordinates{
    size_t i;
    size_t j;
};

template <typename tValue>
tValue mazePath(const tValue N, const std::vector<coordinates> wallsList, typename GrafoP<tValue>::vertice &origin, typename GrafoP<tValue>::vertice &destiny){
    GrafoP<tValue> grafo(N);
    std::vector<typename GrafoP<tValue>::vertice> vectorAux;

    for(size_t i = 0; i < grafo.numVert(); i++) for(size_t j = 0; j < grafo.numVert(); j++) grafo[i][j] = 1;
    for(size_t i = 0; i < wallsList.size(); i++) grafo[wallsList[i].i][wallsList[i].j] = grafo.INFINITO;

    std::vector<tValue> vectorPaths = Dijkstra(grafo,origin,vectorAux);

    return vectorPaths[destiny];
}
#pragma endregion


#pragma region Ejercicio3
template <typename tValue>
struct solution{
    tValue cost;
    std::vector<unsigned int> stored;
};

template <typename tValue>
solution<tValue> companyMoves(const GrafoP<tValue> &grafo, const typename GrafoP<tValue>::vertice origin,
                                std::vector<tValue> subsidy, std::vector<unsigned int> capacity, unsigned int stock){
    tValue minimunCost;
    size_t aux;
    unsigned int cityCapacity;
    solution<tValue> path;
    GrafoP<tValue> grafoAux = grafo;
    std::vector<bool> stored(grafo.numVert(),false);
    std::vector<typename GrafoP<tValue>::vertice> vec;
    std::vector<tValue> costs = Dijkstra(grafoAux,origin,vec);

    for(size_t i = 0; i < subsidy.size(); i++){
        if(origin != i) costs[i] -= ((subsidy[i]/100) * costs[i]);
    }

    while(stock > 0){
        minimunCost = grafo.INFINITO;
        
        for(size_t i = 0; i < costs.size(); i++){
            if(i != origin && minimunCost > costs[i] && capacity[i] > 0 && !stored[i]){
                minimunCost = costs[i];
                cityCapacity = capacity[i];
                aux = i;
            }
        }

        path.cost += minimunCost;
        path.stored[aux] = (cityCapacity <= stock) ? cityCapacity : stock;
        stored[aux] = (cityCapacity <= stock) ? true : false;
        stock -= cityCapacity;
    }

    return path;
}

#pragma endregion


#pragma region Ejercicio4
template <typename tValue>
std::map<tValue,bool> kmExtra(const GrafoP<tValue> &grafo, const typename GrafoP<tValue>::vertice capital, const std::vector<tValue> kmTruck){
    std::map<int,bool> trucks(grafo.numVert());
    GrafoP<tValue> grafoAux = grafo;
    std::vector<typename GrafoP<tValue>::vertice> vec;
    matriz<tValue> costs = Floyd(grafoAux,vec);

    for(size_t i = 0; i < grafo.numVert(); i++){
        if(i != capital && costs[i][capital] != grafo.INFINITO && costs[capital][i] != grafo.INFINITO && kmTruck[i] > 0){
            if((costs[i][capital] + costs[i][capital]) < kmTruck[i]) trucks.insert(std::make_pair(kmTruck[i],true));
            else trucks.insert(std::make_pair(kmTruck[i],false));
        }
    }

    return trucks;
}
#pragma endregion


#pragma region Ejercicio5
enum Allergic{ROAD,TRAIN,PLANE};

template <typename tValue>
std::vector<bool> traveler(const GrafoP<tValue> &road, const GrafoP<tValue> &train, const GrafoP<tValue> &plane, 
                            tValue cost, typename GrafoP<tValue>::vertice origin, const Allergic allergic){
    size_t size = road.numVert();
    GrafoP grafo(size);
    tValue minCost;
    size_t aux;

    for(size_t i = 0; i < size; i++){
        for(size_t j = 0; j < size; j++){
            if(allergic != ROAD) grafo[i][j] = road[i][j];
            if(allergic != TRAIN) grafo[i][j] = std::min(grafo[i][j],train[i][j]);
            if(allergic != PLANE) grafo[i][j] = std::min(grafo[i][j],plane[i][j]);
        }
    }

    std::vector<bool> visits(train.numVert(),false);
    std::vector<typename GrafoP<tValue>::vertice> vert;

    std::vector<tValue> costG = Dijkstra(grafo,origin,vert);

    while(cost > 0){
        minCost = train.INFINITO;
        for(size_t i = 0; i < costG.size(); i++){
            if(i != origin && minCost > costG[i] && !visits[i]){
                minCost = costG[i];
                aux = i;
            }
        }
        
        if(minCost < cost){
            visits[aux] = true;
            cost -= minCost;
        }
    }
    return visits;
}
#pragma endregion


#pragma region Ejercicio6
template <typename tValue>
matriz<tValue> minimumRate(const GrafoP<tValue> &train, const GrafoP<tValue> &bus, const typename GrafoP<tValue>::vertice nexo){
    matriz<tValue> minimunCost(train.numVert());
    matriz<typename GrafoP<tValue>::vertice> vecTrain;
    matriz<typename GrafoP<tValue>::vertice> vecBus;

    matriz<tValue> matrizTrain = Floyd(train,vecTrain);
    matriz<tValue> matrizBus = Floyd(bus,vecBus);

    for(size_t i = 0; i < matrizTrain.dimension(); i++){
        for(size_t j = 0; j < matrizTrain.dimension(); j++){
            minimunCost[i][j] = std::min(matrizTrain[i][j],matrizBus[i][j]);
            minimunCost[i][j] = std::min(minimunCost[i][j],matrizTrain[i][nexo] + matrizBus[nexo][j]);
        }
    }

    return minimunCost;
}
#pragma endregion


#pragma region Ejercicio7
template <typename tValue>
struct travel{
    size_t cost;
    typename GrafoP<tValue>::tCamino path;
};

template <typename tValue>
GrafoP<tValue> superGrafo(const GrafoP<tValue> &train, const GrafoP<tValue> &bus){
    GrafoP<tValue> superGrafo(train.numVert()+bus.numVert());

    for(size_t i = 0; i < train.numVert(); i++){
        for(size_t j = 0; j < train.numVert(); j++){
            superGrafo[i][j] = train[i][j];
        }
    }

    for(size_t i = 0; i < bus.numVert(); i++){
        for(size_t j = 0; j < bus.numVert(); j++){
            superGrafo[i+bus.numVert()][j+bus.numVert()] = bus[i][j];
        }
    }

    return superGrafo;
}

template <typename tValue>
travel<tValue> tripRoute(const GrafoP<tValue> &train, const GrafoP<tValue> &bus, typename GrafoP<tValue>::vertice change1, 
                        typename GrafoP<tValue>::vertice change2,typename GrafoP<tValue>::vertice origin, 
                        typename GrafoP<tValue>::vertice destiny){
    travel<tValue> trip;
    trip.cost = train.INFINITO;
    size_t size = train.numVert() + bus.numVert();
    GrafoP<tValue> superGrafo(size);
    superGrafo = superGrafo(train,bus);

    superGrafo[change1][change1+train.numVert()] = 0;
    superGrafo[change1+train.numVert()][change1] = 0;
    superGrafo[change2][change2+train.numVert()] = 0;
    superGrafo[change2+train.numVert()][change2] = 0;

    matriz<typename GrafoP<tValue>::vertice> vect;
    matriz<tValue> matrizCost = Floyd(superGrafo,vect);

    typename GrafoP<tValue>::vertice intermedio = vect[origin][destiny];
    typename GrafoP<tValue>::vertice vertice = origin;
    typename GrafoP<tValue>::vertice aux;

    // needed for std::vector< typename GrafoP<tValue>::vertice> vector of the path
    /*trip.path.push_back(origin);
    while(intermedio != destiny){
        trip.path.push_back(intermedio);

        aux = intermedio;
        intermedio = vect[vertice][intermedio];
        vertice = aux;
    }    
    trip.path.push_back(destiny);*/
    
    // This with tCamino, is a list of the path
    trip.path = camino(origin,destiny,matrizCost);
    trip.cost = matrizCost[origin][destiny];
    return trip;
}

//ver2

template <typename tValue>
travel<tValue> tripRoute2(const GrafoP<tValue> &train, const GrafoP<tValue> &bus, typename GrafoP<tValue>::vertice change1, 
                        typename GrafoP<tValue>::vertice change2,typename GrafoP<tValue>::vertice origin, 
                        typename GrafoP<tValue>::vertice destiny){
    travel<tValue> trip;
    trip.cost = train.INFINITO;
    std::vector<typename GrafoP<tValue>::vertice> vectTo;
    std::vector<typename GrafoP<tValue>::vertice> vectTd;
    std::vector<typename GrafoP<tValue>::vertice> vectBo;
    std::vector<typename GrafoP<tValue>::vertice> vectBd;
    std::vector<tValue> trainOriginCost = Dijkstra(train,origin,vectTo);
    std::vector<tValue> trainDestinyCost = DijkstraInvFloyd(train,destiny); // or DijkstraInverso(train,destiny,vectTd);
    std::vector<tValue> busOriginCost = Dijkstra(bus,origin,vectBo);
    std::vector<tValue> busDestinyCost = DijkstraInvFloyd(bus,destiny); // or DijkstraInverso(bus,destiny,vectBd);

    trip.cost = busOriginCost[change1] + trainOriginCost[change1];
    trip.path.insertar(origin,trip.path.fin());

    if(busOriginCost[change2] + trainOriginCost[change2] < trip.cost){
        trip.path.insertar(change2,trip.path.fin());
    }else trip.path.insertar(change1,trip.path.fin());

    return trip;
}
#pragma endregion


#pragma region Ejercicio8
template <typename tValue>
std::vector<tValue> DijkstraInvFloyd(const GrafoP<tValue> &grafo, const typename GrafoP<tValue>::vertice destiny){
    std::vector<tValue> cost(grafo.numVert());
    matriz<typename GrafoP<tValue>::vertice> aux;
    matriz<tValue> matrizCost = Floyd(grafo,aux);

    for(size_t i = 0; i < matrizCost.dimension(); i++){
        cost[i] = matrizCost[i][destiny];
    }

    return cost;
}

template <typename tValue>
tValue minimunPath(const GrafoP<tValue> &train, const GrafoP<tValue> &bus, const typename GrafoP<tValue>::vertice origin,
                    const typename GrafoP<tValue>::vertice destiny){
    std::vector<typename GrafoP<tValue>::vertice> vectTo;
    std::vector<typename GrafoP<tValue>::vertice> vectTd;
    std::vector<typename GrafoP<tValue>::vertice> vectBo;
    std::vector<typename GrafoP<tValue>::vertice> vectBd;
    std::vector<tValue> trainOriginCost = Dijkstra(train,origin,vectTo);
    std::vector<tValue> trainDestinyCost = DijkstraInvFloyd(train,destiny); // or DijkstraInverso(train,destiny,vectTd);
    std::vector<tValue> busOriginCost = Dijkstra(bus,origin,vectBo);
    std::vector<tValue> busDestinyCost = DijkstraInvFloyd(bus,destiny); // or DijkstraInverso(bus,destiny,vectBd);
    tValue pathCost = train.INFINITO;

    pathCost = std::min(busOriginCost[destiny],trainOriginCost[destiny]);
    
    for(size_t i = 0; i < train.numVert(); i++){
        if(i != origin && i != destiny) pathCost = std::min(pathCost,trainOriginCost[i] + busDestinyCost[i]);
    }

    for(size_t i = 0; i < train.numVert(); i++){
        if(i != origin && i != destiny) pathCost = std::min(pathCost,busOriginCost[i] + trainDestinyCost[i]);
    }

    return pathCost;
}
#pragma endregion


#pragma region Ejercicio9
template <typename tValue>
tValue trainBusTravel(std::vector<typename GrafoP<tValue>::vertice> &paths, const GrafoP<tValue> &train, const GrafoP<tValue> &bus,
                        const typename GrafoP<tValue>::vertice origin, const typename GrafoP<tValue>::vertice destiny, const tValue taxiCost){
    std::vector<typename GrafoP<tValue>::vertice> vectTo;
    std::vector<typename GrafoP<tValue>::vertice> vectTd;
    std::vector<typename GrafoP<tValue>::vertice> vectBo;
    std::vector<typename GrafoP<tValue>::vertice> vectBd;
    std::vector<tValue> trainOriginCost = Dijkstra(train,origin,vectTo);
    std::vector<tValue> trainDestinyCost = DijkstraInvFloyd(train,destiny); // or DijkstraInverso(train,destiny,vectTd);
    std::vector<tValue> busOriginCost = Dijkstra(bus,origin,vectBo);
    std::vector<tValue> busDestinyCost = DijkstraInvFloyd(bus,destiny); // or DijkstraInverso(bus,destiny,vectBd);
    tValue minimunCost = train.INFINITO;
    size_t aux;

    paths.push_back(origin);

    for(size_t i = 0; i < train.numVert(); i++){
        if(i != origin && i != destiny){
            if(minimunCost > trainOriginCost[i] + busDestinyCost[i] + taxiCost){
                minimunCost = trainOriginCost[i] + busDestinyCost[i] + taxiCost;
                aux = i;
            }
            if(minimunCost > busOriginCost[i] + trainDestinyCost[i] + taxiCost){
                minimunCost = trainOriginCost[i] + busDestinyCost[i] + taxiCost;
                aux = i;
            }
        }
    }

    if(minimunCost > trainOriginCost[destiny] || minimunCost > busOriginCost[destiny]){
        minimunCost = std::min(trainOriginCost[destiny],busOriginCost[destiny]);
    }else paths.push_back(aux);
    paths.push_back(destiny);

    return minimunCost;
}

template <typename tValue>
travel<tValue> trainBusTravel2(const GrafoP<tValue> &train, const GrafoP<tValue> &bus,
                        const typename GrafoP<tValue>::vertice origin, const typename GrafoP<tValue>::vertice destiny, const tValue taxiCost){
    travel<tValue> trav;
    size_t superGrafoSize(train.numVert() + bus.numVert());
    GrafoP<tValue> superGrafo = superGrafo(train,bus);

    for(size_t i = 0; i < train.numVert(); i++){
        for(size_t j = 0; j < train.numVert(); j++){
            superGrafo[train.numVert()+i][j] = taxiCost;
            superGrafo[i][train.numVert()+j] = taxiCost;
        }
    }

    matriz<typename GrafoP<tValue>::vertice> route;
    matriz<tValue> matrizCost = Floyd(superGrafo,route);
    
    if(matrizCost[origin][destiny] < matrizCost[origin][destiny+train.numVert()] && matrizCost[origin][destiny] < matrizCost[origin+train.numVert()][destiny]){
        trav.cost = matrizCost[origin][destiny];
        trav.path = camino(origin,destiny,route);
    }else if(matrizCost[origin][destiny+train.numVert()] < matrizCost[origin+train.numVert()][destiny]){
        trav.cost = matrizCost[origin][destiny+train.numVert()];
        trav.path = camino(origin,destiny+train.numVert(),route);
    }else{
        trav.cost = matrizCost[origin+train.numVert()][destiny];
        trav.path = camino(origin+train.numVert(),destiny,route);
    }
    
    return trav;
}
#pragma endregion


#pragma region Ejercicio10
typedef enum originTravel{BUS,TRAIN,PLANE};

template <typename tValue>
travel<tValue> triTravelFloyd(const GrafoP<tValue> &plane, const GrafoP<tValue> &train, const GrafoP<tValue> &bus, 
                        const typename GrafoP<tValue>::vertice origin, const typename GrafoP<tValue>::vertice destiny, 
                        const tValue taxi, const tValue taxiPlane){
    travel<tValue> trav;
    trav.cost = plane.INFINITO;
    tValue auxI,auxJ;

    GrafoP<tValue> grafoAux = superGrafo(bus,train);
    for(size_t i = 0; i < bus.numVert(); i++){
        for(size_t j = 0; j < bus.numVert(); j++){
            grafoAux[bus.numVert()+i][j] = taxi;
            grafoAux[i][bus.numVert()+j] = taxi;
        }
    }
    GrafoP<tValue> superGrafo = superGrafo(grafoAux,plane);
    for(size_t i = 0; i < bus.numVert(); i++){
        for(size_t j = 0; j < bus.numVert(); j++){
            superGrafo[grafoAux.numVert()+i][j] = taxiPlane;
            superGrafo[i][grafoAux.numVert()+j] = taxiPlane;
        }
    }

    matriz<typename GrafoP<tValue>::vertice> matrizAux;
    matriz<tValue> matrizCost = Floyd(superGrafo,matrizAux);
    for(size_t i = 0; i < 3; i++){
        for(size_t j = 0; j < 3; j++){
            if(trav.cost > matrizCost[origin+(plane.numVert()*i)][destiny+(plane.numVert()*i)]){
                trav.cost = matrizCost[origin+(plane.numVert()*i)][destiny+(plane.numVert()*i)];
                auxI = i;
                auxJ = j;
            }
        }
    }
    trav.path = camino(origin+(plane.numVert()*auxI),destiny+(plane.numVert()*auxJ),matrizAux);

    return trav;
}

template <typename tValue>
travel<tValue> triTravelDijkstra(const GrafoP<tValue> &plane, const GrafoP<tValue> &train, const GrafoP<tValue> &bus, 
                        const typename GrafoP<tValue>::vertice origin, const typename GrafoP<tValue>::vertice destiny, 
                        const tValue taxi, const tValue taxiPlane){
    originTravel ori;
    travel<tValue> trav;
    trav.cost = plane.INFINITO;
    tValue aux;
    std::vector<typename GrafoP<tValue>::vertice> vectBus;
    std::vector<typename GrafoP<tValue>::vertice> vectTrain;
    std::vector<typename GrafoP<tValue>::vertice> vectPlane;
    std::vector<tValue> costPathBus;
    std::vector<tValue> costPathTrain;
    std::vector<tValue> costPathPlane;
    GrafoP<tValue> grafoAux = superGrafo(bus,train);
    for(size_t i = 0; i < bus.numVert(); i++){
        for(size_t j = 0; j < bus.numVert(); j++){
            grafoAux[bus.numVert()+i][j] = taxi;
            grafoAux[i][bus.numVert()+j] = taxi;
        }
    }

    GrafoP<tValue> superGrafo = superGrafo(grafoAux,plane);
    for(size_t i = 0; i < bus.numVert(); i++){
        for(size_t j = 0; j < bus.numVert(); j++){
            superGrafo[grafoAux.numVert()+i][j] = taxiPlane;
            superGrafo[i][grafoAux.numVert()+j] = taxiPlane;
        }
    }

    costPathBus = Dijkstra(superGrafo,origin,vectBus);
    costPathTrain = Dijkstra(superGrafo,origin,vectTrain);
    costPathPlane = Dijkstra(superGrafo,origin,vectPlane);

    for(size_t i = 0; i < 3; i++){
        if(trav.cost > costPathBus[destiny+(costPathBus.size()*i)]){
            trav.cost = costPathBus[destiny+(costPathBus.size()*i)];
            ori = BUS;
            aux = i;
        }
        if(trav.cost > costPathTrain[destiny+(costPathTrain.size()*i)]){
            trav.cost = costPathTrain[destiny+(costPathTrain.size()*i)];
            ori = TRAIN;
            aux = i;
        }
        if(trav.cost > costPathPlane[destiny+(costPathPlane.size()*i)]){
            trav.cost = costPathPlane[destiny+(costPathPlane.size()*i)];
            ori = PLANE;
            aux = i;
        }
    }

    if(ori == BUS) trav.path = camino(origin,destiny+(costPathBus.size()*aux),vectBus);
    if(ori == TRAIN) trav.path = camino(origin+(costPathTrain.size()*1),destiny+(costPathTrain.size()*aux),vectTrain);
    if(ori == PLANE) trav.path = camino(origin+(costPathPlane.size()*2),destiny+(costPathPlane.size()*aux),vectPlane);

    return trav;
}
#pragma endregion


#pragma region Ejercicio11
template <typename tValue>
struct vertices{
    typename GrafoP<tValue>::vertice v1;
    typename GrafoP<tValue>::vertice v2;
};

template <typename tValue>
matriz<tValue> huries(const GrafoP<tValue> &island1, const GrafoP<tValue> &island2, const GrafoP<tValue> &island3, const std::vector<vertices<tValue>> &bridge){
    GrafoP<tValue> grafo = superGrafo(island1,island2);
    GrafoP<tValue> superGrafo = superGrafo(grafo,island3);
    
    for(size_t i = 0; i < bridge.size(); i++){
        superGrafo[bridge[i].v1][bridge[i].v2] = 0;
        superGrafo[bridge[i].v2][bridge[i].v1] = 0;
    }

    matriz<typename GrafoP<tValue>::vertice> aux;
    matriz<tValue> cost = Floyd(superGrafo,aux);
    
    return cost;
}
#pragma endregion