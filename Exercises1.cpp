#include <vector>
#include <algorithm>
#include "grafoPMC.h"
#include "./LibGrafos/alg_grafoPMC.h"
#include "./LibGrafos/matriz.h"
#include "./LibGrafos/particion.h"

/*Create a generic function, called DijkstraInv, to solve the inverse of Dijkstra's problem, with the same
parameter and result types than the function already included for it. The new function, therefore, must find
the minimum cost path to a destination from each vertex of the graph and its corresponding cost.*/

/*Due to the implementation of the weighted graph, it has an infinite value therefore if one of the two elements is
infinite, that value will be returned instead of adding them, otherwise these two elements will be added.*/ 
#pragma region Ejercicio1
#pragma region firstTry
template <typename tValue>
vector<tValue> dijkstraInv(const GrafoP<tValue> &grafo, typename GrafoP<tValue>::vertice destino, vector<typename GrafoP<tValue>::vertice>& aux){
    matriz<typename GrafoP<tValue>::vertice> vecAux;
    matriz<tValue> matrizFloyd = Floyd(grafo,vecAux);
    vector<tValue> sVector;

    for(size_t i = 0; i < grafoFloyd.dimension(); i++) sVector[i] = matrizFloyd[destino][i];

    return sVector;
}
#pragma endregion

/*###################################################################################################################*/

#pragma region secondTry
template <typename tCoste>
vector<tCoste> DijkstraInv(const GrafoP<tCoste>& G, typename GrafoP<tCoste>::vertice destino, vector<typename GrafoP<tCoste>::vertice>& P){
    typedef typename GrafoP<tCoste>::vertice v, w;
    const size_t n = G.numVert();
    vector<bool> S(n, false);
    vector<tCoste> D;

    for(size_t i = 1; i < n; i++) D[i] = G[i][destino];

    D[destino] = 0;
    S[destino] = true;
    P = vector<vertice>(n, destino);

    for(size_t i = 1; i <= n-2; i++){
        tCoste costeMin = GrafoP<tCoste>::INFINITO;

        for(v = 0; v < n; v++){
            if(!S[v] && D[v] <= costeMin){
                costeMin = D[v];
                w = v;
            }
        }

        S[w] = true;

        for(v = 0; v < n; v++){
            if(!S[v]){
                tCoste Owv = suma(D[w], G[w][v]);
                
                if(Owv < D[v]){
                    D[v] = Owv;
                    P[v] = w;
                }
            }
        }
    }
    return D;
}
#pragma endregion
#pragma endregion


/*We will define the pseudocenter of a connected graph as the node that minimizes the sum of the minimum distances to its two furthest nodes.
We will define the diameter of the graph as the sum of the minimum distances to the two nodes furthest from the pseudocenter of the graph.
Given a connected graph represented by a cost matrix, implement a subprogram that returns the length of its diameter.*/
#pragma region Ejercicio2
template <typename tValue>
tValue distance(matriz<tValue> minimunPaths, size_t i, tValue INFINITO){
    tValue maximun1 = -INFINITO;
    tValue maximun2 = -INFINITO;
    
    for(size_t j = 0; j < minimunPaths.dimension(); j++){
        if(j != i && minimunPaths[i][j] >= maximun1){
            maximun2 = maximun1,
            maximun1 = minimunPaths[i][j];
        }else if(minimunPaths[i][j] >= maximun2) maximun2 = minimunPaths[i][j];
    }
    
    if(maximun2 != INFINITO && maximun1 != INFINITO) return suma(maximun1,maximun2);
    else return INFINITO;
}

template <typename tValue>
tValue pseudocenter(const GrafoP<tValue> &grafo){
    tValue minimun = grafo.INFINITO;
    tValue distance = grafo.INFINITO;
    matriz<typename GrafoP<tValue>::vertice> auxMatriz;
    matriz<tValue> minimunPaths = Floyd(grafo,auxMatriz);

    for(size_t i = 0; i < minimunPaths.dimension(); i++){
        distance = distance(minimunPaths,i,grafo.INFINITO);
        if(minimun >= distance) minimun = distance;
    }

    return minimun;
}
#pragma endregion


/*Check if a Cost Graph is Acyclic*/
#pragma region Ejercicio3
#pragma region PartitionOption
template <typename tValue>
bool acyclic(const GrafoP<tValue> &grafo){
    Particion particion(grafo.numVert());
    typename GrafoP<tValue>::vertice i,j;

    for(i = 0; i < grafo.numVert(); i++){
        for(j = 0; j < grafo.numVert(); j++){
            if(grafo[i][j] != grafo.INFINITO){
                int vertI = particion.encontrar(i);
                int vertJ = particion.encontrar(j);

                if(vertI == vertJ) return false;
                else particion.unir(vertI,vertJ);
            }
        }
    }

    return true;
}
#pragma endregion

/*###################################################################################################################*/

#pragma region vectors
template <typename tValue>
bool ciclo(const GrafoP<tValue> &grafo, std::vector<bool> &visitados, std::vector<bool> &apilados, typename GrafoP<tValue>::vertice indice){
    typename GrafoP<tValue>::vertice i;
    visitados[indice] = true; 
    apilados[indice] = true;

    for(i = 0; i<grafo.numVert(); i++){
        if(grafo[indice][i] != grafo.INFINITO && !visitados[i]){
            if(ciclico(grafo,visitados,apilados,i)) return true;
        }else if(grafo[indice][i] != grafo.INFINITO && apilados[i]) return true;
    }

    apilados[indice] = false;
    return apilados[indice];
}

template <typename tValue>
bool grafoAciclico(const GrafoP<tValue> &grafo){
    bool aciclico = true;
    std::vector<bool> visitados(grafo.numVert(),false),apilados(grafo.numVert(),false);;
    typename GrafoP<tValue>::vertice indice;

    for(indice = 0; indice<grafo.numVert(); indice++){
        if(!visitados[indice] && ciclo(grafo,visitados,apilados,indice)){
            aciclico = false;
        }
    }

    return aciclico;
}
#pragma endregion
#pragma endregion


/*A study needs to be done on the minimum distances necessary to travel between
any two cities in a country called Zueland. The problem is simple but
A few small details must be taken into account:
 a) The orography of Zuelandia is a bit special, the roads are very narrow
 and therefore only allow one direction of circulation.

 b) Currently Zueland is a country at war. And in fact there are a series of
 cities in the country that have been taken by the rebels, so they cannot
 be used for traveling.

 c) The rebels have not only taken over certain cities in the country, but
 They have also cut off certain roads, (so these roads cannot be
 used).

 d) But the government cannot remain impassive in the face of the situation and has demanded
 that absolutely all trips made through the country pass through the capital
 of the same, where the pertinent security controls will be carried out.

Given these four conditions, it is requested to implement a subprogram that, given
• the network (cost matrix) of Zuelandia in a normal situation,
• the list of the cities taken by the rebels,
• the list of roads cut by the rebels
• and the capital of Zueland,
Calculate the minimum cost matrix for traveling between any two cities
Zuelandes in this situation.*/
#pragma region Ejercicio4
template <typename tValue>
matriz<tValue> zuelandia(const GrafoP<tValue> &grafo, const std::vector<typename GrafoP<tValue>::vertice> &ciudadesTomadas, 
                        const std::vector<std::vector<bool>> &carreterasCortadas, const typename GrafoP<tValue>::vertice capital){
    typename GrafoP<tValue>::vertice i,j;
    GrafoP<tValue> grafoAux = grafo;
    matriz<tValue> mapaFinal;
    matriz<typename GrafoP<tValue>::vertice> aux;

    for(i = 0; i<ciudadesTomadas.size(), i++){
        for(j = 0; j<grafoAux.numVert(), j++){
            grafoAux[ciudadesTomadas[i]][j] = grafoAux.INFINITO;
            grafoAux[j][ciudadesTomadas[i]] = grafoAux.INFINITO;
        }   
    }

    for(i = 0; i<carreterasCortadas.size(), i++){
        for(j = 0; j<carreterasCortadas.size(), j++){
            if(carreterasCortadas[i][j]) grafoAux[i][j] = grafoAux.INFINITO;
        }   
    }

    mapaFinal = Floyd(grafoAux,aux);

    for(i = 0; i<mapaFinal.dimension(), i++){
        for(j = 0; j<mapaFinal.dimension(), j++){
            if(mapaFinal[i][capital] == grafoAux.INFINITO || mapaFinal[capital][j] == grafoAux.INFINITO) mapaFinal[i][j] = grafoAux.INFINITO;
            else if(i != j) mapaFinal[i][j] = mapaFinal[i][capital] + mapaFinal[capital][j];
        }   
    }

    return mapaFinal;
}

#pragma endregion