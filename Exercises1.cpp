#include <vector>
#include <algorithm>
#include "grafoPMC.h"
#include "./LibGrafos/alg_grafoPMC.h"
#include "./LibGrafos/matriz.h"
#include "./LibGrafos/particion.h"

 
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