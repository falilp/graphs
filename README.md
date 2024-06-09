# Graphs Statements
## Exercises1
 - <details open>
    <summary>1</summary>
    <br>
        Create a generic function, called DijkstraInv, to solve the inverse of Dijkstra's problem, with the same
    parameter and result types than the function already included for it. The new function, therefore, must find
    the minimum cost path to a destination from each vertex of the graph and its corresponding cost.<p>

    Due to the implementation of the weighted graph, it has an infinite value therefore if one of the two elements is
    infinite, that value will be returned instead of adding them, otherwise these two elements will be added.
    </details>
 - <details open>
    <summary>2</summary>
    <br>
        We will define the pseudocenter of a connected graph as the node that minimizes the sum of the minimum distances to its two furthest nodes.<p><p>
        We will define the diameter of the graph as the sum of the minimum distances to the two nodes furthest from the pseudocenter of the graph.<p>
        Given a connected graph represented by a cost matrix, implement a subprogram that returns the length of its diameter.
    </details>
 - <details open>
    <summary>3</summary>
    <br>
        Check if a Cost Graph is Acyclic
    </details>
 - <details open>
    <summary>4</summary>
    <br>
        A study needs to be done on the minimum distances necessary to travel between any two cities in a country called Zueland. The problem is simple but A few small details must be taken into account:<p>
        a) The orography of Zuelandia is a bit special, the roads are very narrow
        and therefore only allow one direction of circulation.<p>
        b) Currently Zueland is a country at war. And in fact there are a series of
        cities in the country that have been taken by the rebels, so they cannot
        be used for traveling.<p>
        c) The rebels have not only taken over certain cities in the country, but
        They have also cut off certain roads, (so these roads cannot be
        used).<p>
        d) But the government cannot remain impassive in the face of the situation and has demanded that absolutely all trips made through the country pass through the capital of the same, where the pertinent security controls will be carried out.<p>
        Given these four conditions, it is requested to implement a subprogram that, given:<p>
        • the network (cost matrix) of Zuelandia in a normal situation,<p>
        • the list of the cities taken by the rebels,<p>
        • the list of roads cut by the rebels<p>
        • and the capital of Zueland,<p>
        Calculate the minimum cost matrix for traveling between any two cities
        Zuelandes in this situation.
    </details>
 - <details open>
    <summary>3</summary>
    <br>
        Three graphs representing the Cost Matrix for Travel in a
        certain country but by different means of transport, of course all graphs
        They will have the same number of nodes. The first graph represents the costs of going for
        road, the second by train and the third by plane. Given a traveler who has a
        certain amount of money, which is allergic to one of the three means of transport, and
        which leaves a given city, implement a subprogram that determines the
        cities that our target traveler could reach.
    </details>
## Exercises2
 - <details open>
    <summary>1</summary>
    <br>
        Your travel agency “OTRAVEZUNGRAFO S.A.” He faces a curious client.<p>
        He is an amazing character, he doesn't care about money and he wants to make the trip more expensive
        possible between the cities you offer. 
        <p>Your goal is to spend as much money as possible. Possible money and does not care about the origin or destination of the trip.<p>
        Knowing that it is impossible to pass through the same city twice, since by chance
        Your travel agency's graph turned out to be acyclical, returning the cost, origin and destination
        of such a curious journey. It is based on the matrix of direct costs between the cities of the graph.
    </details>
 - <details open>
    <summary>2</summary>
    <br>
        There is a maze of NxN squares of which the entry squares are known.
        and exit from it. If you are in one square you can only move in the following ones
        Four directions (up, down, right, left). On the other hand, among some of the
        squares there is a wall that prevents movement between the two squares that separate said wall
        (Otherwise it wouldn't be a real labyrinth).<p>
        Implements an applet that given
            N (maze dimension),
            the maze wall list,
            the entry box and the exit box,
        Calculate the shortest path to go from the entrance to the exit and its length.
    </details>
 - <details open>
    <summary>3</summary>
    <br>
        You are the proud owner of a distribution company. Your mission lies in
        distribute all your stock between the different cities in which your company has
        store.<p>
        You have a graph represented by the cost matrix, in which the cost appears
        (per unit of product) to transport the products between the different cities of the
        graph.<p>
        But it also turns out that the City Councils of the different cities in which
        You have a warehouse, they are very interested in you storing your products in them, so
        who are willing to subsidize you with a percentage of the minimum costs of
        transportation to the city. To make the problem easier, we will consider negligible the
        costs of returning the truck to its base (production center).<p>
        Here is your problem. You have<p>
        -the production center, origin node in which you have your product (it does not have
        store),
        -a number of units of product (quantity),<p>
        -the cost matrix of the distribution graph with N cities,<p>
        -the storage capacity of each of them,<p>
        -the percentage of subsidy (over the minimum expenses) that each offers you
        City hall.<p>
        Different cities (warehouses) may have different capacities, and also the
        total capacity may be greater than the available quantity of product, so
        You must decide how many units of product you store in each of the cities.<p>
        You must also take into account the subsidies that you will receive from the different
        City Councils, which may be different in each one and will be between 0% and
        100% of minimum costs.<p>
        The solution to the problem must include the quantities to be stored in each city under
        these conditions and the minimum total cost of the distribution operation for your
        company
    </details>
 - <details open>
    <summary>4</summary>
    <br>
        You are the proud owner of the company "Cementos de Zuelandia S.A". Company
        Dedicated to the manufacture and distribution of cement, located in the capital of Zuelandia. For
        The distribution of cement among your different clients (cities of Zuelandia)
        You have a fleet of trucks and a template of Zuelandes.<p>
        The problem to solve has to do with the character of the Zuelandés. The Zuelandés is
        a person who takes too many "kmextra" in his work, in fact, you have
        Founded suspicions that your drivers use company trucks to
        particular uses (that is, improper, and at your coast) so you want to control the
        kilometers that travel your trucks.<p>
        Every day the work part is generated, which includes the number of
        cement loads (1 load = 1 truck full of cement) that you must send to each
        client (client = city of Zuelandia). It is unnecessary to indicate that not every day there are
        that send loads to all customers, and also, you can reasonably assume that you
        Truck fleet is able to do the daily work.<p>
        For problem solving, it is perhaps interesting to remember that Zuelandia is a
        country whose special orography only allows roads to have a sense of
        circulation.<p>
        Implements a function that given the graph with the direct distances between the
        different Zuelandesas cities, the daily work part, and the capital of Zueland,
        return the total distance in kilometers that must travel your trucks in the day, to
        that you can discover whether or not they use your trucks in activities outside the
        company.
    </details>
 - <details open>
    <summary>5</summary>
    <br>
        Three graphs representing the Cost Matrix for Travel in a
        certain country but by different means of transport, of course all graphs
        They will have the same number of nodes. The first graph represents the costs of going for
        road, the second by train and the third by plane. Given a traveler who has a
        certain amount of money, which is allergic to one of the three means of transport, and
        which leaves a given city, implement a subprogram that determines the
        cities that our target traveler could reach.
    </details>
 - <details open>
    <summary>6</summary>
    <br>
        The following situation is raised to the owner of a transport agency. The
        Travel Agency offers different combined trajectories between n Spanish cities
        using train and bus. Two graphs representing the costs are available (matrix of
        costs) to travel between different cities, on the one hand by train, and on the other by bus
        (Of course among cities that have a direct line between them). It also coincides
        that taxis throughout Spain are currently in general strike, which
        implies that only transport can be changed in a given city in which,
        By chance, train and bus stations are united.
        Implements a function that calculates the minimum rate (minimum cost matrix)
        Travel between any of the n cities having the cost graph by bus,
        of the graph of costs by train, and of the city that has the United Stations.
    </details>
 - <details open>
    <summary>7</summary>
    <br>
        There are two graphs (cost matrix) that represent the costs of traveling between<p>
        N Spanish cities using the train (first graph) and the bus (second graph).
        Both graphs represent trips between the same N cities.
        Our objective is to find the minimum cost path to travel between two cities
        specific characteristics of the network, origin and destination, under the following conditions:<p>
        -The city of origin only has transportation by train.<p>
        -The destination city only has bus transportation.<p>
        -The taxi sector, quite conflictive in our problems, is still on strike,<p>
        so it is only possible to change transportation in two cities in the country.
        graph, change1 and change2, where the train and bus stations are
        united.<p>
        Implements an applet that calculates the route and the minimum cost to travel between
        the cities of Origin and Destination under these conditions.
    </details>
 - <details open>
    <summary>8</summary>
    <br>
        “ONE TRANSFER, PLEASE.” This is the title that reads in your
        brand new travel company. Your advertising explains, of course, that you offer trips
        TRAIN and/or BUS combinations (i.e. trips by train, bus, or using
        both), between N cities in the country, which offer unbeatable service, very
        competitive, and that you guarantee before a notary something that none of your
        competitors: that on all your trips AT MAXIMUM there will be only one transfer
        (change of means of transport).<p><p>
        Well, today is July 1st and the travel season begins.
        Lucky! A client just showed up at your office. He explains to you that he wants to travel
        between two cities, Origin and Destination, and he wants to know how much it will cost him.
        To answer this question you have two direct cost graphs (cost matrix).
        costs of traveling between the N cities in the country, a graph with the costs of traveling by train and
        another by bus.<p>
        Implement an applet that calculates the minimum rate under these conditions.
    </details>
 - <details open>
    <summary>9</summary>
    <br>
        There are two graphs that represent the cost matrix for trips in a
        certain country, but by different means of transport (train and bus, for example).
        Of course both graphs will have the same number of nodes, N. Given
        both graphs, a city of origin, a city of destination and the cost of the taxi to
        change from one station to another within any city (assumed constant and equal
        for all cities), implements an applet that calculates the path and the cost
        minimum to go from the city of origin to the city of destination.
    </details>
 - <details open>
    <summary>10</summary>
    <br>
        There are three graphs that represent the cost matrix for trips in a
        certain country, but by different means of transport (train, bus and plane). By
        Of course the three graphs will have the same number of nodes, N.
        Given the following data:<p>
        -the three graphs,<p>
        -a city of origin,<p>
        -a city of destination,<p>
        -the cost of the taxi to change, within a city, from the train station to the
        bus or vice versa (taxi-train-bus) and<p>
        -the cost of the taxi from the airport to the train or bus station, or
        vice versa (taxi-airport-train/bus) and assuming that both taxi costs (different from each other, they are two different costs) are
        constant and equal for all cities, implements a subprogram that calculates the
        path and the minimum cost to go from the origin city to the destination city.
    </details>
 - <details open>
    <summary>11</summary>
    <br>
        We have three graphs (cost matrix) that represent the direct costs of
        travel between the cities of three of the islands of the Hurí archipelago (Zuelandia).
        In order to travel from one island to another, there are a series of bridges that connect
        cities of the different islands at a truly affordable price (by decision of the
        Prefect of the Houris, the use of the bridges is absolutely free).<p>
        If the student wants to simplify the problem, he can number the N1 cities on the island
        1, from 0 to N1-1, the N2 cities of island 2, from N1 to N1+N2-1, and the N3 of the last one, from
        N1+N2 to N1+N2+N3-1.<p>
        Having the three direct cost matrices of traveling within each of
        the islands, and the list of bridges between their cities, calculate the minimum costs
        to travel between any two cities on these three islands.
    </details>
## Exercises3