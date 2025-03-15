import networkx as nx
import matplotlib.pyplot as plt


# Read data and construct the graph
# Create an empty dictionary to store the mapping of node ids to city/airport names
cities = {}

# open global-cities.dat
with open("global-cities.dat", "r", encoding="utf-8") as f:
    # Read the contents of a file line by line
    for line in f:
        # Remove whitespace at the beginning and end and separate each row by “|”.
        parts = line.strip().split("|")
        # Ensure that each row has three data items (airport short name, node id, city where airport is located)
        if len(parts) < 3:
            continue  # If NULL exist, pass to next line

        airport_short = parts[0]  # Airport abbreviation
        node_id = parts[1]  # Node id
        airport_city = parts[2]  # City name

        # Combine airport short name and city name into one string
        full_name = f"{airport_short} ({airport_city})"

        # Maps the node id to full_name and stores it in the dictionary
        cities[node_id] = full_name

# Create an empty undirected graph
G = nx.Graph()

# open global-net.dat
with open("global-net.dat", "r", encoding="utf-8") as f:
    # Read the contents of a file line by line
    for line in f:
        parts = line.strip().split()  # Split each line by space
        if len(parts) < 2:
            continue  # If there are fewer than 2 elements after splitting, skip this line.
        node1 = parts[0]  # The first node of the edge
        node2 = parts[1]  # The second node of the edge
        # Adding edges to the graph
        G.add_edge(node1, node2)


# Question 1: Basic graph attributes
print("Question 1: Basic graph attributes")

# Get the number of nodes in the graph
num_nodes = G.number_of_nodes()
# Get the number of edges in the graph
num_edges = G.number_of_edges()

print("Number of nodes:", num_nodes)
print("Number of undirected edges:", num_edges)


# Question 2: Connected components and the largest component
print("\nQuestion 2: Connected components")

# Use networkx's connected_components method to get all the connected components in the graph
# Returned is a list containing several sets, each of which holds the node ids of the same connectivity component
components = list(nx.connected_components(G))

# Get the number of connected components
num_components = len(components)
print("Total number of connected components:", num_components)

# Iterate to compare the size of each collection
largest_cc = components[0]  # Set the first connected component to maximum
for comp in components:
    if len(comp) > len(largest_cc):
        largest_cc = comp

# Create a subgraph based on the largest connected component
G_largest = G.subgraph(largest_cc).copy()

print("Number of nodes in the largest component:", G_largest.number_of_nodes())
print("Number of edges in the largest component:", G_largest.number_of_edges())



# Question 3: Top 10 highest-degree nodes in the largest component (displaying city/airport names)
print("\nQuestion 3: Top 10 highest-degree nodes in the largest component")

# Compute the degree of each node in the maximally connected component, get dictionary {node: degree}
degree_dict = dict(G_largest.degree())

# Turn it into a list
degree_items = list(degree_dict.items())

# Define a function to return the number of degrees in the tuple
def get_degree(item):
    return item[1]

# Use the sort() method to sort by degree in descending order.
degree_items.sort(key=get_degree, reverse=True)

# Print the first 10 nodes after sorting
count = 1
for item in degree_items[:10]:
    node = item[0] # Node id
    degree = item[1] # Corresponding degree
    # Get the name of the airport/city that corresponds to this node from the cities dictionary
    print(str(count) + ". " + str(cities.get(node, node)) + " -- " + str(degree) + " connections")
    count = count + 1



# Question 4: Plot the degree distribution of the largest component
print("\nQuestion 4: Plotting degree distribution")

# Use a list to store the degree of each node in the maximally connected component
degrees = []
for n, d in G_largest.degree():
    degrees.append(d)

# Eliminate duplicates and sort
unique_degrees = list(set(degrees))
unique_degrees.sort()

# Calculate the number of occurrences corresponding to each unique degree
counts = []
for x in unique_degrees:
    count_x = degrees.count(x)
    counts.append(count_x)

# Calculate the proportion of occurrences of each degree
fractions = []
total_nodes = len(degrees)
for count_x in counts:
    fractions.append(count_x / float(total_nodes))

# Plotting degree distributions at linear scales
plt.figure()
plt.plot(unique_degrees, fractions, 'bo')
plt.xlabel("Degree (x)")
plt.ylabel("Fraction of nodes (y)")
plt.title("Degree Distribution (Linear Scale)")
plt.xlim(min(unique_degrees), max(unique_degrees))
plt.savefig("degree_distribution_linear.png", dpi=300)
plt.show()

# Plotting the degree distribution on a log-log scale
plt.figure()
plt.loglog(unique_degrees, fractions, 'bo')
plt.xlabel("Degree (log10)")
plt.ylabel("Fraction of nodes (log10)")
plt.title("Degree Distribution (Log-Log Scale)")
plt.xlim(min(unique_degrees), max(unique_degrees))
plt.savefig("degree_distribution_loglog.png", dpi=300)
plt.show()


# Question 5: Diameter and one longest shortest path (using city/airport names)
print("\nQuestion 5: Diameter and longest shortest path of the largest component")

# Calculated Diameter
diameter = nx.diameter(G_largest)
print("Diameter of the largest component:", diameter)

# Calculate the eccentricity of each node
ecc = nx.eccentricity(G_largest)

# Find a node with centroid equal to the diameter as a starting point
u = None
for node in ecc:
    if ecc[node] == diameter:
        u = node
        break

# Iterate through all the nodes and find the node v whose shortest path length to u is equal to the diameter
v = None
for node in G_largest.nodes():
    if nx.shortest_path_length(G_largest, source=u, target=node) == diameter:
        v = node
        break

# Find the shortest path between u and v. The length of this path is the diameter.
path = nx.shortest_path(G_largest, source=u, target=v)
print("One longest shortest path (sequence of cities/airports):")
for node in path:
    print(" -", cities.get(node, node))


# Question 6: Shortest route from Canberra to Cape Town
print("\nQuestion 6: Shortest route from Canberra to Cape Town")

start_node = None
end_node = None

# Find nodes with “Canberra” and “Cape Town” in their name.
for node in cities:
    # If the airport/city name contains “Canberra”, the corresponding node id is recorded
    if "Canberra" in cities[node]:
        start_node = node
    # If the airport/city name contains “Cape Town”, record the corresponding node id
    if "Cape Town" in cities[node]:
        end_node = node

# Check if all the corresponding nodes have been found
if start_node is None or end_node is None:
    print("Cannot find the nodes for Canberra and/or Cape Town in the dataset.")
else:
    print("Canberra:", cities.get(start_node))
    print("Cape Town:", cities.get(end_node))
    try:
        # Find the shortest path from Canberra to Cape Town using the functions included in the networkx library.
        shortest_route = nx.shortest_path(G_largest, source=start_node, target=end_node)
        print("Shortest path from Canberra to Cape Town:")
        for node in shortest_route:
            print(" -", cities.get(node, node))
        # The number of edges in the path is the number of flights required
        print("Number of flights required:", len(shortest_route) - 1)
    except nx.NetworkXNoPath:
        print("No path found between Canberra and Cape Town in the largest component.")


# Question 7: Top 10 cities/airports by betweenness centrality
print("\nQuestion 7: Top 10 cities/airports by betweenness centrality")

# Calculate the betweenness centrality of all nodes in the maximum connectivity component
bc = nx.betweenness_centrality(G_largest)
bc_items = list(bc.items())

# Returns the centrality value in the tuple
def get_centrality(item):
    return item[1]

# Sort by centrality from largest to smallest
bc_items.sort(key=get_centrality, reverse=True)

# Print the first 10 nodes of the sorted list with a counter.
count = 1
for item in bc_items[:10]:
    node = item[0] # Node id
    centrality = item[1] # Corresponding centrality
    # Print the name of the airport/city and the centrality of the node
    print(str(count) + ". " + str(cities.get(node, node)) + " -- Betweenness centrality: " + str(round(centrality, 3)))
    count = count + 1
