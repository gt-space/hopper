from general_fluid_network import Node, Ambient, Connection, Network

source = Ambient("Krypton", 101325*6000/14.7, 293.15, "Supply")
print(source.phase)
fill = Connection(5e-7, location=1)
tank = Node("Krypton", 0.02, 5, 293, "vehicle")
cart = {fill: (source, tank)}
network = Network(cart)

network.sim(600, 0.1)
network.plot_nodes_overlay([tank], units="SI")
network.plot_connections_overlay([fill])
