class Atom:
    def __init__(self, id, type, diameter, density, x, y, z):
        self.id = int(id)
        self.type = int(type)
        self.diameter = float(diameter)
        self.density = float(density)
        self.position = [float(x), float(y), float(z)]

