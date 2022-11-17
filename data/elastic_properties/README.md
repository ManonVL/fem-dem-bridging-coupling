## Elastic properties of DEM samples

To determine the elastic properties of each box, we first generate the DEM boxes, then we relax the system, and finally we compute the elastic properties.

To generate a DEM box, the user must launch:

- python3 box_0.006/seed_1/create_random.py

To relax a DEM box, the user must launch the script relax.in using LAMMPS:

- lmp -in box_0.006/seed_1/pressure_5.0/relax.in

To generate the elastic properties, the user must launch the script in.elastic using LAMMPS:

- lmp -in box_0.006/seed_1/pressure_5.0/elastic_properties/in.elastic

Once the data are generated, the user can plot the results using the script main.ipynb located in the notebooks. 

