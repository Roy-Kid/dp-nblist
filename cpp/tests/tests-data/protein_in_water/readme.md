This is a model of protein in water.
### Creating a model using packmol and moltemplate:
1. Prepare the PDB files for the [protein](https://www.rcsb.org/structure/7MYZ) and water.
2. Create a packmol input file named gen.inp:
```
# gen.inp
tolerance 2.0
filetype pdb
output 7myz-water.pdb

structure 7myz.pdb
    number 6
    inside cube 0. 0. 0. 300
end structure

structure water.pdb
    number 1000
    inside cube 0. 0. 0. 300
end structure
```
Execute the following command: `packmol < gen.inp`

3. Use Ovito to export the protein and water PDB files into LAMMPS format files, selecting the "atomic" type (containing only atomic coordinate information).
4. Use the `ltemplify.py` script to convert the LAMMPS files into LT files:  
`ltemplify.py -name water -molid "1" water.lmp > water.lt`
5. Create a `system.lt` file with the following content:
```
# system.lt
import "7myz.lt"
import "water.lt"

w7myz = new w7myz[6]  # Corresponding to the order used in packmol
Water = new water[1000]  # Corresponding to the order used in packmol

write_once("Data Boundary") {
   0.0   300.0  xlo xhi
   0.0   300.0  ylo yhi
   0.0   300.0  zlo zhi
}
```
Run the following command:  
`moltemplate.sh -nocheck -atomstyle "atomic" -pdb 7myz-water.pdb system.lt`

Reference Link:
https://www.lammps.org.cn/zh/tools/moltemplate/packmol.html#%E4%B8%AD%E6%80%A7-%E6%BA%B6%E5%89%82%E5%8C%96
