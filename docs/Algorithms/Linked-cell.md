# Linked-List Cell NBL
Consider a particle system on a two-dimensional plane, where numerous particles are scattered in space. By establishing neighbor relationships, one needs to traverse the neighboring particles, which in this case are all particles, and then examine the characteristic distances to determine if they meet the conditions for neighbor relationships, such as Euclidean distance. However, this approach has a complexity of $O(n^2)$, and as the number of particles increases, the computational cost will sharply rise. Building upon this, the method based on a grid-based Cell List is proposed. Its core idea involves partitioning space into grid cells and placing each particle in the corresponding grid cell. When traversing neighboring particles, it is only necessary to traverse particles in the neighboring cells of the particle's cell. The subsequent process of examining characteristic distances remains the same. This method reduces the complexity of searching for neighboring particles to $O(n)$, achieving efficiency improvement.
![verlet-list](https://github.com/okihane/dp-nblist/assets/30775452/ce59ddbb-100b-44e2-b66f-f04786dc730a)
The specific implementation can be roughly divided into three steps: grid partitioning, establishing particle-grid relationships, and constructing a grid-based adjacency list for neighboring particles.
![step_new](https://github.com/okihane/dp-nblist/assets/30775452/0cf8daa2-7296-4ac1-9aa8-ae0dfabe5b27)

## Meshing
First divide the simulation box into small cells of equal size. The edge lengths of each cell, $(r_{cx}, r_{cy}, r_{cz})$, must be at least $r_c$ ; We use $L_{cα}(α = x, y, z)$ to represent the cell length in each direction, where   
$$L_{cα} = \frac{L_α}{r_c}$$
and $L_α$ is the simulation box length in the α direction.  The number of cells to accommodate all these atoms is $L_{cx}L_{cy}L_{cz}$. We identify a cell with a vector cell index, $\vec{c} = (c_x, c_y, c_z) (0 \leq c_x ≤ L_{cx}−1; 0 ≤ c_y ≤ L_{cy}−1; 0 ≤ c_z ≤ L_{cz}−1)$, and a serial cell index (see the figure below),   
$$c = c_xL_{cy}L_{cz} + c_yL_{cz} + c_z$$
or  
$$c_x = c/(L_{cy}L_{cz})$$
$$c_y = (c/L_{cz}) \ mod \ L_{cy}$$
$$c_z = c \ mod \ L_{cz}$$
An atom with coordinate $\vec{r}$ belongs to a cell with the vector cell index,  
$$c_α = \frac{\vec{r}α}{\vec{r}_{cα}} (α = x, y, z)$$
An atom within a cell interacts only with other atoms in the same unit and those in the eight adjacent units in two-dimensional space or 26 adjacent units in three-dimensional space.
![2d-cell](https://github.com/okihane/dp-nblist/assets/30775452/057a956c-58ec-45b4-a3dd-79cfda9301df)

## Linked list
The atoms belonging to a cell is organized using linked lists. The following describes the data structures and algorithms to construct the linked lists and compute the interatomic interaction using them.
![linked-list](https://github.com/okihane/dp-nblist/assets/30775452/09a405e5-271b-4852-bbd0-b5475fffab43)
 ```
size_t n_atoms = xyz.size();
_natoms += n_atoms;
_lscl.resize(_natoms, EMPTY);
Vec3<int> xyz_cell_index;
xyz = _box->wrap(xyz);

for (size_t i = 0; i < n_atoms; i++)
{
    xyz_cell_index = Vec3<int>(xyz[i] / _r_cutoff);
    for (int j = 0; j < 3; j++){
        if (xyz_cell_index[j] == _cell_length[j])
            xyz_cell_index[j] = xyz_cell_index[j] - 1;
    }
    size_t cell_index = get_cell_index(xyz_cell_index);
    _lscl[i] = _head[cell_index];
    _head[cell_index] = i;
}
 ```
`_lscl[NMAX]`: An array implementation of the linked lists. `_lscl[i]` holds the atom index to which the ith atom points.  
`_head[NCLMAX]`: `_head[c]` holds the index of the first atom in the c-th cell, or `_head[c] = EMPTY (= −1)` if there is no atom in the cell.

## Neighbor Cell Searching
In the assumption of three-dimensional space, the box is divided into nx, ny, nz cells along the x, y, z directions, respectively. The vector of the i-th cell is ($x_i$ , $y_i$ , $z_i$), and its 26 neighboring cells are respectively ($x_i-1$, $y_i-1$, $z_i-1$), ($x_i$ , $y_i-1$, $z_i-1$), ($x_i+1$, $y_i-1$, $z_i-1$), ... , ($x_i+1$, $y_i+1$, $z_i+1$). If periodic boundary conditions are considered, wrapping needs to be applied to the vectors of neighboring cells. We assume that the vector ($x_i$ , $y_i$ , $z_i$) of each cell can be expressed as:
$$\textbf{r}=\sum_i{f_i a_i}$$
where $a_i$ represents the three basis vectors of the box, $f_i$ represents the projection of the cell position on the basis vectors. This can also be represented in matrix form as:
$$\textbf{r}=\textbf{A}\textbf{f}$$
where A is composed of the basis vectors of the box. Then, we define $\textbf{B}=\textbf{A}^{-1}$, so $f$ can be expressed as:
$$\textbf{f}=\textbf{B}\textbf{r}$$
To ensure periodic boundary conditions, we only need to ensure $0 \leq f_i \leq 1$, Therefore, we operate on $f_i$ i.e., $f_{i-wrap}=f_i-floor(f_i)$. Then, we obtain the new cell coordinates: 
$$r_{wrap}=\textbf{A}f_{wrap}$$
```
std::vector<size_t> neighbors;
Vec3<int> cell_vector = get_cell_vector(cell_index);
// define a offset matrix 26*3
std::vector<Vec3<int>> matrix = {
    {-1, -1, -1},{-1, -1, 0},{-1, -1, 1},{-1, 0, -1},{-1, 0, 0},{-1, 0, 1},{-1, 1, -1},{-1, 1, 0},{-1, 1, 1},
    {0, -1, -1},{0, -1, 0},{0, -1, 1},{0, 0, -1},{0, 0, 1},{0, 1, -1},{0, 1, 0},{0, 1, 1},
    {1, -1, -1},{1, -1, 0},{1, -1, 1},{1, 0, -1},{1, 0, 0},{1, 0, 1},{1, 1, -1},{1, 1, 0},{1, 1, 1} 
};
Mat3<double> A = {_cell_length[0], 0, 0, 0, _cell_length[1], 0, 0, 0, _cell_length[2]};
Mat3<double> inv_A = A.invert();
Vec3<double> f;
Vec3<double> wrapped_f;
Vec3<int> warp_neb_vector;
for (int i = 0; i < 26; i++){
    Vec3<int> neb_vector = cell_vector + matrix[i];
    f = inv_A * neb_vector;
    wrapped_f = f - floor(f);
    warp_neb_vector = A * wrapped_f;
    neighbors.push_back(get_cell_index(warp_neb_vector));
}
```

## Interaction computation algorithm
The entire computation process is as follows: the first layer iterates over all units, the second layer iterates over the particles within the unit, the third layer iterates over neighboring units, and finally calculates the distance. Subsequently, particles with a neighbor relationship are mutually added to each other's neighbor lists.
```
size_t n_cells = _cell_list.get_ncells();
for (size_t cell_index = 0; cell_index < n_cells; ++cell_index)
{
    for (auto neighbor_cell : _cell_list.get_neighbors(cell_index))
    {
        for (size_t i : _cell_list.get_atoms_in_cell(cell_index))
        {
            for (size_t j : _cell_list.get_atoms_in_cell(neighbor_cell))
            {
                if (i < j) // Avoid double counting
                {
                    Vec3<double> pos_i = xyz[i];
                    Vec3<double> pos_j = xyz[j];
                    double r = _box->calc_distance(pos_i, pos_j);
                    if (r < _r_cutoff)
                    {
                        _neighborListArray[i].push_back(j);
                        _neighborListArray[j].push_back(i);
                    }
                }
            }
        }
    }
}
```  
