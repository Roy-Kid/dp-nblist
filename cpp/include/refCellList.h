#pragma once

#include <vector>
#include <array>
#include <cmath>

#include "vec3.hpp"
#include "box.h"
#include "baseNBL.h"

namespace dpnblist
{

    class CellList
    {
        
    public:

        CellList(Box *box, double cell_width);

        size_t get_cell_index(Vec3<int>&) const;

        Vec3<int> get_cell_vector(size_t) const;

        void build(std::vector<Vec3<double>> &xyz);

        void update(std::vector<Vec3<double>> &xyz);

        void reset();

        size_t get_ncells();

        std::vector<size_t> get_atoms_in_cell(size_t cell_index) const;

        std::vector<size_t> get_neighbors(size_t cell_index) const;

    private:
        Box *_box;
        Vec3<double> _cell_length;
        size_t _natoms;
        double _r_cutoff;
        std::vector<size_t> _head;
        std::vector<size_t> _lscl;
        const size_t EMPTY = std::numeric_limits<size_t>::max();
    };

    class NeighborList : public BaseNBL
    {
        public:
            using NeighborListArray = std::vector<std::vector<size_t>>;

            NeighborList(Box* box, double r_cutoff);

            ~NeighborList();

            void reset();

            void build(std::vector<Vec3<double>> &xyz);

            void update(std::vector<Vec3<double>> &xyz);

            NeighborListArray get_listArray();

        private:
            double _r_cutoff;
            CellList _cell_list;
            NeighborListArray _neighborListArray;
    };

}
