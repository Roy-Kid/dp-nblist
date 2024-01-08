#include <vector>
#include "box.h"

namespace dpnblist
{
    // base class of neighborlist
    class base_NBL{
    public:
        virtual void build(std::vector<std::vector<float>> &xyz) = 0;
        virtual void update(std::vector<std::vector<float>> &xyz) = 0;
        virtual void out() = 0;
        virtual size_t gettime() = 0;
        virtual std::vector<std::vector<size_t>> get_listArray() = 0;
        ~base_NBL(){}
    };
    // cell list neighborlist
    class CudaCellList: public base_NBL
    {

    public:
        //CudaCellList(std::vector<float> &box_length, float r_cutoff, float skin);
        CudaCellList(Box *box, float r_cutoff, float skin);
        void build(std::vector<std::vector<float>> &xyz);
        void update(std::vector<std::vector<float>> &xyz);
        void out();
        size_t gettime();
        std::vector<std::vector<size_t>> get_listArray();
        ~CudaCellList();

    private:
        /* data */
        Box *_box;
        float _r_cutoff;
        float *d_r_cutoff;
        float _skin;
        float *d_skin;
        float _box_len[3];
        float *d_box_len;
        int _cell_len[3];
        int *d_cell_len;
        size_t _natoms;
        size_t *d_natoms;
        size_t _ncells;
        size_t *d_ncells;

        // std::vector<std::vector<int>> off_set_vec;
        // int *d_off_set_vec_1d;
        int *d_nebcell_list;

        int *d_head, *d_lscl, *d_atom_cellindex, *d_cell_atoms_count;
        int *_cell_atoms_count;
        // int *d_neighborListArray;
        int *_neighborListArray;
        size_t time_cost;
    };

    // create a neighborlist class as a factory mode to create different neighborlist
    class NeighborList
    {
    public:
        static base_NBL* createNeighborList(std::string type, Box* box, float r_cutoff, float skin);
        ~NeighborList();
    };

} // namespace dpnblist
