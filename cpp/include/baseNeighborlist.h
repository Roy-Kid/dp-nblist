#include "box.h"

namespace dpnblist
{

class BaseNeighborList
{
  public:
    BaseNeighborList(Box *box) : _box(box)
    {
    }
    ~BaseNeighborList()
    {
        // delete _box;
    }

    void set_lengths_and_angles(Vec3<double> lengths, Vec3<double> angles)
    {
        _box->set_lengths_and_angles(lengths, angles);
    }

    virtual void build(std::vector<Vec3<double>> &xyz) = 0;

    virtual void update(std::vector<Vec3<double>> &xyz) = 0;

  protected:
    Box *_box;
};

} // namespace dpnblist
