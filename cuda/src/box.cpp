#include "box.h"
#include "vec3.cpp"
#include <stdexcept>

namespace dpnblist
{

    // Sinus and Cosine for degree values
    constexpr float pi = 3.141592653589793238463;
    static float deg2rad(float x)
    {
        return x * pi / 180.0;
    }

    static float rad2deg(float x)
    {
        return x * 180.0 / pi;
    }

    static float dacos(float alpha)
    {
        return rad2deg(acos(alpha));
    }

    // static float dasin(float alpha)
    // {
    //     return rad2deg(asin(alpha));
    // }

    static float cosd(float theta)
    {
        return cos(deg2rad(theta));
    }

    // static float sind(float theta)
    // {
    //     return sin(deg2rad(theta));
    // }

    Box::Box() : _matrix(Mat3<float>::unit())
    {
        set_periodic(PBC::P, PBC::P, PBC::P);
    }

    Box::Box(std::vector<float> lengths, std::vector<float> angles){
        set_periodic(PBC::P, PBC::P, PBC::P);
        if (lengths.size() != 3 || angles.size() != 3){
            throw std::runtime_error("lengths and angles must be 3D vectors");
        }
        _lengths = {lengths[0], lengths[1], lengths[2]};
        _angles = {angles[0], angles[1], angles[2]};
        set_lengths_and_angles(_lengths, _angles);
    }

    Box::Box(Vec3<float> lengths, Vec3<float> angles) : Box()
    {
        _lengths = lengths;
        _angles = angles;
        set_lengths_and_angles(lengths, angles);
    }

    void Box::set_lengths_and_angles(Vec3<float> lengths, Vec3<float> angles)
    {
        auto a = lengths[0];
        auto b = lengths[1];
        auto c = lengths[2];

        auto alpha = angles[0];
        auto beta  = angles[1];
        auto gamma = angles[2];

        auto lx = a;
        auto xy = b * cosd(gamma);
        auto xz = c * cosd(beta);
        auto ly = sqrt(b * b - xy * xy);
        auto yz = (b * c * cosd(alpha) - xy * xz) / ly;
        auto lz = sqrt(c * c - xz * xz - yz * yz);

        _matrix[0][0] = lx;
        _matrix[0][1] = xy;
        _matrix[0][2] = xz;
        _matrix[1][0] = 0;
        _matrix[1][1] = ly;
        _matrix[1][2] = yz;
        _matrix[2][0] = 0;
        _matrix[2][1] = 0;
        _matrix[2][2] = lz;

        _inv_matrix = _matrix.invert();
    }

    const Mat3<float> Box::get_matrix() const
    {
        return _matrix;
    }

    const Mat3<float> Box::get_inverse() const
    {
        return _matrix.invert();
    }

    const Vec3<float> Box::get_angles() const
    {
        auto tilts = get_tilts();
        float xy = tilts[0];
        float xz = tilts[1];
        float yz = tilts[2];

        auto lengths = get_lengths();
        auto b = lengths[1];
        auto c = lengths[2];

        float ly = _matrix[1][1];

        auto cos_alpha = (xy * xz + ly * yz) / (b * c);
        auto cos_beta = xz / c;
        auto cos_gamma = xy / b;

        return {dacos(cos_alpha), dacos(cos_beta), dacos(cos_gamma)};    
    }

    void Box::set_periodic(PBC x, PBC y, PBC z)
    {
        _pbc = {x, y, z};
    }

    const std::array<Box::PBC, 3> Box::get_periodic() const
    {
        return _pbc;
    }

    const Vec3<float> Box::get_lengths() const
    {
        auto tilts = get_tilts();
        float xy = tilts[0];
        float xz = tilts[1];
        float yz = tilts[2];
        float ly = _matrix[1][1];
        float lz = _matrix[2][2];

        float a = _matrix[0][0];
        float b = sqrt(ly * ly + xy * xy);
        float c = sqrt(lz*lz + xz*xz + yz*yz);
        return {a, b, c};
    }

    const Vec3<float> Box::get_tilts() const
    {
        float xy = _matrix[0][1];
        float xz = _matrix[0][2];
        float yz = _matrix[1][2];
        return {xy, xz, yz};
    }

    const float Box::get_volume() const
    {
        return _matrix.determinant();
    }

    std::vector<Vec3<float>> Box::wrap(std::vector<Vec3<float>>& positions)
    {
        std::vector<Vec3<float>> wrapped_positions(positions.size());
        for (int i = 0; i < positions.size(); i++)
        {
            wrapped_positions[i] = wrap(positions[i]);
        }
        return wrapped_positions;
    }

    Vec3<float> Box::wrap(Vec3<float>& position)
    {
        // PBC are all P
        if (_pbc[0] == P && _pbc[1] == P && _pbc[2] == P)
        {
            // auto _inv_mat = get_inverse();
            Vec3<float> reci_vec = _inv_matrix * position;
            Vec3<float> wrapped_reci_vec = reci_vec - floor(reci_vec + Vec3<float>(0.000001, 0.000001, 0.000001));
            return _matrix * wrapped_reci_vec;
        }
        else
            throw std::runtime_error("Only PBC = P is implemented");
    }

    float Box::calc_distance(Vec3<float>& r1, Vec3<float>& r2)
    {
        Vec3<float> difference(0,0,0);
        float diff = 0;

        // Vec3<float> wrapped_dr;
        // for (int i = 0; i < 3; i++){
        //     wrapped_dr[i] = r1[i] - r2[i] + _lengths[i] / 2;
        //     diff = std::fmod(wrapped_dr[i], _lengths[i]);
        //     if (diff < 0) diff += _lengths[i];
        //     diff -= _lengths[i] / 2;
        //     difference[i] = diff;
        // }

        for (int i = 0; i < 3; i++){
            diff = r1[i] - r2[i] + 0.5 * _lengths[i];
            if (diff < 0) {
                difference[i] = diff + 0.5 * _lengths[i];
            }
            else if (diff > _lengths[i]) {
                difference[i] = diff - 1.5 * _lengths[i];
            }
            else {
                difference[i] = diff - 0.5 * _lengths[i];
            }
        }

        return norm(difference);
    }

    Vec3<float> Box::calc_vector(Vec3<float>& r1, Vec3<float>& r2)
    {
        auto dr = r1 - r2;
        auto wrapped_dr = wrap(dr);
        return wrapped_dr;
    }

}
