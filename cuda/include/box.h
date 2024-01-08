#pragma once

#include <vector>
#include <array>
#include "vec3.cpp"
#include "mat3.cpp"


namespace dpnblist
{
    class Box
    {

        public:

            enum PBC {
                F = 0,
                P = 1,
            };

            Box();
            Box(std::vector<float> lengths, std::vector<float> angles = {90, 90, 90});
            Box(Vec3<float> lengths, Vec3<float> angles = {90, 90, 90});
            void set_periodic(PBC, PBC, PBC);
            const std::array<Box::PBC, 3> get_periodic() const;
            void set_lengths_and_angles(Vec3<float> lengths, Vec3<float> angles);
            const Mat3<float> get_matrix() const;
            const Mat3<float> get_inverse() const;
            const Vec3<float> get_lengths() const;
            const Vec3<float> get_angles() const;
            const Vec3<float> get_tilts() const;
            const float get_volume() const;
            std::vector<Vec3<float>> wrap(std::vector<Vec3<float>>& points);
            Vec3<float> wrap(Vec3<float>& point);
            float calc_distance(Vec3<float>& point1, Vec3<float>& point2);
            Vec3<float> calc_vector(Vec3<float>&, Vec3<float>&);

        private:
            Vec3<float> _lengths;
            Vec3<float> _angles;
            Mat3<float> _matrix;
            Mat3<float> _inv_matrix;
            std::array<PBC, 3> _pbc;

    };

}
