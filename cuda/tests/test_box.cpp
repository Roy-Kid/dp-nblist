#include "approx.cpp"
#include "box.h"
#include "mat3.cpp"
#include "vec3.cpp"
#include <doctest/doctest.h>

namespace dpnblist {
TEST_SUITE("TestBox") {
  TEST_CASE("test_init") {
    Box box1;
    CHECK_EQ(box1.get_lengths(), Vec3<float>(1, 1, 1));
    CHECK_EQ(box1.get_tilts(), Vec3<float>(0, 0, 0));
    CHECK_EQ(box1.get_angles(), Vec3<float>(90, 90, 90));
    CHECK_EQ(box1.get_volume(), 1);

    Box box2(Vec3<float>(1, 2, 3), Vec3<float>(90, 90, 90));
    CHECK_EQ(box2.get_lengths(), Vec3<float>(1, 2, 3));
    CHECK(approx_eq(box2.get_tilts(), Vec3<float>(0, 0, 0)));
    CHECK_EQ(box2.get_angles(), Vec3<float>(90, 90, 90));
    CHECK_EQ(box2.get_volume(), 6);

    Box box3(Vec3<float>(1, 1, 1), Vec3<float>(45, 45, 45));
    CHECK(approx_eq(box3.get_lengths(), Vec3<float>(1, 1, 1)));
    CHECK(approx_eq(box3.get_tilts(), Vec3<float>(0.707, 0.707, 0.293), 1e-03));
    CHECK(approx_eq(box3.get_angles(), Vec3<float>(45, 45, 45)));
    CHECK(approx_eq(box3.get_volume(), 0.45509, 1e-05));
  }

  TEST_CASE("test_orth_wrap") {

    Box box(Vec3<float>(2, 2, 2));
    std::vector<Vec3<float>> points = {{0.5, 0.5, 0.5},
                                        {1.5, 1.5, 1.5},
                                        {-0.5, -0.5, -0.5},
                                        {-1.5, -1.5, -1.5}};
    std::vector<Vec3<float>> expected = {
        {0.5, 0.5, 0.5}, {1.5, 1.5, 1.5}, {1.5, 1.5, 1.5}, {0.5, 0.5, 0.5}};
    auto wrapped_points = box.wrap(points);
    for (int i = 0; i < wrapped_points.size(); ++i) {
      CHECK(approx_eq(wrapped_points[i], expected[i]));
    }
    //CHECK(std::equal(wrapped_points.begin(), wrapped_points.end(),expected.begin()));
  }

  TEST_CASE("test_tric_wrap") {

    Box box(Vec3<float>(2, 2, 2), Vec3<float>(45, 45, 45));
    std::vector<Vec3<float>> points = {{2, 1, 0}};
    std::vector<Vec3<float>> expected = {{2, 1, 0}}; // TODO: need more tests
    auto wrapped_points = box.wrap(points);
    for (int i = 0; i < wrapped_points.size(); ++i) {
      CHECK(approx_eq(wrapped_points[i], expected[i]));
    }
    // CHECK(std::equal(wrapped_points.begin(), wrapped_points.end(),expected.begin()));
  }
}
} // namespace dpnblist