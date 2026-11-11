#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include "libslic3r/ClipperUtils.hpp"
#include "libslic3r/ElephantFootCompensation.hpp"
#include "libslic3r/ExPolygon.hpp"
#include "libslic3r/Flow.hpp"
#include "libslic3r/Geometry.hpp"
#include "libslic3r/Point.hpp"
#include "libslic3r/Polygon.hpp"
#include <cmath>

using namespace Slic3r;

TEST_CASE("Shrinkage transform includes XY and Z compensation once", "[VoxelCompensation]")
{
    Geometry::Transformation transform;
    transform.set_rotation(Vec3d(0., 0., std::acos(-1.) / 2.));
    transform.set_scaling_factor(Vec3d(1.5, 0.75, 2.));
    transform.set_offset(Vec3d(10., 20., 1.));
    const Vec3d shrinkage(100. / 95., 100. / 95., 100. / 80.);
    const Vec3d point = transform.get_matrix_with_applied_shrinkage_compensation(shrinkage) * Vec3d(2., 4., 6.);
    CHECK(std::abs(point.x() - (10. - 3. * shrinkage.x())) < 1e-9);
    CHECK(std::abs(point.y() - (20. + 3. * shrinkage.y())) < 1e-9);
    CHECK(std::abs(point.z() - 13. * shrinkage.z()) < 1e-9);
}

TEST_CASE("Native elephant foot trimming preserves interior material interfaces", "[VoxelCompensation]")
{
    const ExPolygons solid{ExPolygon(Polygon({Point::new_scale(0., 0.), Point::new_scale(10., 0.),
        Point::new_scale(10., 10.), Point::new_scale(0., 10.)}))};
    const ExPolygons island{ExPolygon(Polygon({Point::new_scale(4., 4.), Point::new_scale(6., 4.),
        Point::new_scale(6., 6.), Point::new_scale(4., 6.)}))};
    const double compensation = GENERATE(0., 0.05, 0.1);
    const ExPolygons trimming = elephant_foot_compensation(solid, Flow(0.2f, 0.1f, 0.2f), compensation);
    const ExPolygons a = intersection_ex(diff_ex(solid, island), trimming);
    const ExPolygons b = intersection_ex(island, trimming);
    CHECK(diff_ex(island, b).empty());
    CHECK(diff_ex(b, island).empty());
    CHECK(intersection_ex(a, b).empty());
    CHECK(diff_ex(trimming, union_ex(a, b)).empty());
    CHECK(diff_ex(union_ex(a, b), trimming).empty());
    if (compensation > 0.)
        CHECK_FALSE(diff_ex(solid, trimming).empty());
    else
        CHECK(diff_ex(solid, trimming).empty());
}
