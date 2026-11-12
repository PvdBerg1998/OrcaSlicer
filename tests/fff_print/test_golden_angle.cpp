#include <catch2/catch_test_macros.hpp>
#include "libslic3r/Fill/FillRectilinear.hpp"
#include "libslic3r/Surface.hpp"
#include "libslic3r/Polygon.hpp"
#include "libslic3r/ExPolygon.hpp"
#include "libslic3r/Point.hpp"
#include "libslic3r/Polyline.hpp"
#include "libslic3r/PrintConfig.hpp"
#include "libslic3r/BoundingBox.hpp"
#include "libslic3r/ClipperUtils.hpp"
#include <cmath>
#include <algorithm>
#include <limits>
#include <string>

using namespace Slic3r;

TEST_CASE("Golden angle solid lines rotate without depending on call order", "[GoldenAngle]")
{
    Surface surface(stInternalSolid, ExPolygon(Polygon({Point::new_scale(0., 0.),
        Point::new_scale(20., 0.), Point::new_scale(20., 20.), Point::new_scale(0., 20.)})));
    FillGoldenAngle fill;
    fill.bounding_box = BoundingBox(Point::new_scale(0., 0.), Point::new_scale(20., 20.));
    fill.spacing = 0.4;
    fill.overlap = 0.;
    fill.angle = 0.;
    FillParams params;
    params.density = 1.;
    params.dont_adjust = false;
    params.anchor_length_max = 0.;
    const double pi = std::acos(-1.);
    for (size_t layer : {0u, 1u, 5u, 29u}) {
        fill.layer_id = layer;
        const auto paths = fill.fill_surface(&surface, params);
        REQUIRE_FALSE(paths.empty());
        const double expected = layer * pi * (3. - std::sqrt(5.)) / 2.;
        size_t long_segments = 0;
        for (const auto &path : paths)
            for (size_t i = 1; i < path.points.size(); ++i) {
                const Vec2d delta = (path.points[i] - path.points[i - 1]).cast<double>();
                if (delta.norm() > scale_(5.)) {
                    // Rectilinear sweeps start vertical, then rotate by angle + pi/2.
                    CHECK(std::abs(std::sin(std::atan2(delta.y(), delta.x()) - expected)) < 1e-5);
                    ++long_segments;
                }
            }
        CHECK(long_segments > 10);
        fill.layer_id = layer + 100;
        fill.fill_surface(&surface, params);
        fill.layer_id = layer;
        const auto repeated = fill.fill_surface(&surface, params);
        REQUIRE(repeated.size() == paths.size());
        for (size_t i = 0; i < paths.size(); ++i)
            CHECK(paths[i].points == repeated[i].points);
        CHECK(std::abs(fill.spacing - 0.4) < 1e-12);
    }
}

TEST_CASE("Golden angle shifts the grid even when solid spacing adjustment is requested", "[GoldenAngle]")
{
    Surface surface(stInternalSolid, ExPolygon(Polygon({Point::new_scale(0., 0.),
        Point::new_scale(20., 0.), Point::new_scale(20., 20.), Point::new_scale(0., 20.)})));
    FillGoldenAngle fill;
    fill.bounding_box = BoundingBox(Point::new_scale(0., 0.), Point::new_scale(20., 20.));
    fill.spacing = 0.4;
    fill.angle = 0.;
    fill.layer_id = 0;
    FillParams params;
    params.density = 1.;
    params.dont_adjust = false;
    params.anchor_length_max = 0.;
    const auto a = fill.fill_surface(&surface, params);
    fill.layer_id = 1;
    // Cancel the rotation to isolate the translation of the grid.
    fill.angle = float(-std::acos(-1.) * (3. - std::sqrt(5.)) / 2.);
    const auto b = fill.fill_surface(&surface, params);
    REQUIRE_FALSE(a.empty());
    REQUIRE_FALSE(b.empty());
    coord_t ay = std::numeric_limits<coord_t>::max(), by = ay;
    for (const auto &path : a)
        for (const auto &p : path.points) ay = std::min(ay, p.y());
    for (const auto &path : b)
        for (const auto &p : path.points) by = std::min(by, p.y());
    CHECK(std::abs(double(ay) - double(by)) > scale_(0.01));
    CHECK(std::abs(double(ay) - double(by)) < scale_(fill.spacing));
}

TEST_CASE("Golden angle is offered only for internal solid infill", "[GoldenAngle]")
{
    CHECK(print_config_def.get("internal_solid_infill_pattern")->has_enum_value("goldenangle"));
    CHECK_FALSE(print_config_def.get("top_surface_pattern")->has_enum_value("goldenangle"));
    CHECK_FALSE(print_config_def.get("bottom_surface_pattern")->has_enum_value("goldenangle"));
    CHECK_FALSE(print_config_def.get("sparse_infill_pattern")->has_enum_value("goldenangle"));
    ConfigOptionEnum<InfillPattern> option(ipGoldenAngle);
    CHECK(option.serialize() == std::string("goldenangle"));
    ConfigOptionEnum<InfillPattern> loaded(ipRectilinear);
    REQUIRE(loaded.deserialize("goldenangle"));
    CHECK(loaded.value == ipGoldenAngle);
}

TEST_CASE("Golden angle lines remain clipped around air cavities", "[GoldenAngle]")
{
    const Polygon outer({Point::new_scale(0., 0.), Point::new_scale(20., 0.),
        Point::new_scale(20., 20.), Point::new_scale(0., 20.)});
    const Polygon hole({Point::new_scale(6., 6.), Point::new_scale(14., 6.),
        Point::new_scale(14., 14.), Point::new_scale(6., 14.)});
    ExPolygon area(outer);
    area.holes.push_back(hole);
    area.holes.back().reverse();
    Surface surface(stInternalSolid, area);
    FillGoldenAngle fill;
    fill.bounding_box = BoundingBox(Point::new_scale(0., 0.), Point::new_scale(20., 20.));
    fill.spacing = 0.4;
    fill.angle = 0.;
    FillParams params;
    params.density = 1.;
    for (size_t layer = 0; layer < 8; ++layer) {
        fill.layer_id = layer;
        const auto paths = fill.fill_surface(&surface, params);
        REQUIRE_FALSE(paths.empty());
        CHECK(intersection_pl(paths, hole).empty());
        CHECK(diff_pl(paths, Polygons{outer}).empty());
    }
}

TEST_CASE("Golden angle bridge output matches ordinary rectilinear", "[GoldenAngle]")
{
    Surface surface(stInternalBridge, ExPolygon(Polygon({Point::new_scale(0., 0.),
        Point::new_scale(10., 0.), Point::new_scale(10., 10.), Point::new_scale(0., 10.)})));
    surface.bridge_angle = 0.3;
    FillGoldenAngle golden;
    golden.bounding_box = BoundingBox(Point::new_scale(0., 0.), Point::new_scale(10., 10.));
    golden.spacing = 0.4;
    golden.angle = 0.;
    golden.layer_id = 8;
    FillRectilinear ordinary = golden;
    FillParams params;
    params.density = 1.;
    const auto a = golden.fill_surface(&surface, params);
    const auto b = ordinary.fill_surface(&surface, params);
    REQUIRE_FALSE(a.empty());
    REQUIRE(a.size() == b.size());
    for (size_t i = 0; i < a.size(); ++i)
        CHECK(a[i].points == b[i].points);
}
