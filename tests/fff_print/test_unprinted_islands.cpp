#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include "libslic3r/BoundingBox.hpp"
#include "libslic3r/ClipperUtils.hpp"
#include "libslic3r/ExPolygon.hpp"
#include "libslic3r/ExtrusionEntityCollection.hpp"
#include "libslic3r/Flow.hpp"
#include "libslic3r/Fill/FillBase.hpp"
#include "libslic3r/Fill/FillRectilinear.hpp"
#include "libslic3r/Layer.hpp"
#include "libslic3r/PerimeterGenerator.hpp"
#include "libslic3r/Point.hpp"
#include "libslic3r/Polygon.hpp"
#include "libslic3r/Print.hpp"
#include "libslic3r/PrintConfig.hpp"
#include "libslic3r/Surface.hpp"
#include "libslic3r/SurfaceCollection.hpp"
#include "libslic3r/libslic3r.h"
#include <cstddef>
#include <string>

using namespace Slic3r;

namespace {
Polygon square(double x, double y, double side)
{
    return Polygon({Point::new_scale(x, y), Point::new_scale(x + side, y),
                    Point::new_scale(x + side, y + side), Point::new_scale(x, y + side)});
}

struct IslandLayer : Layer {
    IslandLayer() : Layer(10, nullptr, 0.1, 1.1, 1.05) {}
};

// Direct polygon inputs exercise the same perimeter generator used after a
// posSlice plugin replaces region.slices, without proprietary meshes or a GUI.
struct IslandFixture {
    PrintConfig print_config;
    PrintObjectConfig object_config;
    PrintRegionConfig region_config;
    PrintRegion material_a, material_b;
    IslandLayer layer;
    const ExPolygons solid{ExPolygon(square(0., 0., 10.))};
    LayerRegion *a, *b;
    double width = 0.2;
    double height = 0.1;

    explicit IslandFixture(double side)
    {
        print_config.nozzle_diameter.values = {width, width};
        print_config.resolution.value = 0.001;
        object_config.wall_generator.value = PerimeterGeneratorType::Arachne;
        object_config.min_feature_size.value = 25.;
        object_config.min_bead_width.value = 85.;
        object_config.min_length_factor.value = 0.5;
        object_config.interface_shells.value = false;
        region_config.wall_loops.value = 2;
        region_config.sparse_infill_density.value = 100.;
        region_config.only_one_wall_top.value = false;
        region_config.only_one_wall_first_layer.value = false;
        region_config.precise_outer_wall.value = false;
        region_config.fuzzy_skin.value = FuzzySkinType::None;
        region_config.outer_wall_filament_id.value = 1;
        material_a.set_config(region_config);
        region_config.outer_wall_filament_id.value = 2;
        material_b.set_config(region_config);
        a = layer.add_region(&material_a);
        b = layer.add_region(&material_b);
        b->slices.set(ExPolygons{ExPolygon(square(5., 5., side))}, stInternal);
        a->slices.set(diff_ex(solid, to_expolygons(b->slices.surfaces)), stInternal);
        layer.make_slices();
        layer.backup_untyped_slices();
    }

    void generate()
    {
        for (LayerRegion *region : layer.regions()) {
            region->perimeters.clear();
            region->thin_fills.clear();
            region->fill_surfaces.clear();
            region->fill_no_overlap_expolygons.clear();
            region->unprinted_islands.clear();
            const LayerRegionPtrs compatible{region};
            const Flow flow{float(width), float(height), float(width)};
            PerimeterGenerator g(&region->slices, &compatible, height, layer.slice_z, flow,
                &region->region().config(), &object_config, &print_config, false, 0.,
                &region->perimeters, &region->thin_fills, &region->fill_surfaces,
                &region->fill_no_overlap_expolygons);
            g.layer_id = int(layer.id());
            g.upper_slices = &solid;
            g.lower_slices = &solid;
            g.upper_slices_same_region = &region->slices;
            g.unprinted_islands = &region->unprinted_islands;
            if (object_config.wall_generator.value == PerimeterGeneratorType::Arachne)
                g.process_arachne();
            else
                g.process_classic();
        }
    }

    ExPolygons uncovered(const ExPolygons &area) const
    {
        Polygons coverage;
        for (const LayerRegion *region : layer.regions()) {
            region->perimeters.polygons_covered_by_width(coverage, 0.f);
            region->thin_fills.polygons_covered_by_width(coverage, 0.f);
            polygons_append(coverage, to_polygons(region->fill_surfaces.surfaces));
        }
        return diff_ex(area, coverage);
    }

    ExPolygons unextruded(const ExPolygons &area) const
    {
        Polygons coverage;
        for (const LayerRegion *region : layer.regions()) {
            region->perimeters.polygons_covered_by_width(coverage, 0.f);
            region->thin_fills.polygons_covered_by_width(coverage, 0.f);
            FillRectilinear fill;
            FillParams params;
            params.density = 1.f;
            params.flow = Flow{float(width), float(height), float(width)};
            params.config = &region->region().config();
            fill.spacing = params.flow.spacing();
            fill.angle = 0.f;
            fill.layer_id = layer.id();
            fill.z = layer.print_z;
            fill.bounding_box = get_extents(solid);
            for (const Surface &surface : region->fill_surfaces.surfaces)
                polygons_append(coverage, offset(fill.fill_surface(&surface, params), float(scale_(width / 2.))));
        }
        return diff_ex(area, coverage);
    }
};

bool same_area(const ExPolygons &a, const ExPolygons &b)
{
    return diff_ex(a, b).empty() && diff_ex(b, a).empty();
}
} // namespace

TEST_CASE("Wholly discarded internal material islands return to the solid surrounding region", "[UnprintedIslands]")
{
    const double side = GENERATE(0.04, 0.08, 0.10);
    IslandFixture f(side);
    f.object_config.wall_generator.value = GENERATE(PerimeterGeneratorType::Arachne, PerimeterGeneratorType::Classic);
    const ExPolygons island = f.b->raw_slices;
    const ExPolygons merged = f.layer.lslices;
    f.generate();
    REQUIRE(f.b->perimeters.empty());
    REQUIRE(f.b->fill_surfaces.empty());
    REQUIRE(f.b->unprinted_islands.size() == 1);
    REQUIRE_FALSE(f.uncovered(island).empty());
    REQUIRE_FALSE(f.unextruded(island).empty());
    REQUIRE(f.layer.reassign_unprinted_islands());
    CHECK(f.b->slices.empty());
    CHECK(same_area(to_expolygons(f.a->slices.surfaces), f.solid));
    CHECK(same_area(f.layer.lslices, merged));
    f.generate();
    CHECK(f.uncovered(island).empty());
    CHECK(f.unextruded(island).empty());
    CHECK_FALSE(f.layer.reassign_unprinted_islands());
}

TEST_CASE("Printable material islands retain their material assignment", "[UnprintedIslands]")
{
    IslandFixture f(GENERATE(0.12, 0.2, 1.));
    const ExPolygons before_a = f.a->raw_slices, before_b = f.b->raw_slices;
    f.generate();
    REQUIRE_FALSE(f.b->perimeters.empty());
    REQUIRE(f.b->unprinted_islands.empty());
    CHECK_FALSE(f.layer.reassign_unprinted_islands());
    CHECK(same_area(to_expolygons(f.a->slices.surfaces), before_a));
    CHECK(same_area(to_expolygons(f.b->slices.surfaces), before_b));
}

TEST_CASE("Island recovery preserves air gaps and external boundaries", "[UnprintedIslands]")
{
    IslandFixture f(0.08);
    f.object_config.wall_generator.value = GENERATE(PerimeterGeneratorType::Arachne, PerimeterGeneratorType::Classic);
    const int geometry = GENERATE(0, 1, 2, 3);
    if (geometry == 0) // A real air hole, without any material B.
        f.b->slices.clear();
    else if (geometry == 1) // Air surrounds B inside a larger hole.
        f.a->slices.set(diff_ex(f.solid, ExPolygons{ExPolygon(square(4.99, 4.99, 0.1))}), stInternal);
    else if (geometry == 2) { // A disconnected external B island.
        f.a->slices.set(f.solid, stInternal);
        f.b->slices.set(ExPolygons{ExPolygon(square(12., 12., 0.08))}, stInternal);
    } else { // B itself contains an air cavity.
        f.b->slices.set(diff_ex(to_expolygons(f.b->slices.surfaces),
            ExPolygons{ExPolygon(square(5.03, 5.03, 0.02))}), stInternal);
    }
    const ExPolygons before_a = to_expolygons(f.a->slices.surfaces);
    const ExPolygons before_b = to_expolygons(f.b->slices.surfaces);
    f.generate();
    CHECK_FALSE(f.layer.reassign_unprinted_islands());
    CHECK(same_area(to_expolygons(f.a->slices.surfaces), before_a));
    CHECK(same_area(to_expolygons(f.b->slices.surfaces), before_b));
}

TEST_CASE("Recovering one island preserves printable islands in the same material region", "[UnprintedIslands]")
{
    IslandFixture f(0.1);
    const ExPolygons printable{ExPolygon(square(7., 7., 1.))};
    f.b->slices.append(printable, stInternal);
    f.a->slices.set(diff_ex(f.solid, to_expolygons(f.b->slices.surfaces)), stInternal);
    f.generate();
    REQUIRE(f.b->unprinted_islands.size() == 1);
    REQUIRE_FALSE(f.b->perimeters.empty());
    REQUIRE(f.layer.reassign_unprinted_islands());
    CHECK(same_area(to_expolygons(f.b->slices.surfaces), printable));
    CHECK(same_area(to_expolygons(f.a->slices.surfaces), diff_ex(f.solid, printable)));
}

TEST_CASE("Changing print settings restores the original material ownership", "[UnprintedIslands]")
{
    IslandFixture f(0.1);
    const ExPolygons before_a = f.a->raw_slices, before_b = f.b->raw_slices;
    f.generate();
    REQUIRE(f.layer.reassign_unprinted_islands());
    // Surface-classification reruns need the recovered geometry.
    f.layer.restore_untyped_slices();
    CHECK(f.b->slices.empty());
    // Perimeter reruns need the original geometry, so a now-printable B can return.
    f.layer.restore_island_ownership();
    CHECK(same_area(to_expolygons(f.a->slices.surfaces), before_a));
    CHECK(same_area(to_expolygons(f.b->slices.surfaces), before_b));
    f.width = 0.05;
    f.height = 0.02;
    f.print_config.nozzle_diameter.values = {f.width, f.width};
    f.generate();
    CHECK_FALSE(f.b->perimeters.empty());
    CHECK(f.b->unprinted_islands.empty());
    CHECK_FALSE(f.layer.reassign_unprinted_islands());
}

TEST_CASE("Material boundary walls can precede the first visible island extrusion", "[UnprintedIslands]")
{
    IslandFixture before(0.1), tiny(0.1), visible(0.12);
    before.a->slices.set(before.solid, stInternal);
    before.b->slices.clear();
    before.generate();
    tiny.generate();
    visible.generate();
    // Every layer has the same solid union above and below: this requires no
    // cross-layer projection. B exists geometrically before B paths survive.
    CHECK(before.b->perimeters.empty());
    CHECK(tiny.b->perimeters.empty());
    CHECK_FALSE(visible.b->perimeters.empty());
    CHECK(tiny.a->perimeters.items_count() > before.a->perimeters.items_count());
    REQUIRE(tiny.layer.reassign_unprinted_islands());
    tiny.generate();
    CHECK(tiny.a->perimeters.items_count() == before.a->perimeters.items_count());
}

TEST_CASE("Island recovery does not override deliberately sparse surrounding material", "[UnprintedIslands]")
{
    IslandFixture f(0.1);
    PrintRegionConfig sparse = f.material_a.config();
    sparse.sparse_infill_density.value = 15.;
    f.material_a.set_config(sparse);
    f.generate();
    REQUIRE(f.b->unprinted_islands.size() == 1);
    CHECK_FALSE(f.layer.reassign_unprinted_islands());
    CHECK(same_area(to_expolygons(f.b->slices.surfaces), f.b->raw_slices));
}

TEST_CASE("A fresh slice replaces the material ownership saved by island recovery", "[UnprintedIslands]")
{
    IslandFixture f(0.1);
    f.generate();
    REQUIRE(f.layer.reassign_unprinted_islands());
    const ExPolygons replacement{ExPolygon(square(6., 6., 1.))};
    f.b->slices.set(replacement, stInternal);
    f.a->slices.set(diff_ex(f.solid, replacement), stInternal);
    f.layer.backup_untyped_slices();
    f.layer.restore_island_ownership();
    CHECK(same_area(to_expolygons(f.b->slices.surfaces), replacement));
    CHECK(same_area(to_expolygons(f.a->slices.surfaces), diff_ex(f.solid, replacement)));
}

TEST_CASE("Island recovery uses the enclosing material rather than the first region", "[UnprintedIslands]")
{
    IslandFixture f(0.04);
    f.object_config.wall_generator.value = GENERATE(PerimeterGeneratorType::Arachne, PerimeterGeneratorType::Classic);
    PrintRegion material_c;
    PrintRegionConfig config = f.material_a.config();
    config.outer_wall_filament_id.value = 3;
    material_c.set_config(config);
    f.print_config.nozzle_diameter.values = {f.width, f.width, f.width};
    LayerRegion *c = f.layer.add_region(&material_c);
    const ExPolygons enclosing{ExPolygon(square(2., 2., 6.))};
    const ExPolygons island = f.b->raw_slices;
    f.a->slices.set(diff_ex(f.solid, enclosing), stInternal);
    c->slices.set(diff_ex(enclosing, island), stInternal);
    const ExPolygons before_a = to_expolygons(f.a->slices.surfaces);
    f.generate();
    REQUIRE(f.b->unprinted_islands.size() == 1);
    REQUIRE(f.layer.reassign_unprinted_islands());
    CHECK(f.b->slices.empty());
    CHECK(same_area(to_expolygons(f.a->slices.surfaces), before_a));
    CHECK(same_area(to_expolygons(c->slices.surfaces), enclosing));
    CHECK(c->region().config().outer_wall_filament_id.value == 3);
    f.generate();
    CHECK(f.unextruded(island).empty());
}

TEST_CASE("Island recovery is disabled by default and survives config serialization", "[UnprintedIslands]")
{
    PrintObjectConfig config;
    CHECK_FALSE(config.fill_unprintable_islands.value);
    PrintObjectConfig stored;
    REQUIRE(stored.option("fill_unprintable_islands")->deserialize("1"));
    config.apply(stored);
    CHECK(config.fill_unprintable_islands.value);
    CHECK(config.fill_unprintable_islands.serialize() == std::string("1"));
    REQUIRE(stored.option("fill_unprintable_islands")->deserialize("0"));
    config.apply(stored);
    CHECK_FALSE(config.fill_unprintable_islands.value);
}

TEST_CASE("Classic recovers wholly discarded islands with thin wall detection on or off", "[UnprintedIslands]")
{
    IslandFixture f(0.04);
    f.object_config.wall_generator.value = PerimeterGeneratorType::Classic;
    PrintRegionConfig config = f.material_b.config();
    config.detect_thin_wall.value = GENERATE(false, true);
    config.filter_out_gap_fill.value = 0.;
    f.material_b.set_config(config);
    const ExPolygons island = f.b->raw_slices;
    f.generate();
    REQUIRE(f.b->perimeters.empty());
    REQUIRE(f.b->thin_fills.empty());
    REQUIRE(f.b->fill_surfaces.empty());
    REQUIRE(f.b->unprinted_islands.size() == 1);
    REQUIRE(f.layer.reassign_unprinted_islands());
    f.generate();
    CHECK(f.unextruded(island).empty());
}

TEST_CASE("Classic preserves printable thin material islands", "[UnprintedIslands]")
{
    IslandFixture f(1.);
    f.object_config.wall_generator.value = PerimeterGeneratorType::Classic;
    PrintRegionConfig config = f.material_b.config();
    config.detect_thin_wall.value = true;
    config.filter_out_gap_fill.value = 0.;
    f.material_b.set_config(config);
    const double width = GENERATE(0.1, 0.2, 1.);
    const ExPolygons island{ExPolygon(Polygon({Point::new_scale(5., 5.), Point::new_scale(5. + width, 5.),
        Point::new_scale(5. + width, 6.), Point::new_scale(5., 6.)}))};
    f.b->slices.set(island, stInternal);
    f.a->slices.set(diff_ex(f.solid, island), stInternal);
    f.generate();
    REQUIRE_FALSE(f.b->perimeters.empty());
    CHECK(f.b->unprinted_islands.empty());
    CHECK_FALSE(f.layer.reassign_unprinted_islands());
    CHECK(same_area(to_expolygons(f.b->slices.surfaces), island));
}
