/*
 * Copyright (c) 2014-2024 Key4hep-Project.
 *
 * This file is part of Key4hep.
 * See https://key4hep.github.io/key4hep-doc/ for further info.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
#include "catch2/catch_test_macros.hpp"
#include "catch2/matchers/catch_matchers_floating_point.hpp"

#include "../../k4ActsTracking/src/components/GeantinoMaterialAssigner.h"
#include <k4ActsTracking/Helpers.hxx>

#include <Acts/Geometry/Blueprint.hpp>
#include <Acts/Geometry/BlueprintOptions.hpp>
#include <Acts/Geometry/CuboidVolumeBounds.hpp>
#include <Acts/Geometry/StaticBlueprintNode.hpp>
#include <Acts/Geometry/TrackingGeometry.hpp>
#include <Acts/Geometry/TrackingVolume.hpp>
#include <Acts/MagneticField/ConstantBField.hpp>
#include <Acts/Material/HomogeneousSurfaceMaterial.hpp>
#include <Acts/Surfaces/PlaneSurface.hpp>
#include <Acts/Surfaces/RectangleBounds.hpp>

#include <memory>
#include <optional>
#include <vector>

namespace {
struct PropagationFixture {
  Acts::GeometryContext gctx = Acts::GeometryContext::dangerouslyDefaultConstruct();
  Acts::MagneticFieldContext mctx;
  std::shared_ptr<Acts::PlaneSurface> face = Acts::Surface::makeShared<Acts::PlaneSurface>(
      Acts::Transform3(Acts::Translation3(0., 0., 50.)), std::make_shared<Acts::RectangleBounds>(90., 90.));
  std::shared_ptr<const Acts::TrackingGeometry> geometry;

  PropagationFixture() {
    face->assignSurfaceMaterial(std::make_shared<Acts::HomogeneousSurfaceMaterial>(
        Acts::MaterialSlab(Acts::Material::fromMolarDensity(100., 200., 28., 14., 0.1), 1.)));
    auto volume = std::make_unique<Acts::TrackingVolume>(
        Acts::Transform3::Identity(), std::make_shared<Acts::CuboidVolumeBounds>(100., 100., 100.), "box");
    volume->addSurface(face);
    Acts::Blueprint::Config config;
    config.envelope = Acts::ExtentEnvelope{{.x = {10., 10.}, .y = {10., 10.}, .z = {10., 10.}}};
    Acts::Blueprint root{config};
    root.addChild(std::make_shared<Acts::StaticBlueprintNode>(std::move(volume)));
    geometry = root.construct(Acts::BlueprintOptions{}, gctx, *Acts::getDefaultLogger("test", Acts::Logging::WARNING));
  }

  ACTSTracking::CaloFacePropagator propagator() const {
    Acts::Navigator::Config cfg{geometry};
    cfg.resolvePassive = true;
    cfg.resolveMaterial = true;
    cfg.resolveSensitive = true;
    return ACTSTracking::CaloFacePropagator(
        Acts::EigenStepper<>(std::make_shared<Acts::ConstantBField>(Acts::Vector3::Zero())), Acts::Navigator(cfg));
  }

  Acts::BoundTrackParameters start() const {
    return Acts::BoundTrackParameters::createCurvilinear(Acts::Vector4::Zero(), Acts::Vector3::UnitZ(), 1.,
                                                         std::nullopt, Acts::ParticleHypothesis::geantino());
  }
};
} // namespace

TEST_CASE("Calorimeter extrapolation stops at the requested face or world boundary", "[propagation]") {
  PropagationFixture fixture;
  auto propagator = fixture.propagator();
  std::vector<Acts::GeometryIdentifier> ids{fixture.face->geometryId()};
  auto expected = ACTSTracking::CaloExtrapolationStatus::Ok;
  std::size_t maxSteps = 100;

  SECTION("Reach the face") {}
  SECTION("An absent face stops at the world boundary") {
    ids = {Acts::GeometryIdentifier{}.withVolume(255)};
    expected = ACTSTracking::CaloExtrapolationStatus::NotReached;
  }
  SECTION("No configured faces") {
    ids.clear();
    expected = ACTSTracking::CaloExtrapolationStatus::NoSurfaces;
  }
  SECTION("Step limit reports a propagation failure") {
    maxSteps = 0;
    expected = ACTSTracking::CaloExtrapolationStatus::PropagationError;
  }

  const auto result =
      ACTSTracking::extrapolateToCaloFace(propagator, fixture.start(), ids, fixture.gctx, fixture.mctx, maxSteps);
  REQUIRE(result.status == expected);
  REQUIRE(result.params.has_value() == (expected == ACTSTracking::CaloExtrapolationStatus::Ok));
  if (result.params) {
    CHECK(result.params->referenceSurface().geometryId() == fixture.face->geometryId());
    CHECK_THAT(result.params->position(fixture.gctx).z(), Catch::Matchers::WithinAbs(50., 1e-6));
  }
}

TEST_CASE("Material assignment reads actor state and counts failed propagations", "[propagation][material]") {
  PropagationFixture fixture;
  unsigned int maxSteps = 100;
  bool fail = false;
  SECTION("Collect material through the world boundary") {}
  SECTION("Discard incomplete assignments on failure") {
    maxSteps = 0;
    fail = true;
  }
  ACTSTracking::GeantinoMaterialAssigner assigner(fixture.propagator(), maxSteps,
                                                  Acts::getDefaultLogger("test", Acts::Logging::WARNING));
  const auto [surfaces, volumes] =
      assigner.assignmentCandidates(fixture.gctx, fixture.mctx, Acts::Vector3::Zero(), Acts::Vector3::UnitZ());
  INFO(assigner.stats()->lastError);
  CHECK(assigner.stats()->nTracks == 1);
  CHECK(assigner.stats()->nFailed == (fail ? 1 : 0));
  if (fail) {
    CHECK(surfaces.empty());
    CHECK(volumes.empty());
    CHECK_FALSE(assigner.stats()->lastError.empty());
  } else {
    REQUIRE(surfaces.size() == 1);
    CHECK(surfaces.front().surface == fixture.face.get());
    CHECK(assigner.stats()->maxStepsObserved > 0);
  }
}
