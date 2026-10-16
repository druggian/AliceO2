// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file   FDPrototype.h
/// \brief  Geometry of the FD prototype (one FV0-like Ring-1 cell read out by bent fibers and one PMT in a black box),
///         installed at the end of data taking ~17 m from the IP, in front of FDD-A. Simulated as part of FV0.

#ifndef ALICEO2_FV0_FDPROTOTYPE_H_
#define ALICEO2_FV0_FDPROTOTYPE_H_

class TGeoVolume;

namespace o2
{
namespace fv0
{
namespace fdprototype
{
/// Name of the sensitive (scintillator) volume
constexpr const char* sensitiveVolumeName = "FDPROTO_SCINT";

/// Builds the prototype and returns it as an assembly in the "box frame":
/// origin at the box centre, X = box length (43 cm), Y = box height (16.1 cm, scintillator at -Y, PMT at +Y), Z = box width (9 cm).
/// Uses the FV0 media created in o2::fv0::Detector::createMaterials().
TGeoVolume* build();
} // namespace fdprototype
} // namespace fv0
} // namespace o2

#endif
