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

/// \file   FDPrototypeParam.h
/// \brief  Configurable parameters of the FD prototype placement, e.g.
///         o2-sim ... --configKeyValues "FV0FDProto.enable=true;FV0FDProto.z=1681.95"

#ifndef ALICEO2_FV0_FDPROTOTYPEPARAM_H_
#define ALICEO2_FV0_FDPROTOTYPEPARAM_H_

#include "CommonUtils/ConfigurableParam.h"
#include "CommonUtils/ConfigurableParamHelper.h"

namespace o2
{
namespace fv0
{
struct FDPrototypeParam : public o2::conf::ConfigurableParamHelper<FDPrototypeParam> {
  bool enable = false; ///< build the FD prototype (only present at the end of data taking)
  // Position of the box centre in the global ALICE frame (cm). The box is placed in caveRB24, like FDD-A.
  // From the installation drawing: FDD-A face towards the IP at z = 1692.0 cm, 2.0 cm gap, so the box face towards
  // FDD-A is at z = 1690.0 cm and the box centre at 1690.0 - 16.1 / 2 = 1681.95 cm (the drawing gives a 17.0 cm deep box).
  // The box lies in the transverse plane, its long side along a radius 45 deg above the +x axis, with its inner end
  // 5.2 cm from the beam axis: box centre at r = 5.2 + 43 / 2 = 26.7 cm, i.e. x = y = 26.7 / sqrt(2) = 18.88 cm.
  float x = 18.88f;
  float y = 18.88f;
  float z = 1681.95f;
  // Rotation of the box frame (deg), applied in the order X, Y, Z (around the global axes).
  // rotX = -90 maps box +Y onto -Z: the PMT side faces the IP and the scintillator is on the FDD-A side, as installed
  // (rotX = 90 would put the scintillator on the IP side).
  // rotZ then turns the box around the beam axis: box +X (the scintillator end) points along rotZ, i.e. 225 deg = towards
  // the beam from the upper +x side.
  float rotX = -90.f;
  float rotY = 0.f;
  float rotZ = 225.f;
  int cellId = 48; ///< cell ID written in the FD hits (FV0 cells are 0..47)

  O2ParamDef(FDPrototypeParam, "FV0FDProto");
};
} // namespace fv0
} // namespace o2

#endif
