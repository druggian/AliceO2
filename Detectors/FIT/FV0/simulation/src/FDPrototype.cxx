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

/// \file   FDPrototype.cxx
/// \brief  Geometry of the FD prototype, ported from the standalone macro buildFDGeometry.C.
///
/// One FV0-like scintillator cell (Ring 1, sector 45-90 deg) wrapped in Teflon, glued to a plastic support, read out by
/// 120 bent 1 mm fibers that end on a single Hamamatsu H11934-100 PMT, everything inside a closed black plastic box.
///
/// Units: cm. Frames:
///   * "cell frame": frame of the FV0 cell shapes, z = thickness axis, sector between 45 and 90 deg in x-y.
///   * "box frame" : origin at the box centre, X = box length, Y = box height, Z = box width.
///
/// LAYOUT (box frame, seen from the front):
///        +Y
///         ^   [PMT] <==== horizontal fibers (narrowing bundle) ====.
///         |   (top left, window faces +X)                          |   90 deg bend (bendRadius)
///         |                                          [SCINTILLATOR]   fibers leave the support vertically (+Y)
///         +-----------------------------------------------------------> +X
///
/// All numbers marked "TO BE VERIFIED" are placeholders: check them against datasheets / CAD.

// \author: Daniela Ruggiano, WUT, poland
//\email: daniela.ruggiano@cern.ch

#include "FV0Simulation/FDPrototype.h"

#include <TGeoManager.h>
#include <TGeoMedium.h>
#include <TGeoVolume.h>
#include <TGeoMatrix.h>
#include <TGeoBBox.h>
#include <TGeoTube.h>
#include <TGeoCompositeShape.h>
#include <TMath.h>
#include <algorithm>
#include <cmath>
#include <map>
#include <string>
#include <utility>
#include <vector>
#include "Framework/Logger.h"

namespace
{
constexpr double eps = 0.005;                                       // numerical tolerance (cm)
constexpr double dzScint = 4.0;                                     // scintillator thickness (cm)
constexpr double dzPlastic = 1.0;                                   // plastic support thickness (cm)
constexpr double dzGlue = 0.1;                                      // optical glue thickness (cm)
constexpr double teflonThickness = 0.18;                            // Teflon wrap thickness (cm)
constexpr double ringR[6] = {4.01, 7.3, 12.9, 21.25, 38.7, 72.115}; // FV0 ring radii (cm)
constexpr double drSep = 0.07;                                      // separation between rings (cm)
constexpr double xShiftInner = -0.15;                               // inner shift along X (cm)
constexpr double dxHoleExt = 0.2;                                   // external hole shift along X (cm)
constexpr double drHoleSmall = 0.265;                               // small hole radius (cm)
constexpr double drHoleLarge = 0.415;                               // large hole radius (cm)

// Outer edge of the Ring 1 sector: straight (chord between the 45 and 90 deg corners) instead of the FV0 arc.
constexpr bool straightOuterEdge = true;
constexpr double chordNormalDeg = 67.5; // direction of the chord normal (middle of the 45-90 sector)
constexpr double outerEdgeOffset = 0.0; // cm; >0 moves the straight edge outwards (0 = chord, ~0.55 = tangent to the arc)

// Fibers: they leave the support perpendicularly (along the cell z axis).
constexpr double fiberDiameter = 0.10;   // fiber diameter (cm)
constexpr double fiberPitch = 0.30;      // maximum pitch (cm); reduced if 120 sites do not fit
constexpr double fiberPitchMin = 0.20;   // minimum allowed pitch (cm)
constexpr int nFibers = 120;             // total number of fibers
constexpr double dzFiber = 35.0;         // fiber length along the bundle axis (cm)
constexpr double fiberEdgeMargin = 0.20; // minimum distance of the fiber centres from the support edges / holes (cm)

// PMT Hamamatsu H11934-100. DIMENSIONS TO BE VERIFIED on the datasheet.
constexpr double pmtSize = 3.0;        // side of the square module (cm)
constexpr double pmtBodyDz = 3.0;      // body length behind the window (cm)
constexpr double pmtWindowR = 0.9;     // glass window radius (cm)
constexpr double pmtWindowDz = 0.1;    // window thickness (cm)
constexpr double pmtActiveR = 0.8;     // photocathode radius (cm)
constexpr double pmtPcDz = 0.005;      // photocathode thickness (cm)
constexpr double fiberStartGap = 0.01; // gap between fibers and support / PMT window (cm)

// Black light-tight plastic box, OUTER dimensions. Wall thickness TO BE VERIFIED.
constexpr double boxLength = 43.0; // along X (cm)
constexpr double boxHeight = 16.1; // along Y (cm)
constexpr double boxWidth = 9.0;   // along Z (cm)
constexpr double boxWall = 0.643;  // wall thickness (cm)
constexpr double boxGap = 0.2;     // minimum clearance between components and inner walls (cm)
// "U" = coordinate along the sector symmetry axis (67.5 deg in the cell frame), from the ring centre.
constexpr double stackCenterU = 5.2; // approximate centre of the sector along U (cm)
constexpr double stackHalfU = 1.65;  // half-extent of the sector (with Teflon) along U (cm)
constexpr double bendRadius = 6.0;   // bend radius of the bundle axis (cm), TO BE VERIFIED vs fiber minimum bending radius
constexpr int nBendSegments = 24;    // each fiber bend is approximated by this many straight tubes

TGeoMedium* medium(const char* name)
{
  auto* med = gGeoManager->GetMedium(Form("FV0_%s", name));
  if (!med) {
    LOG(fatal) << "FD prototype: medium FV0_" << name << " not found";
  }
  return med;
}

// Fiber centres in the cell-frame xy plane on a grid of the given pitch, kept only inside the support
// (at least fiberEdgeMargin from the edges, the 45/90 deg cuts, the inner cut and the screw holes).
std::vector<std::pair<double, double>> gridSites(double pitch)
{
  const double rMin = ringR[0], rMax = ringR[1];
  const double rMinS = rMin + drSep, rMaxS = rMax - drSep;
  const double xHole = drSep + dxHoleExt;
  const double s45 = TMath::Sin(TMath::PiOver4()), c45 = TMath::Cos(TMath::PiOver4());
  const double holes[2][2] = {{xHole, TMath::Cos(TMath::ASin(xHole / rMax)) * rMax}, {s45 * rMax, c45 * rMax}};
  const double m = fiberEdgeMargin;

  std::vector<std::pair<double, double>> pts;
  const int n = static_cast<int>(rMax / pitch) + 1;
  for (int i = -n; i <= n; ++i) {
    for (int j = -n; j <= n; ++j) {
      const double x = i * pitch, y = j * pitch;
      if (TMath::Hypot(x - xShiftInner, y) < rMinS + m || TMath::Hypot(x, y) > rMaxS - m) {
        continue;
      }
      if (x < drSep + m) { // 90 deg edge
        continue;
      }
      if ((y - x) / TMath::Sqrt2() < drSep + m) { // 45 deg edge
        continue;
      }
      if (straightOuterEdge && x * TMath::Cos(chordNormalDeg * TMath::DegToRad()) + y * TMath::Sin(chordNormalDeg * TMath::DegToRad()) >
                                 rMaxS * TMath::Cos(22.5 * TMath::DegToRad()) + outerEdgeOffset - m) { // straight outer edge
        continue;
      }
      bool nearHole = false;
      for (const auto& h : holes) {
        if (TMath::Hypot(x - h[0], y - h[1]) < drHoleSmall + fiberDiameter / 2 + 0.1) {
          nearHole = true;
        }
      }
      if (!nearHole) {
        pts.emplace_back(x, y);
      }
    }
  }
  return pts;
}

// Reduces the pitch in 0.1 mm steps until there are at least nFibers sites, then keeps the nFibers closest to the centroid.
std::vector<std::pair<double, double>> fiberPositions(double& pitchUsed)
{
  std::vector<std::pair<double, double>> pts;
  for (double pitch = fiberPitch; pitch >= fiberPitchMin - 1e-9; pitch -= 0.01) {
    pts = gridSites(pitch);
    pitchUsed = pitch;
    if (static_cast<int>(pts.size()) >= nFibers) {
      break;
    }
  }
  if (static_cast<int>(pts.size()) < nFibers) {
    LOG(warning) << "FD prototype: only " << pts.size() << " fiber sites for " << nFibers << " fibers with pitch " << 10 * pitchUsed << " mm";
    return pts;
  }
  if (pitchUsed < fiberPitch - 1e-9) {
    LOG(info) << "FD prototype: fiber pitch reduced to " << 10 * pitchUsed << " mm to fit " << nFibers << " fibers";
  }
  double cx = 0, cy = 0;
  for (const auto& p : pts) {
    cx += p.first;
    cy += p.second;
  }
  cx /= pts.size();
  cy /= pts.size();
  std::sort(pts.begin(), pts.end(), [&](const auto& a, const auto& b) {
    return TMath::Hypot(a.first - cx, a.second - cy) < TMath::Hypot(b.first - cx, b.second - cy);
  });
  pts.resize(nFibers);
  return pts;
}

// Creates and registers a named translation, so that it can be referred to by name inside a composite-shape formula.
void trans(const std::string& name, double x, double y, double z = 0)
{
  (new TGeoTranslation(name.c_str(), x, y, z))->RegisterYourself();
}

// Builds the shape of one FV0-style scintillator cell (ring sector minus separators and screw holes).
//   tag     : prefix for all shape names (unique per kind of cell)
//   type    : 'a' = sector between 45 and 90 deg, 'b' = sector between 0 and 45 deg
//   ir      : ring index (0 = Ring 1)
//   dz      : thickness along z; extraR / extraDz : enlarge the cell radially / in thickness (Teflon wrap)
TGeoCompositeShape* makeCellShape(const std::string& tag, char type, int ir, double dz, double extraR, double extraDz)
{
  const double rMin = ringR[ir] - extraR, rMax = ringR[ir + 1] + extraR, rMid = 0.5 * (rMin + rMax);
  const double rMinS = rMin + drSep, rMaxS = rMax - drSep;
  const double xHole = drSep + dxHoleExt;
  const double d2r = TMath::DegToRad();
  const std::string name = tag + "CELL" + type + std::to_string(ir);

  const std::string sepName = tag + "SecSep";
  const std::string rot45 = tag + "SecSepRot45", rot90 = tag + "SecSepRot90";
  if (!gGeoManager->GetListOfShapes()->FindObject(sepName.c_str())) {
    new TGeoBBox(sepName.c_str(), ringR[5] - drSep + eps, drSep, (dz + extraDz) / 2);
    (new TGeoRotation(rot45.c_str(), 45, 0, 0))->RegisterYourself();
    (new TGeoRotation(rot90.c_str(), 90, 0, 0))->RegisterYourself();
    new TGeoTube((tag + "HoleSmall").c_str(), 0, drHoleSmall, (dz + extraDz) / 2);
    new TGeoTube((tag + "HoleLarge").c_str(), 0, drHoleLarge, (dz + extraDz) / 2);
    new TGeoBBox((tag + "HoleSmallCut").c_str(), dxHoleExt, drHoleSmall, (dz + extraDz) / 2);
    new TGeoBBox((tag + "HoleLargeCut").c_str(), dxHoleExt, drHoleLarge, (dz + extraDz) / 2);
  }
  auto hole = [&](bool small) { return tag + (small ? "HoleSmall" : "HoleLarge"); };
  auto holeCut = [&](bool small) { return tag + (small ? "HoleSmallCut" : "HoleLargeCut"); };

  const std::string baseName = name + "Shape";
  const double phi1 = (type == 'a') ? 45 : 0, phi2 = (type == 'a') ? 90 : 45;
  if (ir == 0) {
    // with a straight edge the outer radius must not round the corners: the chord cut defines the edge
    const double rFull = (straightOuterEdge && type == 'a') ? 1.5 * rMax : rMaxS - eps;
    new TGeoTubeSeg((baseName + "Full").c_str(), 0, rFull, (dz + extraDz) / 2 - eps, phi1, phi2);
    new TGeoTube((baseName + "Cut").c_str(), 0, rMinS, dz + extraDz);
    trans(baseName + "CutTrans", xShiftInner, 0);
    std::string baseFormula = baseName + "Full-" + baseName + "Cut:" + baseName + "CutTrans";
    if (straightOuterEdge && type == 'a') {
      // cut away everything beyond the chord: box with its inner face on the chord
      const double dChord = (rMaxS - eps) * TMath::Cos(22.5 * d2r) + outerEdgeOffset, hb = rMax;
      const double nx = TMath::Cos(chordNormalDeg * d2r), ny = TMath::Sin(chordNormalDeg * d2r);
      new TGeoBBox((baseName + "OutCut").c_str(), hb, rMax + 1, (dz + extraDz));
      (new TGeoCombiTrans((baseName + "OutCutTrans").c_str(), (dChord + hb) * nx, (dChord + hb) * ny, 0,
                          new TGeoRotation((baseName + "OutCutRot").c_str(), chordNormalDeg, 0, 0)))
        ->RegisterYourself();
      baseFormula += "-" + baseName + "OutCut:" + baseName + "OutCutTrans";
    }
    new TGeoCompositeShape(baseName.c_str(), baseFormula.c_str());
  } else {
    new TGeoTubeSeg(baseName.c_str(), rMinS + eps, rMaxS - eps, (dz + extraDz) / 2 - eps, phi1, phi2);
  }

  std::string f = baseName;
  auto sub = [&](const std::string& shape, const std::string& id, double x, double y) {
    trans(name + "Tr" + id, x, y);
    f += "-" + shape + ":" + name + "Tr" + id;
  };
  const double s45 = TMath::Sin(45 * d2r), c45 = TMath::Cos(45 * d2r);
  const double s22 = TMath::Sin(22.5 * d2r), c22 = TMath::Cos(22.5 * d2r);

  if (type == 'a') {
    f += "-" + sepName + ":" + rot45 + "-" + sepName + ":" + rot90;
    sub(hole(ir < 2), "H1", xHole, TMath::Cos(TMath::ASin(xHole / rMax)) * rMax);
    sub(holeCut(ir < 2), "H1c", 0, TMath::Cos(TMath::ASin(xHole / rMax)) * rMax);
    sub(hole(ir < 2), "H3", s45 * rMax, c45 * rMax);
    if (ir > 0) {
      sub(hole(ir < 3), "H2", xHole, TMath::Cos(TMath::ASin(xHole / rMin)) * rMin);
      sub(holeCut(ir < 3), "H2c", 0, TMath::Cos(TMath::ASin(xHole / rMin)) * rMin);
      sub(hole(ir < 3), "H4", s45 * rMin, c45 * rMin);
    }
    if (ir > 1) {
      sub(hole(false), "H5", s22 * rMax, c22 * rMax);
    }
    if (ir > 2) {
      sub(hole(false), "H6", s22 * rMin, c22 * rMin);
    }
    if (ir == 4) {
      sub(hole(false), "H7", xHole, TMath::Cos(TMath::ASin(xHole / rMid)) * rMid);
      sub(holeCut(false), "H7c", 0, TMath::Cos(TMath::ASin(xHole / rMid)) * rMid);
      sub(hole(false), "H8", s45 * rMid, c45 * rMid);
    }
  } else {
    f += "-" + sepName + "-" + sepName + ":" + rot45;
    sub(hole(ir < 2), "H1", s45 * rMax, c45 * rMax);
    sub(hole(ir < 2), "H3", rMax, 0);
    if (ir > 0) {
      sub(hole(ir < 3), "H2", s45 * rMin, c45 * rMin);
      sub(hole(ir < 3), "H4", rMin, 0);
    }
    if (ir > 1) {
      sub(hole(false), "H5", c22 * rMax, s22 * rMax);
    }
    if (ir > 2) {
      sub(hole(false), "H6", c22 * rMin, s22 * rMin);
    }
    if (ir == 4) {
      sub(hole(false), "H7", s45 * rMid, c45 * rMid);
      sub(hole(false), "H8", rMid, 0);
    }
  }
  return new TGeoCompositeShape((name + "CS" + std::string(extraR ? "ext" : "base")).c_str(), f.c_str());
}
} // namespace

TGeoVolume* o2::fv0::fdprototype::build()
{
  LOG(info) << "FV0: Building FD prototype geometry";

  TGeoMedium* mScint = medium("Scintillator$");
  TGeoMedium* mPlast = medium("Plastic$");
  TGeoMedium* mGlue = medium("FDGlue$");
  TGeoMedium* mTeflon = medium("FDTeflon$");
  TGeoMedium* mFiber = medium("FDFiber$");
  TGeoMedium* mGlass = medium("FDPMTGlass$");
  TGeoMedium* mPhotocathode = medium("FDPhotocathode$");
  TGeoMedium* mAluminium = medium("Aluminium$");

  auto* assembly = new TGeoVolumeAssembly("FDPROTO");

  // 1. Scintillator and Teflon wrap. The wrap is the same cell, slightly larger, made hollow (scintillator subtracted) and
  //    open on the +z face, where the glue and the support sit, so that no volumes overlap.
  auto* scintShape = makeCellShape("FDSCINT", 'a', 0, dzScint, 0.0, 0.0);
  auto* teflonFullShape = makeCellShape("FDTEFLON", 'a', 0, dzScint, teflonThickness, 2 * teflonThickness);
  new TGeoBBox("FDTeflonTopCut", 2 * ringR[1], 2 * ringR[1], 1.);
  trans("FDTeflonTopCutTrans", 0, 0, dzScint / 2 + 1.);
  auto* teflonShape = new TGeoCompositeShape("FDTeflonShellShape",
                                             (std::string(teflonFullShape->GetName()) + "-" + scintShape->GetName() +
                                              "-FDTeflonTopCut:FDTeflonTopCutTrans")
                                               .c_str());
  auto* scint = new TGeoVolume(sensitiveVolumeName, scintShape, mScint);
  auto* teflonWrap = new TGeoVolume("FDPROTO_TEFLON", teflonShape, mTeflon);
  scint->SetLineColor(kAzure + 1);
  scint->SetTransparency(30);
  teflonWrap->SetLineColor(kWhite);
  teflonWrap->SetTransparency(50);

  // 2. Optical glue layer and plastic support (same cell shape, different thickness/material)
  auto* glueLayer = new TGeoVolume("FDPROTO_GLUE", makeCellShape("FDGLUE", 'a', 0, dzGlue, 0.0, 0.0), mGlue);
  auto* plastHolder = new TGeoVolume("FDPROTO_SUPPORT", makeCellShape("FDPLAST", 'a', 0, dzPlastic, 0.0, 0.0), mPlast);
  glueLayer->SetLineColor(kCyan - 10);
  glueLayer->SetTransparency(70);
  plastHolder->SetLineColor(kCyan - 6);
  plastHolder->SetTransparency(40);

  // 3. Black plastic box: outer box minus inner box (closed shell).
  const double innerHalfX = boxLength / 2 - boxWall, innerHalfY = boxHeight / 2 - boxWall;
  new TGeoBBox("FDBoxOuter", boxLength / 2, boxHeight / 2, boxWidth / 2);
  new TGeoBBox("FDBoxInner", innerHalfX, innerHalfY, boxWidth / 2 - boxWall);
  auto* box = new TGeoVolume("FDPROTO_BOX", new TGeoCompositeShape("FDBoxShape", "FDBoxOuter-FDBoxInner"), mPlast);
  box->SetLineColor(kBlack);
  box->SetTransparency(80);

  // Stack = Teflon + scintillator + glue + plastic support, in the cell frame (z = thickness, support on the +z side).
  // Placed in the box with its cell z axis along +Y; the sector symmetry axis points to -X (wide outer edge towards the PMT side).
  auto* stack = new TGeoVolumeAssembly("FDPROTO_STACK");
  stack->AddNode(teflonWrap, 1, new TGeoTranslation(0, 0, 0));
  stack->AddNode(scint, 1, new TGeoTranslation(0, 0, 0));
  const double zGluePos = dzScint / 2 + dzGlue / 2;
  const double zPlastPos = zGluePos + dzGlue / 2 + dzPlastic / 2;
  stack->AddNode(glueLayer, 1, new TGeoTranslation(0, 0, zGluePos));
  stack->AddNode(plastHolder, 1, new TGeoTranslation(0, 0, zPlastPos));

  const double xStackC = innerHalfX - boxGap - stackHalfU;                    // box X of the sector centre (pivot)
  const double yStack = -innerHalfY + boxGap + teflonThickness + dzScint / 2; // box Y of the stack centre
  // gam = angle (from +X towards +Z) of the cell-frame x axis. 247.5 deg makes the sector symmetry axis point to -X.
  const double gam = (180 + 67.5) * TMath::DegToRad(), cg = TMath::Cos(gam), sg = TMath::Sin(gam);
  const double cxl = stackCenterU * TMath::Cos(67.5 * TMath::DegToRad()), cyl = stackCenterU * TMath::Sin(67.5 * TMath::DegToRad());
  // Maps a point (px, py) of the cell xy plane to the box (X, Z), with the pivot going to (xStackC, 0).
  auto toBox = [&](double px, double py, double& X, double& Z) {
    X = xStackC + cg * (px - cxl) + sg * (py - cyl);
    Z = sg * (px - cxl) - cg * (py - cyl);
  };
  const double rotStackM[9] = {cg, sg, 0, 0, 0, 1, sg, -cg, 0}; // columns = cell-frame axes in the box frame
  auto* rotStack = new TGeoRotation("FDPROTO_STACKROT");
  rotStack->SetMatrix(rotStackM);

  // 4. 120 fibers, each made of (a) a vertical leg, (b) a 90 deg bend towards -X, (c) a horizontal leg converging on the PMT.
  const double zPlastTop = dzScint / 2 + dzGlue + dzPlastic; // z (cell frame) of the top face of the plastic support
  const double yFibStart = yStack + zPlastTop + fiberStartGap;

  double fiberPitchUsed = fiberPitch;
  const auto fibers = fiberPositions(fiberPitchUsed);
  double bx = 0, by = 0;
  for (const auto& p : fibers) {
    bx += p.first;
    by += p.second;
  }
  bx /= fibers.size();
  by /= fibers.size();
  double rMaxFiber = 0;
  for (const auto& p : fibers) {
    rMaxFiber = std::max(rMaxFiber, TMath::Hypot(p.first - bx, p.second - by));
  }

  // Convergence factor: fiber-to-fiber distance on the PMT window = 1.01 x fiber diameter
  const double k = 1.01 * fiberDiameter / fiberPitchUsed;
  double xc, zc; // box X, Z of the bundle axis at the support
  toBox(bx, by, xc, zc);

  const double pmtTotDz = pmtWindowDz + pmtPcDz + pmtBodyDz;
  double oxMax = 0;
  for (const auto& p : fibers) {
    double X, Z;
    toBox(p.first, p.second, X, Z);
    oxMax = std::max(oxMax, X - xc);
  }
  const double yPmt = innerHalfY - boxGap - std::max(pmtSize / 2, oxMax + fiberDiameter / 2);

  const double R = bendRadius;
  const double yArc = yPmt - R;                          // Y where the bend starts
  const double L1 = yArc - yFibStart;                    // vertical leg length
  const double L2 = dzFiber - TMath::PiOver2() * R - L1; // horizontal leg length
  const double xWin = xc - R - L2;                       // box X of the outer face of the PMT window
  if (L1 < 0 || L2 < 0) {
    LOG(fatal) << "FD prototype: fiber path not feasible (L1 = " << L1 << " cm, L2 = " << L2 << " cm): change bendRadius";
  }
  if (xWin - pmtTotDz < -innerHalfX) {
    LOG(fatal) << "FD prototype: PMT outside the box (back at x = " << xWin - pmtTotDz << " cm, inner wall at " << -innerHalfX << " cm)";
  }

  auto* fiberMatrix = new TGeoVolumeAssembly("FDPROTO_FIBERS");
  // Straight segments with the same length (to 1 um) share one volume, to keep the number of volumes small.
  std::map<long, TGeoVolume*> segVolumes;
  int nSeg = 0;
  // Adds one fiber segment: a straight tube from point p to point q (box frame). Returns its length.
  auto addTube = [&](const double* p, const double* q) {
    const double dx = q[0] - p[0], dy = q[1] - p[1], dz = q[2] - p[2];
    const double len = TMath::Sqrt(dx * dx + dy * dy + dz * dz);
    const double zx = dx / len, zy = dy / len, zz = dz / len;
    // x' = ref x z' (normalised), ref = Y, or Z if z' is almost parallel to Y; y' = z' x x'
    const double ry = TMath::Abs(zy) > 0.9 ? 0 : 1, rz = 1 - ry;
    double xx = ry * zz - rz * zy, xy = rz * zx, xz = -ry * zx;
    const double xn = TMath::Sqrt(xx * xx + xy * xy + xz * xz);
    xx /= xn;
    xy /= xn;
    xz /= xn;
    const double yx = zy * xz - zz * xy, yy = zz * xx - zx * xz, yz = zx * xy - zy * xx;
    const double rot[9] = {xx, yx, zx, xy, yy, zy, xz, yz, zz};
    const long key = std::lround(len * 1e4);
    auto& seg = segVolumes[key];
    if (!seg) {
      seg = new TGeoVolume(Form("FDPROTO_FIBSEG%zu", segVolumes.size()), new TGeoTube(0, fiberDiameter / 2, 0.5e-4 * key), mFiber);
      seg->SetLineColor(kOrange + 1);
    }
    auto* r = new TGeoRotation();
    r->SetMatrix(rot);
    fiberMatrix->AddNode(seg, ++nSeg, new TGeoCombiTrans(0.5 * (p[0] + q[0]), 0.5 * (p[1] + q[1]), 0.5 * (p[2] + q[2]), r));
    return len;
  };

  double fibLenMin = 1e9, fibLenMax = 0;
  for (const auto& p : fibers) {
    double X, Z;
    toBox(p.first, p.second, X, Z);
    const double ox = X - xc, oz = Z - zc;             // offset of this fiber from the bundle axis
    const double A[3] = {xc + ox, yFibStart, zc + oz}; // start of the fiber, just above the support
    const double rho = R + ox;                         // bend radius of this fiber
    double prev[3] = {xc + ox, yArc, zc + oz};         // end of the vertical leg = start of the bend
    double total = addTube(A, prev);
    for (int i = 1; i <= nBendSegments; ++i) {
      const double th = TMath::PiOver2() * i / nBendSegments;
      const double cur[3] = {xc - R + rho * TMath::Cos(th), yArc + rho * TMath::Sin(th), zc + oz};
      total += addTube(prev, cur);
      std::copy(cur, cur + 3, prev);
    }
    const double end[3] = {xWin + fiberStartGap, yPmt + k * ox, zc + k * oz};
    total += addTube(prev, end);
    fibLenMin = std::min(fibLenMin, total);
    fibLenMax = std::max(fibLenMax, total);
  }
  const double rMaxBundle = k * rMaxFiber; // radius of the narrowed bundle on the PMT window

  // 5. PMT H11934-100: glass window, photocathode, metal body; its axis is on the bundle axis.
  auto* pmt = new TGeoVolumeAssembly("FDPROTO_PMT");
  auto* pmtWindow = new TGeoVolume("FDPROTO_PMTWINDOW", new TGeoTube("FDPMTWindowShape", 0, pmtWindowR, pmtWindowDz / 2), mGlass);
  auto* pmtCathode = new TGeoVolume("FDPROTO_PMTCATHODE", new TGeoTube("FDPMTCathodeShape", 0, pmtActiveR, pmtPcDz / 2), mPhotocathode);
  auto* pmtBody = new TGeoVolume("FDPROTO_PMTBODY", new TGeoBBox("FDPMTBodyShape", pmtSize / 2, pmtSize / 2, pmtBodyDz / 2), mAluminium);
  pmtWindow->SetLineColor(kGray + 1);
  pmtWindow->SetTransparency(60);
  pmtCathode->SetLineColor(kRed);
  pmtBody->SetLineColor(kGray + 3);
  pmtBody->SetTransparency(40);

  // PMT local frame: z = 0 is the outer face of the window; the body extends towards +z.
  pmt->AddNode(pmtWindow, 1, new TGeoTranslation(0, 0, pmtWindowDz / 2));
  pmt->AddNode(pmtCathode, 1, new TGeoTranslation(0, 0, pmtWindowDz + pmtPcDz / 2));
  pmt->AddNode(pmtBody, 1, new TGeoTranslation(0, 0, pmtWindowDz + pmtPcDz + pmtBodyDz / 2));
  const double rotPmtM[9] = {0, 0, -1, 1, 0, 0, 0, -1, 0}; // local z -> -X, x -> +Y, y -> -Z
  auto* rotPmt = new TGeoRotation("FDPROTO_PMTROT");
  rotPmt->SetMatrix(rotPmtM);

  // Placement in the box frame
  assembly->AddNode(box, 1, new TGeoTranslation(0, 0, 0));
  assembly->AddNode(stack, 1, new TGeoCombiTrans(xStackC - (cg * cxl + sg * cyl), yStack, -(sg * cxl - cg * cyl), rotStack));
  assembly->AddNode(fiberMatrix, 1, new TGeoTranslation(0, 0, 0));
  assembly->AddNode(pmt, 1, new TGeoCombiTrans(xWin, yPmt, zc, rotPmt));

  LOG(info) << "FD prototype: box " << boxLength << " x " << boxHeight << " x " << boxWidth << " cm; " << fibers.size()
            << " fibers (" << segVolumes.size() << " segment volumes), length " << fibLenMin << "-" << fibLenMax
            << " cm; bundle radius on PMT " << rMaxBundle + fiberDiameter / 2 << " cm (photocathode " << pmtActiveR << " cm)";
  if (rMaxBundle + fiberDiameter / 2 > pmtActiveR) {
    LOG(warning) << "FD prototype: the fiber bundle is not fully covered by the photocathode";
  }
  return assembly;
}
