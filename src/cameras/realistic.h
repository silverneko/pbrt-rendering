#if defined(_MSC_VER)
#pragma once
#endif

#ifndef PBRT_CAMERAS_REALISTIC_H
#define PBRT_CAMERAS_REALISTIC_H

#include "camera.h"
#include "paramset.h"
#include "film.h"

class Lens;
class RealisticCamera;

// RealisticCamera Declarations
class RealisticCamera : public Camera {
public:
  RealisticCamera(const AnimatedTransform &cam2world,
                  float hither, float yon, float sopen,
                  float sclose, float filmdistance, float aperture_diameter,
                  string specfile, float filmdiag, Film *film);
  float GenerateRay(const CameraSample &sample, Ray *) const;

private:
  vector<Lens> Lenses;
  float exitPupilR, exitPupilZ;

  Transform RasterToCamera;
  Transform CameraToLens, LensToCamera;
};

class Lens {
public:
  float radius;
  float posZ;
  float nd;
  float aperture;

  Lens(float radius, float posZ, float nd, float aperture)
    : radius(radius), posZ(posZ), nd(nd), aperture(aperture) {}

  bool getRefractedRay(const Ray&, float n0, Ray *) const;
};


RealisticCamera *CreateRealisticCamera(const ParamSet &params,
                                       const AnimatedTransform &cam2world,
                                       Film *film);

#endif	// PBRT_CAMERAS_REALISTIC_H
