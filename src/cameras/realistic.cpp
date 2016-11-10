#include "cameras/realistic.h"
#include "sampler.h"
#include "stdafx.h"
#include "transform.h"

#include <cstdio>
#include <sstream>

using namespace std;

bool Lens::getRefractedRay(const Ray &r, float eta, Ray *refracted) const {
    if (radius == 0) {
      float t0;
      t0 = (posZ - r.o.z) / r.d.z;
      if (t0 < 0) {
        return false;
      }
      Point phit(r(t0));
      if (phit.x * phit.x + phit.y * phit.y > aperture) {
        return false;
      }
      *refracted = r;
      return true;
    }

    // Transform _Ray_ to object space
    Ray ray;
    Translate(Vector(0, 0, -posZ))(r, &ray);

    // Compute quadratic sphere coefficients
    float A = ray.d.x*ray.d.x + ray.d.y*ray.d.y + ray.d.z*ray.d.z;
    float B = 2 * (ray.d.x*ray.o.x + ray.d.y*ray.o.y + ray.d.z*ray.o.z);
    float C = ray.o.x*ray.o.x + ray.o.y*ray.o.y +
              ray.o.z*ray.o.z - radius*radius;

    // Solve quadratic equation for _t_ values
    float t0, t1;
    if (!Quadratic(A, B, C, &t0, &t1))
        return false;

    // Compute intersection distance along ray
    if (t0 > ray.maxt || t1 < ray.mint)
        return false;
    float thit = t0;
    if (t0 < ray.mint) {
        thit = t1;
        if (thit > ray.maxt) return false;
    }

    Point phit(r(thit));

    if (phit.x * phit.x + phit.y * phit.y > aperture) {
      return false;
    }

    Vector N(phit - Vector(0, 0, posZ)), I(r.d), T;
    N /= radius;
    if (N.z > 0) N *= -1;
    float c1 = Dot(N, -I);
    assert(c1 > 0);
    float t = 1 - eta * eta * (1 - c1 * c1);
    if (t < 0) {
      return false;
    }
    float c2 = sqrt(t);
    T = eta * I + (eta * c1 - c2) * N;
    *refracted = Ray(phit, T, 0, INFINITY);
    return true;
}

void parseLensSpec(const string &specfile, float zpos, vector<Lens> &Lenses) {
  char buffer[512];
  FILE *Specfile = fopen(specfile.c_str(), "r");
  while (fgets(buffer, sizeof(buffer), Specfile)) {
    string buf(buffer);
    if (buf.find('#') != -1) {
      continue;
    }
    istringstream ssin(buf);
    float rad, thick, n, aper;
    ssin >> rad >> thick >> n >> aper;
    if (n == 0) {
      n = 1;
    }
    aper = aper * aper * 0.25;
    Lenses.push_back(Lens(rad, thick, n, aper));
  }
  reverse(Lenses.begin(), Lenses.end());
  for (int i = 0; i < Lenses.size(); ++i) {
    Lens &lens = Lenses[i];
    zpos += lens.posZ;
    lens.posZ = zpos - lens.radius;
    if (lens.radius < 0) {
      lens.radius *= -1;
    }
    if (i+1 < Lenses.size()) {
      lens.nd /= Lenses[i+1].nd;
    }
  }
  fclose(Specfile);
}

RealisticCamera::RealisticCamera(const AnimatedTransform &cam2world,
                                 float hither, float yon,
                                 float sopen, float sclose,
                                 float filmdistance, float aperture_diameter,
                                 string specfile, float filmdiag, Film *film)
  : Camera(cam2world, sopen, sclose, film) // pbrt-v2 doesnot specify hither and yon
{
  // YOUR CODE HERE -- build and store datastructures representing the given lens
  // and film placement.
  parseLensSpec(specfile, filmdistance, Lenses);
  // exitPupilR = sqrt(Lenses[0].aperture);
  // exitPupilZ = Lenses[0].posZ;
  exitPupilR = aperture_diameter / 2;
  exitPupilZ = filmdistance;

  float aspectRatio = (float)film->yResolution / film->xResolution;
  float x = sqrt(filmdiag * filmdiag / (1 + aspectRatio * aspectRatio));
  float y = x * aspectRatio;
  RasterToCamera = Scale(x / film->xResolution,
                         y / film->yResolution, 1) *
    Scale(-1, 1, 1) *
    Translate(Vector(-film->xResolution, 0, 0));
  CameraToLens = Translate(Vector(x * -0.5, y * -0.5, 0));
  LensToCamera = Inverse(CameraToLens);
}

float RealisticCamera::GenerateRay(const CameraSample &sample,
                                   Ray *generated_ray) const {
  // YOUR CODE HERE -- make that ray!
  // use sample->imageX and sample->imageY to get raster-space coordinates
  // of the sample point on the film.
  // use sample->lensU and sample->lensV to get a sample position on the lens
  Point ro(sample.imageX, sample.imageY, 0);
  ro = RasterToCamera(ro);
  ro = CameraToLens(ro);
  float r = exitPupilR * sqrt(sample.lensU), theta = 2 * M_PI * sample.lensV;
  Point lensSample(r * cos(theta), r * sin(theta), exitPupilZ);
  Vector rd(Normalize(lensSample - ro));
  Ray ray(ro, rd, 0, INFINITY);
  float n0 = 1;
  for (int i = 0; i < Lenses.size(); ++i) {
    const Lens &lens = Lenses[i];
    if (lens.getRefractedRay(ray, lens.nd, &ray) == false) {
      return 0;
    }
    n0 = lens.nd;
  }
  ray = LensToCamera(ray);
  // ray.o.z = 0;
  CameraToWorld(ray, generated_ray);
  float Weight = rd.z * rd.z * rd.z * rd.z / exitPupilZ / exitPupilZ * M_PI * exitPupilR * exitPupilR;
  return Weight;
}


RealisticCamera *CreateRealisticCamera(const ParamSet &params,
                                       const AnimatedTransform &cam2world,
                                       Film *film) {
  // Extract common camera parameters from \use{ParamSet}
  float hither = params.FindOneFloat("hither", -1);
  float yon = params.FindOneFloat("yon", -1);
  float shutteropen = params.FindOneFloat("shutteropen", -1);
  float shutterclose = params.FindOneFloat("shutterclose", -1);

  // Realistic camera-specific parameters
  string specfile = params.FindOneString("specfile", "");
  // about 70 mm default to film
  float filmdistance = params.FindOneFloat("filmdistance", 70.0);
  float fstop = params.FindOneFloat("aperture_diameter", 1.0);	
  float filmdiag = params.FindOneFloat("filmdiag", 35.0);

  Assert(hither != -1 && yon != -1 && shutteropen != -1 &&
         shutterclose != -1 && filmdistance != -1);
  if (specfile == "") {
    Severe( "No lens spec file supplied!\n" );
  }
  return new RealisticCamera(cam2world, hither, yon,
                             shutteropen, shutterclose, filmdistance, fstop,
                             specfile, filmdiag, film);
}
