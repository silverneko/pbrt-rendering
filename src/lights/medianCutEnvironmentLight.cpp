
/*
    pbrt source code Copyright(c) 1998-2012 Matt Pharr and Greg Humphreys.

    This file is part of pbrt.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are
    met:

    - Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
    IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
    PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
    HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 */


// lights/infinite.cpp*
#include <queue>

#include "stdafx.h"
#include "lights/medianCutEnvironmentLight.h"
#include "sh.h"
#include "montecarlo.h"
#include "paramset.h"
#include "imageio.h"

using std::queue;


class P4 {
public:
    int x1, y1, x2, y2;
    P4 (int a, int b, int c, int d) : x1(a), y1(b), x2(c), y2(d) {}
    void get(int &a, int &b, int &c, int &d) {
      a = x1;
      b = y1;
      c = x2;
      d = y2;
    }
};

// MedianCutEnvironmentLight Utility Classes
struct MedianCutEnvironmentCube {
    // InfiniteAreaCube Public Methods
    MedianCutEnvironmentCube(const MedianCutEnvironmentLight *l, const Scene *s,
                     float t, bool cv, float pe)
        : light(l), scene(s), time(t), pEpsilon(pe), computeVis(cv) { }
    Spectrum operator()(int, int, const Point &p, const Vector &w) {
        Ray ray(p, w, pEpsilon, INFINITY, time);
        if (!computeVis || !scene->IntersectP(ray))
            return light->Le(RayDifferential(ray));
        return 0.f;
    }
    const MedianCutEnvironmentLight *light;
    const Scene *scene;
    float time, pEpsilon;
    bool computeVis;
};


// MedianCutEnvironmentLight Method Definitions
MedianCutEnvironmentLight::~MedianCutEnvironmentLight() {
    delete [] Intensity;
    delete [] LightSourceU;
    delete [] LightSourceV;
    delete radianceMap;
}


MedianCutEnvironmentLight::MedianCutEnvironmentLight(const Transform &light2world,
        const Spectrum &L, int ns, const string &texmap)
    : Light(light2world, ns), nLightSources(ns) {
    int width = 0, height = 0;
    RGBSpectrum *texels = NULL;
    // Read texel data from _texmap_ into _texels_
    if (texmap != "") {
        texels = ReadImage(texmap, &width, &height);
        if (texels)
            for (int i = 0; i < width * height; ++i)
                texels[i] *= L.ToRGBSpectrum();
    }
    if (!texels) {
        width = height = 1;
        texels = new RGBSpectrum[1];
        texels[0] = L.ToRGBSpectrum();
    }
    radianceMap = new MIPMap<RGBSpectrum>(width, height, texels);

    // Initialize sum area table
    float *SumAreaTable = new float[(width+1)*(height+1)];
    for (int v = 0; v < height+1; ++v) {
        SumAreaTable[v*(width+1)] = 0.f;
    }
    for (int u = 1; u < width+1; ++u) {
        SumAreaTable[u] = 0.f;
    }
#define Sum(_x, _y) (SumAreaTable[ ((_x) + 1) + ((_y) + 1) * (width+1) ])
    for (int v = 0; v < height; ++v) {
        for (int u = 0; u < width; ++u) {
            float PixelIntensity = texels[u+v*width].y();
            Sum(u, v) = Sum(u-1, v) + Sum(u, v-1) - Sum(u-1, v-1)
              + PixelIntensity * sinf((float)v / height * M_PI);
        }
    }
    // Do median cut
#define Rec(_x1, _y1, _x2, _y2) \
        (Sum(_x2-1, _y2-1) - Sum(_x1-1, _y2-1) - Sum(_x2-1, _y1-1) + Sum(_x1-1, _y1-1))
    queue<P4> WorkList;
    WorkList.push(P4(0, 0, width, height));
    while (WorkList.size() < nLightSources) {
      int x1, y1, x2, y2;
      WorkList.front().get(x1, y1, x2, y2);
      WorkList.pop();
      int h, w;
      h = y2 - y1;
      w = (x2 - x1) * sinf((float)(y2 + y1) / 2.f / height * M_PI);
      if (h > w) {
        int lb = y1, hb = y2;
        while (lb != hb-1) {
          int mid = (lb + hb) / 2;
          if (Rec(x1, y1, x2, mid) > Rec(x1, mid, x2, y2)) {
            hb = mid;
          } else {
            lb = mid;
          }
        }
        if (y1 != lb) WorkList.push(P4(x1, y1, x2, lb));
        if (lb != y2) WorkList.push(P4(x1, lb, x2, y2));
      } else {
        int lb = x1, hb = x2;
        while (lb != hb-1) {
          int mid = (lb + hb) / 2;
          if (Rec(x1, y1, mid, y2) > Rec(mid, y1, x2, y2)) {
            hb = mid;
          } else {
            lb = mid;
          }
        }
        if (x1 != lb) WorkList.push(P4(x1, y1, lb, y2));
        if (lb != x2) WorkList.push(P4(lb, y1, x2, y2));
      }
    }

    // Compute centroid and place light source at each centroid
    invnLightSources = 1.f / nLightSources;
    Intensity = new RGBSpectrum[nLightSources];
    LightSourceU = new float[nLightSources];
    LightSourceV = new float[nLightSources];
    for (int i = 0; i < nLightSources; ++i) {
      int x1, y1, x2, y2;
      WorkList.front().get(x1, y1, x2, y2);
      WorkList.pop();
      float sx = 0.f, sy = 0.f;
      for (int u = x1; u < x2; ++u) {
        sx += Rec(u, y1, u+1, y2) * u;
      }
      sx /= Rec(x1, y1, x2, y2);
      for (int v = y1; v < y2; ++v) {
        sy += Rec(x1, v, x2, v+1) * v;
      }
      sy /= Rec(x1, y1, x2, y2);
      // Place light source at (sx, sy)
      RGBSpectrum RegionIntensity(0.f);
      for (int v = y1; v < y2; ++v) {
        for (int u = x1; u < x2; ++u) {
          RegionIntensity += texels[u+v*width] * sinf((float)v / height * M_PI);
          assert(sinf((float)v / height * M_PI) >= 0.f);
        }
      }
      LightSourceU[i] = sx / width * 2.f * M_PI;
      LightSourceV[i] = sy / height * M_PI;
      Intensity[i] = RegionIntensity * (M_PI / height * 2.f * M_PI / width);
    }
#undef Rec
#undef Sum
    delete[] SumAreaTable;
    delete[] texels;
}


Spectrum MedianCutEnvironmentLight::Power(const Scene *scene) const {
    Point worldCenter;
    float worldRadius;
    scene->WorldBound().BoundingSphere(&worldCenter, &worldRadius);
    return M_PI * worldRadius * worldRadius *
        Spectrum(radianceMap->Lookup(.5f, .5f, .5f), SPECTRUM_ILLUMINANT);
}


Spectrum MedianCutEnvironmentLight::Le(const RayDifferential &r) const {
    Vector wh = Normalize(WorldToLight(r.d));
    float s = SphericalPhi(wh) * INV_TWOPI;
    float t = SphericalTheta(wh) * INV_PI;
    return Spectrum(radianceMap->Lookup(s, t), SPECTRUM_ILLUMINANT);
}


void MedianCutEnvironmentLight::SHProject(const Point &p, float pEpsilon,
        int lmax, const Scene *scene, bool computeLightVis,
        float time, RNG &rng, Spectrum *coeffs) const {
    // Project _MedianCutEnvironmentLight_ to SH using Monte Carlo if visibility needed
    if (computeLightVis) {
        Light::SHProject(p, pEpsilon, lmax, scene, computeLightVis,
                         time, rng, coeffs);
        return;
    }
    for (int i = 0; i < SHTerms(lmax); ++i)
        coeffs[i] = 0.f;
    int ntheta = radianceMap->Height(), nphi = radianceMap->Width();
    if (min(ntheta, nphi) > 50) {
        // Project _MedianCutEnvironmentLight_ to SH from lat-long representation

        // Precompute $\theta$ and $\phi$ values for lat-long map projection
        float *buf = new float[2*ntheta + 2*nphi];
        float *bufp = buf;
        float *sintheta = bufp;  bufp += ntheta;
        float *costheta = bufp;  bufp += ntheta;
        float *sinphi = bufp;    bufp += nphi;
        float *cosphi = bufp;
        for (int theta = 0; theta < ntheta; ++theta) {
            sintheta[theta] = sinf((theta + .5f)/ntheta * M_PI);
            costheta[theta] = cosf((theta + .5f)/ntheta * M_PI);
        }
        for (int phi = 0; phi < nphi; ++phi) {
            sinphi[phi] = sinf((phi + .5f)/nphi * 2.f * M_PI);
            cosphi[phi] = cosf((phi + .5f)/nphi * 2.f * M_PI);
        }
        float *Ylm = ALLOCA(float, SHTerms(lmax));
        for (int theta = 0; theta < ntheta; ++theta) {
            for (int phi = 0; phi < nphi; ++phi) {
                // Add _MedianCutEnvironmentLight_ texel's contribution to SH coefficients
                Vector w = Vector(sintheta[theta] * cosphi[phi],
                                  sintheta[theta] * sinphi[phi],
                                  costheta[theta]);
                w = Normalize(LightToWorld(w));
                Spectrum Le = Spectrum(radianceMap->Texel(0, phi, theta),
                                       SPECTRUM_ILLUMINANT);
                SHEvaluate(w, lmax, Ylm);
                for (int i = 0; i < SHTerms(lmax); ++i)
                    coeffs[i] += Le * Ylm[i] * sintheta[theta] *
                        (M_PI / ntheta) * (2.f * M_PI / nphi);
            }
        }

        // Free memory used for lat-long theta and phi values
        delete[] buf;
    }
    else {
        // Project _MedianCutEnvironmentLight_ to SH from cube map sampling
        SHProjectCube(MedianCutEnvironmentCube(this, scene, time, computeLightVis,
                                       pEpsilon),
                      p, 200, lmax, coeffs);
    }
}


MedianCutEnvironmentLight *
CreateMedianCutEnvironmentLight(const Transform &light2world,
                                const ParamSet &paramSet) {
    Spectrum L = paramSet.FindOneSpectrum("L", Spectrum(1.0));
    Spectrum sc = paramSet.FindOneSpectrum("scale", Spectrum(1.0));
    string texmap = paramSet.FindOneFilename("mapname", "");
    int nSamples = paramSet.FindOneInt("nsamples", 1);
    if (PbrtOptions.quickRender) nSamples = max(1, nSamples / 4);
    return new MedianCutEnvironmentLight(light2world, L * sc, nSamples, texmap);
}


Spectrum MedianCutEnvironmentLight::Sample_L(const Point &p, float pEpsilon,
        const LightSample &ls, float time, Vector *wi, float *pdf,
        VisibilityTester *visibility) const {
    int li = Float2Int(ls.uPos[0] * nLightSources);
    float phi = LightSourceU[li];
    float theta = LightSourceV[li];
    RGBSpectrum Ls = Intensity[li];

    float costheta = cosf(theta), sintheta = sinf(theta);
    float sinphi = sinf(phi), cosphi = cosf(phi);
    *wi = LightToWorld(Vector(sintheta * cosphi, sintheta * sinphi,
                              costheta));
    *pdf = invnLightSources;
    visibility->SetRay(p, pEpsilon, *wi, time);
    return Ls;
}


float MedianCutEnvironmentLight::Pdf(const Point &, const Vector &w) const {
    return 0.f;
}


Spectrum MedianCutEnvironmentLight::Sample_L(const Scene *scene,
        const LightSample &ls, float u1, float u2, float time,
        Ray *ray, Normal *Ns, float *pdf) const {
    int li = Float2Int(ls.uPos[0] * nLightSources);
    float phi = LightSourceU[li];
    float theta = LightSourceV[li];
    RGBSpectrum Ls = Intensity[li];

    // Compute direction for infinite light sample ray
    float costheta = cosf(theta), sintheta = sinf(theta);
    float sinphi = sinf(phi), cosphi = cosf(phi);
    Vector d = -LightToWorld(Vector(sintheta * cosphi, sintheta * sinphi,
                                    costheta));
    *Ns = (Normal)d;

    // Compute origin for infinite light sample ray
    Point worldCenter;
    float worldRadius;
    scene->WorldBound().BoundingSphere(&worldCenter, &worldRadius);
    Vector v1, v2;
    CoordinateSystem(-d, &v1, &v2);
    float d1, d2;
    ConcentricSampleDisk(u1, u2, &d1, &d2);
    Point Pdisk = worldCenter + worldRadius * (d1 * v1 + d2 * v2);
    *ray = Ray(Pdisk + worldRadius * -d, d, 0., INFINITY, time);

    // Compute _MedianCutEnvironmentLight_ ray PDF
    float areaPdf = 1.f / (M_PI * worldRadius * worldRadius);
    *pdf = invnLightSources * areaPdf;
    return Ls;
}


