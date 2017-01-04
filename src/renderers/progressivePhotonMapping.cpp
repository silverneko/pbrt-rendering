
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


// renderers/progressivePhotonMapping.cpp*
#include "stdafx.h"
#include "renderers/progressivePhotonMapping.h"
#include "scene.h"
#include "film.h"
#include "volume.h"
#include "sampler.h"
#include "samplers/lowdiscrepancy.h"
#include "integrator.h"
#include "montecarlo.h"
#include "progressreporter.h"
#include "camera.h"
#include "intersection.h"

class HitPointUpdateTask : public Task {
public:
    HitPointUpdateTask(Camera *camera, vector<HitPoint> &hitPoints,
                       const KdTree<Photon> *photonMap,
                       ProgressReporter &prog)
      : camera(camera), hitPoints(hitPoints), photonMap(photonMap),
      progress(prog) {}
    void Run();
private:
    Camera &camera;
    vector<HitPoint> &hitPoints;
    const KdTree<Photon> *photonMap;
    ProgressReporter &progress;
};

void HitPointUpdateTask::Run() {
  
}

struct Photon {
    Photon(const Point &pp, const Spectrum &wt, const Vector &w)
        : p(pp), alpha(wt), wi(w) { }
    Photon() { }
    Point p;
    Spectrum alpha;
    Vector wi;
};

class PhotonTracingTask : public Task {
public:
    PhotonTracingTask(int tn, float ti, Mutex &m, int want,
        ProgressReporter &prog, vector<Photon> &photons,
        Distribution1D *distrib, const Scene *sc, const Renderer *sr)
    : taskNum(tn), time(ti), mutex(m), want(want), progress(prog),
    photons(photons), lightDistribution(distrib), scene(sc), renderer (sr) { }
    void Run();
private:
    int taskNum;
    float time;
    Mutex &mutex;
    int want;
    ProgressReporter &progress;
    vector<Photon> &photons;
    const Distribution1D *lightDistribution;
    const Scene *scene;
    const Renderer *renderer;
};

void PhotonTracingTask::Run() {
    MemoryArena arena;
    RNG rng(31 * taskNum);
    vector<Photon> localPhotons;
    uint32_t totalPaths = 0;
    PermutedHalton halton(6, rng);

    for (uint32_t i = 0; i < want; ++i, progress.Update()) {
      float u[6];
      halton.Sample(++totalPaths, u);
      // Choose light to shoot photon from
      float lightPdf;
      int lightNum = lightDistribution->SampleDiscrete(u[0], &lightPdf);
      const Light *light = scene->lights[lightNum];

      // Generate _photonRay_ from light source and initialize _alpha_
      RayDifferential photonRay;
      float pdf;
      LightSample ls(u[1], u[2], u[3]);
      Normal Nl;
      Spectrum Le = light->Sample_L(scene, ls, u[4], u[5],
                                    time, &photonRay, &Nl, &pdf);
      if (pdf == 0.f || Le.IsBlack()) continue;
      Spectrum alpha = (AbsDot(Nl, photonRay.d) * Le) / (pdf * lightPdf);
      if (!alpha.IsBlack()) {
        // Follow photon path through scene and record intersections
        bool specularPath = true;
        Intersection photonIsect;
        int nIntersections = 0;
        while (scene->Intersect(photonRay, &photonIsect)) {
          ++nIntersections;
          // Handle photon/surface intersection
          alpha *= renderer->Transmittance(scene, photonRay, NULL, rng, arena);
          BSDF *photonBSDF = photonIsect.GetBSDF(photonRay, arena);
          BxDFType specularType = BxDFType(BSDF_REFLECTION |
                                           BSDF_TRANSMISSION | BSDF_SPECULAR);
          bool hasNonSpecular = (photonBSDF->NumComponents() >
                                 photonBSDF->NumComponents(specularType));
          Vector wo = -photonRay.d;
          if (hasNonSpecular) {
            // Deposit photon at surface
            Photon photon(photonIsect.dg.p, alpha, wo);
            localPhotons.push_back(photon);
          }

          // Sample new photon ray direction
          Vector wi;
          float pdf;
          BxDFType flags;
          Spectrum fr = photonBSDF->Sample_f(wo, &wi, BSDFSample(rng),
                                             &pdf, BSDF_ALL, &flags);
          if (fr.IsBlack() || pdf == 0.f) break;
          Spectrum anew = alpha * fr *
            AbsDot(wi, photonBSDF->dgShading.nn) / pdf;

          // Possibly terminate photon path with Russian roulette
          float continueProb = min(1.f, anew.y() / alpha.y());
          if (rng.RandomFloat() > continueProb)
            break;
          alpha = anew / continueProb;
          specularPath &= ((flags & BSDF_SPECULAR) != 0);

          if (!specularPath) break;
          photonRay = RayDifferential(photonIsect.dg.p, wi, photonRay,
                                      photonIsect.rayEpsilon);
        }
      }
      arena.FreeAll();
    }

    // Merge local photon data with data in _PhotonIntegrator_
    {
      MutexLock lock(mutex);

      for (uint32_t i = 0; i < localPhotons.size(); ++i) {
        photons.push_back(localPhotons[i]);
      }
    }
}

struct HitPoint {
    Point p;
    Normal n;
    Vector wo;
    BSDF *bsdf;
    float x, y;
    Spectrum L;
    float r;
    int N;
    Spectrum tau;
};

class RayTracingTask : public Task {
public:
    RayTracingTask(const Scene *sc, Renderer *ren, Camera *c,
                   ProgressReporter &pr, Sampler *ms, Sample *sam,
                   vector<HitPoint> &hitPoints, MemoryArena &arena,
                   int maxDepth,
                   int tn, int tc)
      : reporter(pr), hitPoints(hitPoints), globalArena(arena), maxDepth(maxDepth)
    {
        scene = sc; renderer = ren; camera = c; mainSampler = ms;
        origSample = sam; taskNum = tn; taskCount = tc;
    }
    void Run();
private:
    const Scene *scene;
    const Renderer *renderer;
    Camera *camera;
    Sampler *mainSampler;
    ProgressReporter &reporter;
    Sample *origSample;
    int taskNum, taskCount;
    vector<HitPoint> &hitPoints;
    MemoryArena &globalArena;
    int maxDepth;
};

void RayTracingTask::Run() {
    // Get sub-_Sampler_ for _ProgressivePhotonMappingTask_
    Sampler *sampler = mainSampler->GetSubSampler(taskNum, taskCount);
    if (!sampler) {
        reporter.Update();
        return;
    }

    // Declare local variables used for rendering loop
    MemoryArena arena;
    RNG rng(taskNum);

    // Allocate space for samples and intersections
    int maxSamples = sampler->MaximumSampleCount();
    Sample *samples = origSample->Duplicate(maxSamples);

    // Get samples from _Sampler_ and update image
    int sampleCount;
    while ((sampleCount = sampler->GetMoreSamples(samples, rng)) > 0) {
        // Generate camera rays and compute radiance along rays
        for (int i = 0; i < sampleCount; ++i) {
            // Find camera ray for _sample[i]_
            Sample &sample = samples[i];
            RayDifferential r;
            float rayWeight = camera->GenerateRayDifferential(sample, &r);
            r.ScaleDifferentials(1.f / sqrtf(sampler->samplesPerPixel));

            bool newHitPoint = false;
            Spectrum L(0.0);
            Intersection localIsect;
            Intersection *isectp = &localIsect;
            if (!scene->Intersect(r, isectp)) {
              // Handle ray that doesn't intersect any geometry
              for (uint32_t i = 0; i < scene->lights.size(); ++i)
                  L += scene->lights[i]->Le(r);
              L = rayWeight * L;
            } else {
              Spectrum pathThroughput(rayWeight);
              RayDifferential ray(r);
              bool specularBounce = false;
              for (int bounces = 0; ; ++bounces) {
                // Possibly add emitted light at path vertex
                if (bounces == 0 || specularBounce)
                  L += pathThroughput * isectp->Le(-ray.d);

                // Sample illumination from lights to find path contribution
                BSDF *bsdf = isectp->GetBSDF(ray, arena);
                const Point &p = bsdf->dgShading.p;
                const Normal &n = bsdf->dgShading.nn;
                Vector wo = -ray.d;
                L += pathThroughput *
                  UniformSampleOneLight(scene, renderer, arena, p, n, wo,
                                        isectp->rayEpsilon, ray.time, bsdf, &sample, rng);
                BxDFType nonSpecular = BxDFType(BSDF_REFLECTION |
                                                BSDF_TRANSMISSION |
                                                BSDF_DIFFUSE |
                                                BSDF_GLOSSY);
                if (bsdf->NumComponents(nonSpecular) > 0) {
                  bsdf = isectp->GetBSDF(ray, globalArena);
                  hitPoints.push_back((HitPoint){
                                      .p = p,
                                      .n = n,
                                      .wo = wo,
                                      .bsdf = bsdf,
                                      .x = sample.imageX,
                                      .y = sample.imageY,
                                      .L = L,
                                      .r = 0,
                                      .N = 0,
                                      .tau = 0});
                  newHitPoint = true;
                  break;
                }

                // Sample BSDF to get new path direction
                Vector wi;
                BxDFType flags;
                float pdf;
                BSDFSample outgoingBSDFSample = BSDFSample(rng);
                Spectrum f = bsdf->Sample_f(wo, &wi, outgoingBSDFSample, &pdf,
                                            BSDF_ALL, &flags);
                if (f.IsBlack() || pdf == 0.)
                  break;
                specularBounce = (flags & BSDF_SPECULAR) != 0;
                assert(specularBounce == true);
                pathThroughput *= f * AbsDot(wi, n) / pdf;
                ray = RayDifferential(p, wi, ray, isectp->rayEpsilon);

                // Possibly terminate the path
                if (bounces > 3) {
                  float continueProbability = min(.5f, pathThroughput.y());
                  if (rng.RandomFloat() > continueProbability)
                    break;
                  pathThroughput /= continueProbability;
                }

                // Find next vertex of path
                if (!scene->Intersect(ray, &localIsect)) {
                  if (specularBounce)
                    for (uint32_t i = 0; i < scene->lights.size(); ++i)
                      L += pathThroughput * scene->lights[i]->Le(ray);
                  break;
                }
                pathThroughput *= renderer->Transmittance(scene, ray, NULL, rng, arena);
                isectp = &localIsect;
              }
            }

            if (!newHitPoint) {
              camera->film->AddSample(sample, L);
            }
        }

        arena.FreeAll();
    }

    // Clean up after _ProgressivePhotonMappingTask_ is done with its image region
    delete sampler;
    delete[] samples;
    reporter.Update();
}



// ProgressivePhotonMapping Method Definitions
ProgressivePhotonMapping::ProgressivePhotonMapping(Camera *c,
                                 int perPixelSamples,
                                 int maxDepth,
                                 int nIterations,
                                 int nPhotonsPerIter)
  : nPixelSamples(perPixelSamples), maxDepth(maxDepth)
  , nIterations(nIterations), nPhotonsPerIter(nPhotonsPerIter)
{
    camera = c;
}


ProgressivePhotonMapping::~ProgressivePhotonMapping() {
    delete camera;
}


void ProgressivePhotonMapping::Render(const Scene *scene) {
    PBRT_FINISHED_PARSING();
    // Allow integrators to do preprocessing for the scene
    PBRT_STARTED_PREPROCESSING();
    PBRT_FINISHED_PREPROCESSING();
    PBRT_STARTED_RENDERING();
    // Allocate and initialize _sample_
    int x0, x1, y0, y1;
    camera->film->GetPixelExtent(&x0, &x1, &y0, &y1);
    float t0 = camera->shutterOpen, t1 = camera->shutterClose;
    LDSampler sampler(x0, x1, y0, y1, nPixelSamples, t0, t1);
    Sample *sample = new Sample(&sampler, NULL, NULL, scene);

    // Create and launch _ProgressivePhotonMappingTask_s for rendering image

    // Compute number of _ProgressivePhotonMappingTask_s to create for rendering
    int nPixels = camera->film->xResolution * camera->film->yResolution;
    int nTasks = max(32 * NumSystemCores(), nPixels / (16*16));
    nTasks = RoundUpPow2(nTasks);
    ProgressReporter reporter(nTasks, "Ray Tracing Pass");
    vector<Task *> rayTracingTasks;
    vector<vector<HitPoint> > hitPoints(nTasks);
    MemoryArena *memoryArenas = new MemoryArena[nTasks];
    for (int i = 0; i < nTasks; ++i)
        rayTracingTasks.push_back(new RayTracingTask(scene, this, camera,
                                                     reporter, &sampler, sample,
                                                     hitPoints[i], memoryArenas[i],
                                                     maxDepth,
                                                     nTasks-1-i, nTasks));
    EnqueueTasks(rayTracingTasks);
    WaitForAllTasks();
    for (uint32_t i = 0; i < rayTracingTasks.size(); ++i)
        delete rayTracingTasks[i];
    reporter.Done();

    // Compute light power CDF for photon shooting
    Distribution1D *lightDistribution = ComputeLightSamplingCDF(scene);
    for (int iter = 0; iter < nIterations; ++iter) {
      if (scene->lights.size() == 0) break;
      ProgressReporter progress(nPhotonsPerIter, "Photon Tracing Pass");
      Mutex *mutex = Mutex::Create();
      vector<Photon> photons;
      photons.reserve(nPhotonsPerIter);

      vector<Task *> photonTracingTasks;
      int nTasksPhoton = NumSystemCores();
      int photonPerCore = nPhotonsPerIter / nTasksPhoton;
      int photonRes = nPhotonsPerIter - photonPerCore * nTasksPhoton;
      for (int i = 0; i < nTasksPhoton; ++i) {
        int want = photonPerCore;
        if (i == nTasksPhoton - 1) want += photonRes;
        photonTracingTasks.push_back(new PhotonTracingTask(
            i, camera->shutterOpen, *mutex, want, progress,
            photons, lightDistribution, scene, this));
      }
      EnqueueTasks(photonTracingTasks);
      WaitForAllTasks();
      for (uint32_t i = 0; i < photonTracingTasks.size(); ++i)
        delete photonTracingTasks[i];
      progress.Update();
      Mutex::Destroy(mutex);

      KdTree<Photon> *photonMap = NULL;
      if (photons.size() > 0)
        photonMap = new KdTree<Photon>(photons);
    }


    PBRT_FINISHED_RENDERING();
    // Clean up after rendering and store final image
    delete []memoryArenas;
    delete sample;
    camera->film->WriteImage();
}


Spectrum ProgressivePhotonMapping::Li(const Scene *scene,
        const RayDifferential &ray, const Sample *sample, RNG &rng,
        MemoryArena &arena, Intersection *isect, Spectrum *T) const {
    Severe("Unimplemented ProgressivePhotonMapping::Li() method called");
    return Spectrum(0.0);
}


Spectrum ProgressivePhotonMapping::Transmittance(const Scene *scene,
        const RayDifferential &ray, const Sample *sample, RNG &rng,
        MemoryArena &arena) const {
  return 1.0;
}


ProgressivePhotonMapping *
CreateProgressivePhotonMappingRenderer(Camera *camera,
                                       const ParamSet &params) {
    int perPixelSamples = params.FindOneInt("samplesperpixel", 64);
    int maxDepth = params.FindOneInt("maxdepth", 5);
    int nIterations = params.FindOneInt("iterations", 32);
    int x0, x1, y0, y1;
    camera->film->GetPixelExtent(&x0, &x1, &y0, &y1);
    int area = abs((x1-x0) * (y1-y0));
    int photonsPerIter = params.FindOneInt("photonsperiteration", area);
    return new ProgressivePhotonMapping(camera, perPixelSamples, maxDepth,
                                        nIterations, photonsPerIter);
}

