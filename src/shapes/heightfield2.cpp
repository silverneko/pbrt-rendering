
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


// shapes/heightfield2.cpp*
#include "stdafx.h"
#include "shapes/heightfield2.h"
#include "shapes/trianglemesh.h"
#include "paramset.h"

// Heightfield Method Definitions
Heightfield2::Heightfield2(const Transform *o2w, const Transform *w2o,
        bool ro, int x, int y, const float *zs)
    : Shape(o2w, w2o, ro) {
    nx = x;
    ny = y;
    nVoxels[0] = nx - 1;
    nVoxels[1] = ny - 1;
    for (int i = 0; i < 2; ++i) {
      width[i] = (nVoxels[i] == 0 ? 0.f : 1.0 / nVoxels[i]);
    }
    z = new float[nx*ny];
    memcpy(z, zs, nx*ny*sizeof(float));

    float minz = z[0], maxz = z[0];
    for (int i = 1; i < nx*ny; ++i) {
        if (z[i] < minz) minz = z[i];
        if (z[i] > maxz) maxz = z[i];
    }
    objectBounds = BBox(Point(0,0,minz), Point(1,1,maxz));

    {
      const int ntris = 2*(nx-1)*(ny-1);
      int *verts = new int[3*ntris];
      Point *P = new Point[nx*ny];
      Normal *N = new Normal[nx*ny];
      float *uvs = new float[2*nx*ny];
      int nverts = nx*ny;
      int x, y;
      // Compute heightfield vertex positions
      int pos = 0;
      for (y = 0; y < ny; ++y) {
          for (x = 0; x < nx; ++x) {
              P[pos].x = uvs[2*pos]   = (float)x / (float)(nx-1);
              P[pos].y = uvs[2*pos+1] = (float)y / (float)(ny-1);
              P[pos].z = z[pos];
              ++pos;
          }
      }

      // Compute heightfield vertex normals
      for (y = 0; y < ny-1; ++y) {
        for (x = 0; x < nx-1; ++x) {
#define VERT(x,y) ((x)+(y)*nx)
          const int x1 = x + 1, y1 = y;
          const int x2 = x + 1, y2 = y + 1;
          const int x3 = x, y3 = y + 1;
          Vector e1 = P[VERT(x1, y1)] - P[VERT(x, y)];
          Vector e2 = P[VERT(x2, y2)] - P[VERT(x, y)];
          Vector e3 = P[VERT(x3, y3)] - P[VERT(x, y)];
          Normal n1(Cross(e1, e2));
          Normal n2(Cross(e2, e3));
          N[VERT(x, y)] += n1 + n2;
          N[VERT(x1, y1)] += n1;
          N[VERT(x2, y2)] += n1 + n2;
          N[VERT(x3, y3)] += n2;
#undef VERT
        }
      }

      // Fill in heightfield vertex offset array
      int *vp = verts;
      for (y = 0; y < ny-1; ++y) {
          for (x = 0; x < nx-1; ++x) {
#define VERT(x,y) ((x)+(y)*nx)
              *vp++ = VERT(x, y);
              *vp++ = VERT(x+1, y);
              *vp++ = VERT(x+1, y+1);

              *vp++ = VERT(x, y);
              *vp++ = VERT(x+1, y+1);
              *vp++ = VERT(x, y+1);
          }
#undef VERT
      }
      ParamSet paramSet;
      paramSet.AddInt("indices", verts, 3*ntris);
      paramSet.AddFloat("uv", uvs, 2 * nverts);
      paramSet.AddPoint("P", P, nverts);
      paramSet.AddNormal("N", N, nverts);
      Reference<TriangleMesh> TM =
        CreateTriangleMeshShape(ObjectToWorld, WorldToObject,
                                ReverseOrientation, paramSet);
      triangle.reserve(ntris);
      TM->Refine(triangle);
      delete[] P;
      delete[] uvs;
      delete[] verts;
      delete[] N;
    }
}


Heightfield2::~Heightfield2() {
    delete[] z;
}


BBox Heightfield2::ObjectBound() const {
    return objectBounds;
}


bool Heightfield2::Intersect(const Ray &rayWorld, float *tHit,
               float *rayEpsilon, DifferentialGeometry *dg) const {
    // Check ray against overall grid bounds
    Ray ray;
    (*WorldToObject)(rayWorld, &ray);
    float rayT;
    if (objectBounds.Inside(ray(ray.mint))) {
        rayT = ray.mint;
    } else if (!objectBounds.IntersectP(ray, &rayT)) {
        return false;
    }
    Point gridIntersect = ray(rayT);

    // Set up 3D DDA for ray
    float NextCrossingT[2], DeltaT[2];
    int Step[2], Out[2], Pos[2];
    for (int axis = 0; axis < 2; ++axis) {
        // Compute current voxel for axis
        Pos[axis] = posToVoxel(gridIntersect, axis);
        if (ray.d[axis] >= 0) {
            // Handle ray with positive direction for voxel stepping
            NextCrossingT[axis] = rayT +
                (voxelToPos(Pos[axis]+1, axis) - gridIntersect[axis]) / ray.d[axis];
            DeltaT[axis] = width[axis] / ray.d[axis];
            Step[axis] = 1;
            Out[axis] = nVoxels[axis];
        }
        else {
            // Handle ray with negative direction for voxel stepping
            NextCrossingT[axis] = rayT +
                (voxelToPos(Pos[axis], axis) - gridIntersect[axis]) / ray.d[axis];
            DeltaT[axis] = -width[axis] / ray.d[axis];
            Step[axis] = -1;
            Out[axis] = -1;
        }
    }

    // Walk ray through voxel grid
    bool hitSomething = false;
    Ray rayW(rayWorld);
    for (;;) {
        // Check for intersection in current voxel and advance to next
        // Remember that triangle's xyz are in the world space
        int Offset = 2 * offset(Pos[0], Pos[1]);
        if (triangle[Offset]->Intersect(rayW, tHit, rayEpsilon, dg)) {
          hitSomething |= true;
          rayW.maxt = *tHit;
        }
        if (triangle[Offset+1]->Intersect(rayW, tHit, rayEpsilon, dg)) {
          hitSomething |= true;
          rayW.maxt = *tHit;
        }
        if (hitSomething) break;

        // Advance to next voxel

        // Find _stepAxis_ for stepping to next voxel
        int stepAxis = NextCrossingT[0] > NextCrossingT[1];
        if (rayW.maxt < NextCrossingT[stepAxis])
            break;
        Pos[stepAxis] += Step[stepAxis];
        if (Pos[stepAxis] == Out[stepAxis])
            break;
        NextCrossingT[stepAxis] += DeltaT[stepAxis];
    }
    return hitSomething;
}


bool Heightfield2::IntersectP(const Ray &rayWorld) const {
    // Check ray against overall grid bounds
    Ray ray;
    (*WorldToObject)(rayWorld, &ray);
    float rayT;
    if (objectBounds.Inside(ray(ray.mint))) {
        rayT = ray.mint;
    } else if (!objectBounds.IntersectP(ray, &rayT)) {
        return false;
    }
    Point gridIntersect = ray(rayT);

    // Set up 3D DDA for ray
    float NextCrossingT[2], DeltaT[2];
    int Step[2], Out[2], Pos[2];
    for (int axis = 0; axis < 2; ++axis) {
        // Compute current voxel for axis
        Pos[axis] = posToVoxel(gridIntersect, axis);
        if (ray.d[axis] >= 0) {
            // Handle ray with positive direction for voxel stepping
            NextCrossingT[axis] = rayT +
                (voxelToPos(Pos[axis]+1, axis) - gridIntersect[axis]) / ray.d[axis];
            DeltaT[axis] = width[axis] / ray.d[axis];
            Step[axis] = 1;
            Out[axis] = nVoxels[axis];
        }
        else {
            // Handle ray with negative direction for voxel stepping
            NextCrossingT[axis] = rayT +
                (voxelToPos(Pos[axis], axis) - gridIntersect[axis]) / ray.d[axis];
            DeltaT[axis] = -width[axis] / ray.d[axis];
            Step[axis] = -1;
            Out[axis] = -1;
        }
    }

    // Walk ray through voxel grid
    for (;;) {
        // Check for intersection in current voxel and advance to next
        // Remember that triangle's xyz are in the world space
        int Offset = 2 * offset(Pos[0], Pos[1]);
        if (triangle[Offset]->IntersectP(rayWorld)) {
          return true;
        }
        if (triangle[Offset+1]->IntersectP(rayWorld)) {
          return true;
        }

        // Advance to next voxel

        // Find _stepAxis_ for stepping to next voxel
        int stepAxis = NextCrossingT[0] > NextCrossingT[1];
        if (rayWorld.maxt < NextCrossingT[stepAxis])
            break;
        Pos[stepAxis] += Step[stepAxis];
        if (Pos[stepAxis] == Out[stepAxis])
            break;
        NextCrossingT[stepAxis] += DeltaT[stepAxis];
    }
    return false;
}


void Heightfield2::GetShadingGeometry(const Transform &obj2world,
        const DifferentialGeometry &dg,
        DifferentialGeometry *dgShading) const {
    dg.shape->GetShadingGeometry(obj2world, dg, dgShading);
    return;
}


Heightfield2 *CreateHeightfield2Shape(const Transform *o2w,
        const Transform *w2o, bool reverseOrientation, const ParamSet &params) {
    int nu = params.FindOneInt("nu", -1);
    int nv = params.FindOneInt("nv", -1);
    int nitems;
    const float *Pz = params.FindFloat("Pz", &nitems);
    Assert(nitems == nu*nv);
    Assert(nu != -1 && nv != -1 && Pz != NULL);
    return new Heightfield2(o2w, w2o, reverseOrientation, nu, nv, Pz);
}
