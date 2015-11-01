/* =======================================================================================
   Copyright 2014-2015 Texas Advanced Computing Center, The University of Texas
   at Austin
   All rights reserved.

   Licensed under the BSD 3-Clause License, (the "License"); you may not use
   this file
   except in compliance with the License.
   A copy of the License is included with this software in the file LICENSE.
   If your copy does not contain the License, you may obtain a copy of the
   License at:

       http://opensource.org/licenses/BSD-3-Clause

   Unless required by applicable law or agreed to in writing, software
   distributed under
   the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR
   CONDITIONS OF ANY
   KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under
   limitations under the License.

   pvOSPRay is derived from VTK/ParaView Los Alamos National Laboratory Modules
   (PVLANL)
   Copyright (c) 2007, Los Alamos National Security, LLC
   =======================================================================================
   */

#include "ospray/ospray.h"
#include "ospray/common/OSPCommon.h"

#include "vtkOSPRay.h"
#include "vtkOSPRayActor.h"
#include "vtkOSPRayManager.h"
#include "vtkOSPRayPolyDataMapper.h"
#include "vtkOSPRayProperty.h"
#include "vtkOSPRayRenderer.h"
#include "vtkOSPRayTexture.h"

#include "vtkAppendPolyData.h"
#include "vtkCellArray.h"
#include "vtkCellData.h"
#include "vtkCommand.h"
#include "vtkDataArray.h"
#include "vtkDoubleArray.h"
#include "vtkFloatArray.h"
#include "vtkGenericCell.h"
#include "vtkGlyph3D.h"
#include "vtkImageData.h"
#include "vtkMatrix4x4.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkPoints.h"
#include "vtkPolyData.h"
#include "vtkProperty.h"
#include "vtkScalarsToColors.h"
#include "vtkSphereSource.h"
#include "vtkTransform.h"
#include "vtkTubeFilter.h"
#include "vtkUnsignedCharArray.h"

#include "vtkInformation.h"
#include "vtkInformationVector.h"

#include <vector>

#include <math.h>
#include <algorithm>


vtkStandardNewMacro(vtkOSPRayPolyDataMapper);

namespace vtkosp {
class Vec3 {
  public:
  Vec3(float x, float y, float z) {
    vals[0] = x;
    vals[1] = y;
    vals[2] = z;
  }
  float operator[](unsigned int i) const { return vals[i]; }
  float x() { return vals[0]; }
  float y() { return vals[1]; }
  float z() { return vals[2]; }

  float vals[3];
};
class Vec4 {
  public:
  Vec4(float x, float y, float z, float w) {
    vals[0] = x;
    vals[1] = y;
    vals[2] = z;
    vals[3] = w;
  }
  float operator[](unsigned int i) const { return vals[i]; }
  float x() { return vals[0]; }
  float y() { return vals[1]; }
  float z() { return vals[2]; }
  float w() { return vals[3]; }

  float vals[4];
};
class Vec2 {
  public:
  Vec2(float x, float y) {
    vals[0] = x;
    vals[1] = y;
  }
  float operator[](unsigned int i) const { return vals[i]; }
  float x() { return vals[0]; }
  float y() { return vals[1]; }

  float vals[2];
};
class Mesh {
  public:
  size_t size() { return vertex_indices.size() / 3; }
  std::vector<size_t> vertex_indices;
  std::vector<Vec3> vertices;
  std::vector<Vec3> vertexNormals;
  std::vector<Vec2> texCoords;
  std::vector<size_t> texture_indices;
  std::vector<size_t> normal_indices;
  std::vector<Vec4> colors;
  std::vector<ospray::vec3fa> wireframe_vertex;
  std::vector<int> wireframe_index;
};
}

int vtkOSPRayPolyDataMapper::timestep = 0;  // HACK!

//----------------------------------------------------------------------------
// Construct empty object.
vtkOSPRayPolyDataMapper::vtkOSPRayPolyDataMapper() {
#if 1
  // cerr << "MM(" << this << ") CREATE" << endl;
  this->InternalColorTexture = NULL;
  this->OSPRayManager = NULL;
  this->PointSize = 1.0;
  this->LineWidth = 1.0;
  this->Representation = VTK_SURFACE;
#endif
}

//----------------------------------------------------------------------------
// Destructor (don't call ReleaseGraphicsResources() since it is virtual
vtkOSPRayPolyDataMapper::~vtkOSPRayPolyDataMapper() {
#if 1
  // cerr << "MM(" << this << ") DESTROY" << endl;
  if (this->InternalColorTexture) {
    this->InternalColorTexture->Delete();
  }

  if (this->OSPRayManager) {
    // cerr << "MM(" << this << ") DESTROY "
    //     << this->OSPRayManager << " "
    //     << this->OSPRayManager->GetReferenceCount() << endl;
    this->OSPRayManager->Delete();
  }

#endif
}

//----------------------------------------------------------------------------
// Release the graphics resources used by this mapper.  In this case, release
// the display list if any.
void vtkOSPRayPolyDataMapper::ReleaseGraphicsResources(vtkWindow *win) {
#if 1
  // cerr << "MM(" << this << ") RELEASE GRAPHICS RESOURCES" << endl;
  this->Superclass::ReleaseGraphicsResources(win);

  if (this->InternalColorTexture) {
    this->InternalColorTexture->Delete();
  }
  this->InternalColorTexture = NULL;
#endif
}

//----------------------------------------------------------------------------
// Receives from Actor -> maps data to primitives
// called by Mapper->Render() (which is called by Actor->Render())
void vtkOSPRayPolyDataMapper::RenderPiece(vtkRenderer *ren, vtkActor *act) {
#if 1
  vtkOSPRayRenderer *OSPRayRenderer = vtkOSPRayRenderer::SafeDownCast(ren);
  if (!OSPRayRenderer) {
    return;
  }
  if (!this->OSPRayManager) {
    this->OSPRayManager = OSPRayRenderer->GetOSPRayManager();
    // cerr << "MM(" << this << ") REGISTER " << this->OSPRayManager << " "
    //     << this->OSPRayManager->GetReferenceCount() << endl;
    this->OSPRayManager->Register(this);
  }

  // write geometry, first ask the pipeline to update data
  vtkPolyData *input = this->GetInput();
  // std::cerr << "polydata input: " << input << std::endl;
  if (input == NULL) {
    vtkErrorMacro(<< "No input to vtkOSPRayPolyDataMapper!");
    return;
  } else {
    this->InvokeEvent(vtkCommand::StartEvent, NULL);

    // Static = 1:  this mapper does NOT need to propagate updates to other
    // mappers
    // down the pipeline and therefore saves the time that would be otherwise
    // taken
    if (!this->Static) {
      this->Update();
    }

    this->InvokeEvent(vtkCommand::EndEvent, NULL);

    vtkIdType numPts = input->GetNumberOfPoints();
    if (numPts == 0) {
      vtkDebugMacro(<< "No points from the input to vtkOSPRayPolyDataMapper!");
      input = NULL;
      return;
    }
  }

  if (this->LookupTable == NULL) {
    this->CreateDefaultLookupTable();
  }

  // TODO: vtkOpenGLPolyDataMapper uses OpenGL clip planes here

  // For vertex coloring, this sets this->Colors as side effect.
  // For texture map coloring, this sets ColorCoordinates
  // and ColorTextureMap as a side effect.
  this->MapScalars(act->GetProperty()->GetOpacity());

  if (this->ColorTextureMap) {
    if (!this->InternalColorTexture) {
      this->InternalColorTexture = vtkOSPRayTexture::New();
      this->InternalColorTexture->RepeatOff();
    }
    this->InternalColorTexture->SetInputData(this->ColorTextureMap);
  }

  // if something has changed, regenerate OSPRay primitives if required
  if (this->GetMTime() > this->BuildTime ||
      input->GetMTime() > this->BuildTime ||
      // act->GetMTime()   > this->BuildTime ||
      act->GetProperty()->GetMTime() > this->BuildTime ||
      act->GetMatrix()->GetMTime() > this->BuildTime) {
    // this->ReleaseGraphicsResources( ren->GetRenderWindow() );

    // If we are coloring by texture, then load the texture map.
    // Use Map as indicator, because texture hangs around.
    if (this->ColorTextureMap) {
      this->InternalColorTexture->Load(ren);
    }

    this->Draw(ren, act);
    this->BuildTime.Modified();
  }

// input->Modified();
this->Update();
  input = NULL;
// TODO: deal with timer ??
#endif
}

//----------------------------------------------------------------------------
void vtkOSPRayPolyDataMapper::DrawPolygons(vtkPolyData *polys,
                                           vtkPoints *ptarray,
                                           vtkosp::Mesh *mesh
                                           /*,
                                          OSPRay::Group *points,
                                          OSPRay::Group *lines*/) {
  // std::cout << __PRETTY_FUNCTION__ << std::endl;

  int total_triangles = 0;
  vtkCellArray *cells = polys->GetPolys();
  vtkIdType npts = 0, *index = 0, cellNum = 0;

  // cerr << this->Representation << endl;
  switch (this->Representation) {
    case VTK_POINTS: {
      std::cout << "VTK_POINTS\n";
      for (cells->InitTraversal(); cells->GetNextCell(npts, index); cellNum++) {
        total_triangles++;
      }
      // cerr << "polygons: # of triangles = " << total_triangles << endl;
    }  // VTK_POINTS;
    break;
    case VTK_WIREFRAME: {
      // std::cout << "VTK_WIREFRAME\n";
      double coord0[3];
      for (cells->InitTraversal(); cells->GetNextCell(npts, index); cellNum++) {
        ptarray->GetPoint(index[0], coord0);
        mesh->wireframe_vertex.push_back(
            ospray::vec3fa(coord0[0], coord0[1], coord0[2]));
        for (vtkIdType i = 1; i < npts; i++) {

          // std::cout << "points in cell: " << npts << "\n";
          mesh->wireframe_index.push_back(mesh->wireframe_vertex.size() - 1);
          ptarray->GetPoint(index[i], coord0);
          mesh->wireframe_vertex.push_back(
              ospray::vec3fa(coord0[0], coord0[1], coord0[2]));
        }
      }
    }  // VTK_WIREFRAME:
    break;
    case VTK_SURFACE: {
      // std::cout << "VTK_SURFACE\n";
      // write polygons with on the fly triangulation, assuming polygons are
      // simple and
      // can be triangulated into "fans"
      for (cells->InitTraversal(); cells->GetNextCell(npts, index); cellNum++) {
        int triangle[3];

        // the first triangle
        triangle[0] = index[0];
        triangle[1] = index[1];
        triangle[2] = index[2];
        mesh->vertex_indices.push_back(triangle[0]);
        mesh->vertex_indices.push_back(triangle[1]);
        mesh->vertex_indices.push_back(triangle[2]);
        // mesh->face_material.push_back(0);

        if (!mesh->vertexNormals.empty()) {
          mesh->normal_indices.push_back(triangle[0]);
          mesh->normal_indices.push_back(triangle[1]);
          mesh->normal_indices.push_back(triangle[2]);
        }

        if (!mesh->texCoords.empty()) {
          if (this->CellScalarColor) {
            mesh->texture_indices.push_back(cellNum);
            mesh->texture_indices.push_back(cellNum);
            mesh->texture_indices.push_back(cellNum);
          } else {
            mesh->texture_indices.push_back(triangle[0]);
            mesh->texture_indices.push_back(triangle[1]);
            mesh->texture_indices.push_back(triangle[2]);
          }
        }
        total_triangles++;

        // the remaining triangles, of which
        // each introduces a triangle after extraction
        for (int i = 3; i < npts; i++) {
          triangle[1] = triangle[2];
          triangle[2] = index[i];
          mesh->vertex_indices.push_back(triangle[0]);
          mesh->vertex_indices.push_back(triangle[1]);
          mesh->vertex_indices.push_back(triangle[2]);
          // mesh->face_material.push_back(0);

          if (!mesh->vertexNormals.empty()) {
            mesh->normal_indices.push_back(triangle[0]);
            mesh->normal_indices.push_back(triangle[1]);
            mesh->normal_indices.push_back(triangle[2]);
          }

          if (!mesh->texCoords.empty()) {
            if (this->CellScalarColor) {
              mesh->texture_indices.push_back(cellNum);
              mesh->texture_indices.push_back(cellNum);
              mesh->texture_indices.push_back(cellNum);
            } else {
              mesh->texture_indices.push_back(triangle[0]);
              mesh->texture_indices.push_back(triangle[1]);
              mesh->texture_indices.push_back(triangle[2]);
            }
          }
          total_triangles++;
        }
      }
      // cerr << "polygons: # of triangles = " << total_triangles << endl;

      for (cells->InitTraversal();
           this->Edges && cells->GetNextCell(npts, index); cellNum++) {
        double coord0[3];
        ptarray->GetPoint(index[0], coord0);
        mesh->wireframe_vertex.push_back(
            ospray::vec3fa(coord0[0], coord0[1], coord0[2]));
        for (vtkIdType i = 1; i < npts; i++) {
          mesh->wireframe_index.push_back(mesh->wireframe_vertex.size() - 1);
          ptarray->GetPoint(index[i], coord0);
          mesh->wireframe_vertex.push_back(
              ospray::vec3fa(coord0[0], coord0[1], coord0[2]));
        }
      }

    }  // VTK_SURFACE
    break;
    default:
      std::cerr << "unknwon representation type\n";
      break;
  }
}

//----------------------------------------------------------------------------
void vtkOSPRayPolyDataMapper::DrawTStrips(vtkPolyData *polys,
                                          vtkPoints *ptarray,
                                          vtkosp::Mesh *mesh)
    // OSPRay::Mesh *mesh,
    // OSPRay::Group *points,
    // OSPRay::Group *lines)
{
  // std::cout << __PRETTY_FUNCTION__ << std::endl;

  // total number of triangles
  int total_triangles = 0;

  vtkCellArray *cells = polys->GetStrips();
  vtkIdType npts = 0, *index = 0, cellNum = 0;
  ;

  // cerr << this->Representation << endl;
  switch (this->Representation) {
    case VTK_POINTS: {
      // std::cout << "VTK_POINTS\n";
      for (cells->InitTraversal(); cells->GetNextCell(npts, index); cellNum++) {
        total_triangles++;
      }
      // cerr << "polygons: # of triangles = " << total_triangles << endl;
    }  // VTK_POINTS;
    break;
    case VTK_WIREFRAME: {
      // std::cout << "VTK_WIREFRAME\n";
      double coord0[3];
      double coord1[3];
      double coord2[3];
      // OSPRay::Vector noTC(0.0,0.0,0.0);
      // OSPRay::TextureCoordinateCylinder *segment;
      for (cells->InitTraversal(); cells->GetNextCell(npts, index); cellNum++) {
        ptarray->GetPoint(index[0], coord0);
        mesh->wireframe_vertex.push_back(
            ospray::vec3fa(coord0[0], coord0[1], coord0[2]));
        for (vtkIdType i = 2; i < npts; i++) {
          mesh->wireframe_index.push_back(mesh->wireframe_vertex.size() - 1);
          ptarray->GetPoint(index[0], coord0);
          mesh->wireframe_vertex.push_back(
              ospray::vec3fa(coord0[0], coord0[1], coord0[2]));
        }
      }
    }  // VTK_WIREFRAME:
    break;
    case VTK_SURFACE: {
      // std::cout << "VTK_SURFACE\n";
      for (cells->InitTraversal(); cells->GetNextCell(npts, index); cellNum++) {
        // count of the i-th triangle in a strip
        int numtriangles2 = 0;

        int triangle[3];
        // the first triangle
        triangle[0] = index[0];
        triangle[1] = index[1];
        triangle[2] = index[2];
        mesh->vertex_indices.push_back(triangle[0]);
        mesh->vertex_indices.push_back(triangle[1]);
        mesh->vertex_indices.push_back(triangle[2]);
        // mesh->face_material.push_back(0);

        if (!mesh->vertexNormals.empty()) {
          mesh->normal_indices.push_back(triangle[0]);
          mesh->normal_indices.push_back(triangle[1]);
          mesh->normal_indices.push_back(triangle[2]);
        }

        if (!mesh->texCoords.empty()) {
          if (this->CellScalarColor) {
            mesh->texture_indices.push_back(cellNum);
            mesh->texture_indices.push_back(cellNum);
            mesh->texture_indices.push_back(cellNum);
          } else {
            mesh->texture_indices.push_back(triangle[0]);
            mesh->texture_indices.push_back(triangle[1]);
            mesh->texture_indices.push_back(triangle[2]);
          }
        }

        total_triangles++;
        numtriangles2++;

        // the rest of triangles
        for (int i = 3; i < npts; i++) {
          int tmp[3];
          if (numtriangles2 % 2 == 1) {
            // an odd triangle
            tmp[0] = triangle[1];
            tmp[1] = triangle[2];
            tmp[2] = index[i];

            triangle[0] = tmp[0];
            triangle[1] = tmp[2];
            triangle[2] = tmp[1];
          } else {
            // an even triangle
            tmp[0] = triangle[1];
            tmp[1] = triangle[2];
            tmp[2] = index[i];

            triangle[0] = tmp[1];
            triangle[1] = tmp[0];
            triangle[2] = tmp[2];
          }

          mesh->vertex_indices.push_back(triangle[0]);
          mesh->vertex_indices.push_back(triangle[1]);
          mesh->vertex_indices.push_back(triangle[2]);
          // mesh->face_material.push_back(0);

          if (!mesh->vertexNormals.empty()) {
            mesh->normal_indices.push_back(triangle[0]);
            mesh->normal_indices.push_back(triangle[1]);
            mesh->normal_indices.push_back(triangle[2]);
          }

          if (!mesh->texCoords.empty()) {
            if (this->CellScalarColor) {
              mesh->texture_indices.push_back(cellNum);
              mesh->texture_indices.push_back(cellNum);
              mesh->texture_indices.push_back(cellNum);
            } else {
              mesh->texture_indices.push_back(triangle[0]);
              mesh->texture_indices.push_back(triangle[1]);
              mesh->texture_indices.push_back(triangle[2]);
            }
          }

          total_triangles++;
          numtriangles2++;
        }
      }

    } break;

    default:
      std::cerr << "unkown representation type\n";
  }
}

//----------------------------------------------------------------------------
// Draw method for OSPRay.
void vtkOSPRayPolyDataMapper::Draw(vtkRenderer *renderer, vtkActor *actor) {
#if 1
  // printf("ospPolyDataMapper::Draw\n");
  vtkOSPRayActor *OSPRayActor = vtkOSPRayActor::SafeDownCast(actor);
  if (!OSPRayActor) {
    return;
  }
  vtkOSPRayProperty *OSPRayProperty =
      vtkOSPRayProperty::SafeDownCast(OSPRayActor->GetProperty());
  if (!OSPRayProperty) {
    return;
  }
  vtkPolyData *input = this->GetInput();

  vtkInformation *inputInfo = this->GetInput()->GetInformation();
  // // vtkInformation* outputInfo = outputVector->GetInformationObject(0);

  // std::cerr << "ospPDM Actor: " << actor << std::endl;
  if (inputInfo && inputInfo->Has(vtkDataObject::DATA_TIME_STEP())) {
    double time = inputInfo->Get(vtkDataObject::DATA_TIME_STEP());
    // cerr << "MA time: " << time << std::endl;
    timestep = time;
    if (OSPRayActor->cache[time] != NULL) {
      // std::cerr << "using cache at time " << time << "\n";
      // this->OSPRayModel = cache[time];

      OSPRayActor->OSPRayModel = OSPRayActor->cache[time];
      return;

      // this->MeshMTime.Modified();
      // UpdateObjects(ren);
    }
    // return;

  } else if (!inputInfo) {
    // cerr << "MA time: didn't have info\n";
  } else {
    // cerr << "MA time: didn't have time\n";
  }
  OSPRayActor->MeshMTime.Modified();

  // Compute we need to for color
  this->Representation = OSPRayProperty->GetRepresentation();
  this->Edges = OSPRayProperty->GetEdgeVisibility();

  this->CellScalarColor = false;
  if ((this->ScalarMode == VTK_SCALAR_MODE_USE_CELL_DATA ||
       this->ScalarMode == VTK_SCALAR_MODE_USE_CELL_FIELD_DATA ||
       this->ScalarMode == VTK_SCALAR_MODE_USE_FIELD_DATA ||
       !input->GetPointData()->GetScalars()) &&
      this->ScalarMode != VTK_SCALAR_MODE_USE_POINT_FIELD_DATA) {
    this->CellScalarColor = true;
  }

  OSPMaterial ospMaterial = NULL;
  vtkosp::Mesh *mesh = new vtkosp::Mesh();

  // force create a new material every time this is called in case the renderer
  // has changed
  // osp::Material* osmat = OSPRayProperty->GetOSPRayMaterial();
  osp::Material *osmat = 0;
  if (!osmat) {
    OSPRayProperty->CreateOSPRayProperty();
    ospMaterial = ((OSPMaterial)OSPRayProperty->GetOSPRayMaterial());
  } else
    ospMaterial = ((OSPMaterial)osmat);

  if (!this->ScalarVisibility || (!this->Colors && !this->ColorCoordinates)) {
    // cerr << "poly colors: Solid color from actor's property" << endl;
  } else if (this->Colors) {
    // cerr << "poly colors: Color scalar values directly (interpolation in color "
            // "space)" << endl;
    // OSPRay::Texture<OSPRay::Color> *texture = new OSPRay::TexCoordTexture();

    // this table has one RGBA for every point (or cell) in object
    for (int i = 0; i < this->Colors->GetNumberOfTuples(); i++) {
      unsigned char *color = this->Colors->GetPointer(4 * i);
      // texCoords.push_back
      // (OSPRay::Vector(color[0]/255.0, color[1]/255.0, color[2]/255.0) );
      // mesh->texCoords.push_back(vtkosp::Vec3(color[0]/255.0, color[1]/255.0,
      // color[2]/255.0));
      mesh->colors.push_back(vtkosp::Vec4(color[0] / 255.0, color[1] / 255.0,
                                          color[2] / 255.0, 1));
    }
    // printf("texture coords: using rgba every point\n");

  } else if (this->ColorCoordinates) {
    // printf(
        // "poly colors: texture coords: using color coordinates for a texture\n");
    // cerr << "interpolate in data space, then color map each pixel" << endl;
    // OSPRay::Texture<OSPRay::Color> *texture =
    //   this->InternalColorTexture->GetOSPRayTexture();
    osp::Texture2D *texture = this->InternalColorTexture->GetOSPRayTexture();
    // PRINT((OSPTexture2D)texture);
    Assert(texture);
    ospSetParam(ospMaterial, "map_Kd", ((OSPTexture2D)(texture)));
    ospCommit(ospMaterial);

    // material = new OSPRay::Lambertian(texture);

    // //this table is a color transfer function with colors that cover the
    // scalar range
    // //I think
    for (int i = 0; i < this->ColorCoordinates->GetNumberOfTuples(); i++) {
      double *tcoord = this->ColorCoordinates->GetTuple(i);
      //   texCoords.push_back( OSPRay::Vector(tcoord[0], 0, 0) );
      // mesh->texCoords.push_back(vtkosp::Vec2(tcoord[0], 0));
			if (tcoord[0] >= 1.0) tcoord[0] = 0.99999;	// avoid sampling texture at 1
      mesh->texCoords.push_back(vtkosp::Vec2(tcoord[0], tcoord[1]));
      // texCoords.push_back(vtkosp::Vec2(tcoord[0],0));
      // mesh->colors.push_back(vtkosp::Vec4(color[0]/255.0,color[1]/255.0,color[2]/255.0,1));
      // printf("texCoord: %f %f\n", tcoord[0], 0);
    }

    // printf("NEED TO IMPLEMENT COLORCOORDINATES\n");
  } else if (input->GetPointData()->GetTCoords() && actor->GetTexture()) {
    // printf("poly colors: texture coords: using texture\n");
#if 1
    // cerr << "color using actor's texture" << endl;
    vtkOSPRayTexture *osprayTexture =
        vtkOSPRayTexture::SafeDownCast(actor->GetTexture());
    if (osprayTexture) {
      //   OSPRay::Texture<OSPRay::Color> *texture =
      //     OSPRayTexture->GetOSPRayTexture();
      //   material = new OSPRay::Lambertian(texture);
      ospSetParam(ospMaterial, "map_Kd",
                  ((OSPTexture2D)(osprayTexture->GetOSPRayTexture())));
      ospCommit(ospMaterial);
    }

    // // convert texture coordinates to OSPRay format
    vtkDataArray *tcoords = input->GetPointData()->GetTCoords();
    for (int i = 0; i < tcoords->GetNumberOfTuples(); i++) {
      double *tcoord = tcoords->GetTuple(i);
			if (tcoord[0] >= 1.0) tcoord[0] = 0.99999;	// avoid sampling texture at 1
      mesh->texCoords.push_back(vtkosp::Vec2(tcoord[0], tcoord[1]));
      //     ( OSPRay::Vector(tcoord[0], tcoord[1], tcoord[2]) );
    }
#endif
    // printf("NEED TO IMPLEMENT TEXTURES\n");
  }

  // transform point coordinates according to actor's transformation matrix
  // TODO: Use OSPRay instancing to transform instead of doing it brute force
  // here
  // to reduce number of copies
  vtkTransform *transform = vtkTransform::New();
  transform->SetMatrix(actor->GetMatrix());
  vtkPoints *points = vtkPoints::New();
  transform->TransformPoints(input->GetPoints(), points);

  // obtain the OpenGL-based point size and line width
  // that are specified through vtkProperty
  this->PointSize = OSPRayProperty->GetPointSize();
  this->LineWidth = OSPRayProperty->GetLineWidth();
  if (this->PointSize < 0.0) {
    this->PointSize = 1.0;
  }
  if (this->LineWidth < 0.0) {
    this->LineWidth = 1.0;
  }
  this->PointSize = sqrt(this->PointSize) * 0.010;
  this->LineWidth = sqrt(this->LineWidth) * 0.005;

  std::vector<ospray::vec3fa> slVertex;
  std::vector<ospray::vec3fa> slColors;
  std::vector<int> slIndex;
  float slRadius;

  // convert VTK_LINE type cells to OSPRay cylinders
  if (input->GetNumberOfLines() > 0) {
    vtkCellArray *ca = input->GetLines();
    ca->InitTraversal();
    vtkIdType npts;
    vtkIdType *pts;
    vtkPoints *ptarray = points;
    std::vector<ospray::vec3fa> tmpColors;
    int scalarSize = ptarray->GetNumberOfPoints();
    double coord0[3];
    vtkIdType cell;
    // output = new unsigned char[scalarSize * 4];

    // this->ColorTextureMap->GetLoop

    vtkScalarsToColors *vstc = GetLookupTable();
    vtkDataArray *scalar = input->GetPointData()->GetScalars(NULL);
    if (scalar) {
      int vectorSize = (scalar) ? scalar->GetNumberOfComponents() : 0;
      unsigned char *output = new unsigned char[scalarSize * 4];
      vstc->SetVectorModeToMagnitude();
      vstc->MapVectorsThroughTable(scalar->GetVoidPointer(0), output,
                                   scalar->GetDataType(), scalarSize,
                                   vectorSize, VTK_RGBA);
      for (int ii = 0; ii < scalarSize; ii++) {
        double color[3];
        for (int jj = 0; jj < 3; jj++) {
          color[jj] = float(output[ii * 4 + jj]) / 255.0;
        }
        tmpColors.push_back(ospray::vec3fa(color[0], color[1], color[2]));
      }
    } else if (vstc) {
      std::cout << "I have a lookup table" << std::endl;

      if (this->ColorCoordinates) {
        std::cout << "Tex coords " << this->ColorCoordinates->GetSize() << std::endl;
        double* minmax = vstc->GetRange();
        std::cout << "m: " << minmax[0] << " M:" << minmax[1] << std::endl;
        //vstc->SetRange(0,1);
        
        double scale = minmax[1] - minmax[0];

        for (int i = 0; i < scalarSize; i++) {
          double *tcoord = this->ColorCoordinates->GetTuple(i);
          double *color = vstc->GetColor((tcoord[0] * scale) + minmax[0]);
          tmpColors.push_back(ospray::vec3fa(color[0], color[1], color[2]));
        }

      } else {
        double solidColor[3];
        OSPRayProperty->GetDiffuseColor(solidColor);
        for (int i = 0; i < scalarSize; i++) {
          tmpColors.push_back(
              ospray::vec3fa(solidColor[0], solidColor[1], solidColor[2]));
        }
      }
    }

    std::vector<ospray::vec3fa> tmpPoints;
    for (int ii = 0; ii < ptarray->GetNumberOfPoints(); ii++) {
      ptarray->GetPoint(ii, coord0);
      tmpPoints.push_back(ospray::vec3fa(coord0[0], coord0[1], coord0[2]));
    }

    slRadius = this->LineWidth / 0.005;

    while ((cell = ca->GetNextCell(npts, pts))) {
      if (npts <= 2) continue;
      slVertex.push_back(tmpPoints[pts[0]]);
      slColors.push_back(tmpColors[pts[0]]);
      for (vtkIdType i = 1; i < npts; i++) {
        slIndex.push_back(slVertex.size() - 1);
        slVertex.push_back(tmpPoints[pts[i]]);
        slColors.push_back(tmpColors[pts[i]]);
      }
    }
  }

  // convert coordinates to OSPRay format
  // TODO: eliminate the copy
  for (int i = 0; i < points->GetNumberOfPoints(); i++) {
    double *pos = points->GetPoint(i);
    bool wasNan = false;
    int fixIndex = i - 1;
    do {
      wasNan = false;
      for (int j = 0; j < 3; j++) {
        if (std::isnan(pos[j])) {
          wasNan = true;
        }
      }
      if (wasNan && fixIndex >= 0) pos = points->GetPoint(fixIndex--);
    } while (wasNan == true && fixIndex >= 0);
    mesh->vertices.push_back(vtkosp::Vec3(pos[0], pos[1], pos[2]));
  }

  // Do flat shading by not supplying vertex normals to OSPRay
  if (OSPRayProperty->GetInterpolation() != VTK_FLAT) {
    vtkPointData *pointData = input->GetPointData();
    if (pointData->GetNormals()) {
      vtkDataArray *normals = vtkFloatArray::New();
      normals->SetNumberOfComponents(3);
      transform->TransformNormals(pointData->GetNormals(), normals);
      for (int i = 0; i < normals->GetNumberOfTuples(); i++) {
        double *normal = normals->GetTuple(i);
        mesh->vertexNormals.push_back(
            vtkosp::Vec3(normal[0], normal[1], normal[2]));
      }
      normals->Delete();
    }
  }

  // convert polygons to OSPRay format
  if (input->GetNumberOfPolys() > 0) {
    this->DrawPolygons(input, points, mesh /*, sphereGroup, tubeGroup*/);
  }

  // convert triangle strips to OSPRay format
  if (input->GetNumberOfStrips() > 0) {
    this->DrawTStrips(input, points, mesh /*, sphereGroup, tubeGroup*/);
  }

  // delete transformed point coordinates
  transform->Delete();
  points->Delete();

  if (mesh->size() || mesh->wireframe_vertex.size() || slVertex.size()) {

#if USE_OSPRAY
    OSPRenderer renderer = ((OSPRenderer) this->OSPRayManager->OSPRayRenderer);
    OSPRayActor->OSPRayModel = ospNewModel();

    if (mesh->size()) {

      size_t numNormals = mesh->vertexNormals.size();
      size_t numTexCoords = mesh->texCoords.size();
      size_t numPositions = mesh->vertices.size();
      size_t numTriangles = mesh->vertex_indices.size() / 3;

      // std::cout << "normals: " << numNormals
      //           << " normal indices: " << mesh->normal_indices.size()
      //           << " numPositions: " << numPositions
      //           << " numTriangles: " << numTriangles << std::endl;

      ospray::vec3fa *vertices = (ospray::vec3fa *)embree::alignedMalloc(
          sizeof(ospray::vec3fa) * numPositions);
      ospray::vec3i *triangles = (ospray::vec3i *)embree::alignedMalloc(
          sizeof(ospray::vec3i) * numTriangles);
      ospray::vec3fa *normals = (ospray::vec3fa *)embree::alignedMalloc(
          sizeof(ospray::vec3fa) * numNormals);

      for (size_t i = 0; i < numPositions; i++) {
        vertices[i] =
            ospray::vec3fa(mesh->vertices[i].x(), mesh->vertices[i].y(),
                           mesh->vertices[i].z());
        // vertices[i] = ospray::vec3fa(float(i)*.01,float(i)*.01,float(i)*.01);
        // printf("vert: %f %f %f\n",mesh->vertices[i].x(),
        // mesh->vertices[i].y(), mesh->vertices[i].z());
      }
      for (size_t i = 0, mi = 0; i < numTriangles; i++, mi += 3) {
        triangles[i] = embree::Vec3i(mesh->vertex_indices[mi + 0],
                                     mesh->vertex_indices[mi + 1],
                                     mesh->vertex_indices[mi + 2]);
      }

      for (size_t i = 0; i < numNormals; i++) {
        normals[i] = ospray::vec3fa(mesh->vertexNormals[i].x(),
                                    mesh->vertexNormals[i].y(),
                                    mesh->vertexNormals[i].z());
      }

      OSPGeometry ospMesh = ospNewTriangleMesh();
      OSPData position = ospNewData(numPositions, OSP_FLOAT3A, &vertices[0]);
      ospSetData(ospMesh, "position", position);

      if (!mesh->normal_indices.empty()) {
        OSPData normal =
            ospNewData(mesh->vertexNormals.size(), OSP_FLOAT3A, &normals[0]);
        ospSetData(ospMesh, "vertex.normal", normal);
      }

      OSPData index = ospNewData(numTriangles, OSP_INT3, &triangles[0]);
      ospSetData(ospMesh, "index", index);

      if (!mesh->texCoords.empty()) {
        OSPData texcoord =
            ospNewData(mesh->texCoords.size(), OSP_FLOAT2, &mesh->texCoords[0]);
        assert(mesh->texCoords.size() > 0);
        ospSetData(ospMesh, "vertex.texcoord", texcoord);
      }
      if (!mesh->colors.empty()) {
        std::cerr << "using color coordinates\n";
        // note: to share data use OSP_DATA_SHARED_BUFFER
        OSPData colors =
            ospNewData(mesh->colors.size(), OSP_FLOAT4, &mesh->colors[0]);
        ospSetData(ospMesh, "vertex.color", colors);
      }

      if (!ospMaterial) {
        // printf("no material specification used, using default\n");
        OSPRayProperty->CreateOSPRayProperty();
        ospMaterial = ((OSPMaterial)OSPRayProperty->GetOSPRayMaterial());
      }
      // PRINT(ospMaterial);

      ospSetMaterial(ospMesh, ospMaterial);
      ospCommit(ospMesh);

      // static OSPModel ospModel;
      // OSPModel ospModel = ospNewModel();
      ospAddGeometry(OSPRayActor->OSPRayModel, ospMesh);
    }

    if (mesh->wireframe_vertex.size()) {
      double edgeColor[3];
      OSPRayProperty->GetEdgeColor(edgeColor);
      OSPMaterial wireMat = ospNewMaterial(renderer, "default");
      if (wireMat) {
        ospSet3f(wireMat, "kd", edgeColor[0], edgeColor[1], edgeColor[2]);
        ospCommit(wireMat);
      }
      OSPGeometry wireGeometry = ospNewGeometry("streamlines");
      Assert(wireGeometry);
      OSPData vertex = ospNewData(mesh->wireframe_vertex.size(), OSP_FLOAT3A,
                                  &mesh->wireframe_vertex[0]);
      OSPData index = ospNewData(mesh->wireframe_index.size(), OSP_INT,
                                 &mesh->wireframe_index[0]);
      ospSetObject(wireGeometry, "vertex", vertex);
      ospSetObject(wireGeometry, "index", index);
      ospSet1f(wireGeometry, "radius", this->LineWidth);

      if (wireMat) ospSetMaterial(wireGeometry, wireMat);

      ospCommit(wireGeometry);
      ospAddGeometry(OSPRayActor->OSPRayModel, wireGeometry);
    }

    if (slVertex.size()) {
      double solidColor[3];
      OSPGeometry slGeometry = ospNewGeometry("streamlines");
      Assert(slGeometry);
      OSPData vertex = ospNewData(slVertex.size(), OSP_FLOAT3A, &slVertex[0]);
      OSPData color = ospNewData(slColors.size(), OSP_FLOAT3A, &slColors[0]);
      OSPData index = ospNewData(slIndex.size(), OSP_INT, &slIndex[0]);
      ospSetObject(slGeometry, "vertex", vertex);
      ospSetObject(slGeometry, "vertex.color", color);
      ospSetObject(slGeometry, "index", index);
      ospSet1f(slGeometry, "radius", slRadius);
      ospCommit(slGeometry);
      ospAddGeometry(OSPRayActor->OSPRayModel, slGeometry);

      // std::cerr << " Commited geometry correctly " << std::endl;
    }

    // std::cerr << "Trying to commit model" << timestep << "\n";
    ospCommit(OSPRayActor->OSPRayModel);
    if (inputInfo && inputInfo->Has(vtkDataObject::DATA_TIME_STEP())) {
      double time = inputInfo->Get(vtkDataObject::DATA_TIME_STEP());
      OSPRayActor->cache[time] = OSPRayActor->OSPRayModel;
    } else {
      OSPRayActor->cache[timestep] = OSPRayActor->OSPRayModel;
      // std::cerr << "added nontime actor at timestep" << timestep << "\n";
    }

#endif

  } else {
    delete mesh;
  }

#endif
}
