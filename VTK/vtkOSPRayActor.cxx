/* ======================================================================================= 
   Copyright 2014-2015 Texas Advanced Computing Center, The University of Texas at Austin  
   All rights reserved.
                                                                                           
   Licensed under the BSD 3-Clause License, (the "License"); you may not use this file     
   except in compliance with the License.                                                  
   A copy of the License is included with this software in the file LICENSE.               
   If your copy does not contain the License, you may obtain a copy of the License at:     
                                                                                           
       http://opensource.org/licenses/BSD-3-Clause                                         
                                                                                           
   Unless required by applicable law or agreed to in writing, software distributed under   
   the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY 
   KIND, either express or implied.                                                        
   See the License for the specific language governing permissions and limitations under   
   limitations under the License.

   pvOSPRay is derived from VTK/ParaView Los Alamos National Laboratory Modules (PVLANL)
   Copyright (c) 2007, Los Alamos National Security, LLC
   ======================================================================================= */

#include "ospray/ospray.h"

#define GL_GLEXT_PROTOTYPES

#include "vtkOSPRay.h"
#include "vtkOSPRayActor.h"
#include "vtkOSPRayManager.h"
#include "vtkOSPRayProperty.h"
#include "vtkOSPRayRenderer.h"
#include "vtkMapper.h"

#include "vtkDataSet.h"
#include "vtkObjectFactory.h"
#include "vtkRendererCollection.h"
#include "vtkTimerLog.h"

//VBO includes
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xatom.h>
#include <GL/glx.h>
#ifndef __APPLE__
#include <GL/glu.h>
#else
#include <OpenGL/glu.h>
#endif


#include <GL/glext.h>
#include "vtkOpenGL.h"
#include "vtkOpenGLError.h"
#include <map>
#include <algorithm>


#include "vtkInformation.h"
#include "vtkInformationVector.h"

//===========================================================================

vtkStandardNewMacro(vtkOSPRayActor);

//----------------------------------------------------------------------------
vtkOSPRayActor::vtkOSPRayActor()
{
  // std::cout << __PRETTY_FUNCTION__ << " " << this << std::endl;
  //cerr << "MA(" << this << ") CREATE" << endl;
  this->OSPRayManager = NULL;
  this->OSPRayModel = NULL;
}

//----------------------------------------------------------------------------
// now some OSPRay resources, ignored previously, can be de-allocated safely
//
vtkOSPRayActor::~vtkOSPRayActor()
{
  //cerr << "MA(" << this << ") DESTROY" << endl;
  if (this->OSPRayManager)
  {
    this->ReleaseGraphicsResources(NULL);
    //cerr << "MA(" << this << " DESTROY " << this->OSPRayManager << " "
    //     << this->OSPRayManager->GetReferenceCount() << endl;
    this->OSPRayManager->Delete();
  }
}

//----------------------------------------------------------------------------
void vtkOSPRayActor::PrintSelf( ostream & os, vtkIndent indent )
{
  this->Superclass::PrintSelf( os, indent );
}

//----------------------------------------------------------------------------
vtkProperty *vtkOSPRayActor::MakeProperty()
{
  return vtkOSPRayProperty::New();
}

//----------------------------------------------------------------------------
void vtkOSPRayActor::ReleaseGraphicsResources( vtkWindow * win )
{
}

//----------------------------------------------------------------------------
void vtkOSPRayActor::Render( vtkRenderer * ren, vtkMapper * mapper )
{
  // cerr << __PRETTY_FUNCTION__ << " " << this << endl;
  if ( vtkOSPRayRenderer * OSPRayRenderer = vtkOSPRayRenderer::SafeDownCast( ren ) )
  {
    if (!this->OSPRayManager)
    {
      this->OSPRayManager = OSPRayRenderer->GetOSPRayManager();
      //cerr << "MA(" << this << " REGISTER " << this->OSPRayManager << " "
      //     << this->OSPRayManager->GetReferenceCount() << endl;
      this->OSPRayManager->Register(this);
    }

    {
      mapper->Render(ren, this);
    }
      UpdateObjects(ren);
  }
}

//----------------------------------------------------------------------------
void vtkOSPRayActor::SetVisibility(int newval)
{
  //cerr << "MA(" << this << ") SET VISIBILITY " << newval << endl;
  if (newval == this->GetVisibility())
  {
    return;
  }
  this->Superclass::SetVisibility(newval);
}

//----------------------------------------------------------------------------
void vtkOSPRayActor::RemoveObjects()
{
}

//----------------------------------------------------------------------------
void vtkOSPRayActor::UpdateObjects( vtkRenderer * ren )
{
  #if 1
  // cerr << "MA(" << this << ") UPDATE" << endl;
  vtkOSPRayRenderer * OSPRayRenderer =
  vtkOSPRayRenderer::SafeDownCast( ren );
  if (!OSPRayRenderer)
  {
    return;
  }

  //Remove whatever we used to show in the scene
  if (!this->OSPRayManager)
  {
    return;
  }

  if (!this->OSPRayModel)
    return;

  //Remove what was shown.
  this->RemoveObjects();

  if (!this->GetVisibility())
    return;

  //Add what we are now supposed to show.
  #if 1
  {
    vtkTimerLog::MarkStartEvent("Execute AccelStructBuild ");
    OSPGeometry inst = ospNewInstance((((OSPModel)this->OSPRayModel)), osp::affine3f(embree::one));
      /*ospray::affine3f(embree::LinearSpace3f(mt(0,0), mt(0,1), mt(0,2), mt(1,0), mt(1,1), mt(1,2), mt(2,0),mt(2,1),mt(2,2)), embree::Vec3fa(mt(0,3),mt(1,3),mt(2,3))*/
    ospAddGeometry((((OSPModel)this->OSPRayManager->OSPRayModel)),inst);

    vtkTimerLog::MarkEndEvent("Execute AccelStructBuild ");
  }
  #endif
  #endif
}
