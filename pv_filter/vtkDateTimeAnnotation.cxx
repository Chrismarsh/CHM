/*=========================================================================

  Program:   ParaView
  Module:    $RCSfile: vtkJulianDayToTextConvertor.cxx,v $

  Copyright (c) Kitware, Inc.
  All rights reserved.
  See Copyright.txt or http://www.paraview.org/HTML/Copyright.html for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkDateTimeAnnotation.h"

#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkStringArray.h"
#include "vtkTable.h"

vtkStandardNewMacro(vtkDateTimeAnnotation);

//----------------------------------------------------------------------------
vtkDateTimeAnnotation::vtkDateTimeAnnotation()
{
  this->Format = 0;
  this->Shift  = 0.0;
  this->Scale  = 1.0;

  this->SetFormat("Time: %04i/%02i/%02i %02i:%02i:%02i");  // YYYY/MM/DD HH:MM:SS
}

//----------------------------------------------------------------------------
vtkDateTimeAnnotation::~vtkDateTimeAnnotation()
{
  this->SetFormat(0);
}

//----------------------------------------------------------------------------
int vtkDateTimeAnnotation::FillInputPortInformation(
  int vtkNotUsed(port), vtkInformation* info)
{
  info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataObject");
  info->Set(vtkAlgorithm::INPUT_IS_OPTIONAL(), 1);
  return 1;
}

//----------------------------------------------------------------------------
int vtkDateTimeAnnotation::RequestInformation(
  vtkInformation *request,
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector)
{
  if (!this->Superclass::RequestInformation(
      request, inputVector, outputVector))
    {
    return 0;
    }
  double timeRange[2];
  timeRange[0] = VTK_DOUBLE_MIN;
  timeRange[1] = VTK_DOUBLE_MAX;

  vtkInformation* outInfo = outputVector->GetInformationObject(0);
  outInfo->Set(vtkStreamingDemandDrivenPipeline::TIME_RANGE(), timeRange, 2);
  return 1;
}
//----------------------------------------------------------------------------
inline double vtkDateTimeAnnotation_ForwardConvert(double T0, double shift, double scale)
{
  return T0*scale + shift;
}
//----------------------------------------------------------------------------
int vtkDateTimeAnnotation::RequestData(
  vtkInformation* vtkNotUsed(request),
  vtkInformationVector** inputVector,
  vtkInformationVector* outputVector)
{
  vtkDataObject* input = vtkDataObject::GetData(inputVector[0]);
  vtkTable* output = vtkTable::GetData(outputVector);

  char *buffer = new char[strlen(this->Format)+1024];
  strcpy(buffer, "?");

  vtkInformation* inputInfo = input? input->GetInformation() : 0;
  vtkInformation* outputInfo = outputVector->GetInformationObject(0);

  double time = 0;
  if (inputInfo && inputInfo->Has(vtkDataObject::DATA_TIME_STEP())
      && this->Format)
  {
    time = inputInfo->Get(vtkDataObject::DATA_TIME_STEP());//[0];

  }
  else if (outputInfo &&
           outputInfo->Has(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP())
           && this->Format)
  {
    time = outputInfo->Get(
            vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP());//[0];

  }

  time = vtkDateTimeAnnotation_ForwardConvert(time, this->Shift, this->Scale);


  const ptime epoch = boost::posix_time::from_time_t(0);
  ptime datetime = epoch + seconds(static_cast<long>(time));


  std::tm tm = boost::posix_time::to_tm(datetime);
  int year =  tm.tm_year + 1900; //convert from epoch
  int month =  tm.tm_mon + 1;//conert jan == 0
  int day =   tm.tm_mday; //starts at 1, ok
  int hour = tm.tm_hour; // 0 = midnight, ok
  int min = tm.tm_min; // 0, ok
  int sec = tm.tm_sec; // [0,60] in c++11, ok http://en.cppreference.com/w/cpp/chrono/c/tm


  sprintf(buffer, this->Format, year, month, day,hour,min,sec);

  vtkStringArray* data = vtkStringArray::New();
  data->SetName("Text");
  data->SetNumberOfComponents(1);
  data->InsertNextValue(buffer);
  output->AddColumn(data);
  data->Delete();

  delete[] buffer;
  return 1;
}

//----------------------------------------------------------------------------
void vtkDateTimeAnnotation::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
  os << indent << "Format: " 
    << (this->Format? this->Format : "(none)") << endl;
}


